#' FUNCTIONS TO FIT NON-LINEAR MODELS TO ESTIMATE MAXB

#' Wrapper function to fit non-linear growth model per species.
#'
#' The function will create sensible parameter ranges for A, k and p parameters
#'   of the Chapman-Richards and Logistic growth curves and attempt to run a
#'   forward step-AIC procedure to add covariates to the linear component of the
#'   model (on the A parameter -- the asymptote). The maximum number of covariates
#'   to add is determined by `maxNoCoefs`.
#'
#' @param sp species name -- only used for messaging.
#' @param predictorVarsData a `data.table` of predictor variables including
#'    those in `predictorVars` and `age`, as well as `pixelIndex`.
#'    Note that `age` should be in the original scale (e.g., not logged).
#' @param sppVarsB s `data.table` of species biomass (`B`) and `pixelIndex`.
#' @param predictorVars character vector of predictor variables to be included
#'   in the linear component of the model affecting the asymptote (need to correspond to
#'   `names(predictorVarsData[[sp]])`) the same predictors will be considered for all
#'   species.
#' @param predictorVarsCombos a list of sets of covariates in `predictorVars` to add to
#'    the fitted models. If this list has several entries with sets of covariates,
#'    each will be fitted as part of the model selection process.
#' @param doFwdSelection should covariates be added one at a time to the
#'    linear component of the model? If `TRUE`, and `is.null(predictorVarsCombos)`,
#'    then each entry in `predictorVarsCombos` is used as the set of covariates to
#'    test. Otherwise, `predictorVarsCombos` will be created from combinations of
#'    `predictorVars`, with `maxNoCoefs` determining the maximum number of covariates
#'     to add. If `FALSE` the full model is fitted.
#' @param maxNoCoefs how many covariates from `predictorVars` should be added to
#'    the linear component of the model affecting the asymptote? Note that the
#'    more covariates are added the longer the model takes to fit, as all
#'    combinations are attempted. For 2 or more covariates, only combinations with
#'    "cover" are attempted.
#' @param sampleSize how many data points should be randomly sampled to fit the model?
#'   If `NA` the full dataset will be used. Note that this may result in long
#'   computation times. Biomass data will be binned into 10 regular bins before sampling points in number
#'   equal to `sampleSize`.
#' @param Ntries how many times should the models be fit with new randomly
#'   generated starting values? Only used if `randomStarts == TRUE`.
#' @param maxCover numeric. Value indicating maximum cover/dominance.
#' @param models character vector of models to fit. Only Chapman-Richards ('CR')
#'   and 'Logistic' can be chosen at the moment.
#' @param modelOutputsPrev previous outputs of `fitNLMmodels`. The model will try
#'   refitting and comparing AIC with the last results.
#' @param randomStarts logical. Should random starting values of A, k and p non-linear
#'   parameters be picked from a range sensible values, or should all combinations
#'   of values within this range be used? If FALSE, the default, the starting values
#'   are spaced at regular intervals within an acceptable range for each parameter -- 20
#'   values for A, 10 for k and p -- and all combinations are used (2000 starting
#'   values in total). Parameter ranges are estimated from data following Fekedulegn et al. (1999)
#'   as follows:
#'   \itemize{
#'     \item{range of `A` starting values (B0 parameter in Fekedulegn et al. 1999)
#'     varies between \eqn{ObsMaxB \times 0.3} and \eqn{ObsMaxB \times 0.9}, where ObsMaxB is the
#'     maximum observed B across the full dataset (not the sampled data for fitting)}
#'     \item{`k` (CR model) and `p` (Logistic model; both are B2 parameter in
#'     Fekedulegn et al. 1999) are estimated as a constant rate to get to ObsMaxB,
#'     calculated as \eqn{\frac{\frac{Bobs2 - Bobs1}{age2 - age1}}{ObsMaxB}}, where *B1/2* and
#'     *age1/2* are are observed values at two points in time. We draw 100 samples
#'     of two *age* values, and corresponding *B*, to calculate a sample of rates.
#'     After excluding rates <= 0, we take the minimum and maximum values as
#'     the range of `k` (CR model) or `p` (Logistic model) parameters of the growth model}
#'     \item{the `p` parameter (CR model) (related to B3 parameter in Fekedulegn
#'     et al. 1999) should be > 1. Here we use a range of values between 1.1 and
#'     80 which provided suitable fitting using data from the Northwest Territories,
#'     Canada}
#'     \item{the `k` parameter (Logistic model; B1 parameter in Fekedulegn et al. 1999)
#'     is estimated as \eqn{B_0 = \frac{ObsMaxB}{1 + k}}, using a small positive
#'     number for \eqn{B_0}, e.g. 2. Here, we estimate `k` values for \eqn{B_0}
#'     values 1 to 5, and use the minimum and maximum to determine the range from
#'     where to draw starting values.}
#'   }
#' @param lowerBounds a named vector of lower parameter boundaries. If FALSE, no lower
#'   boundaries are applied. If TRUE, coefficients of the linear model on the A
#'   parameter (intercept, cover, k and p) are bound (intercept = observed maximum B * 0.5,
#'   cover = 0, k = 0.05 and p = 1). Alternatively, pass a named vector of parameter
#'   boundaries.
#' @param upperBounds a named vector of upper parameter boundaries. If FALSE, no lower
#'   boundaries are applied. If TRUE, coefficient of the linear model on the A parameter
#'   (intercept and k) are bound (intercept = observed maximum B * 1.5, k = 0.2). Alternatively, pass
#'   a named vector of parameter boundaries.
#' @param nbWorkers integer. If > 1, the number of workers to use in `future.apply::future_apply`, otherwise
#'   no parallelisation is done.
#'
#' @importFrom crayon blue magenta

fitNLMModels <- function(sp = NULL, predictorVarsData, sppVarsB, predictorVars,
                         predictorVarsCombos = NULL, maxNoCoefs = 4, doFwdSelection = FALSE,
                         sampleSize = 3000, Ntries = 2000, maxCover = 1L, models = c("CR", "Logistic"),
                         modelOutputsPrev = NULL, randomStarts = FALSE, lowerBounds = TRUE,
                         upperBounds = TRUE, nbWorkers = 1L) {

  ## Checks ------------------------

  if (any(!models %in% c("CR", "Logistic"))) {
    stop("only 'CR' and/or 'Logistic' models are accepted")
  }

  if (isFALSE(is(predictorVarsData, "data.table"))) {
    stop("'predictorVarsData' must be a data.table")
  }

  if (isFALSE(is(predictorVars, "character"))) {
    stop("'predictorVars' must be a character vector of column names")
  }

  predsNeeded <- c(predictorVars, "age")
  notFound <- predsNeeded[!predsNeeded %in% colnames(predictorVarsData)]
  if (length(notFound)) {
    stop("The following predictors were not found in predictorVarsData$", sp, ": ",
         paste(notFound, collapse = ", "))
  }

  if (isFALSE(is(sppVarsB, "data.table"))) {
    stop("'sppVarsB' must be a data.table")
  }

  if (isFALSE(is(randomStarts, "logical"))) {
    stop("'randomStarts' must be TRUE/FALSE")
  }

  if (isFALSE(is(lowerBounds, c("vector")))) {
    stop("'lowerBounds' must be a vector")
  }

  if (isFALSE(is(upperBounds, c("vector")))) {
    stop("'upperBounds' must be a vector")
  }

  if (!is.null(modelOutputsPrev)) {
    if (isFALSE(is(modelOutputsPrev, "list")) | is.null(names(modelOutputsPrev))) {
      stop("modelOutputsPrev must be a list names after 'models'")
    }

    missingModels <- setdiff(models, names(modelOutputsPrev))
    if (length(missingModels)) {
      stop(paste("No 'modelOutputsPrev' models found for", missingModels))
    }
    modelOutputsPrev <- modelOutputsPrev[models]
  } else {
    modelOutputsPrev <- sapply(models, function(x) list(), USE.NAMES = TRUE, simplify = FALSE)
  }

  if (!"pixelIndex" %in% names(predictorVarsData)) {
    stop("'predictorVarsData' is missing 'pixelIndex' column")
  }

  if (!"pixelIndex" %in% names(sppVarsB)) {
    stop("'predictorVarsData' is missing 'pixelIndex' column")
  }

  ## Join predictors and response variable -------------------

  ## drop coordinates
  colsPred <- setdiff(names(predictorVarsData), c("x", "y"))
  colsResp <- setdiff(names(sppVarsB), c("x", "y"))
  specDat2 <- predictorVarsData[, ..colsPred][sppVarsB[, ..colsResp], on = c("pixelIndex")]

  if (any(is.na(specDat2)))
    message(magenta("Found NA's in species biomass or maxB model predictors.",
                    "These lines will be removed"))
  specDat2 <- specDat2[complete.cases(specDat2)]
  # specDat2 <- split(specDat, by = "speciesCode")[[1]]
  # specDat2[, age := exp(logAge)]
  # folds <- dismo::kfold(specDat2, k = 10)
  # trainData <- specDat2[folds == 1]
  # testData <- specDat2[folds != 1]

  ## (Re-)sample data -----------------------
  if (is.na(sampleSize) || sampleSize > nrow(specDat2)) {
    sampleIDs <- 1:nrow(specDat2)
  } else {
    ## this is no longer necessary as the stratification should be done elsewhere
    ## stratify sampling by B to ensure good coverage
    # specDat2[, Bbins := cut(B, breaks = 10)]
    # specDat2[, Bsamples := round((.N/nrow(specDat2))*sampleSize), by = Bbins]
    # specDat2[Bsamples < 1, Bsamples := 1]
    # sampleIDs <- specDat2[, sample(.I, unique(Bsamples), replace = FALSE), by = Bbins]$V1
    # set(specDat2, NULL, c("Bbins", "Bsamples"), NULL)

    sampleIDs <- sample(1:nrow(specDat2), size = sampleSize)
  }
  trainData <- specDat2[sampleIDs]

  ## calculate log-age
  trainData[, logAge := log(age)]
  trainData[age == 0, logAge := 0]

  ## Data at origin -----------------------

  ## Add a B(0) age(0) point at the origin as per Fekedulegn et al. (1999)
  if (nrow(trainData[age == 0 & B == 0]) == 0) {
    cols <- names(which(sapply(trainData, is.numeric)))
    newRow <- trainData[, lapply(.SD, mean), .SDcols = cols]
    cols <- setdiff(names(trainData), cols)
    newRowFact <- trainData[, lapply(.SD, function(x)  {
      xTemp <- as.numeric(as.factor(x))
      meanX <- round(mean(xTemp))
      unique(x[xTemp == meanX])
    }), .SDcols = cols]

    newRow <- cbind(newRow, newRowFact)
    newRow[, `:=`(logAge = 0, cover = 0, age = 0, B = 0)]
    trainData <- rbind(trainData, newRow, use.names = TRUE)
  }

  ## Variable combinations to try -----------------------

  if (is.null(predictorVarsCombos)) {
    if (doFwdSelection) {
      predictorVarsCombos <- sapply(1:maxNoCoefs, function(x) {
        combos <- combn(predictorVars, x, simplify = FALSE)
        if (x > 1) {
          toKeep <- sapply(combos, function(x) any(x == "cover"))
          combos <- combos[toKeep]
        }
        # apply(combos, 2, function(x) x, simplify = FALSE)
        combos
      }, simplify = TRUE)

      predictorVarsCombos <- unlist(predictorVarsCombos, recursive = FALSE)
    } else {
      predictorVarsCombos <- list(predictorVars)
    }
  }

  ## Initial parameter values ---------------------------

  ## Fekedulegn et al. (1999) suggest the following:
  ## To estimate A starting values (=B0 parameter), use the maximum observed B
  ## To estimate k (CR) and p (Log) (=B2 parameter) use the following to calculate a constant rate to get to maxB
  ## ((Bobs2 - Bobs1) / (age2 - age1)) / maxB  -- where B and age are observed values at two points in time and maxB is the maximum biomass
  ## p (CR) (related to B3 parameter -- p = 1/(1-B3) and 0 < B3 < 1) should be p > 1
  ## To estimate k (Log) (B1 parameter) use the following:
  ## B(0) = maxB/(1 + k), using a small positive number for B(0), e.g. 2

  ## smaller k lead to reaching asymptote slower
  ## smaller p lead to reaching the asymptote faster

  ## maximum observed biomass across full dataset
  ObsMaxB <- max(specDat2$B, na.rm = TRUE)   ## maximum across whole data for the focal species

  ## this is probably a bad way of estimating rates, because its a space-for-time across a potentially very large area
  if (FALSE) {
    ## initial parameter ranges for k (CR model) and p (Logistic model) parameters:
    ## sample a wide range of ages
    ageSamples <- lapply(1:100, function(x) sort(sample(specDat2[cover > maxCover * 0.8, age], size = 2, replace = FALSE)))
    ageSamples <- lapply(ageSamples, function(x) if (length(unique(x)) > 1) x)

    rateEstimates <- sapply(ageSamples, function(ages) {
      if (!is.null(ages)) {
        Bobs1 <- specDat2[age %in% ages[1] & cover > maxCover * 0.8, B]
        Bobs2 <- specDat2[age %in% ages[2] & cover > maxCover * 0.8, B]

        sampSize <- min(100, length(Bobs1), length(Bobs2))
        Bobs1 <- sample(Bobs1, sampSize)
        Bobs2 <- sample(Bobs2, sampSize)

        ((Bobs2 - Bobs1) / (ages[2] - ages[1])) / ObsMaxB
      }
    })
    rateEstimates <- unlist(rateEstimates)
    rateEstimates <- rateEstimates[rateEstimates > 0]
  }


  ## initial parameter ranges for k (Logistic model) parameter:
  B0s <- seq(1, 5, 0.1)
  kLogEstimates <- (ObsMaxB - B0s)/B0s

  ## Parameter boundaries ---------------------------

  ## parameter bounds -- k and p restricted below, as they depend on the type of model
  vars <- c("A.(Intercept)", paste0("A.", predictorVars), "k", "p")
  if (!is.list(lowerBounds)) {
    lowerLims <- rep(-Inf, length(vars))
    names(lowerLims) <- vars
    if (isTRUE(lowerBounds)) {
      lowerLims["A.(Intercept)"] <- ObsMaxB * 0.5
      lowerLims["A.cover"] <- 0 ## restrict cover effect to positive values
      # lowerLims["k"] <- min(rateEstimates)
      lowerLims["k"] <- 0.05
      lowerLims["p"] <- 1
    }
  } else {
    lowerLims <- lowerBounds
    if (!all(vars %in% names(lowerLims))) {
      ## make defaults and complete provided list/vector necessary
      defaultLowerLims <- rep(-Inf, length(vars))
      names(defaultLowerLims) <- vars
      defaultLowerLims["A.(Intercept)"] <- ObsMaxB * 0.5
      defaultLowerLims["A.cover"] <- 0 ## restrict cover effect to positive values
      # lowerLims["k"] <- min(rateEstimates)
      defaultLowerLims["k"] <- 0.05
      defaultLowerLims["p"] <- 1

      missingBounds <- setdiff(names(defaultLowerLims), names(lowerLims))
      lowerLims <- c(lowerLims, defaultLowerLims[missingBounds])
    }
  }
  if (!is.list(upperBounds)) {
    upperLims <- rep(Inf, length(vars))
    names(upperLims) <- vars
    if (isTRUE(upperBounds)) {
      upperLims["A.(Intercept)"] <- ObsMaxB * 1.5
      # upperLims["k"] <- max(rateEstimates)  ## too high
      upperLims["k"] <- 0.2
    }
  } else {
    upperLims <- upperBounds
    if (!all(vars %in% names(upperLims))) {
      ## make defaults and complete provided list/vector necessary
      defaultUpperLims <- rep(Inf, length(vars))
      defaultUpperLims["A.(Intercept)"] <- ObsMaxB * 1.5
      defaultUpperLims["k"] <- 0.2

      missingBounds <- setdiff(names(defaultUpperLims), names(upperLims))
      upperLims <- c(upperLims, defaultLowerLims[missingBounds])
    }
  }

  ## Model setup ----------------------------------------
  ## linear equation for A, the asymptote (same for CR and Logistic models)
  # lnEqn <- paste("A ~", paste(predictorVars, collapse = " + "))
  if ("CR" %in% models) {
    message(blue("Fitting CR model for", sp))
    modelParams <- list(
      ## dpois represents the stochastic component generating B, for any set of conditions.
      nonLinEqn = quote(B ~ dpois(lambda = A * (1 - exp(-k * age))^p)),
      # nonLinEqn = quote(B ~ dpois(lambda = A * (1 - exp(-0.07 * age))^p)),
      # nonLinEqn = quote(B ~ dexp(rate = A * (1 - exp(-k * age))^p)),
      # linEqn = list(parse(text = lnEqn)),  ## below
      # plim = c(min = 1.1, max = 80),   ## Eliot's
      plim = c(min = 1, max = 20), ## higher values lead to reaching A slower
      Alim = c(min = ObsMaxB * 0.3, max = ObsMaxB * 0.9),  ## Eliot's
      # klim = c(min = 0.0001, max = 0.13)  ## Eliot's
      klim = c(min = 0.05, max = 0.2)    ## higher values lead to reaching A faster
      # klim = c(min = min(rateEstimates), max = max(rateEstimates))
    )

    paramRanges <- modelParams[c("Alim"
                                 , "klim"
                                 , "plim")]

    ## starting values
    if (randomStarts) {
      starts <- data.table(A = runif(Ntries, paramRanges$Alim[["min"]], paramRanges$Alim[["max"]]),
                           k = runif(Ntries, paramRanges$klim[["min"]], paramRanges$klim[["max"]]),
                           p = runif(Ntries, paramRanges$plim[["min"]], paramRanges$plim[["max"]]))
    } else {
      ## instead of drawing randomly, explore parameter space
      message(blue("Using starting values exploring parameter space. 'Ntries' will be ignored"))
      starts <- list(A = seq(paramRanges$Alim[["min"]], paramRanges$Alim[["max"]], length.out = 10),
                     k = seq(paramRanges$klim[["min"]], paramRanges$klim[["max"]], length.out = 10),
                     p = seq(paramRanges$plim[["min"]], paramRanges$plim[["max"]], length.out = 10))
      starts <- as.data.table(expand.grid(starts))
    }

    # this function can be run again and again
    #  to keep trying new start values, i.e., "iterate". To iterate,
    #  don't reset modelOutputsPrev$CR
    for (pred in predictorVarsCombos) {
      lnEqn <- paste("A ~", paste(pred, collapse = " + "))
      # lnEqn <- paste("A ~ 1") ## testing
      modelParams$linEqn <- list(lnEqn)  ## don't parse/eval: mem leak! keep as character - .fitNLMwCovariates converts to formula

      message(blue("using ", lnEqn))

      ## tests with theoretical data
      ## generate B curves for changing age, k and p
      # trainData2 <- expand.grid(list(age = unique(trainData$age),
      #                                A = ObsMaxB,
      #                                k = round(seq(0.01, 1, length.out = 5), 3),
      #                                # p = c(0.01, 0.1, 1, 10)
      #                                p = 10
      #                           ))
      # trainData2 <- as.data.table(trainData2)
      # trainData2[, lambda := A * (1 - exp(-k * age))^p]
      # trainData2[, B := sapply(lambda, function(x) rpois(1, x))]
      # ggplot(trainData2, aes(y = B, x = age, colour = as.factor(k))) +
      #   geom_line(size = 1) +
      #   theme_pubr() +
      #   facet_grid(~ p)
      #
      # ggplot(trainData2, aes(y = B, x = age, colour = as.factor(k))) +
      #   geom_line(size = 1) +
      #   theme_pubr() +
      #   facet_grid(~ p)

      ## generate theoretical data
      # trainData2 <- unique(trainData[, .(age)])
      # trainData2[, lambda := mean(starts$A) * (1 - exp(-0.2 * age))^200]  ## for k = 0.05 and p = 100 the model coefs were very far
      # trainData2[, B := sapply(lambda, function(x) rpois(1, x))]
      # preds <- trainData[, lapply(.SD, mean), .SDcols = pred]
      # trainData2 <- cbind(trainData2, preds)
      # browser()
      modelOutputsPrev$CR$mllsOuter <- .fitNLMwCovariates(data = trainData,
                                                          # data = trainData2,  ## for testing
                                                          maxCover = maxCover,
                                                          nonLinModelQuoted = modelParams$nonLinEqn,
                                                          linModelQuoted = modelParams$linEqn,
                                                          starts = starts,
                                                          lower = lowerLims,
                                                          upper = upperLims,
                                                          mllsOuterPrev = modelOutputsPrev$CR$mllsOuter,
                                                          nbWorkers = nbWorkers)
    }
  }

  if ("Logistic" %in% models) {
    message(blue("Fitting Logistic model for", sp))

    modelParams <- list(
      ## dpois represents the stochastic component generating B, for any set of conditions.
      nonLinEqn = quote(B ~ dpois(lambda = A / (1 + k * exp(-p * age)))),
      # linEqn = list(parse(text = lnEqn)),  ## below
      # plim = c(min = 0.001, max = 1), ## Eliot's
      plim = c(min = min(rateEstimates), max = max(rateEstimates)),
      Alim = c(min = ObsMaxB * 0.3, max = ObsMaxB * 0.9),
      klim <- c(min = min(kLogEstimates), max = max(kLogEstimates))
    )
    # modelParams$klim <- c(min = 10, max = max(modelParams$Alim)) ## Eliot's

    paramRanges <- modelParams[c("plim", "Alim", "klim")]

    ## starting values
    if (randomStarts) {
      starts <- data.table(A = runif(Ntries, paramRanges$Alim[["min"]], paramRanges$Alim[["max"]]),
                           k = runif(Ntries, paramRanges$klim[["min"]], paramRanges$klim[["max"]]),
                           p = runif(Ntries, paramRanges$plim[["min"]], paramRanges$plim[["max"]]))
    } else {
      message(blue("Using starting values exploring parameter space. 'Ntries' will be ignored"))
      ## instead of drawing randomly, explore parameter space
      starts <- list(A = seq(paramRanges$Alim[["min"]], paramRanges$Alim[["max"]], length.out = 20),
                     k = seq(paramRanges$klim[["min"]], paramRanges$klim[["max"]], length.out = 10),
                     p = seq(paramRanges$plim[["min"]], paramRanges$plim[["max"]], length.out = 10))
      starts <- as.data.table(expand.grid(starts))
    }

    ## restrict k and p -- this did not result in good models.
    # lowerLims["k"] <- min(paramRanges$klim) * 0.1
    # upperLims["k"] <- max(paramRanges$klim) * 2
    # lowerLims["p"] <- min(paramRanges$plim) * 0.1
    # upperLims["p"] <- min(paramRanges$plim) * 2


    for (pred in predictorVarsCombos) {
      lnEqn <- paste("A ~", paste(pred, collapse = " + "))
      # lnEqn <- paste("A ~ 1") ## testing
      modelParams$linEqn <- list(parse(text = lnEqn))

      message(blue("using", lnEqn))

      modelOutputsPrev$Logistic$mllsOuter <- .fitNLMwCovariates(data = trainData,
                                                                maxCover = maxCover,
                                                                nonLinModelQuoted = modelParams$nonLinEqn,
                                                                linModelQuoted = modelParams$linEqn,
                                                                starts = starts,
                                                                lower = lowerLims,
                                                                upper = upperLims,
                                                                mllsOuterPrev = modelOutputsPrev$Logistic$mllsOuter,
                                                                nbWorkers = nbWorkers)
    }
  }
  return(modelOutputsPrev)
}


#' Fit non-linear growth model under various starting conditions
#'
#' Uses likelihood parameter estimation to fit non linear models
#' while attempting several starting values.
#'
#' @param data a \code{data.table} or \code{data.frame} with all covariates
#'   and the response variable. Note that incomplete lines are removed.
#' @param nonLinModelQuoted The non-linear equation as a \code{call}
#'   (quoted expression) passed to \code{mle2(minuslog1)}. See \code{?mle}.
#'   Accepts equations with three parameters 'A', 'p' and 'k'.
#' @param linModelQuoted A list of linear equations/modes relating each
#'   parameter ('A', 'p' and 'k') with a set of covariates. A \code{call}
#'   (quoted expression) passed to \code{mle2(..., parameters)}. Note that for the
#'   purpose of tree growth, the linear equation determining 'A' should include a
#'   'cover' predictor indicating the tree cover or dominance in the stand. Should be
#'   scaled between 0 and \code{maxCover}.
#' @param mllsOuterPrev the output of a previous \code{fitNLMwCovariates} run which
#'   is used to extract last best AIC and maximum biomass estimate and judge if new
#'   iterations are better.
#' @param model character. Non-linear model form used to estimate average maximum
#'   biomass. One of "CR" (Chapman-Richards) or "Logistic". In both cases, maximum biomass
#'   is equivalent to the 'A' asymptote parameter, which is estimated using observed mean
#'   values of predictors entering its linear equation and \code{cover == maxCover}, if this
#'   predictor is included (as it should). Passed to \code{extractMaxB}
#' @param maxCover numeric. Value indicating maximum cover/dominance.
#' @param starts `data.table` or `data.frame` of parameter starting values. Will be coerced to named list
#'   with names being parameter names.
#' @param lower passed to `mle2`
#' @param upper passed to `mle2`
#'
#' @return a \code{list} with entries 'mll' (the maximum likelihood-estimated
#' coefficients) and 'AICbest' (the AIC of the best models generating these
#' coefficients)
#'
#' @importFrom crayon cyan red

.fitNLMwCovariates <- function(data, nonLinModelQuoted, linModelQuoted, mllsOuterPrev,
                               model = c("CR", "Logistic"), maxCover = 1L, cores = 1L,
                               starts = NULL, lower = NULL, upper = NULL, nbWorkers = 1L) {
  if (requireNamespace("bbmle", quietly = TRUE)) {
    linModelQuoted <- lapply(linModelQuoted, as.formula, env = .GlobalEnv) # .GlobalEnv keeps object small. don't eval/parse!
    nonLinModelQuoted <- as.formula(nonLinModelQuoted, env = .GlobalEnv)
    data <- data[complete.cases(data),]

    if (!is.null(starts)) {
      if (isFALSE(is(starts, "data.frame"))) {
        stop("starts must be a data.table or data.frame")
      }

      Ntries <- nrow(starts)
    } else {
      Ntries <- 1L
    }

    if (any(grepl("cover", linModelQuoted))) {
      if (!"cover" %in% names(data)) {
        stop("Please provide a 'cover' variable as a linear predictor to, e.g., 'A'")
      }

      if (max(data$cover, na.rm = TRUE) > maxCover |
          min(data$cover, na.rm = TRUE) < 0) {
        stop("data$cover must be scaled between 0 and 'maxCover'")
      }
    }

    if (length(mllsOuterPrev)) {
      aicPrev <- mllsOuterPrev$AICbest
      mllKeep <- mllsOuterPrev$mll

      if (!inherits(mllsOuterPrev$mll, c("error", "try-error"))) {
        newdata <- copy(data)
        newdata[, cover := maxCover]
        maxBEst <- extractMaxB(mllsOuterPrev$mll, newdata = newdata, average = TRUE, model = model)
      } else {
        maxBEst <- NA
      }
      message(cyan("    best AIC so far:", round(aicPrev, 0), "; maxB est = ", maxBEst))
    }

    if (!exists("aicPrev", inherits = FALSE)) {
      aicPrev <- 1e13
    }

    if (nbWorkers > 1) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package MASS not installed. Install using `install.packages('MASS')`.")
      }

      if (!requireNamespace("future", quietly = TRUE)) {
        stop("Package future not installed. Install using `install.packages('future')`.")
      }
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        stop("Package future.apply not installed. Install using `install.packages('future.apply')`.")
      }

      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Package parallel not installed. Install using `install.packages('parallel')`.")
      }

      if (.Platform$OS.type == "unix") {
        cl <- parallel::makeForkCluster(nbWorkers)
      } else {

        if (!requireNamespace("parallelly")) {
          stop("Package parallelly not installed. Install using `install.packages('parallelly')`.")
        }
        cl <- parallelly::makeClusterPSOCK(nbWorkers, rscript_libs = .libPaths())
      }

      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, {
        library("MASS")

        ## limit spawning of additional threads from workers
        data.table::setDTthreads(1)
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
      })

      future::plan(future::cluster, workers = cl)
      applyFUN <- future.apply::future_apply
    } else {
      applyFUN <- apply
    }

    mle2Args <- list(
      "minuslogl" = nonLinModelQuoted
      , "data" = data
      , "parameters" = linModelQuoted
      # , "method" = "L-BFGS-B"
      # , skip.hessian = TRUE
      # , control = list(maxit = 10000)
    )

    if (!is.null(lower)) {
      mle2Args$lower <- lower
    }
    if (!is.null(upper)) {
      mle2Args$upper <- upper
    }

    mll <- applyFUN(
      X = starts,
      MARGIN = 1,
      mle2Args = mle2Args,
      simplify = FALSE,
      FUN = .mllWrapper)

    mllAICs <- sapply(mll, FUN = function(mll) {
      if (inherits(mll, c("error", "try-error"))) {
        return(NA)
      } else {
        return(bbmle::AIC(mll))
      }
    })

    if (any(!is.na(mllAICs))) {
      bestModel <- last(which.min(mllAICs))  ## takes simpler model in case there are more than 1 min(AIC) (first to be fitted)
      aicTry <- mllAICs[bestModel]

      if (aicTry < aicPrev) {
        newdata <- copy(data)
        newdata[, cover := maxCover]

        maxBEst <- extractMaxB(mll[[bestModel]], newdata = newdata, average = TRUE, model = model)
        message(cyan("    best AIC so far:", round(aicTry, 0), "; maxB est = ", maxBEst))
        aicPrev <- aicTry
        mllKeep <- mll[[bestModel]]
      }

    } else {
      warning("Could not converge. Try increasing Ntries or changing paramRanges values.",
              "See output for last error message")
      if (!exists("mllKeep", inherits = FALSE)) {
        mllKeep <- mll[[1]]
      }
    }

    ## old non parallelised code:
    # else {
    #   ## TODO: use plan sequential?
    #   for (ii in 1:Ntries) {
    #     if (ii %% (Ntries/10) == 0) {
    #       message(cyan("    ... at attempt", ii, "of", Ntries))
    #     }
    #
    #     mle2Args <- list("minuslogl" = nonLinModelQuoted
    #                      , "data" = data
    #                      , "parameters" = linModelQuoted
    #                      # , "method" = "L-BFGS-B"
    #                      # , skip.hessian = TRUE
    #                      # , control = list(maxit = 10000)
    #     )
    #     if (!is.null(starts)) {
    #       mle2Args$start <- as.list(starts[ii,])
    #     }
    #
    #     if (!is.null(lower)) {
    #       mle2Args$lower <- lower
    #     }
    #     if (!is.null(upper)) {
    #       mle2Args$upper <- upper
    #     }
    #
    #     mll <- tryCatch({
    #       R.utils::withTimeout({
    #         suppressWarnings(
    #           do.call(mle2, mle2Args)
    #         )
    #       }, timeout = 5)
    #     }, TimeoutException = function(ex) {
    #       NULL
    #     }, error = function(e) e)
    #     ## put new tryCatch outside of TimeoutException function for readability
    #     if (is.null(mll)) {
    #       message(red("Timeout. Trying again without boundaries"))
    #
    #       mle2Args$lower <- NULL
    #       mle2Args$upper <- NULL
    #       mll <- tryCatch({
    #         R.utils::withTimeout({
    #           suppressWarnings(
    #             do.call(mle2, mle2Args)
    #           )
    #         }, timeout = 5)}
    #         , TimeoutException = function(ex) {
    #           stop("Convergence timed-out even without boundaries.")
    #         }, error = function(e) e)
    #     }
    #
    #     if (inherits(mll, c("error", "try-error"))) {
    #       next
    #     }
    #     aicTry <- AIC(mll)
    #
    #     if (aicTry < aicPrev) {
    #       newdata <- copy(data)
    #       newdata[, cover := maxCover]
    #
    #       maxBEst <- extractMaxB(mll, newdata = newdata, average = TRUE, model = model)
    #       message(cyan("    best AIC so far:", round(aicTry, 0), "; maxB est = ", maxBEst))
    #       aicPrev <- aicTry
    #       mllKeep <- mll
    #     }
    #   }
    #   if (!exists("mllKeep", inherits = FALSE)) {   ## may still not have converged after all tries
    #     warning("Could not converge. Try increasing Ntries or changing paramRanges values.",
    #             "See output for last error message")
    #     mllKeep <- mll
    #   }
    # }

    list(mll = mllKeep, AICbest = aicPrev)
  } else {
    stop("Package bbmle not installed. Install using `install.packages('bbmle')`.")
  }
}

.mllWrapper <- function(X, mle2Args) {
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package R.utils not installed. Install using `install.packages('R.utils')`.")
  }
  if (!is.null(X)) {
    mle2Args$start <- as.list(X)
  }
  mll <- tryCatch({
    R.utils::withTimeout({
      suppressWarnings({
        mll <- do.call(bbmle::mle2, mle2Args)
      }
      )
    }, timeout = 5)
  }, TimeoutException = function(ex) {
    stop("Convergence timed-out.")
  }, error = function(e) e)

  # put new tryCatch outside of TimeoutException function for readability
  # if (inherits(mll, c("error", "try-error"))) {
  #   message(red("Timeout. Trying again without boundaries"))
  #
  #   mle2Args$lower <- NULL
  #   mle2Args$upper <- NULL
  #   mll <- tryCatch({
  #     R.utils::withTimeout({
  #       suppressWarnings(
  #         do.call(mle2, mle2Args)
  #       )
  #     }, timeout = 5)}
  #     , TimeoutException = function(ex) {
  #       stop("Convergence timed-out even without boundaries.")
  #     }, error = function(e) e)
  # }
  return(mll)
}


#' Maximum biomass estimator
#'
#' Estimation of maximum biomass as the A parameter in the Chapman-Richards and
#'   Logistic growth equations. Since A is modelled as linear term, it is a matrix
#'   product of its linear coefficients. It is assumed that all coefficients are
#'   related additively.
#'
#' @param mll the output of an \code{bbmle::mle2} call (the fitted non-linear model),
#'   from which coefficient values will be extracted
#' @param newdata data for estimation of 'A'
#' @param average should 'A' be estimated for average values of its predictors.
#' @param model character. Non-linear model form used to estimate average maximum
#'   biomass. One of "CR" (Chapman-Richards) or "Logistic".
#'
#' @export
extractMaxB <- function(mll, newdata, average = FALSE, model = c("CR", "Logistic"))  {
  if (requireNamespace("bbmle", quietly = TRUE)) {
    if (is.null(dim(newdata))) {
      stop("newdata should have at least two dimensions and be coercible to a data.table")
    }

    a <- bbmle::coef(mll)
    paramNames <- .getMaxBCoefs(mll, model = model)

    newdata <- as.data.table(newdata)
    cols <- paramNames[["origCoefNames"]]
    newdata <- newdata[, ..cols]
    newdata <- cbind(rep(1L, nrow(newdata)), newdata)
    setnames(newdata, "V1", grep("Intercept", paramNames[["mllCoefNames"]],value = TRUE))

    if (average) {
      newdata <- suppressWarnings(sapply(newdata, mean))
    }

    preds <- mapply(`*`, x = newdata, y = a[paramNames[["mllCoefNames"]]])
    preds <- if (is.null(dim(preds))) {
      round(sum(preds), 0)
    } else {
      round(rowSums(preds), 0)
    }

    return(preds)
  } else {
    stop("Package bbmle not installed. Install using `install.packages('bbmle')`.")
  }
}

#' Plot estimated maximum biomass by age
#'
#' Plots a maximum biomass estimated at maximum 'cover'
#'   (or dominance) levels as a function of age
#'
#' @param mll a named list with outputs of an \code{bbmle::mle2} call (the fitted non-linear
#'  model), from which coefficient values will be extracted. If several model outputs
#'  are provided all fitted models will be plotted, with plot labels corresponding to list names.
#' @param data data for estimation of maximum biomass. Should contain at least
#'  an 'age' column. Note that other covariates will be averaged and 'cover' values
#'  will be replaced with the maximum cover value (\code{maxCover}). If `mll` is a list
#'  data is assumed to be the same for the two models.
#' @param maxCover numeric. Value indicating maximum cover/dominance.
#' @param xCovar the variable shown in the x axis. Defaults to `age`.
#' @param plotTitle character. Passed to ggplot2::labs(title)
#' @param nonLinModelQuoted a named list of non-linear equations as a \code{call}
#'   (quoted expression) passed to \code{mle2(minuslog1)}. See \code{?mle}.
#'   Accepts equations with three parameters 'A', 'p' and 'k'. List names and
#'   length must the same as in `mll`.
#' @param linModelQuoted A named list of lists of linear equations/modes relating each
#'   parameter ('A', 'p' and 'k') with a set of covariates. A \code{call}
#'   (quoted expression) passed to \code{mle2(..., parameters)}. Note that for the
#'   purpose of tree growth, the linear equation determining 'A' should include a
#'   'cover' predictor indicating the tree cover or dominance in the stand. Should be
#'   scaled between 0 and \code{maxCover}. List names and length must the same as in `mll`.
#' @param averageCovariates should covariates other than age/cover be averaged for
#'   biomass predictions? If not, for each age (at maximum cover) there will be as
#'   many predictions as other covariate values. If `observedAge == TRUE` and
#'   `averageCovariates == FALSE` then the original data is used, with cover
#'   changed to `maxCover`.
#' @param observedAge should observed age values be used, or should these be generated
#'   as `round(seq(min(age), max(age)*1.5, length.out = 100), 0)`? If `observedAge == TRUE` and
#'   `averageCovariates == FALSE` then the original data is used, with cover
#'   changed to `maxCover`.
#' @param plotCIs should confidence intervals be calculated and plotted?
#'
#' @importFrom ggplot2 ggplot geom_point labs geom_line geom_ribbon theme_classic scale_color_distiller
#'
ggplotMLL_maxB <- function(mll, data, maxCover = 1L, xCovar = "age",
                           plotTitle = NULL, nonLinModelQuoted, linModelQuoted,
                           averageCovariates = TRUE, observedAge = FALSE, plotCIs = TRUE) {
  ## checks
  if (!(is.list(mll) | is.null(names(mll)))) {
    stop("mll must be a named list")
  }

  if (any(!vapply(mll, FUN = function(x) is(x, "mle2"),
                  FUN.VALUE = logical(1)))) {
    stop("mll must be a list of mle2 outputs or an mle2 output")
  }

  if (is.list(nonLinModelQuoted) & length(nonLinModelQuoted) == length(mll)) {
    if (any(sort(names(nonLinModelQuoted)) != sort(names(mll)))) {
      stop("nonLinModelQuoted names must be identical to mll names")
    }
  } else {
    stop("nonLinModelQuoted must be a list of the same length as mll")
  }

  if (is.list(linModelQuoted) & length(linModelQuoted) == length(mll)) {
    if (any(sort(names(linModelQuoted)) != sort(names(mll)))) {
      stop("linModelQuoted names must be identical to mll names")
    }
  } else {
    stop("linModelQuoted must be a list of the same length as mll")
  }

  ## make sure lists follow the same order
  nonLinModelQuoted <- nonLinModelQuoted[names(mll)]
  linModelQuoted <- linModelQuoted[names(mll)]
  outs <- Map(f = .MLLMaxBplotData,
              mll = mll,
              nonLinModelQuoted = nonLinModelQuoted,
              linModelQuoted = linModelQuoted,
              MoreArgs = list(data = data,
                              maxCover = maxCover,
                              averageCovariates = averageCovariates,
                              observedAge = observedAge,
                              plotCIs = plotCIs))

  allPlotData <- rbindlist(lapply(outs, function(out) out$plotData), idcol = "model")
  allConfInts <- rbindlist(lapply(outs, function(out) out$confIntervals), idcol = "model")
  ## if not averaging covariates then we need to have as many upper/lower values as
  ## combos of covariate values. take the min(lower) and max(lower)

  ## extract average predicted maxB across ages
  maxBmean <- Map(f = function(mll, df) extractMaxB(mll, df, average = TRUE),
                  mll = mll,
                  df = lapply(outs, function(out) out$plotData))

  ## extract quantiles of maxB for maximum age used in plots (the asymptotes of each quantile)
  maxBfittedQuants <- allPlotData[age == max(age),
                                  list(value = quantile(pred1, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)),
                                       quant = c("5%", "25%", "50%", "75%", "95%")),
                                  by = model]

  subtitle <- paste(paste("average maxB =", maxBmean, paste0("(", names(maxBmean), ")")),
                    collapse = "\n")
  gg1 <- ggplot() +
    geom_point(data = data, aes(x = get(xCovar), y = B, color = cover)) +
    scale_color_distiller(palette = "Greens", direction = 1) +
    theme_classic() +
    labs(title = plotTitle, subtitle = subtitle, x = xCovar)

  if (nrow(allConfInts) & averageCovariates) {
    gg1 <- gg1 +
      geom_ribbon(data = allConfInts,
                  aes(x = get(xCovar), ymin = lower, ymax = upper),
                  fill = "grey", alpha = 0.5)
  }

  lineCols <- RColorBrewer::brewer.pal(5, "Blues")
  gg1 <- gg1 +
    stat_summary(data = allPlotData,
                 mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                 fun = mean, geom = "line",
                 linewidth = 1, colour = lineCols[3]) +
    stat_summary(data = allPlotData,
                 mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                 fun = quantile, fun.args = list(probs = 0.05),
                 geom = "line", linewidth = 1, colour = lineCols[1]) +
    stat_summary(data = allPlotData,
                 mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                 fun = quantile, fun.args = list(probs = 0.25),
                 geom = "line", linewidth = 1, colour = lineCols[2]) +
    stat_summary(data = allPlotData,
                 mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                 fun = quantile, fun.args = list(probs = 0.75),
                 geom = "line", linewidth = 1, colour = lineCols[4]) +
    stat_summary(data = allPlotData,
                 mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                 fun = quantile, fun.args = list(probs = 0.95),
                 geom = "line", linewidth = 1, colour = lineCols[5])

  gg1 <- gg1 +
    geom_hline(yintercept = maxBfittedQuants[quant == "5%", value], colour = lineCols[1],
               linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = maxBfittedQuants[quant == "25%", value], colour = lineCols[2],
               linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = maxBfittedQuants[quant == "50%", value], colour = lineCols[3],
               linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = maxBfittedQuants[quant == "75%", value], colour = lineCols[4],
               linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = maxBfittedQuants[quant == "95%", value], colour = lineCols[5],
               linetype = "dashed", linewidth = 1)


  gg1
}

#' Prepare data for model plotting
#'
#' @param mll outputs of an \code{bbmle::mle2} call (the fitted non-linear
#'  model), from which coefficient values will be extracted.
#' @param data data for estimation of maximum biomass. Should contain at least
#'  an 'age' column. Note that other covariates will be averaged and 'cover' values
#'  will be replaced with the maximum cover value (\code{maxCover}).
#' @param maxCover numeric. Value indicating maximum cover/dominance.
#' @param nonLinModelQuoted The non-linear equation as a \code{call}
#'   (quoted expression) passed to \code{mle2(minuslog1)}. See \code{?mle}.
#'   Accepts equations with three parameters 'A', 'p' and 'k'.
#' @param linModelQuoted A list of linear equations/modes relating each
#'   parameter ('A', 'p' and 'k') with a set of covariates. A \code{call}
#'   (quoted expression) passed to \code{mle2(..., parameters)}. Note that for the
#'   purpose of tree growth, the linear equation determining 'A' should include a
#'   'cover' predictor indicating the tree cover or dominance in the stand. Should be
#'   scaled between 0 and \code{maxCover}.
#' @param averageCovariates should covariates other than age/cover be averaged for
#'   biomass predictions? If not, for each age (at maximum cover) there will be as
#'   many predictions as other covariate values. If `observedAge == TRUE` and
#'   `averageCovariates == FALSE` then the original data is used, with cover
#'   changed to `maxCover`.
#' @param observedAge should observed age values be used, or should these be generated
#'   as `round(seq(min(age), max(age)*1.5, length.out = 100), 0)`? If `observedAge == TRUE` and
#'   `averageCovariates == FALSE` then the original data is used, with cover
#'   changed to `maxCover`.
#' @param plotCIs should confidence intervals be calculated and plotted?
#'
#' @importFrom data.table data.table as.data.table

.MLLMaxBplotData <- function(mll, nonLinModelQuoted, linModelQuoted, maxCover, data,
                             averageCovariates = TRUE, observedAge = FALSE, plotCIs = TRUE) {
  if (requireNamespace("bbmle", quietly = TRUE)) {
    linModelQuoted <- lapply(linModelQuoted, eval)
    nonLinModelQuoted <- eval(nonLinModelQuoted)

    cols <- c("age", .getMaxBCoefs(mll)[[2]])
    missingCols <- setdiff(c("B", cols), names(data))
    if (length(missingCols)) {
      stop("The following colums were not found in data: ",
           paste(missingCols, collapse = ", "))
    }

    df <- data[, ..cols]

    if (!observedAge) {
      ageVals <- round(seq(min(df$age), max(df$age)*1.5, length.out = 100), 0)
      ageCover <- data.table(cover = maxCover, age = ageVals)
    }

    if (averageCovariates) {
      if (exists("ageCover")) {
        df <- df[, as.list(colMeans(.SD)), .SDcols = cols]
        cols <- setdiff(cols, c("cover", "age"))   ## remove cover/age for binding
        df <- cbind(ageCover, df[, ..cols])
      } else {
        cols <- setdiff(cols, c("cover", "age"))   ## remove cover/age for averaging others
        df[, (cols) := as.list(colMeans(.SD)), .SDcols = cols]
      }
    } else {
      ## if ageCover does not exist, then we take the original data with maxCover
      if (exists("ageCover")) {
        ## if not averaging other covariates, then they need to be repeated for each
        ## age/cover combination
        cols <- setdiff(cols, c("cover", "age"))
        df <- apply(ageCover, 1, function(x) {
          cbind(data.table(cover = x["cover"], age = x["age"]),
                df[, ..cols])
        })
        df <- rbindlist(df)
        df <- df[, lapply(.SD, round, digits = 2)]   ## round to reduce size
        df <- unique(df)
      } else {
        df[, cover := maxCover]
      }
    }

    mll@call$parameters <- eval(mll@call$parameters)
    mll@call$minuslogl <- eval(mll@call$minuslogl)

    ## get biomass predictions for maximum cover at each age, given other covariates' values
    df[, pred1 := bbmle::predict(mll, newdata = df)]

    ## estimate 95% population prediction intervals
    if (all(is.na(vcov(mll)))) {
      plotCIs <- FALSE
    }
    if (plotCIs) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package MASS not installed. Install using `install.packages('MASS')`.")
      }

      if (!requireNamespace("R.utils", quietly = TRUE)) {
        stop("Package R.utils not installed. Install using `install.packages('R.utils')`.")
      }
      ## generate new parameter values from a multivariate normal distribution
      ## see https://stats.stackexchange.com/questions/221426/95-confidence-intervals-on-prediction-of-censored-binomial-model-estimated-usin
      ## and Bolker's book Ecological Models and Data in R (here the solve(hessian)) is suggested as an alternative if using optim)
      newparams <- tryCatch({
        R.utils::withTimeout({
          MASS::mvrnorm(10000, mu = bbmle::coef(mll), Sigma = solve(mll@details$hessian))
        }, timeout = 300) ## 5min
      }, TimeoutException = function(ex) {
        stop("Convergence timed-out.")
      }, error = function(e) e) ## may not work if matrix isn't postive definite

      if (is(newparams, "error") & any(grepl("not positive definite", newparams))) {
        if (!requireNamespace("lqmm", quietly = TRUE)) {
          stop("Package lqmm not installed. Install using `install.packages('lqmm')`.")
        }

        ## try coercing to a positive definite matrix
        ## see https://stackoverflow.com/questions/69804459/multivariete-distribution-error-sigma-is-not-positive-definite
        ## note that this may be due to some variables being linear combinations of others, or bad model specification
        # https://stats.stackexchange.com/questions/30465/what-does-a-non-positive-definite-covariance-matrix-tell-me-about-my-data
        SigmaMatrix <- tryCatch(lqmm::make.positive.definite(solve(mll@details$hessian)), error = function(e) e)
        if (isFALSE(is(SigmaMatrix, "error"))) {
          if (!requireNamespace("MASS", quietly = TRUE)) {
            stop("Package MASS not installed. Install using `install.packages('MASS')`.")
          }
          newparams <- tryCatch({
            R.utils::withTimeout({
              MASS::mvrnorm(10000, mu = bbmle::coef(mll), Sigma = SigmaMatrix)
            }, timeout = 300) ## 5min
          }, TimeoutException = function(ex) {
            stop("Convergence timed-out.")
          }, error = function(e) e)
        }
      }

      if (is(newparams, "error")) {
        confIntervals <- NULL
      } else {
        preds <- apply(newparams, 1, function(x) bbmle::predict(mll, newdata = df, newparams = unlist(x))) ## predict biomass values for each age and new parameter combination
        confIntervals <- apply(preds, 1, function(x) quantile(x, c(0.025, 0.975))) ## estimate 95% quantiles for each iteration of age
        confIntervals <- as.data.table(t(confIntervals))
        confIntervals[, age := df$age]
        setnames(confIntervals, c("2.5%", "97.5%"), c("lower", "upper"))
      }
    } else {
      confIntervals <- NULL
    }
    return(list(plotData = df, confIntervals = confIntervals))
  } else {
    stop("Package bbmle not installed. Install using `install.packages('bbmle')`.")
  }
}


#' Get maximum biomass coefficient names
#'
#' Extracts the names of linear coefficients
#'   for the maximum-biomass-equivalent parameter in
#'   the non-linear growth equations
#'
#' @param mll the output of an \code{bbmle::mle2} call (the fitted non-linear model),
#'   from which coefficient values will be extracted
#' @param model character. Non-linear model form used to estimate average maximum
#'   biomass. One of "CR" (Chapman-Richards) or "Logistic".
#'
#' @return a list of two vectors of parameter names one following
#'   coefficient names in \code{mll} ('mllCoefNames'), the other using the original
#'   names as in the data used for model fitting ('origCoefNames')
#'
.getMaxBCoefs <- function(mll, model = c("CR", "Logistic")) {
  if (requireNamespace("bbmle", quietly = TRUE)) {
    model <- match.arg(model)
    a <- bbmle::coef(mll)
    paramNames <- grep("A\\.", names(a), value = TRUE)
    paramNames2 <- grep("Intercept", paramNames, invert = TRUE, value = TRUE)
    paramNames2 <- sub("(.*\\.)(.*)", "\\2", paramNames2)

    list(mllCoefNames = paramNames, origCoefNames = paramNames2)
  } else {
    stop("Package bbmle not installed. Install using `install.packages('bbmle')`.")
  }
}


#' Partial effect plots of maximum biomass estimates by age
#'
#' Plots the maximum biomass estimates along an age gradient as a function
#'  of abother target covariate, with all others held at average values.
#'
#' @param mll a named list with outputs of an \code{bbmle::mle2} call (the fitted non-linear
#'  model), from which coefficient values will be extracted. If several model outputs
#'  are provided all fitted models will be plotted, with plot labels corresponding to list names.
#' @param data data for estimation of maximum biomass. Should contain at least
#'  an 'age' column. Note that other covariates will be averaged and 'cover' values
#'  will be replaced with the maximum cover value (\code{maxCover}). If `mll` is a list
#'  data is assumed to be the same for the two models.
#' @param targetCovar the covariate for which variation in maxB values will be shown per
#'  age value. Defaults to showing how maxB values change with "cover" at any given age.
#'  Age values are generated as `round(seq(min(age), max(age)*1.5, length.out = 100), 0)`.
#'  When `targetCovar != "cover"`, "cover" may be fixed at `maxCover`. See `fixMaxCover`.
#' @param fixMaxCover logical. If `TRUE` and `targetCovar != "cover"`, cover is
#'  not averaged and is fixed to `maxCover`.
#' @param xCovar the variable shown in the x axis. Defaults to "age". When `xCovar == "age"`
#'  the output plots are not true "Partial Effects" plots. Instead, they show the variation
#'  is B values across values of `targetCovar` for each value of age.
#' @param showQuantiles controls whether quantile predictions will be shown. If
#'  "allQuantiles", quantile values (5%, 25%, 50%, 75%, 95%) of B will be
#'  plotted as blue lines, with their respective asymptote values (quantile values
#'  at maximum age) as dashed lines. If "maximum" the 100% quantile will be plotted.
#'  If "none", only the 50% quantile line (average prediction) is plotted. Quantiles
#'  are always calculated at `max(age)`.
#' @param maxCover numeric. Value indicating maximum cover/dominance.
#' @param plotTitle character. Passed to ggplot2::labs(title)
#' @param nonLinModelQuoted a named list of non-linear equations as a \code{call}
#'   (quoted expression) passed to \code{mle2(minuslog1)}. See \code{?mle}.
#'   Accepts equations with three parameters 'A', 'p' and 'k'. List names and
#'   length must the same as in `mll`.
#' @param linModelQuoted A named list of lists of linear equations/modes relating each
#'   parameter ('A', 'p' and 'k') with a set of covariates. A \code{call}
#'   (quoted expression) passed to \code{mle2(..., parameters)}. Note that for the
#'   purpose of tree growth, the linear equation determining 'A' should include a
#'   'cover' predictor indicating the tree cover or dominance in the stand. Should be
#'   scaled between 0 and \code{maxCover}. List names and length must the same as in `mll`.
#' @param fun passed to `.MLLMaxBPartialPlotData`.
#' @param plotCIs should confidence intervals be calculated and plotted?
#'
#' @details Note that the original data, not the predicted values is shown.
#'
#' @importFrom ggplot2 ggplot geom_point labs geom_line geom_ribbon theme_classic scale_color_distiller
#'
partialggplotMLL_maxB <- function(mll, data, targetCovar = "cover", maxCover = 1, fixMaxCover = TRUE,
                                  xCovar = "age", showQuantiles = "allQuantiles", plotTitle = NULL,
                                  nonLinModelQuoted, linModelQuoted, fun = "mean", plotCIs = TRUE) {
  ## checks
  if (!(is.list(mll) | is.null(names(mll)))) {
    stop("mll must be a named list")
  }

  if (any(!vapply(mll, FUN = function(x) is(x, "mle2"),
                  FUN.VALUE = logical(1)))) {
    stop("mll must be a list of mle2 outputs or an mle2 output")
  }

  if (is.list(nonLinModelQuoted) & length(nonLinModelQuoted) == length(mll)) {
    if (any(sort(names(nonLinModelQuoted)) != sort(names(mll)))) {
      stop("nonLinModelQuoted names must be identical to mll names")
    }
  } else {
    stop("nonLinModelQuoted must be a list of the same length as mll")
  }

  if (is.list(linModelQuoted) & length(linModelQuoted) == length(mll)) {
    if (any(sort(names(linModelQuoted)) != sort(names(mll)))) {
      stop("linModelQuoted names must be identical to mll names")
    }
  } else {
    stop("linModelQuoted must be a list of the same length as mll")
  }

  ## make sure lists follow the same order
  nonLinModelQuoted <- nonLinModelQuoted[names(mll)]
  linModelQuoted <- linModelQuoted[names(mll)]
  outs <- Map(f = .MLLMaxBPartialPlotData,
              mll = mll,
              nonLinModelQuoted = nonLinModelQuoted,
              linModelQuoted = linModelQuoted,
              MoreArgs = list(data = data,
                              targetCovar = targetCovar,
                              maxCover = maxCover,
                              fun = fun,
                              fixMaxCover = fixMaxCover,
                              plotCIs = plotCIs))

  allPlotData <- rbindlist(lapply(outs, function(out) out$plotData), idcol = "model")
  allConfInts <- rbindlist(lapply(outs, function(out) out$confIntervals), idcol = "model")

  if (showQuantiles != "none") {
    ## extract quantiles of maxB for maximum age used in plots (the asymptotes of each quantile)
    maxBfittedQuants <- allPlotData[age == max(age),
                                    list(value = quantile(pred1, probs = c(0.05, 0.25, 0.5, 0.75, 0.95, 1.0)),
                                         quant = c("5%", "25%", "50%", "75%", "95%", "100%")),
                                    by = model]
  }

  gg1 <- ggplot() +
    geom_point(data = data,
               aes(x = get(xCovar), y = B, color = get(targetCovar))) +
    scale_color_distiller(palette = "Greens", direction = 1) +
    theme_classic() +
    labs(title = plotTitle, color = targetCovar, x = xCovar)

  if (nrow(allConfInts)) {
    gg1 <- gg1 +
      geom_ribbon(data = allConfInts,
                  aes(x = get(xCovar), ymin = lower, ymax = upper),
                  fill = "grey", alpha = 0.5)
  }

  lineCols <- RColorBrewer::brewer.pal(5, "Blues")

  if (showQuantiles == "none") {
    gg1 <- gg1 +
      stat_summary(data = allPlotData,
                   mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                   fun = mean, geom = "line",
                   linewidth = 1, colour = lineCols[3])
  }
  if (showQuantiles == "allQuantiles") {
    gg1 <- gg1 +
      stat_summary(data = allPlotData,
                   mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                   fun = quantile, fun.args = list(probs = 0.05),
                   geom = "line", linewidth = 1, colour = lineCols[1]) +
      stat_summary(data = allPlotData,
                   mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                   fun = quantile, fun.args = list(probs = 0.25),
                   geom = "line", linewidth = 1, colour = lineCols[2]) +
      stat_summary(data = allPlotData,
                   mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                   fun = quantile, fun.args = list(probs = 0.75),
                   geom = "line", linewidth = 1, colour = lineCols[4]) +
      stat_summary(data = allPlotData,
                   mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                   fun = quantile, fun.args = list(probs = 0.95),
                   geom = "line", linewidth = 1, colour = lineCols[5]) +
      geom_hline(yintercept = maxBfittedQuants[quant == "5%", value], colour = lineCols[1],
                 linetype = "dashed", linewidth = 1) +
      geom_hline(yintercept = maxBfittedQuants[quant == "25%", value], colour = lineCols[2],
                 linetype = "dashed", linewidth = 1) +
      geom_hline(yintercept = maxBfittedQuants[quant == "50%", value], colour = lineCols[3],
                 linetype = "dashed", linewidth = 1) +
      geom_hline(yintercept = maxBfittedQuants[quant == "75%", value], colour = lineCols[4],
                 linetype = "dashed", linewidth = 1) +
      geom_hline(yintercept = maxBfittedQuants[quant == "95%", value], colour = lineCols[5],
                 linetype = "dashed", linewidth = 1)
  }

  if (showQuantiles == "maximum") {
    gg1 <- gg1 +
      stat_summary(data = allPlotData,
                   mapping = aes(x = get(xCovar), y = pred1, linetype = model),
                   fun = quantile, fun.args = list(probs = 1.0),
                   geom = "line", linewidth = 1, colour = lineCols[5]) +
      geom_hline(yintercept = maxBfittedQuants[quant == "100%", value], colour = lineCols[5],
                 linetype = "dashed", linewidth = 1)
  }

  gg1
}

#' Prepare data for model plotting
#'
#' @param mll outputs of an \code{bbmle::mle2} call (the fitted non-linear
#'  model), from which coefficient values will be extracted.
#' @param data data for estimation of maximum biomass. Should contain at least
#'  an 'age' column. Note that other covariates will be averaged and 'cover' values
#'  will be replaced with the maximum cover value (\code{maxCover}).
#' @param targetCovar the covariate for which variation in maxB values will be shown.
#'  Defaults to showing how maxB values change with "cover". All other covariates except
#'  "age" are averaged. Age values are generated as `round(seq(min(age), max(age)*1.5, length.out = 100), 0)`.
#'  When `targetCovar != "cover"`, "cover" will be fixed at `maxCover`. See `fixMaxCover`.
#' @param fixMaxCover logical. If `TRUE` and `targetCovar != "cover"`, cover is
#'  not averaged and is fixed to `maxCover`.
#' @param maxCover numeric. Value indicating maximum cover/dominance.
#' @param nonLinModelQuoted The non-linear equation as a \code{call}
#'   (quoted expression) passed to \code{mle2(minuslog1)}. See \code{?mle}.
#'   Accepts equations with three parameters 'A', 'p' and 'k'.
#' @param linModelQuoted A list of linear equations/modes relating each
#'   parameter ('A', 'p' and 'k') with a set of covariates. A \code{call}
#'   (quoted expression) passed to \code{mle2(..., parameters)}. Note that for the
#'   purpose of tree growth, the linear equation determining 'A' should include a
#'   'cover' predictor indicating the tree cover or dominance in the stand. Should be
#'   scaled between 0 and \code{maxCover}.
#' @param fun The function to apply when summarizing other variables. By default,
#'   the all other variables except age are  averaged (`"mean"`). Other options are:
#'   `"median"`, `"min"`, `"max"`.
#' @param plotCIs should confidence intervals be calculated and plotted?
#'
#' @importFrom data.table data.table as.data.table
#' @importFrom crayon magenta blue

.MLLMaxBPartialPlotData <- function(mll, nonLinModelQuoted, linModelQuoted,
                                    targetCovar = "cover", fixMaxCover = TRUE,
                                    maxCover = 1, data, fun = "mean", plotCIs = TRUE) {
  if (requireNamespace("bbmle", quietly = TRUE)) {
    linModelQuoted <- lapply(linModelQuoted, eval)
    nonLinModelQuoted <- eval(nonLinModelQuoted)

    cols <- unique(c("age", targetCovar, .getMaxBCoefs(mll)[[2]]))
    missingCols <- setdiff(c("B", cols), names(data))
    if (length(missingCols)) {
      stop("The following colums were not found in data: ",
           paste(missingCols, collapse = ", "))
    }

    df <- data[, ..cols]

    ## generate regular gradients of age and targetCovar values and combine them
    ageTargetCovar <- expand.grid(age = round(seq(min(df$age), max(df$age)*1.5, length.out = 100), 0),
                                  var = round(seq(min(df[, ..targetCovar]), max(df[, ..targetCovar]), length.out = 100), 2))
    ageTargetCovar <- as.data.table(ageTargetCovar)
    setnames(ageTargetCovar, new = c("age", targetCovar))

    ## average other variables
    cols <- setdiff(cols, c("age", targetCovar))
    fun <- get(fun)
    df <- df[, lapply(.SD, fun), .SDcols = cols]
    df <- cbind(ageTargetCovar, df)

    if (targetCovar != "cover" & fixMaxCover) {
      df[, cover := maxCover]
    }

    mll@call$parameters <- eval(mll@call$parameters)
    mll@call$minuslogl <- eval(mll@call$minuslogl)

    ## get biomass predictions for maximum cover at each age, given other covariates' values
    df[, pred1 := bbmle::predict(mll, newdata = df)]

    if (nrow(df[is.na(pred1)])) {
      warning("The prediction returned NA's or NaNs")
    }

    ## estimate 95% population prediction intervals
    if (all(is.na(vcov(mll)))) {
      message(magenta("Cannot estimate confidence intervals"))
      plotCIs <- FALSE
    }
    if (plotCIs) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package MASS not installed. Install using `install.packages('MASS')`.")
      }

      if (!requireNamespace("R.utils", quietly = TRUE)) {
        stop("Package R.utils not installed. Install using `install.packages('R.utils')`.")
      }
      ## generate new parameter values from a multivariate normal distribution
      ## see https://stats.stackexchange.com/questions/221426/95-confidence-intervals-on-prediction-of-censored-binomial-model-estimated-usin
      ## and Bolker's book Ecological Models and Data in R (here the solve(hessian)) is suggested as an alternative if using optim)
      newparams <- tryCatch({
        R.utils::withTimeout({
          MASS::mvrnorm(10000, mu = bbmle::coef(mll), Sigma = solve(mll@details$hessian))
        }, timeout = 300) ## 5min
      }, TimeoutException = function(ex) {
        stop("Convergence timed-out.")
      }, error = function(e) e) ## may not work if matrix isn't postive definite

      if (is(newparams, "error") & any(grepl("not positive definite", newparams))) {
        if (!requireNamespace("lqmm", quietly = TRUE)) {
          stop("Package lqmm not installed. Install using `install.packages('lqmm')`.")
        }
        ## try coercing to a positive definite matrix
        ## see https://stackoverflow.com/questions/69804459/multivariete-distribution-error-sigma-is-not-positive-definite
        ## note that this may be due to some variables being linear combinations of others, or bad model specification
        # https://stats.stackexchange.com/questions/30465/what-does-a-non-positive-definite-covariance-matrix-tell-me-about-my-data
        SigmaMatrix <- tryCatch({
          R.utils::withTimeout({
            lqmm::make.positive.definite(solve(mll@details$hessian))
          }, timeout = 300) ## 5min
        }, TimeoutException = function(ex) {
          stop("Convergence timed-out.")
        }, error = function(e) e)

        if (isFALSE(is(SigmaMatrix, "error"))) {
          if (!requireNamespace("MASS", quietly = TRUE)) {
            stop("Package MASS not installed. Install using `install.packages('MASS')`.")
          }
          newparams <- tryCatch({
            R.utils::withTimeout({
              MASS::mvrnorm(10000, mu = bbmle::coef(mll), Sigma = SigmaMatrix)
            }, timeout = 300) ## 5min
          }, TimeoutException = function(ex) {
            stop("Convergence timed-out.")
          }, error = function(e) e)
        }
      }

      if (is(newparams, "error")) {
        message(magenta("Cannot estimate confidence intervals"))
        confIntervals <- NULL
      } else {
        message(blue("Estimating confidence intervals..."))
        preds <- apply(newparams, 1, function(x) bbmle::predict(mll, newdata = df, newparams = unlist(x))) ## predict biomass values for each age and new parameter combination
        confIntervals <- apply(preds, 1, function(x) quantile(x, c(0.025, 0.975))) ## estimate 95% quantiles for each iteration of age
        confIntervals <- as.data.table(t(confIntervals))
        confIntervals[, age := df$age]
        setnames(confIntervals, c("2.5%", "97.5%"), c("lower", "upper"))
      }
    } else {
      confIntervals <- NULL
    }
    return(list(plotData = df, confIntervals = confIntervals))
  } else {
    stop("Package bbmle not installed. Install using `install.packages('bbmle')`.")
  }
}



## OLD CODE FROM ELIOT:
# library(bbmle)
# specDat[, age := exp(logAge)]
# ll <- split(specDat, by = c("speciesCode"))
#
# models <- list()
# nlsoutInner <- list()
#
# # Chapman Richards
# # https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs092/gtr_srs092-068-coble.pdf
# models <- list(
#   nonLinEqn = quote(B ~ dpois(A * (1 - exp(-k * age))^p)),
#   linEqn = quote(A ~ cover + MAT),  ## Ceres: A is the asymptote
#   plim = c(min = 1, max = 80),
#   Alim = c(max(specDat$B) * c(min = 0.3, max = 0.9)), ## (specDat in B_borealDP)
#   klim = c(min = 0.0001, max = 0.13))
#
# # models <- list(
# #   nonLinEqn = "B ~ dpois(A / (1 + k * exp(-p * age)))",
# #   plim = c(min = 0.001, max = 1),
# #   Alim = c(max(specDat$B) * c(min = 0.3, max = 0.9)) )
# # models$klim <- c(min = 10, max = max(models$Alim))
#
#
# # Two levels:
# # 1. Model type -> e.g., logistic, chapman richards
# # 2. Species
# mllsOuter <- list() # this is here as a placeholder; this function can be run again and again
# #  to keep trying new start values, i.e., "interate". To iterate,
# #  don't reset this object
#
#
# mllsOuter <- fitNLMwCovariates(ll = ll,
#                                nonLinModelQuoted = models$nonLinEqn,
#                                linModelQuoted = models$linEqn, paramRanges = models[3:5],
#                                modelType = "CR", Ntries = 300, mllsOuterPrev = mllsOuter)
#
# maxBs <- extractMaxB(mllsOuter$mlls)
# ggs <- ggplotMLLs(mllsOuter$mlls, ll, linModelQuoted = models$linEqn,
#                   nonLinModelQuoted = models$nonLinEqn)
#
# library(ggpubr)
# ggarrange(plotlist = ggs)

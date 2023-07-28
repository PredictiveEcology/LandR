## BIOMOD FUNCTION WRAPPERS USED TO ESTIMATE SEP
## for caching, looping and parallelizing

#' BIOMOD_Modeling wrapper
#' @param sp character of species name to subset `responseVarData` table
#' @param responseVar character vector, or list of character vectors, of response
#'   variable (column) to subset `responseVarData` table. If a list, it should
#'   be named according to `sp`, for subsetting.
#' @param responseVarData a data.table or list of data.tables with response data
#'   (species presences/absence data). If a list, it should be named according to
#'    `sp`, for subsetting. If the response variable is not binary, values > 0 will be
#'   converted to 1s.
#' @param predictorVars character vector, or list of character vectors, of environmental
#'   variables (columns) to subset `predictorVarsData` table. If a list, it should
#'   be named according to `sp`, for subsetting.
#' @param predictorVarsData environmental data.
#' @param responseVarData a data.table or list of data.tables with environmental data.
#'   If a list, it should be named according to `sp`, for subsetting.#'
#' @param dir.name passed to `biomod2::BIOMOD_FormatingData`
#' @param BIOMOD_ModelingArgs a named list of arguments passed to `biomod2::BIOMOD_Modeling`
#' @param ... further arguments passed to `reproducible::Cache`
#'
#' @importFrom data.table as.data.table
#' @importFrom stats complete.cases
#' @importFrom crayon blue
#' @export
biomodModelingWrapper <- function(sp, responseVar, responseVarData, predictorVars, predictorVarsData,
                                  dir.name, BIOMOD_ModelingArgs = list(
                                    models = c("GLM", "MARS"),
                                    bm.options = NULL,
                                    CV.k = 5, CV.perc = 100,
                                    metric.eval = c("TSS", "ROC"),
                                    modeling.id = "test"),
                                  ...) {
  if (requireNamespace("biomod2", quietly = TRUE)) {
    if (is.null(bm.options)) {
      bm.options <- biomod2::BIOMOD_ModelingOptions(GLM = list(type = "simple"),
                                                    GAM = list(k = 2))
    }
    ## Checks ------------------
    if (is(responseVar, "list") & is.null(names(responseVar))) {
      stop("'responseVar' must be a named list with 'sp' as names, or a character vector")
    }

    if (is(responseVar, "list")) {
      if (!sp %in% names(responseVar)) {
        stop(paste("No 'responseVar' values found for", sp))
      }
      responseVar <- responseVar[[sp]]
    }

    if (is(responseVarData, "list") & is.null(names(responseVarData))) {
      stop("'responseVarData' must be a named list with 'sp' as names, or a single data.table")
    }

    if (is(responseVarData, "list")) {
      if (!sp %in% names(responseVarData)) {
        stop(paste("No 'responseVarData' data found for", sp))
      }
      responseVarData <- responseVarData[[sp]]
    }

    if (is(predictorVars, "list") & is.null(names(predictorVars))) {
      stop("'predictorVars' must be a named list with 'sp' as names, or a character vector")
    }

    if (is(predictorVars, "list")) {
      if (!sp %in% names(predictorVars)) {
        stop(paste("No 'predictorVars' values found for", sp))
      }
      predictorVars <- predictorVars[[sp]]
    }

    if (is(predictorVarsData, "list") & is.null(names(predictorVarsData))) {
      stop("'predictorVarsData' must be a named list with 'sp' as names, or a single data.table")
    }

    if (is(predictorVarsData, "list")) {
      if (!sp %in% names(predictorVarsData)) {
        stop(paste("No 'predictorVarsData' values found for", sp))
      }
      predictorVarsData <- predictorVarsData[[sp]]
    }

    notFound <- predictorVars[!c("x", "y", predictorVars) %in% colnames(predictorVarsData)]
    if (length(notFound)) {
      stop("The following predictors were not found in 'predictorVarsData': ",
           paste(notFound, collapse = ", "))
    }

    ## create outputs dir
    suppressWarnings({dir.create(dir.name, recursive = TRUE)})

    ## Merge response and predictor variables --------------------
    envCols <- c("pixelIndex", "x", "y", predictorVars)
    predictorVarsData <- predictorVarsData[, ..envCols]

    ## x and y may have extremely small differences between the two data.tables (see test below)
    # round coordinates to solve issue
    ## tests:
    # test <- responseVarData[predictorVarsData, on = "pixelIndex", nomatch = 0]
    # head(test$x - test$i.x)
    # [1] -1.862645e-09 -1.862645e-09 -1.862645e-09 -1.862645e-09 -1.862645e-09 -1.862645e-09
    responseVarData[, `:=`(x = round(x, 6), y = round(y, 6))]
    predictorVarsData[, `:=`(x = round(x, 6), y = round(y, 6))]
    # get complete.cases now for a faster join
    responseVarData <- responseVarData[complete.cases(responseVarData)]
    predictorVarsData <- predictorVarsData[complete.cases(predictorVarsData)]

    allData <- merge(responseVarData, predictorVarsData,
                     on = c("pixelIndex", "x", "y"), all = TRUE)
    allData <- allData[complete.cases(allData)]

    ## BIOMOD operations -----------------------------------------
    bm.format <- biomod2::BIOMOD_FormatingData(resp.var = as.data.frame(allData[, ..responseVar]),
                                               resp.xy = as.data.frame(allData[, .(x, y)]),
                                               expl.var = as.data.frame(allData[, ..predictorVars]),
                                               resp.name = sp,
                                               dir.name = dir.name)
    BIOMOD_ModelingArgs <- append(list(bm.format = bm.format,
                                       do.full.models	= FALSE),
                                  BIOMOD_ModelingArgs)

    bm.mod <- Cache(.BIOMOD_ModelingRetry,
                    BIOMOD_ModelingArgs = BIOMOD_ModelingArgs,
                    ...)

    return(bm.mod)
  } else {
    stop("Package biomod2 not installed. Install using: `remotes::install_github('CeresBarros/biomod2@dev_currentCeres')`.")
  }
}

.BIOMOD_ModelingRetry <- function(BIOMOD_ModelingArgs) {
  i <- 1   ## try 5 times if models fail
  while (i < 6) {
    bm.mod <- do.call(what = biomod2::BIOMOD_Modeling,
                      args = BIOMOD_ModelingArgs)
    if (all(bm.mod@models.failed == "none")) {
      break()
    } else {
      i <- i + 1
      message(blue("Some SEP models failed, retrying... attempt", i))
    }
  }
  if (i == 6 & !all(bm.mod@models.failed == "none")) {
    sp <- sub("SEP_model_", "", BIOMOD_ModelingArgs$modeling.id)
    warning(red("Some/all SEP models could not be computed for", sp,
                "\nConsider rerunning '.BIOMOD_ModelingRetry', simplifying the models,",
                "or choosing another algorithm.\nTo clear cached results and run again try:\n",
                paste0("clearCache(..., userTags = c('.BIOMOD_ModelingRetry', '", sp, "'))")))

  }
  bm.mod
}

#' BIOMOD_EnsembleModeling wrapper
#' @param bm.mod output of `biomod2::BIOMOD_Modeling`
#' @param metric.select.thresh passed `biomod2::BIOMOD_EnsembleModeling`.
#'  By default no thresholding is applied.
#' @param ... passed to Cache
#'
#' @export
biomodEnsembleWrapper <- function(bm.mod, metric.select.thresh = NULL, ...) {
  if (requireNamespace("biomod2", quietly = TRUE)) {
    # ## provide thresholds
    metric.select <- unique(bm.mod@models.evaluation@val$metric.eval)
    metric.evalEnsemble <- metric.select   ## use same metrics to evaluate ensemble

    ## exclude metrics wihout defined thresholds.
    ## e.g. RMSE can't be used, since it can't be used with a minimum (but a maximum)
    metric.select[!metric.select %in% c("TSS", "ROC", "R2")] <- NULL

    bm.em <- Cache(biomod2::BIOMOD_EnsembleModeling,
                   bm.mod = bm.mod,
                   models.chosen = 'all',
                   em.by = 'all',
                   metric.select = metric.select,
                   metric.select.thresh = metric.select.thresh,
                   metric.eval = metric.evalEnsemble,
                   prob.mean = FALSE,
                   prob.mean.weight = TRUE,
                   ...)

    return(bm.em)
  } else {
    stop("Package biomod2 not installed. Install using: `remotes::install_github('CeresBarros/biomod2@dev_currentCeres')`.")
  }
}



#' BIOMOD_Projection wrapper
#'
#' @param bm.mod output of `biomod2::BIOMOD_Modeling`
#' @param bm.mod passed to `biomod2::BIOMOD_Projection`. If not supplied
#'   the data used to fit the model will be used.
#' @param proj.name passed to `biomod2::BIOMOD_Projection`
#' @param new.env passed to `biomod2::BIOMOD_Projection`
#' @param new.env.xy passed to `biomod2::BIOMOD_Projection`
#' @param ... passed to Cache
#'
#' @export
biomodProjWrapper <- function(bm.mod, proj.name = "testProj", new.env = NULL, new.env.xy = NULL, ...) {
  if (requireNamespace("biomod2", quietly = TRUE)) {
    if (is.null(new.env)) {
      new.env <- biomod2::get_formal_data(bm.mod)@data.env.var
    }

    if (inherits(new.env, "data.table")) {
      new.env <- as.data.frame(new.env)
    }

    if (!is.null(new.env.xy) & inherits(new.env.xy, "data.table")) {
      new.env.xy <- as.data.frame(new.env.xy)
    }

    ## there may be additional columns in the new data, biomod doesn't like this
    ##  so make sure only the predictors are here.
    new.env <- new.env[, bm.mod@expl.var.names]

    bm.proj <- Cache(biomod2::BIOMOD_Projection,
                     bm.mod = bm.mod,
                     new.env = new.env,
                     new.env.xy = new.env.xy,
                     proj.name = proj.name,
                     models.chosen = "all",
                     compress = FALSE,
                     build.clamping.mask = FALSE,
                     output.format = ".RData",
                     do.stack = TRUE,
                     ...)

    return(bm.proj)
  } else {
    stop("Package biomod2 not installed. Install using: `remotes::install_github('CeresBarros/biomod2@dev_currentCeres')`.")
  }
}


#' BIOMOD_EnsembleForecasting wrapper
#'
#' @param bm.em output of `biomod2::BIOMOD_EnsembleModeling`
#' @param bm.proj output of `biomod2::BIOMOD_Projection`
#' @param new.env passed to `biomod2::BIOMOD_EnsembleForecasting`
#' @param new.env.xy passed to `biomod2::BIOMOD_EnsembleForecasting`
#' @param keep.in.memory passed to `biomod2::BIOMOD_EnsembleForecasting`
#' @param ... passed to `reproducible::Cache`
#'
#' @export
biomodEnsembleFrcstWrapper <- function(bm.em, bm.proj = NULL, proj.name = NULL, new.env = NULL, new.env.xy = NULL,
                                       keep.in.memory = TRUE, ...) {
  if (requireNamespace("biomod2", quietly = TRUE)) {
    if (!is.null("new.env") & inherits(new.env, "data.table")) {
      new.env <- as.data.frame(new.env)
    }

    if (!is.null("new.env.xy") & inherits(new.env.xy, "data.table")) {
      new.env.xy <- as.data.frame(new.env.xy)
    }
    bm.em.proj <- Cache(biomod2::BIOMOD_EnsembleForecasting,
                        bm.em = bm.em,
                        bm.proj = bm.proj,
                        proj.name = proj.name,
                        new.env = new.env,
                        new.env.xy = new.env.xy,
                        models.chosen = "all",
                        compress = "zip",
                        keep.in.memory = keep.in.memory,
                        ...)

    return(bm.em.proj)
  } else {
    stop("Package biomod2 not installed. Install using: `remotes::install_github('CeresBarros/biomod2@dev_currentCeres')`.")
  }
}


#' Make maps from BIOMOD_EnsembleForecasting results
#'
#' @param bm.em.proj output of `biomod2::BIOMOD_EnsembleForecasting`. Note that
#'   x, y coordinates extracted from `bm.em.proj@coord` must be in the same projection
#'   as `rasTemplate`
#' @param predModel character. Which model predictions should be used. Chose one of
#'   `bm.em.proj@models.projected`
#' @param rasTemplate a template RasterLayer or SpatRaster that can be used to map the projections
#'
#' @return a SpatRaster
#'
#' @importFrom terra rast extract vect crs
#' @importFrom data.table as.data.table
#' @export
biomodEnsembleProjMaps <- function(bm.em.proj, predModel, rasTemplate, rasterToMatch) {
  if (requireNamespace("biomod2", quietly = TRUE)) {
    if (!predModel %in% bm.em.proj@models.projected) {
      stop("Can't find model ", predModel,
           ".\n  Choose one of '", paste(bm.em.proj@models.projected, collapse = "', '"), "'.")
    }

    if (!is(rasTemplate, "SpatRaster")) {
      rasTemplate <- rast(rasTemplate)
    }

    ## need to get RTM cell ID _actually_ used in the models
    ## as the data is joined and potentially reduced by biomodModelingWrapper
    predDataCoords <- bm.em.proj@coord
    predDataCoords <- vect(predDataCoords, geom = c("x", "y"), crs = crs(rasterToMatch, proj = TRUE))
    predDataCoords <- as.data.table(terra::extract(rasTemplate, predDataCoords, cells = TRUE))  ## needs pkg name here
    setnames(predDataCoords, "cell", "pixelIndex")

    predDT <- as.data.table(biomod2::get_predictions(bm.em.proj))
    setnames(predDT, c("points", "pred"), c("ID", "SEPpred"))

    ## make sure values are wiped from the template raster
    rasTemplate <- rasTemplate
    rasTemplate[] <- NA

    if (!is.null(predModel)) {
      predDT <- predDT[full.name == predModel]
    }

    if (any(duplicated(predDT$ID))) {
      stop("There seem to be several SEP predictions for the same cell",
           ".\n  Try passing 'predModel'.")
    }

    predDT <- predDT[predDataCoords, on = "ID", nomatch = NA]

    rasTemplate[predDT$pixelIndex] <- predDT[["SEPpred"]]/1000   ## to save disk space BIOMOD saves probs multiplied by 1000
    names(rasTemplate) <- predModel

    return(rasTemplate)
  } else {
    stop("Package biomod2 not installed. Install using: `remotes::install_github('CeresBarros/biomod2@dev_currentCeres')`.")
  }
}


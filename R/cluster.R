#' High level wrapper around parallelly::makeClusterPSOCK for many multi-core machines
#'
#' This is create to only work with machines that are successfully linked with ssh keys. It
#' uses reverse tunnelling, so the keys are only needed in one direction (e.g., from the
#' main to the workers).
#' One common configuration for High Performance Clusters is several virtual machines (VM),
#' each with multiple cores. Usually, 1/2 of those cores are hyperthreaded, which means
#' that when all cores are utilized on a single VM, the per-core speed will be substantially
#' slower than if only the physical cores are used. Similarly, when R sets up a
#' cluster, it doesn't know that collections of the multiple cores are actually on
#' the same virtual machine. When moving large objects "to the cores", it is actually
#' substantially faster to move the files to one core per machine, write the objects
#' to disk locally on that virtual machine, then spawn all the desired cores per machine,
#' reading locally the objects back into memory in each R session core. Similarly,
#' installing required packages should also be done only once per VM, as each core
#' shares a common disk.
#'
#' This function does these steps:
#'
#' \enumerate{
#'   \item spawn 1 core per worker;
#'   \item copy all needed files to that core,
#'   \item write those objects to disk locally,
#'   \item install any missing packages or package versions. The default is to
#'         use the same library path on the remote machines as the local machine. It
#'         is likely a good idea to make sure the \code{libPath} is not the default
#'         personal library on any of the other \code{workers};
#'   \item optionally, test brute cpu speed on each unique worker to determine the
#'          relative speed of the VMs ;
#'   \item if \code{workers} is of length 1 per ip, or if \code{workers} is not
#'            sufficiently long to create the cluster of \code{numCoresNeeded}, then
#'            it will distribute the \code{numCoresNeeded} across workers, minimizing
#'            the number of hyperthreaded cores, and using more faster cores than slower
#'            cores;
#'   \item shut down the cluster of "single-core-per-VM",
#'   \item start a new possibly multi-core cluster with all desired cores per VM,
#'           either given as a precise
#'           set of \code{workers}, or letting the function determine the set;
#'   \item read all objects into memory in each R-session (i.e., 1 per core);
#'   \item return the multi-core cluster object for use by user.
#' }
#'
#' 1) spawn 1 core per worker, 1a) optionally, test brute cpu speed
#' on each unique worker to determine the relative speed of the VMs 1b) if \code{workers}
#' is of length 1 per ip, or if \code{workers} is not sufficiently long to create the
#' cluster of \code{numCoresNeeded}
#'
#' @author Eliot McIntire
#' @export
#' @param workers A character vector of ips (including possily 'localhost') of length
#'    one each ip, or a character vector of possibly repeated ips, and possibly including
#'    'localhost', such that the length is equal to \code{workers}
#' @param objsToExport A character vector naming the objects that should be copied over
#'    to each worker in the cluster.
#' @param reqdPkgs A character vector of packages (using `Require` formatting, so can
#' specify GitHub packages and minimum version number) that should be (installed if necessary)
#' and loaded (i.e., in a `require` call).
#' @param quotedExtra An optional quoted command to run in single cluster (i.e., once per unique
#'        worker) after `.libPaths(libPath)` is run, e.g., forcing an `install.packages`. See
#'        example.
#' @param libPaths The path in which the R packages should be installed and loaded from.
#' @param doSpeedTest Numeric. This should be a proportion of "all real cores". Thus,
#' if set to 1, it will run a test on 28 cores on a 56 core machine. If not 0,
#' a small (e.g., 10 second) test will be run \code{doSpeedTest * numberRealCores} times on
#' each worker to test the raw cpu speed using loops and \code{rnorm}.
#' If this number is large, then the test may take a >60 seconds simply because of
#' the cost of setting up the PSOCK cluster.
#' @param envir The environment in which to find the \code{objsToExport}, defaults to
#'   \code{parent.frame()}
#' @param numCoresNeeded A numeric of length 1, indicating the desired number of workers.
#' @param adjustments A vector of length \code{workers} that can be used to artificially
#'   adjust (i.e., weight) the workers more or less than others. These numbers are
#'   relative to each other, so \code{c(1.2, 0.8)} will result in 50% more cores
#'   than estimated during \code{doSpeedTest} on the first worker than the second.
#'   Default is 1.
#'
#' @note
#' If working on Ubuntu machines, the default is to install binary packages from the
#' Rstudio binary mirror of CRAN. If `source` packages are required (e.g., for compiling
#' spatial R packages on a system that has recently updated its system spatial libraries,
#' such as GDAL, GEOS etc.), then this will have to be done manually. See example.
#' @importFrom parallelly makeClusterPSOCK
#' @importFrom parallel clusterEvalQ clusterExport
#' @importFrom Require Require
#' @examples
#' \dontrun{
#' largeObj <- rnorm(1e6)
#' ips <- c("localhost", "10.20.0.213")
#' objsToExport <- "largeObj"
#' reqdPkgs <- c("achubaty/amc@development", "raster (>= 3.4-5)")
#' cl <- clusterSetup(workers = ips, objsToExport = objsToExport,
#'                    reqdPkgs = reqdPkgs,
#'                    numCoresNeeded = 4)
#' # Now use the cl with parallel or DEoptim
#' parallel::clusterEvalQ(cl, rnorm(1))
#' parallel::stopCluster(cl)
#'
#' # Installing some packages from source on linux (source is default)
#' cl <- clusterSetup(workers = ips, objsToExport = objsToExport,
#'                    reqdPkgs = reqdPkgs,
#'                    quotedExtra = quote(install.packages(c("rgdal", "rgeos", "sf",
#'                                        "sp", "raster", "terra", "lwgeom"),
#'                                        repos = "https://cran.rstudio.com")),
#'                    numCoresNeeded = 4)
#' parallel::stopCluster(cl)
#'
#' }
clusterSetup <- function(workers, objsToExport, reqdPkgs, quotedExtra,
                         libPaths = .libPaths()[1],
                         doSpeedTest = 0, envir = parent.frame(),
                         # fn = ".allObjs.rda",
                         maxPerWorker = Inf,
                         numCoresNeeded = ceiling(detectCores() * 0.8),
                         adjustments = rep(1, length(workers))) {

  workersInLocalhost <- workers %in% "localhost"
  if (packageVersion("parallelly") == "1.26.0" &&
      ( any(workersInLocalhost) && !all( workersInLocalhost) ) )  {
    stop("package parallelly version 1.26.0 has a bug when mixing 'localhost' with other ips; ",
         "please install an older version: \n",
         "install.packages('https://cran.r-project.org/src/contrib/Archive/parallelly/parallelly_1.25.0.tar.gz', repos = NULL)")
  }
  fn <- reproducible::.suffix(paste0(".allObjs.rda"), paste0("_", amc::rndstr(1)))
  if (doSpeedTest > 1) doSpeedTest <- 1
  clusterIPs <- clusterSetupSingles(workers = workers, objsToExport = objsToExport,
                                    reqdPkgs = reqdPkgs, fn = fn,
                                    libPaths = libPaths, doSpeedTest = doSpeedTest,
                                    maxPerWorker = maxPerWorker,
                                    numCoresNeeded = numCoresNeeded, adjustments = adjustments)

  message("Starting cluster with all cores per machine")
  cl <- parallelly::makeClusterPSOCK(clusterIPs, revtunnel = TRUE)
  # cl <- parallelly::makeClusterPSOCK(workers = clusterIPs, revtunnel = TRUE)
  # on.exit(try(parallel::stopCluster(cl)), add = TRUE)
  clusterExport(cl, varlist = c("fn", "reqdPkgs"), envir = environment())
  st <- system.time(clusterEvalQ(cl, {
    load(file = fn, envir = .GlobalEnv)
    .libPaths(libPaths)
    suppressMessages(
      Require::Require(reqdPkgs, install = FALSE)
    )
  }))
  st <- system.time(clusterEvalQ(cl, {
    try(unlink(fn), silent = TRUE)
  }))

  return(cl)
}

#' @importFrom parallel stopCluster clusterExport
#' @importFrom reproducible messageDF
clusterSetupSingles <- function(workers, objsToExport, reqdPkgs, quotedExtra,
                                libPaths = .libPaths()[1], doSpeedTest = 0, envir = parent.frame(),
                                fn = ".allObjs.rda", maxPerWorker = Inf,
                                numCoresNeeded, adjustments = rep(1, length(workers))) {

  message("Starting cluster with 1 core per machine -- install reqdPkgs; copy objects; write to disk")
  uniqueWorkers <- unique(workers)

  oneWorkerPerIPRequested <- (identical(workers, uniqueWorkers) )
  oneWorkerPerIPIsNotEnough <- oneWorkerPerIPRequested && (length(workers) < numCoresNeeded)

  #if (!oneWorkerPerIPIsNotEnough) # close it if it isn't enough
  on.exit(try(parallel::stopCluster(clSingle), silent = TRUE), add = TRUE)


  clSingle <- parallelly::makeClusterPSOCK(workers = uniqueWorkers, revtunnel = TRUE)

  if (doSpeedTest > 0 || oneWorkerPerIPIsNotEnough )  {
    detCores <- clusterEvalQ(clSingle, parallel::detectCores())
    doSpeedTest <- ceiling(unlist(detCores) * doSpeedTest / 2)
    names(doSpeedTest) <- uniqueWorkers
    clus <- if (all(doSpeedTest == 1)) clSingle else workers
    out2 <- determineClusters(clus = clus, doSpeedTest = doSpeedTest,
                              uniqueWorkers = uniqueWorkers, numCoresNeeded = numCoresNeeded,
                              adjustments = adjustments, maxPerWorker = maxPerWorker)
    reproducible::messageDF(out2)

  }

  if (oneWorkerPerIPIsNotEnough) {
    clusterIPs <- rep(out2$.id, out2$cores)
  } else {
    clusterIPs <- if (length(workers) > numCoresNeeded)
      sample(workers, size = numCoresNeeded)
    else
      workers
    reproducible::messageDF(as.data.frame(table(clusterIPs)))
  }
  # parallel::stopCluster(clSingle)

  # clSingle <- parallelly::makeClusterPSOCK(workers = unique(clusterIPs), revtunnel = TRUE)
  clusterExport(clSingle, varlist = objsToExport, envir = envir)
  clusterExport(clSingle, varlist = c("fn", "reqdPkgs", "libPaths", "objsToExport"), envir = environment())

  clusterEvalQ(clSingle, {
    if (any(!dir.exists(libPaths)))
      dir.create(libPaths[1], recursive = TRUE)
    .libPaths(libPaths)
    ## set CRAN repos; use binary linux packages if on Ubuntu
    local({
      options(Ncpus = parallel::detectCores() / 2)
      options("repos" = c(CRAN = "https://cran.rstudio.com"))

      if (Sys.info()["sysname"] == "Linux" && grepl("Ubuntu", utils::osVersion)) {
        .os.version <- strsplit(system("lsb_release -c", intern = TRUE), ":\t")[[1]][[2]]
        .user.agent <- paste0(
          "R/", getRversion(), " R (",
          paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"]),
          ")"
        )
        options(repos = c(CRAN = paste0("https://packagemanager.rstudio.com/all/__linux__/",
                                        .os.version, "/latest")))
        options(HTTPUserAgent = .user.agent)
      }
    })

    if (!suppressWarnings(require("Require"))) {
      install.packages("Require")
    }
    if (exists("quotedExtra"))
      eval(quotedExtra)
    # install.packages(c("rgdal", "rgeos", "sf", "sp", "raster", "terra", "lwgeom"), repos = "https://cran.rstudio.com")
    suppressMessages(Require::Require(reqdPkgs, install = TRUE, require = FALSE))
    save(list = objsToExport, file = fn)
  })
  # if (!oneWorkerPerIPIsNotEnough)
  parallel::stopCluster(clSingle)

  return(clusterIPs)
}

determineClusters <- function(clus, doSpeedTest, uniqueWorkers, numCoresNeeded, adjustments, maxPerWorker) {
  if (any(doSpeedTest > 1)) {
    actualWorkers <- if (length(doSpeedTest) > 1) {
      rep(names(doSpeedTest), doSpeedTest)
    } else {
      rep(uniqueWorkers, each = doSpeedTest)
    }
    message("Starting cluster with ", sum(doSpeedTest),
            " cores -- to run tests. This may take 1 second per unique worker ",
            " (est: ", length(actualWorkers), " secs.)")
    out <- list()
    groups <- if (length(actualWorkers) > 110) {
      grp <- list()
      i <- 0
      doSpeedTestRunning <- doSpeedTest
      while (length(doSpeedTestRunning)) {
        i <- i + 1
        grp[[i]] <- which(cumsum(doSpeedTestRunning) < 110)
        doSpeedTestRunning <- doSpeedTestRunning[-grp[[i]]]
        grp[[i]][] <- i
      }
      message("Testing more than 110 workers; doing in ", i, " batches")
      whGroup <- do.call(c, grp)
      groups <- rep(whGroup, doSpeedTest)
      split(actualWorkers, groups)
    } else {
      list(actualWorkers)
    }
    for (grp in seq(groups)) {
      clus <- parallelly::makeClusterPSOCK(workers = groups[[grp]], revtunnel = TRUE)
      on.exit(try(stopCluster(clus), silent = TRUE), add = TRUE)
      out[[grp]] <- runSpeedTest(clus)
      names(out[[grp]]) <- groups[[grp]]
      out[[grp]] <- data.table::rbindlist(out[[grp]], idcol = TRUE)
      stopCluster(clus)
    }
    out2 <- data.table::rbindlist(out)
  } else {
    if (doSpeedTest == 1 && length(uniqueWorkers) > 1) {
      message("testing speed on each to estimate number cores to use")
      out <- runSpeedTest(clus)
    } else {
      out <- clusterEvalQ(clus, {
        data.table::data.table(elapsed = 1, trueCores = parallel::detectCores())
      })
    }
    names(out) <- sapply(clus, function(x) x$host)
    out2 <- rbindlist(out, idcol = TRUE)
  }

  # parallel::clusterEvalQ(clus, system("pkill -f workRSOCK"))
  if (length(adjustments) != length(uniqueWorkers))
    adjustments <- rep(adjustments, length(uniqueWorkers))
  out2 <- out2[, list(elapsedMax = max(elapsed),
                      elapsed = mean(elapsed),
                      trueCores = trueCores[1]), by = ".id"]
  #relSpeed <- out2$trueCores/out2$elapsed*numCoresNeeded
  relSpeed <- 1/out2$elapsedMax
  relSpeed <- relSpeed/max(relSpeed)
  relSpeedThresh <- 1/1.4 # this is Hyperthreading penalty -- 40% loss
  nonHTcores <- out2$trueCores/2
  out2[, nonHTcores := nonHTcores]
  sumNonHTcores <- sum(nonHTcores)

  out2[, cores := round(nonHTcores * adjustments / max(adjustments))]

  out2[, maxPerWorker := maxPerWorker]
  out2[, cores := pmin(maxPerWorker, cores)]

  speedGoodEnoughToUse <- relSpeed > relSpeedThresh
  if (any (!speedGoodEnoughToUse)) {
    message("worker(s) ", paste(out2$.id[!speedGoodEnoughToUse], collapse = ", "), " too slow to use at all, given speed test")
    if (sum(pmin(out2$maxPerWorker[speedGoodEnoughToUse], out2$trueCores[speedGoodEnoughToUse])) < numCoresNeeded) {
      message("but, not enough cores without them; using them")
      speedGoodEnoughToUse <- rep(TRUE, length(speedGoodEnoughToUse))
    }
    adjustments <- adjustments[speedGoodEnoughToUse]
  }
  #relSpeed <- relSpeed/numCoresNeeded
  relSpeed[!speedGoodEnoughToUse] <- 0
  out2[, relSpeed := round(relSpeed, 2)]
  # out3 <- copy(out2)
  if (any(relSpeed == 0))
    didntUse <- out2[relSpeed == 0]
  out2 <- out2[relSpeed > 0]


  out2[, estCores := ceiling(cores * relSpeed)]
  # Too many
  if (sum(out2$estCores) > numCoresNeeded) {
    out2[, estCores := ceiling(estCores * numCoresNeeded / sum(out2$estCores) ) ]
  }
  out2[, newRelSpeed := relSpeed]

  while (isTRUE(max(out2[estCores < maxPerWorker]$newRelSpeed) >= min(out2$relSpeed)) && (sum(out2$estCores) < numCoresNeeded))  {
    (toChangeUpNew <- out2[estCores < maxPerWorker][newRelSpeed == max(newRelSpeed)]$.id[1])
    out2[toChangeUpNew == .id, estCores := estCores + 1L]
    out2[toChangeUpNew == .id, newRelSpeed := relSpeed - max(0, estCores - nonHTcores)/nonHTcores * (1 - relSpeedThresh)   ]
  }

  # out4 <- copy(out2)
  # out2 <- copy(out4)
  while (sum(out2$estCores) > numCoresNeeded) {
    (toChangeDownNew <- out2[newRelSpeed == min(newRelSpeed) ]$.id[1])
    out2[toChangeDownNew == .id, estCores := estCores - 1L]
    out2[toChangeDownNew == .id, newRelSpeed := relSpeed - max(0, estCores - nonHTcores)/nonHTcores  ]
    if (any(out2$estCores == 0)) {
      message("the raw speed of ", out2[estCores == 0]$.id, " is too slow to include")
      out2 <- out2[estCores > 0]
    }
  }

  out2[, cores := estCores]
  out2[, estCores := NULL]


  # Increase ncores upwards
  m <- 0
  while (sum(out2$cores) < numCoresNeeded) {
    m <- ((m + 1) - 1) %% NROW(out2) + 1
    out2[m, cores := cores + relSpeed]
  }
  # Decrease down to 1 each
  out2[, cores := ceiling(cores)]
  m <- 0
  set(out2, NULL, "cores", as.numeric(out2$cores))
  while ((sum(floor(out2$cores)) > numCoresNeeded) && any(out2$cores > 1)) {
    m <- ((m + 1) - 1) %% NROW(out2) + 1
    if (out2[m,]$cores > 1) {
      maxC <- max(out2$cores)
      out2[m, cores := cores - 1]
    }
  }
  # set(out2, NULL, "cores", floor(out2$cores))
  # Decrease down to 0 or 1 each
  m <- 0
  while (sum(out2$cores) > numCoresNeeded && any(out2$cores > 0)) {
    m <- ((m + 1) - 1) %% NROW(out2) + 1
    if (out2[m,]$cores > 0) {
      out2[m, cores := cores - 1]
    }
  }

  if (any(relSpeed == 0))
    out2 <- rbindlist(list(out2, didntUse[, cores := 0]), fill = TRUE)

  return(out2)
}


runSpeedTest <- function(clus) {
  out <- clusterEvalQ(clus, {
    ss <- system.time({
      for (i in 1:15000) rnorm(1e4)
    })
    data.table::data.table(elapsed = ss[3], trueCores = parallel::detectCores())
  })
  out
}

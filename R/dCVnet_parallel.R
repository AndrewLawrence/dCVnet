# This file contains an optimised version of a hack applied to mclapply to give
#   parallelisation options which "just work" over multiple platforms.
#
# Behaviour:
#
#   For linux / mac forking will be used (parallel::mclapply)
#
#   For Windows (which cannot support forking) a socket based wrapper for
#     mclapply is created as long as options("mc.cores") is greater than 1.
#
#   This behaviour is motivated by the enormous time+memory overhead of the
#     socket method. In testing it slowed down smaller jobs, and quickly
#     caused memory swapping with larger data.
#
#   There is a sweet spot where dCVnet can run 2-3 times faster on windows
#     as long as the problem involves a high number of repetitions/folds.
#
#   To run single threaded on windows simply set options(mc.cores = 1L),
#     Experiment with setting values greater than one to get jobs
#     with a large compute time to run faster.
#
# The code is based on https://github.com/nathanvan/parallelsugar
#   but optimised for dCVnet to reduce memory load
#   (a particular problem for dCVnet)

# Make a sockets version of mclapply:
#' @import parallel
#' @importFrom utils sessionInfo
mclapply_socket <- function(X,
                            FUN,
                            ...,
                            mc.preschedule = TRUE,
                            mc.set.seed = TRUE,
                            mc.silent = FALSE,
                            mc.cores = NULL,
                            mc.cleanup = TRUE,
                            mc.allow.recursive = TRUE) {
  ## Create a cluster
  if (is.null(mc.cores)) {
    dcores <- detectCores()
    mc.cores <- min(length(X), dcores - 1)
  }
  cl <- parallel::makeCluster(mc.cores)

  tryCatch({
    ## Find out the names of the loaded packages
    loaded.package.names <- c(sessionInfo()$basePkgs, ## Additional packages
                              names(sessionInfo()$otherPkgs))

    ### Ship it to the clusters
    parallel::clusterExport(cl, "loaded.package.names", envir = environment())

    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parallel::parLapply(cl, seq.int(length(cl)), function(xx) {
      lapply(loaded.package.names, function(yy) {
        require(yy, character.only = TRUE)
      })
    })

    clusterExport_function(cl, FUN)

    ## Run the lapply in parallel, with a special case for the ... arguments
    if (length(list(...)) == 0) {
      return(parallel::parLapply(cl = cl, X = X, fun = FUN))
    } else {
      return(parallel::parLapply(
        cl = cl,
        X = X,
        fun = FUN,
        ...
      ))
    }
  }, finally = {
    ## Stop the cluster
    parallel::stopCluster(cl)
  })
}


# Overwrite the serial version of mclapply on Windows if mc.cores > 1:
if (Sys.info()[["sysname"]] == "Windows" &&
      !is.null(getOption("mc.cores")) &&
      getOption("mc.cores") > 1L) {
  mclapply <- mclapply_socket
}

clusterExport_function <- function(cl,
                                   FUN,
                                   object_list = c("outfolds",
                                                   "x",
                                                   "y",
                                                   "pp_fn",
                                                   "family",
                                                   "cl.marcvglm",
                                                   "cutoff",
                                                   "offset")) {
  ## We want the enclosing environment, not the calling environment
  ## (I had tried parent.frame, which was not what we needed)
  ##
  ## Written by Hadley Wickham, off the top of his head, when I asked him
  ##   for help at one of his Advanced R workshops.

  keep_ls <- function(keepers, e) {
    r <- ls(all.names = TRUE, envir = e)
    r[r %in% keepers]
  }

  env <- environment(FUN)
  while (!identical(env, globalenv())) {
    env <- parent.env(env)
    parallel::clusterExport(cl, keep_ls(object_list, env), envir = env)
  }
  parallel::clusterExport(cl, keep_ls(object_list, env), envir = env)
  ## // End Hadley Wickham
}


# Single Alpha Functions: ------------------------------------------------------
#   These functions wrap cv.glmnet to run over multiple specified folds and
#   pass back results along with selected best performing lambdas.
#   as suggested by the name, only one alpha is present.


#' repeated.cv.glmnet
#'
#' Repeatedly runs a \code{\link[glmnet]{cv.glmnet}} and returns averaged
#'     results. *This is intended as a dCVnet internal function*.
#'
#' The code will run for any glmnet family, but folds & lambdas must be
#' correctly specified.
#'
#'
#' @inheritParams glmnet::glmnet
#' @inheritParams glmnet::cv.glmnet
#' @param lambdas use a fixed, user supplied lambda sequence (descending)
#'     see \code{\link[glmnet]{glmnet}}
#' @param opt.lambda.type Method for selecting optimum lambda. One of
#'                         \itemize{
#'                         \item{\code{"minimum"} - returns the lambda with best
#'                         CV score.}
#'                         \item{\code{"se"} - returns the +1 se lambda}
#'                         \item{\code{"percentage"} - returns minimum lambda
#'                         scaled by a factor, e.g. allowing lambda+3pc}
#'                         }
#' @param opt.lambda.type.value determines the se multiplier or percentage
#'                               for \code{opt.lambda.type}.
#     e.g. 'percentage' & type.value = 1.03 gives lambda + 3%.
#     e.g. 'se' & type.value = 1.0 gives the 'standard' lambda+1se.
#' @param folds This is a list where each element is an integer vector
#'     of length *n_cases*. The integer for each case labels it as belonging
#'     to a fold *1:n_folds*. This argument implicitly sets the number of repeats
#'     and the k in repeated k-fold cv.
#' @param ... arguments passed to \code{\link[glmnet]{cv.glmnet}}
#' @param debug if TRUE return models and unaveraged results (default: FALSE).
#'
#' @return a data.frame object of class \code{\link{repeated.cv.glmnet}}
#'     containing averaged metrics. Has the following columns:
#'     \itemize{
#'     \item{lambda - lambda at which performance evaluated}
#'     \item{cvm - average performance metric}
#'     \item{cvsd - average sd of performance metric}
#'     \item{cvup - average cvm + cvsd}
#'     \item{cvlo - average cvm - cvsd}
#'     \item{nzero - average number of nonzero predictors}
#'     \item{lambda.min - logical indicating 'best' performing lambda
#'         (see opt.lambda.type and opt.lambda.type.value)}
#'     }
#'     Also contains attributes for the response family (\code{family})
#'         and cv-type (\code{type.measure})
#' @export
repeated.cv.glmnet <- function(x, y,
                               folds = NULL,
                               lambdas = NULL,
                               alpha = 1,
                               nfolds = NULL,
                               family = c("gaussian", "binomial",
                                          "poisson", "multinomial",
                                          "cox", "mgaussian"),
                               opt.lambda.type = "minimum",
                               opt.lambda.type.value = 1,
                               ...,
                               debug = FALSE) {
  cl <- c(as.list(environment()), list(...))

  # We typically want to use fixed folds and fixed lambda sequence, but for
  #   convenience/generality include fallback modes (with warnings):
  if ( missing(lambdas) ) {
    warning("no lambdas provided: extracting glmnet default lambdas")
    # get elements of call suitable for glmnet:
    cl.gnet <- cl[names(cl) %in% methods::formalArgs(glmnet::glmnet)]
    lambdas <- do.call(glmnet::glmnet, cl.gnet)$lambda # extract lambda list
  }
  if ( missing(folds) ) {  # nolint
    warning("no folds provided: generating 5x10-fold")
    if ( missing(nfolds) ) nfolds <- 10
    folds <- lapply(seq.int(5),
                    function(i) {
                      caret::createFolds(y = y,
                                         k = nfolds,
                                         list = FALSE)
                    })
  }

  # prepare a safe call to cv.glmnet:
  cvgnet_args <- unique(c(methods::formalArgs(glmnet::cv.glmnet),
                          methods::formalArgs(glmnet::glmnet)))
  cvgnet_args <- cvgnet_args[!(cvgnet_args %in% "...")]

  cl.cvgnet <- cl[names(cl) %in% cvgnet_args]
  cl.cvgnet$lambda <- lambdas

  # estimate models over folds:
  models <- lapply(seq_along(folds), function(i) {
    cl.cvgnet$foldid <- folds[[i]]
    do.call(glmnet::cv.glmnet, cl.cvgnet)
  } )
  measure_name <- names(models[[1]]$name)

  # extract results:
  results <- mapply(cv.glmnet.modelsummary,
                    models,
                    alpha,
                    as.character(seq_along(models)),
                    SIMPLIFY = FALSE)
  results <- as.data.frame(data.table::rbindlist(results))

  # mean average results over the repetitions for each lambda:
  av <- aggregate(x = results[, c("cvm", "cvsd", "cvup", "cvlo", "nzero")],
                  by = list(lambda = results$lambda),
                  FUN = mean)

  theLambda <- cvlambdafinder(lambda = av$lambda,
                              cvm = av$cvm,
                              cvsd = av$cvsd,
                              minimise = !("auc" %in% measure_name),
                              type = opt.lambda.type,
                              type.value = opt.lambda.type.value)
  # What optimising function should we use?:

  # Which is the optimum lambda?
  av$lambda.min <- FALSE
  av$lambda.min[match(theLambda, av$lambda)] <- TRUE

  av[order(av$lambda, decreasing = TRUE), ]
  attr(av, "type.measure") <- measure_name
  attr(av, "family") <- family
  attr(av, "class") <- c("repeated.cv.glmnet", "data.frame")

  # models and results could be useful, but we do not usually want these.
  if ( !debug ) {
    return(av)
  } else {
    return(list(averaged = av,
                results = results,
                models = models))
  }
}

#' @export
print.repeated.cv.glmnet <- function(x, ...) {
  alpha <- unique(x$alpha)
  lambda <- range(x$lambda)

  cat(paste("a repeated.cv.glmnet object from package dCVnet\n"))
  cat(paste("Model family:", attr(x, "family"), "\n"))
  cat(paste("Tuning metric (cvm):", attr(x, "type.measure"), "\n"))
  cat(paste0("\talpha:\t\t", prettyNum(alpha), "\n"))
  cat(paste0("\tnlambda:\t", length(unique(x$lambda)), "\n"))
  cat(paste0("\tlambda range:\t",
             paste(prettyNum(lambda), collapse = ", "),
             "\n\n"))
  invisible(x)
}

#' summary.repeated.cv.glmnet
#'
#' a summary of key options and results for a \code{\link{repeated.cv.glmnet}}
#'     object.
#'
#' @param object a a \code{\link{repeated.cv.glmnet}} object.
#' @param ... NULL
#' @export
summary.repeated.cv.glmnet <- function(object, ...) {

  best_up <- object$cvup[object$lambda.min]
  best_lo <- object$cvlo[object$lambda.min]

  close <- vapply(X = unique(object$lambda),
                  FUN = function(x) {
                    probe <- object$cvm[object$lambda == x]
                    return( (probe < best_up) & (probe > best_lo) )
                  },
                  FUN.VALUE = c(TRUE))

  cat("Summary of ")
  print(object)

  cat("Best fit:\n")
  print.data.frame(object[object$lambda.min, ], row.names = FALSE)
  cat(paste0("\n", sum(close), "/", length(close),
             " (", prettyNum(round(mean(close), 2)), "%) ",
             "lambdas with CV fit within 1 S.E. of best fitting\n"))

  invisible(object[object$lambda.min, ])
}


# Multiple Alpha Functions: --------------------------------------------------



#' multialpha.repeated.cv.glmnet
#'
#' Runs a \code{\link{repeated.cv.glmnet}} for a list of alpha values and
#'     returns averaged results, selects the 'best' alpha.
#'     *This is intended to be a dCVnet internal function*
#' @inheritParams repeated.cv.glmnet
#' @param lambdas a list of lambda sequences corresponding to alphalist
#' @param alphalist a vector of alpha values to search.
#' @param k the number of folds for k-fold cross-validation.
#' @param nrep the number of repetitions
#'
#' @param opt.ystratify Boolean.
#'     Outer and inner sampling is stratified by outcome.
#'     This is implemented with \code{\link[caret]{createFolds}}
#' @param opt.uniquefolds Boolean.
#'     In most circumstances folds will be unique. This requests
#'     that random folds are checked for uniqueness in inner and outer loops.
#'     Currently it warns if non-unique values are found.
#'
#' @return an object of class \code{\link{multialpha.repeated.cv.glmnet}}.
#'     This is a 3 item list: \itemize{
#'     \item{inner_results - merged \code{\link{repeated.cv.glmnet}} with
#'         additional columns indicating *alpha* and logical for *best* overall}
#'     \item{inner_best - best selected row from inner_results}
#'     \item{inner_folds - record of folds used}
#'     }
#' @export
multialpha.repeated.cv.glmnet <- function(
  y,
  alphalist = round(seq(0.2, 1, len = 6) ^ exp(1), 2),
  lambdas = NULL,
  k = 10,
  nrep = 5,
  opt.lambda.type = "minimum",
  opt.lambda.type.value = 1,
  opt.ystratify = TRUE,
  opt.uniquefolds = FALSE,
  family,
  ...) {
  cl <- c(as.list(environment()), list(...))

  # We typically want to use a fixed lambda sequence over all folds of the
  #   outer CV, but for convenience/generality include a fallback mode:
  if ( missing(lambdas) ) {
    warning("no lambdas provided: extracting glmnet default lambdas")
    # extract lambda lists
    lambdas <- lapply(alphalist,
                      function(ii) {
                        # get elements of call suitable for glmnet:
                        cl.gnet <- cl[names(cl) %in%
                                        methods::formalArgs(glmnet::glmnet)]
                        cl.gnet$alpha <- ii
                        return(do.call(glmnet::glmnet, cl.gnet)$lambda)
                        })
  }

  # Fold generation:
  ystrat <- y
  if ( identical(opt.ystratify, FALSE) | family %in% c("cox", "mgaussian") ) {
    # caret stratification isn't sensible for cox / mgauss:
    ystrat <- rep("x", NROW(y))
  }
  folds <- lapply(1:nrep, function(i) {
    caret::createFolds(y = ystrat, k = k, list = FALSE, returnTrain = FALSE)
  })

  if ( identical(opt.uniquefolds, TRUE) ) checkForDuplicateCVFolds(folds)

  # prepare the repeated.cv.glmnet call:
  cl.rcvglm <- cl[names(cl) %in% c(methods::formalArgs(repeated.cv.glmnet),
                                   methods::formalArgs(glmnet::glmnet),
                                   methods::formalArgs(glmnet::cv.glmnet))]
  cl.rcvglm$folds <- folds

  # run for each alpha/lambdas pair:
  malist <- lapply(seq_along(alphalist),
                   function(i) {
                     if ( getOption("mc.cores", default = 1) == 1 ) {
                       cat(paste("\tInner Alpha", i, "of",
                                 length(alphalist), Sys.time(), "\n"))
                     }
                     cl.rcvglm$alpha <- alphalist[[i]]
                     cl.rcvglm$lambdas <- lambdas[[i]]

                     repeated <- do.call(repeated.cv.glmnet, cl.rcvglm)
                     repeated$alpha <- as.character(alphalist[[i]])

                     return(repeated)
                   })

  # assemble results:
  tmeas <- attr(malist[[1]], "type.measure")
  tfam  <- attr(malist[[1]], "family")

  malist <- as.data.frame(data.table::rbindlist(malist))
  attr(malist, "type.measure") <- tmeas
  attr(malist, "family") <- tfam

  # pick the optimal alpha:
  bestfun <- ifelse(tmeas == "auc", max, min)
  bestcandidates <- malist[malist$lambda.min, ]
  best <- bestcandidates[bestcandidates$cvm == bestfun(bestcandidates$cvm), ]

  # ties are broken by smaller cvsd followed by sparser solutions (unlikely)
  if ( nrow(best) > 1 ) {
    best <- best[order(best$cvsd, best$nzero), ]
    best <- best[1, ]
  }

  # add column to the multialpha list.
  malist$best <- (malist$alpha == best$alpha) & (malist$lambda == best$lambda)

  return(structure(list(inner_results = malist,
                        inner_best = best,
                        inner_folds = folds),
                   class = "multialpha.repeated.cv.glmnet"))
}

#' @export
print.multialpha.repeated.cv.glmnet <- function(x, ...) {

  best <- x$inner_best
  x <- x$inner_results

  alpha <- unique(x$alpha)
  lambdas <- lapply(alpha, function(i) x$lambda[x$alpha == i]  )

  l_lengths <- vapply(lambdas, length, FUN.VALUE = c(1.0))
  l_min <- prettyNum(vapply(lambdas, min, FUN.VALUE = c(1.0)))
  l_max <- prettyNum(vapply(lambdas, max, FUN.VALUE = c(1.0)))
  l_ranges <- paste(l_min, l_max, "\n")

  rcols <- c("alpha", "lambda", "nzero", "cvm", "cvsd", "cvup", "cvlo")

  R <- x[x$lambda.min, rcols]
  attr(R, "class") <- "data.frame"

  selected <- (R$alpha == best$alpha) & (R$lambda == best$lambda)

  cat("A dCVnet::multialpha.repeated.cv.glmnet object\n")
  cat(paste("Model family:", attr(x, "family"), "\n"))
  cat(paste("Tuning metric (cvm):", attr(x, "type.measure"), "\n\n"))
  cat(paste0("\t", length(alpha), " alpha(s): ",
             paste(alpha, collapse = ", "), "\n\n"))
  cat(paste0("\tLambda Counts:\n"))
  cat(paste0("\t", l_lengths))
  cat(paste0("\n\n\tLambda Ranges:\n"))
  cat(paste0("\t", l_ranges))
  cat("\n")

  R$best <- selected
  invisible(R)
}

#' summary.multialpha.repeated.cv.glmnet
#'
#' a summary of key options and results for a
#'     \code{\link{multialpha.repeated.cv.glmnet}} object.
#'
#' @param object a a \code{\link{multialpha.repeated.cv.glmnet}} object.
#' @param print if FALSE silently returns the summary results table.
#' @param ... NULL
#'
#' @export
summary.multialpha.repeated.cv.glmnet <- function(object, print = TRUE, ...) {
  .get_nsimilar <- function(marc) {
    kk <- lapply(unique(marc$alpha), function(A){

      tdf <- marc[marc$alpha == A, ]
      best_up <- tdf$cvup[tdf$lambda.min]
      best_lo <- tdf$cvlo[tdf$lambda.min]
      close <- vapply(X = unique(tdf$lambda),
                      FUN = function(x) {
                        probe <- tdf$cvm[tdf$lambda == x]
                        return( (probe < best_up) & (probe > best_lo) )
                      },
                      FUN.VALUE = c(TRUE))
      nclose <- sum(close)
      ntotal <- length(close)
      pclose <- nclose / ntotal
      return(list(alpha = A,
                  nclose = nclose,
                  ntotal = ntotal,
                  pclose = pclose))
    })
    kk <- as.data.frame(data.table::rbindlist(kk))
    return(kk)
  }

  if ( print ) {
    cat("Summary of ")
    print(object)
  }

  best <- object$inner_best
  object <- object$inner_results

  rcols <- c("alpha", "lambda", "nzero", "cvm", "cvsd", "cvup", "cvlo")

  R <- object[object$lambda.min, rcols]
  R <- merge(R, .get_nsimilar(object), by = "alpha")

  selected <- (R$alpha == best$alpha) & (R$lambda == best$lambda)

  if ( print ) {

    cat("\nBest fitting alpha:\n")
    print(R[selected, ], row.names = FALSE)

    cat("\n\nAll Alphas:\n")
    print(R, row.names = FALSE)

    R$best <- selected
    invisible(R)

  } else {

    R$best <- selected
    return(R)
  }

}


# Single Alpha Functions: ------------------------------------------------------
#   These functions wrap cv.glmnet to run over multiple specified folds and
#   pass back results along with selected best performing lambdas.
#   as suggested by the name, only one alpha is present.


#' repeated.cv.glmnet
#'
#' Repeatedly runs a \code{\link[glmnet]{cv.glmnet}} and returns averaged
#'     results. *This is intended to be a dCVnet internal function*.
#'
#' @inheritParams glmnet::cv.glmnet
#' @inheritParams glmnet::glmnet
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
#'     containting averaged metrics. Has the following columns:
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
#' @export
repeated.cv.glmnet <- function(x, y,
                               folds,
                               lambdas,
                               alpha,
                               opt.lambda.type = "minimum",
                               opt.lambda.type.value = 1,
                               ...,
                               debug = F) {
  # We need to use fixed folds, fixed lambdas and a single fixed alpha
  #   most arguments are passed straight through to cv.glmnet.

  # estimate models:
  models <- lapply(1:length(folds), function(i) {
    f <- folds[[i]]
    glmnet::cv.glmnet(x = x,
                      y = y,
                      lambda = lambdas,
                      family = "binomial",
                      foldid = f,
                      alpha = alpha,
                      ...)
  } )
  measure_name <- names(models[[1]]$name)
  # extract results:
  results <- mapply(cv.glmnet.modelsummary,
                    models,
                    alpha,
                    as.character(1:length(models)),
                    SIMPLIFY = F)
  results <- do.call(rbind, results)

  # mean average results over the repetitions for each lambda:
  av <- aggregate(
    x = results[, !names(results) %in% c("lambda", "rep",
                                         "lambda.min", "lambda.1se")],
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
  av$lambda.min <- F
  av$lambda.min[match(theLambda, av$lambda)] <- TRUE

  av[order(av$lambda, decreasing = T), ]
  attr(av, "type.measure") <- measure_name
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

  cat(paste("A repeated.cv.glmnet object from package dCVnet\n"))
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

  close <- sapply(unique(object$lambda), function(x) {
    probe <- object$cvm[object$lambda == x]
    return(probe < best_up & probe > best_lo)
  })

  print(object)

  cat("Best fit:\n")
  print.data.frame(object[object$lambda.min, ], row.names = F)
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
#' @param alphalist a vector of alpha values to search.
#' @param k the number of folds for k-fold cross-validation.
#' @param nrep the number of repetitions
#'
#' @param opt.ystratify Boolian.
#'     Outer and inner sampling is stratified by outcome.
#'     This is implemented with \code{\link[caret]{createFolds}}
#' @param opt.uniquefolds Boolian.
#'     In most circumstances folds will be unique. This requests
#'     that random folds are checked for uniqueness in inner and outer loops.
#'     Currently it warns if non-unqiue values are found.
#'
#' @return an object of class \code{\link{multialpha.repeated.cv.glmnet}}.
#'     This is a 3 item list: \itemize{
#'     \item{inner_results - merged \code{\link{repeated.cv.glmnet}} with
#'         additional columns indicating *alpha* and logical for *best* overall}
#'     \item{inner_best - best selected row from inner_results}
#'     \item{inner_folds - record of folds used}
#'     }
#' @export
multialpha.repeated.cv.glmnet <- function(alphalist,
                                          lambdas,
                                          y,
                                          k,
                                          nrep,
                                          opt.lambda.type = "minimum",
                                          opt.lambda.type.value = 1,
                                          opt.ystratify = T,
                                          opt.uniquefolds = F,
                                          ...) {
  # Fold generation:
  ystrat <- y
  if ( identical(opt.ystratify, FALSE) ) ystrat <- rep("x", length(y))

  folds <- lapply(1:nrep, function(i) {
    caret::createFolds(y = ystrat, k = k, list = FALSE, returnTrain = FALSE)
  })

  if ( identical(opt.uniquefolds, TRUE) ) checkForDuplicateCVFolds(folds)

  malist <- lapply(1:length(alphalist),
                   function(i) {
                     if ( getOption("mc.cores", default = 1) == 1 ) {
                       cat(paste("\tInner Alpha", i, "of",
                                 length(alphalist), Sys.time(), "\n"))
                     }
                     a <- alphalist[[i]]

                     repeated <- repeated.cv.glmnet(
                       alpha = a,
                       lambdas = lambdas[[i]],
                       y = y,
                       folds = folds,
                       opt.lambda.type = opt.lambda.type,
                       opt.lambda.type.value = opt.lambda.type.value,
                       ...)
                     repeated$alpha <- as.character(a)

                     return(repeated)
                   })
  tmeas <- attr(malist[[1]], "type.measure")
  malist <- do.call(rbind, malist)
  attr(malist, "type.measure") <- tmeas

  bestfun <- ifelse(tmeas == "auc", max, min)
  bestcandidates <- malist[malist$lambda.min, ]
  best <- bestcandidates[bestcandidates$cvm == bestfun(bestcandidates$cvm), ]

  # ties are broken by smaller cvsd followed by
  #   sparser solutions.
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

  l_lengths <- sapply(lambdas, length)
  l_min <- prettyNum(sapply(lambdas, min))
  l_max <- prettyNum(sapply(lambdas, max))
  l_ranges <- paste(l_min, l_max, "\n")

  rcols <- c("alpha", "lambda", "nzero", "cvm", "cvsd", "cvup", "cvlo")

  R <- x[x$lambda.min, rcols]
  attr(R, "class") <- "data.frame"

  selected <- (R$alpha == best$alpha) & (R$lambda == best$lambda)

  cat("A multialpha.repeated.cv.glmnet object, from the dCVnet Package\n\n")
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
summary.multialpha.repeated.cv.glmnet <- function(object, print = T, ...) {
  .get_nsimilar <- function(marc) {
    kk <- lapply(unique(marc$alpha), function(A){

      tdf <- marc[marc$alpha == A, ]
      best_up <- tdf$cvup[tdf$lambda.min]
      best_lo <- tdf$cvlo[tdf$lambda.min]
      close <- sapply(unique(tdf$lambda), function(x) {
        probe <- tdf$cvm[tdf$lambda == x]
        return(probe < best_up & probe > best_lo)
      })
      nclose <- sum(close)
      ntotal <- length(close)
      pclose <- nclose / ntotal
      return(list(alpha = A,
                  nclose = nclose,
                  ntotal = ntotal,
                  pclose = pclose))
    })
    kk <- as.data.frame(do.call(rbind, kk))
    return(kk)
  }

  best <- object$inner_best
  object <- object$inner_results

  alpha <- unique(object$alpha)

  lambdas <- lapply(alpha, function(i) object$lambda[object$alpha == i]  )

  l_lengths <- sapply(lambdas, length)
  l_min <- prettyNum(sapply(lambdas, min))
  l_max <- prettyNum(sapply(lambdas, max))
  l_ranges <- paste(l_min, l_max, "\n")

  rcols <- c("alpha", "lambda", "nzero", "cvm", "cvsd", "cvup", "cvlo")

  R <- object[object$lambda.min, rcols]
  R <- merge(R, .get_nsimilar(object), by = "alpha")

  selected <- (R$alpha == best$alpha) & (R$lambda == best$lambda)

  if ( print ) {
    cat("Summary of multiple-alpha repeated k-fold cross-validated glmnet\n\n")

    cat(paste0("\t", length(alpha), " alpha(s): ",
               paste(alpha, collapse = ", "), "\n\n"))
    cat(paste0("\tLambda Counts:\n"))
    cat(paste0("\t", l_lengths))
    cat(paste0("\n\n\tLambda Ranges:\n"))
    cat(paste0("\t", l_ranges))

    cat("\nBest fitting alpha:\n")
    print(R[selected, ], row.names = F)

    cat("\n\nAll Alphas:\n")
    print(R, row.names = F)

    R$best <- selected
    invisible(R)

  } else {

    R$best <- selected
    return(R)
  }

}


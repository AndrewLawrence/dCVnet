
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
#' @param lambda use a fixed, user supplied lambda sequence (descending)
#'     see \code{\link[glmnet]{glmnet}}
#' @param folds This is a list where each element is an integer vector
#'     of length *n_cases*. The integer for each case labels it as belonging
#'     to a fold *1:n_folds*. This argument overrides the number of repeats
#'     and the k in repeated k-fold cv.
#' @param nreps if folds are not specified, how many times to repeat k-fold
#'     cross-validation? The default (nreps=NULL) uses 5 repeats.
#' @param nfolds if folds are not specified, how many folds should be used in
#'     cross-validation? The default (nfolds=NULL) uses 10-fold.
#' @param ... arguments passed to \code{\link[glmnet]{cv.glmnet}}
#' @param debug return an unmerged list.
#'
#' @return a \code{\link{cv.glmnet}} object with cv performance averaged.
#'
#' @seealso \code{\link{amalgamate_cv.glmnet}}
#' @examples
#' \dontrun{
#' data("CoxExample", package = "glmnet") # x and y
#' mod <- repeated.cv.glmnet(x = x, y = y, family = "cox")
#' }
#' @export
repeated.cv.glmnet <- function(x, y,
                               folds = NULL,
                               lambda = NULL,
                               alpha = 1,
                               nfolds = NULL,
                               nreps = NULL,
                               family = c("gaussian", "binomial",
                                          "poisson", "multinomial",
                                          "cox", "mgaussian"),
                               ...,
                               debug = FALSE) {
  cl <- as.list(match.call(expand.dots = TRUE))[-1]

  if ( "relax" %in% names(cl) ) {
    if ( identical(cl$relax, TRUE) ) {
      stop("The dCVnet package does not support relaxed fit models.")
    }
  }
  # We typically want to use fixed folds and fixed lambda sequence, but for
  #   convenience/generality include fallback modes (with warnings):
  if ( is.null(lambda) || missing(lambda) ) {
    warning("no lambda sequence provided: extracting glmnet default")
    # get elements of call suitable for glmnet:
    cl.gnet <- cl[names(cl) %in% methods::formalArgs(glmnet::glmnet)]
    # extract lambda list
    lambdaseq <- do.call(glmnet::glmnet, cl.gnet)$lambda # nolint
  } else {
    lambdaseq <- lambda
  }
  cl$lambda <- quote(lambdaseq)
  cl$x <- quote(x)
  cl$y <- quote(y)
  cl$family <- quote(family)

  if ( missing(folds) ) {  # nolint
    if ( missing(nfolds) ) nfolds <- 10
    if ( missing(nreps) ) nreps <- 5
    warning(paste0("no folds provided: generating an unstratified ",
                   nreps, "x", nfolds, "fold scheme"))

    folds <- lapply(seq.int(nreps),
                    function(i) {
                      caret::createFolds(y = rep(".", times = NROW(x)),
                                         k = nfolds,
                                         list = FALSE)
                    })
  }

  # prepare a 'safe' call to cv.glmnet:
  cvgnet_args <- unique(c(methods::formalArgs(glmnet::cv.glmnet),
                          methods::formalArgs(glmnet::glmnet)))
  cvgnet_args <- cvgnet_args[!(cvgnet_args %in% "...")]
  cl.cvgnet <- cl[names(cl) %in% cvgnet_args] # nolint

  # estimate models over folds:
  models <- lapply(seq_along(folds), function(i) {
    cl.cvgnet$foldid <- quote(folds[[i]])
    return(set_glmnet_alpha(
      do.call("cv.glmnet", cl.cvgnet),
      setalpha = alpha))
  } )
  # merge and return:
  if (debug) return(models)
  return(amalgamate_cv.glmnet(models))
}


# Multiple Alpha Functions: --------------------------------------------------

#' multialpha.repeated.cv.glmnet
#'
#' Runs \code{\link{repeated.cv.glmnet}} for a list of alpha values and
#'     returns averaged results, selects the 'best' alpha.
#'     One key difference between (repeated.)cv.glmnet and this function is
#'     that a single 'best' lambda/alpha combination is identified
#'     based on opt.lambda.type.
#'     *This is intended to be a dCVnet internal function*
#' @inheritParams repeated.cv.glmnet
#' @param lambdas a list of lambda sequence lists
#'     (corresponding to alphas given in alphalist)
#' @param alphalist a vector of alpha values to search.
#' @param k the number of folds for k-fold cross-validation.
#' @param nrep the number of repetitions
#'
#' @param opt.lambda.type Method for selecting optimum lambda. One of
#'                         \itemize{
#'                         \item{\code{"min"} - returns the lambda with best
#'                         CV score.}
#'                         \item{\code{"1se"} - returns the +1 se lambda}
#'                         }
#' @param opt.ystratify Boolean.
#'     Outer and inner sampling is stratified by outcome.
#'     This is implemented with \code{\link[caret]{createFolds}}
#' @param opt.uniquefolds Boolean.
#'     In most circumstances folds will be unique. This requests
#'     that random folds are checked for uniqueness in inner and outer loops.
#'     Currently it warns if non-unique values are found.
#' @param opt.keep_models The models take up memory. What should we return?
#'     \itemize{
#'       \item{ best - model with the alpha value selected as optimal. }
#'       \item{ none - no models, just cv results. }
#'       \item{ all - list of models at all alphas. }
#'     }
#'
#' @return an object of class \code{\link{multialpha.repeated.cv.glmnet}}.
#'     Containing: \itemize{
#'     \item{results - merged \code{\link{repeated.cv.glmnet}} with
#'         additional columns indicating *alpha* and logical for *best* overall}
#'     \item{best - best selected row from results}
#'     \item{folds - record of folds used}
#'     \item{models - models requested by opt.keep_models.}
#'     \item{bestmodel - index of the best model such that
#'         \code{models[[bestmodel]]} returns the model selected as optimal.}
#'     }
#' @seealso \code{\link{repeated.cv.glmnet}}
#' @export
multialpha.repeated.cv.glmnet <- function(
  x, y,
  alphalist = round(seq(0.2, 1, len = 6) ^ exp(1), 2),
  lambdas = NULL,
  k = 10,
  nrep = 5,
  opt.lambda.type = c("min", "1se"),
  opt.ystratify = TRUE,
  opt.uniquefolds = FALSE,
  family,
  opt.keep_models = c("best", "none", "all"),
  ...) {

  cl <- as.list(match.call())[-1]
  opt.lambda.type <- match.arg(opt.lambda.type)
  opt.keep_models <- match.arg(opt.keep_models)

  # check & name alphas:
  alphalist <- parse_alphalist(alphalist)

  # We typically use a fixed lambda sequence over all folds of the
  #   outer CV, but also include a fallback mode:
  if ( missing(lambdas) || is.null(lambdas) ) {
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
  cl.rcvglm.argset <- unique(
    c(methods::formalArgs(repeated.cv.glmnet),
      methods::formalArgs(glmnet::glmnet),
      methods::formalArgs(glmnet::cv.glmnet))
  )
  cl.rcvglm <- cl[names(cl) %in% cl.rcvglm.argset]
  cl.rcvglm$folds <- folds

  # run for each alpha/lambdas pair:
  malist <- lapply(seq_along(alphalist),
                   function(i) {
                     if ( getOption("mc.cores", default = 1) == 1 ) {
                       cat(paste("\tInner Alpha", i, "of",
                                 length(alphalist), Sys.time(), "\n"))
                     }

                     cl.rcvglm$alpha <- alphalist[[i]]
                     cl.rcvglm$lambda <- lambdas[[i]]

                     repeated <- do.call("repeated.cv.glmnet", cl.rcvglm)
                     return(repeated)
                   })
  names(malist) <- base::prettyNum(alphalist)

  # assemble results:
  # extract 1 model for each alpha:
  #   (N.b. model is fit to all data passed to inner-loop,
  #     not just the inner-loop training set)
  mods <- lapply(malist, "[[", "glmnet.fit")
  # extract names:
  tmeas <- names(malist[[1]]$name)
  # get the summary of glmnet::cv.glmnet performance and merge into data.frame:
  malist <- mapply(cv.glmnet.modelsummary,
                   mod = malist,
                   alpha = alphalist,
                   SIMPLIFY = FALSE)
  malist <- as.data.frame(data.table::rbindlist(malist),
                          stringsAsFactors = FALSE)

  # Note that mods contains models fit to the complete dataset received by
  #   this function.

  # pick the optimal alpha:
  bestfun <- ifelse(tmeas == "auc", max, min)
  lselector <- ifelse(opt.lambda.type == "min", "lambda.min", "lambda.1se")
  # lselector is the name of a column containing logical values in malist:
  bestcandidates <- malist[malist[[lselector]], ]
  best <- bestcandidates[bestcandidates$cvm == bestfun(bestcandidates$cvm), ]

  # ties in the objective function are broken by smaller cvsd
  #     followed by sparser solutions (unlikely)
  if ( nrow(best) > 1 ) {
    best <- best[order(best$cvsd, best$nzero), ]
    best <- best[1, ]
  }

  # add column to the multialpha list.
  malist$best <- (malist$alpha == best$alpha) & (malist$lambda == best$lambda)

  # initialise the return:
  R <- list(results = malist,
            best = best,
            alphas = alphalist,
            folds = folds)

  if ( opt.keep_models == "all" ) {
    R <- append(R,
                values = list(models = mods,
                              bestmodel = which(alphalist == best$alpha)))
  } else {
    R <- append(R,
                values = list(models = mods[which(alphalist == best$alpha)],
                              bestmodel = 1))
  }

  R <- structure(R,
                 class = "multialpha.repeated.cv.glmnet",
                 type.measure = tmeas,
                 family = family,
                 type.lambda = lselector,
                 opt.keep_models = opt.keep_models)
  return(R)
}

#' @export
print.multialpha.repeated.cv.glmnet <- function(x, ...) {

  type.lambda <- attr(x, "type.lambda")
  opt.keep_models <- attr(x, "opt.keep_models")

  cat("A dCVnet::multialpha.repeated.cv.glmnet object\n\n")
  if ( opt.keep_models == "all" ) {
    cat("(includes fitted models for each alpha)\n\n")
  }
  if ( opt.keep_models == "none" ) {
    cat("(includes no fitted models, only cv results)\n\n")
  }
  if ( opt.keep_models == "best" ) {
    cat("(includes fitted models at optimum lambda/alpha)\n\n")
  }
  cat(paste0("Model family:\t\t", attr(x, "family"), "\n"))
  cat(paste0("Tuning metric (cvm):\t", attr(x, "type.measure"), "\n"))
  cat(paste0("Lambda Selection:\t", type.lambda, "\n"))
  cat("\n")

  best <- x$best
  x <- x$results

  alpha <- unique(x$alpha)
  lambdas <- lapply(alpha, function(i) x$lambda[x$alpha == i]  )

  l_lengths <- vapply(lambdas, length, FUN.VALUE = c(1.0))
  l_min <- prettyNum(vapply(lambdas, min, FUN.VALUE = c(1.0)))
  l_max <- prettyNum(vapply(lambdas, max, FUN.VALUE = c(1.0)))
  l_ranges <- paste(l_min, l_max, "\n")

  rcols <- c("s", "alpha", "lambda", "nzero", "cvm", "cvsd", "cvup", "cvlo")

  R <- x[x[[type.lambda]], rcols]
  attr(R, "class") <- "data.frame"

  selected <- (R$alpha == best$alpha) & (R$lambda == best$lambda)
  R$best <- selected

  cat(paste0("    ", length(alpha), " alpha(s): ",
             paste(alpha, collapse = ", "), "\n\n"))
  cat(paste0("    Lambda Counts:\n"))
  cat(paste0("    ", l_lengths))
  cat(paste0("\n\n    Lambda Ranges:\n"))
  cat(paste0("\t", l_ranges))
  cat("\n")

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

  if ( print ) {
    cat("Summary of ")
    print(object)
  }

  type.lambda <- attr(object, "type.lambda")

  best <- object$best
  object <- object$results

  rcols <- c("alpha", "lambda", "nzero", "cvm", "cvsd", "cvup", "cvlo")

  R <- object[object[[type.lambda]], rcols]

  selected <- (R$alpha == best$alpha) & (R$lambda == best$lambda)

  if ( print ) {

    cat("\nBest fitting alpha:\n")
    print(best[, rcols], row.names = FALSE)

    cat("\n\nBest per-alpha:\n")
    print(R, row.names = FALSE)

    R$best <- selected
    invisible(R)

  } else {

    R$best <- selected
    return(R)
  }

}

#' predict.multialpha.repeated.cv.glmnet
#'
#' obtain predictions from a
#'     \code{\link{multialpha.repeated.cv.glmnet}} object.
#'     Uses the 'best' alpha \& lambda hyperparameters determined by the
#'     internal cross-validation results. For lambda this will be lambda.min
#'     or lambda.1se (determined at model runtime) unless requested otherwise.
#'
#' @param object a a \code{\link{multialpha.repeated.cv.glmnet}} object.
#' @param newx matrix of new values for x at which predictions are required.
#'     Note: no impact when type is "coefficients", "class" or "nonzero".
#' @param alpha the penalty type alpha at which prediction is required.
#'     Leave NULL to use the cv-optimised value.
#' @param s The penalty amount lambda at which prediction is required.
#'     Leave NULL to use the cv-optimised value.
#' @param ... passed to \code{\link[glmnet]{predict.glmnet}}
#'
#' @seealso \code{\link[glmnet]{predict.cv.glmnet}},
#'     \code{\link[glmnet]{predict.glmnet}}
#'
#' @export
predict.multialpha.repeated.cv.glmnet <- function(object,
                                                  newx,
                                                  alpha = NULL,
                                                  s = NULL,
                                                  ...) {
  # function to return "best" predictions from a multi-alpha object
  if ( attr(object, "opt.keep_models") == "none" ) {
    stop(paste0("The object ", deparse(substitute(object)),
                " does not include the fitted models required for predict.\n",
                "Rerun multialpha.repeated.cv.glmnet with ",
                "opt.keep_models = best",
                "or opt.keep_models = all",
                "to make predictions"))
  }
  if ( !missing(alpha) && attr(object, "opt.keep_models") == "best" ) {
    stop(paste0("The object ", deparse(substitute(object)),
                " does not include all fitted models required for predict.\n",
                "Rerun multialpha.repeated.cv.glmnet with ",
                "opt.keep_models = all",
                "in order to predict at the non-best alpha."))
  }

  # if an alpha value was specified:
  if ( !missing(alpha) ) {
    if ( alpha %in% object$alphas ) {
      mod <- object$models[[ which(alpha %in% object$alphas) ]]
      if ( missing(s) ) {
        type.lambda <- attr(object, "type.lambda")
        sel <- (object$results$alpha == alpha) &
          (object$results[[type.lambda]])
        if ( sum(sel) != 1 ) stop("ERROR: this should be unique.")
        s <- object$results$lambda[sel]
      }
      return(predict(mod, s = s, newx = newx, ...))
    } else {
      stop("Requested alpha not present")
    }
  }

  # otherwise we will use best alpha:
  mod <- object$models[[object$bestmodel]]

  # if s missing use best s:
  if ( missing(s) ) {
    s <- object$best$lambda
  } else {
    warning(paste("alpha is selected optimally,",
                  "but lambda is manually specified.",
                  "manual alpha: ", s,
                  "best alpha:", object$best$lambda
                  ))
  }

  # run the prediction:
  predict(mod, s = s, newx = newx, ...)
}



#' coef.multialpha.repeated.cv.glmnet
#'
#' obtain coefficients from a
#'     \code{\link{multialpha.repeated.cv.glmnet}} object.
#'     This is a wrapper for \code{\link{predict.multialpha.repeated.cv.glmnet}}
#'
#' @inheritParams predict.multialpha.repeated.cv.glmnet
#' @seealso \code{\link{predict.multialpha.repeated.cv.glmnet}}
#'     \code{\link[glmnet]{predict.cv.glmnet}},
#'     \code{\link[glmnet]{predict.glmnet}}
#'
#' @export
coef.multialpha.repeated.cv.glmnet <- function(object,
                                               newx = NULL,
                                               alpha = NULL,
                                               s = NULL,
                                               ...) {
  # function to return coefficients from a multi-alpha object
  cl <- as.list(match.call())[-1]
  cl$object <- object # pass object not name
  cl$type <- "coefficients"
  cl$newx <- NA # newx not required for type = coefficients.


  do.call("predict.multialpha.repeated.cv.glmnet", args = cl)
}

# Remove models from a multialpha.repeated.cv.glmnet object, but retain the
#   attributes.
#   This is a workaround of an R quirk: `[` drops most attributes (see ?Extract)
#   Packages like vctrs/sticky attempt to address this more formally.
drop_models.multialpha.repeated.cv.glmnet <- function(object) {
  attr <- attributes(object)
  R <- structure(object[!names(object) %in% c("models", "bestmodel")],
                 class = attr$class,
                 type.measure = attr$type.measure,
                 family = attr$family,
                 type.lambda = attr$type.lambda,
                 opt.keep_models = "none")
  return(R)
}

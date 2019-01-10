
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
  #optfun <- ifelse("auc" %in% measure_name, max, min)

  # Which is the optimum lambda?
  av$lambda.min <- F
  av$lambda.min[match(theLambda, av$lambda)] <- TRUE
#
#   # Deal with ties (multiple lambdas with the same minimum).
#   nties <- sum(av$lambda.min)
#   if ( nties > 1 ) {
#     tiedlambdas <- av$lambda[av$lambda.min]
#     # For odd length vectors the median exists in the data.
#     #   otherwise return the larger lambda of the middle two.
#     if ( length(tiedlambdas) %% 2 == 1 ) {
#       av$lambda.min <- av$lambda == median(tiedlambdas)
#     } else {
#       larger_middle_lambda <- (length(tiedlambdas) / 2) + 1
#       av$lambda.min <- av$lambda == sort(tiedlambdas)[larger_middle_lambda]
#     }
#   }

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
  if ( identical(opt.ystratify, FALSE) ) { ystrat <- rep("x", length(y)) }

  folds <- lapply(1:nrep, function(i) {
    caret::createFolds(y = y, k = k, list = FALSE, returnTrain = FALSE)
  })

  if ( identical(opt.uniquefolds, TRUE) ) checkForDuplicateCVFolds(folds)

  malist <- lapply(1:length(alphalist),
                   function(i) {
                     if ( getOption("mc.cores", default = 1) == 1 ) {
                       cat(paste("\tInner Alpha", i, "of",
                                 length(alphalist), Sys.time(), "\n"))
                     }
                     a <- alphalist[[i]]

                     repeated <- repeated.cv.glmnet(alpha = a,
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

#' summary.repeated.cv.glmnet
#'
#' a summary of key options and results for a
#'     \code{\link{multialpha.repeated.cv.glmnet}} object.
#'
#' @param object a a \code{\link{multialpha.repeated.cv.glmnet}} object.
#' @param print if FALSE silently returns the summary results table.
#' @param ... NULL
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


# Outer Loop functions ----------------------------------------------------

#' Fit a doubly cross-validated elastic-net regularised logistic regression
#'
#' Double (nested) repeated k-fold cross-validation used to:
#'     \itemize{
#'     \item{Produce unbiased estimates of out-of-sample
#'         classification performance (outer CV).}
#'     \item{Select optimal hyperparameters for the elasticnet (inner CV).}}
#'     Elasticnet hyperparameters are
#'     \bold{lambda} (the total regularisation penalty)
#'     and \bold{alpha} (the balance of L1 and L2 regularisation types).
#'
#' @inheritParams parse_dCVnet_input
#' @inheritParams multialpha.repeated.cv.glmnet
#'
#' @param k_inner an integer, the k in the inner k-fold CV.
#' @param k_outer an integer, the k in the outer k-fold CV.
#' @param nrep_inner an integer, the number of repetitions (k-fold inner CV)
#' @param nrep_outer an integer, the number of repetitions (k-fold outer CV)
#' @param alphalist a numeric vector of values in \(0,1\].
#'     This sets the search space for optimising hyperparameter alpha.
#' @param nlambda an integer, number of gradations between
#'     lambda.min and lambda.max to search.
#'     See \code{glmnet} argument \code{nlambda}.
#' @param type.measure passed to \code{\link[glmnet]{cv.glmnet}}.
#'     This sets the metric used for hyperparameter optimisation in the
#'     inner cross-validation. Options: \code{"deviance"}, \code{"class"},
#'     \code{"mse"}, \code{"mae"}
#' @param opt.empirical_cutoff Boolian.
#'     Use the empirical proportion of cases as the cutoff for outer CV
#'     classification (affects outer CV performance only).
#'     Otherwise classify at 50\% probability.
#' @param ... Arguments to pass through to cv.glmnet
#'     (beware, making use of this may break things).
#' @return a dCVnet object.
#' @examples
#' \dontrun{
#' iris_class <- dCVnet(f = Species ~ Sepal.Length + Sepal.Width +
#'                          Petal.Length + Petal.Width,
#'                      data = subset(iris, iris$Species != "versicolor"),
#'                      alphalist = 0.5)
#' #Note: in most circumstances larger values of nrep_inner and nrep_outer
#' #      will be required.
#' summary(classperformance(iris_class))
#' plot(iris_class)
#' }
#' @importFrom stats aggregate as.formula coef glm model.frame model.matrix
#' @importFrom stats predict sd terms var
#' @export
dCVnet <- function(
  f,
  data,

  nrep_outer = 2,
  k_outer = 10,
  nrep_inner = 5,
  k_inner = 10,
  alphalist = c(0.2, 0.5, 0.8),
  nlambda = 100,
  type.measure = "deviance",
  positive = 1,

  opt.lambda.type = c("minimum", "se", "percentage"),
  opt.lambda.type.value = 1.0,
  opt.empirical_cutoff = FALSE,
  opt.uniquefolds = FALSE,
  opt.ystratify = TRUE,

  ...) {

  thecall <- match.call()

  time_start <- Sys.time() # for logging.

  parsed <- parse_dCVnet_input(f = f, data = data, positive = positive)
  x <- parsed$x_mat
  y <- parsed$y

  min_lambda_ratio <- ifelse(ncol(parsed$x_mat) > nrow(parsed$x_mat),
                             0.01, 0.0001)

  if ( type.measure == "auc" ) {
    # magic number for innerloop is:
    #   Ncases * outer train proportion * inner test proportion
    auc_magic <- (nrow(parsed$x_mat) * (1 - (1 / k_outer)) * (1 / k_inner))
    if (auc_magic < 11) {
      warning(
        paste("AUC is not possible due to small sample size!\n
              Estimated inner foldsize =", auc_magic,
              "\nUsing 'class' instead."))
      type.measure <- "class"
    }
  }

  # Check the alpha values:
  alphalist <- parse_alphalist(alphalist)
  # how many alphas did we feed in:
  nalpha <- length(alphalist)

  # Print some general run info to screen.
  startup_message(k_inner = k_inner, nrep_inner = nrep_inner,
                  k_outer = k_outer, nrep_outer = nrep_outer,
                  nalpha = nalpha, nlambda = nlambda,
                  parsed = parsed, time.start = time_start)
  # Main work starts here:

  # Step 1: make repeated outer folds in x & y according to nrep_outer, k_outer.
  #
  # Note: this is default stratified by y, we obtain unstratified sampling
  #         by giving caret::createMultiFolds a single-level factor/char.
  ystrat <- y
  if ( identical(opt.ystratify, FALSE) ) { ystrat <- rep("x", length(y)) }
  outfolds <- caret::createMultiFolds(y = ystrat,
                                      k = k_outer,
                                      times = nrep_outer)
  rm(ystrat)
  # Warn if repeated folds were obtained.
  if ( opt.uniquefolds ) checkForDuplicateCVFolds(outfolds)

  imax <- length(outfolds)
  names(outfolds) <- paste0("Out", names(outfolds)) # give names

  cat("Lambda rangefinding...\n")
  outmaxlambdas <- lambda_rangefinder(
    y = y, x = x,
    alphalist = alphalist,
    prop = 1 * (1 - (1 / k_outer)) * (1 - (1 / k_inner)) )

  # define the complete list of lambdas for each alpha
  lambdas <- lapply(outmaxlambdas,
                    function(maxlambda) {
                      exp(seq(from = log(maxlambda),
                              to =   log(maxlambda * min_lambda_ratio),
                              length.out = nlambda))
                    })

  # are we working with empirical cutoffs?
  cutoff <- 0.5
  if ( opt.empirical_cutoff ) {
    cutoff <- (as.numeric(table(y)[1]) / sum(as.numeric(table(y))))
  }
  cat(paste("Cutoff: ", cutoff, "\n"))

  # Main outer loop - for each outer fold run an inner cv loop
  outers <- parallel::mclapply(
    seq(along = outfolds),
    mc.cores = getOption("mc.cores", 1L),
    function(i) {
      if ( getOption("mc.cores", 1L) > 1 ) {
        cat(paste0("Outerloop:", i, " of ", imax, "\n"))
      } else {
        cat(paste0("\nOuterloop:", i, " of ", imax, "\n"))
      }

      of <- outfolds[[i]]

      # Preprocessing
      # first split the data into train and test data according to the fold.
      sel_train <- 1:nrow(x) %in% of
      sel_test <- !sel_train

      trainx <- subset(x, sel_train)
      trainy <- subset(y, sel_train)

      testx <- subset(x, sel_test)
      testy <- subset(y, sel_test)

      # next apply scaling (based just on the train data,
      #                     but applied to the test data)
      PPx <- caret::preProcess(trainx)
      trainx <- predict(PPx, trainx)
      testx <- predict(PPx, testx)

      # Run the tuning for the training data:
      inners <- multialpha.repeated.cv.glmnet(nrep = nrep_inner,
                                              k = k_inner,
                                              alphalist = alphalist,
                                              lambdas = lambdas,
                                              y = trainy, x = trainx,
                                              type.measure = type.measure,
                                              standardize = F,
                                              opt.lambda.type = opt.lambda.type,
                                              opt.lambda.type.value = opt.lambda.type.value,
                                              opt.ystratify = opt.ystratify,
                                              opt.uniquefolds = opt.uniquefolds,
                                              ...)

      # extract the best alpha/lambda based on out of sample performance:
      fit_alpha <- as.numeric(inners$inner_best$alpha)
      fit_lambda <- inners$inner_best$lambda

      # fit a model to ALL the train data using the selected alpha/lambda:
      model <- glmnet::glmnet(x = trainx,
                              y = trainy,
                              family = "binomial",
                              alpha = fit_alpha,
                              standardize = F,
                              lambda = lambdas[[which(alphalist == fit_alpha)]])

      # how well did it do?:
      newx_prediction <- predict(model,
                                 newx = testx,
                                 type = "response",
                                 s = fit_lambda)
      # Format results into a data.frame:
      newx_performance <- data.frame(
        rowid = rownames(newx_prediction),     # row references
        foldid = names(outfolds)[[i]],         # fold names
        label = strsplit(names(outfolds)[[i]],   # repetition groups
                         split = ".", fixed = T)[[1]][2],
        reference = testy,
        probability = c(newx_prediction),
        classification = factor(c(newx_prediction) > cutoff,
                                levels = c(F, T),
                                labels = model$classnames),
        stringsAsFactors = F)

      rownames(newx_performance) <- rownames(newx_prediction)

      return(list(tuning = inners,
                  model = model,
                  performance = structure(newx_performance,
                                          class = c("classperformance",
                                                    "data.frame"))))
    }
  )
  names(outers) <- names(outfolds)

  # Gather performance from the outers into a single dataframe.
  performance <- do.call(rbind, lapply(outers, "[[", "performance"))
  rownames(performance) <- NULL
  performance <- structure(performance,
                           class = c("classperformance",
                                     "data.frame"))
  # and remove it in the source:
  for (ii in 1:length(outers)) {
    outers[[ii]]$performance <- NULL
  }

  # Final model:
  final_PPx <- caret::preProcess(x)   # preprocessed data (centered & scaled)
  xs <- predict(final_PPx, x)

  # Run the inner loop for the final model:
  cat("\n\nFinal Model\n")

  final_tuning <- multialpha.repeated.cv.glmnet(
    alphalist = alphalist,
    lambdas = lambdas,
    k = k_inner,
    nrep = nrep_inner,
    y = y, x = xs,
    type.measure = type.measure,
    opt.lambda.type = opt.lambda.type,
    opt.lambda.type.value = opt.lambda.type.value,
    opt.ystratify = opt.ystratify,
    opt.uniquefolds = opt.uniquefolds
  )

  final_model <- glmnet(
    x = xs,
    y = y,
    family = "binomial",
    alpha = final_tuning$inner_best$alpha,
    lambda = lambdas[[which(alphalist == final_tuning$inner_best$alpha)]])

  final_prediction <- predict(final_model,
                              newx = xs,
                              type = "response",
                              s = final_tuning$inner_best$lambda)

  final_performance <- data.frame(
    rowid = rownames(final_prediction),
    label = "Final",
    foldid = "None",
    reference = y,
    probability = c(final_prediction),
    classification = factor(c(final_prediction) > cutoff,
                            levels = c(F, T),
                            labels = final_model$classnames))

  rownames(final_performance) <- rownames(final_prediction)

  final_performance <- structure(final_performance,
                                 class = c("classperformance",
                                           "data.frame"))

  final <- list(tuning = final_tuning,
                performance = final_performance,
                model = final_model)

  time_stop <- Sys.time()
  run_time_mins <- as.numeric(round( (time_stop - time_start) / 60, 2))

  # The final object:
  obj <- structure(list(tuning = outers,
                        performance = performance,
                        folds = outfolds,
                        final = final,
                        input = list(call = thecall,
                                     runtime.mins = run_time_mins,
                                     parsed = parsed,
                                     lambdas = lambdas,
                                     alphas = alphalist)),
                   class = c("dCVnet", "list"))

  cat("\n\n")
  print(thecall)
  cat("\n")
  cat(paste0("Finished.\n"))
  cat(paste0(time_stop, "\n"))
  cat(paste0("Runtime: ", run_time_mins, " mins\n"))
  cat("------------")

  return(obj)
}

# Extract logreg coefficients from the outerloop best models:
#' coef.dCVnet
#'
#' Coefficients from a dCVnet object.
#'
#' @param object a dCVnet object
#' @param type how to return coefficients.
#'     \itemize{
#'     \item{\code{"all"} - return separate coefficients for each rep/fold.}
#'     \item{\code{"rep"} - return separate coefficients for each rep
#'         (mean average over folds).}
#'     \item{\code{"mean"} - return mean of \code{"rep"}.}
#'     \item{\code{"median"} - return median of \code{"rep"}.}
#'     }
#' @param ... " "
#'
#' @name coef.dCVnet
#'
#' @export
coef.dCVnet <- function(object, type = "all", ...) {
  # Type can be:
  #   all - coefficients for each rep/fold.
  #   rep - mean average per rep.
  #   mean - mean of per-rep means.
  #   median = median of per-rep means.

  # next handle dCVnet cases:

  # works on the output of the outerloop.
  #   Given the 'best' alpha/lambda from the inner loop,
  #     what were the coefficients in a model selected with the best alpha.
  R <- lapply(1:length(object$tuning), function(ii) {
    tt <- object$tuning[[ii]]
    coefs <- as.matrix(predict(tt$model,
                               type = "coef",
                               s = tt$tuning$inner_best$lambda))
    RR <- data.frame(Predictor = rownames(coefs),
                     Coef = c(coefs),
                     fold = names(object$folds)[[ii]],
                     stringsAsFactors = F)
    return(RR)
  })
  R <- do.call(rbind, R)

  R$Predictor <- factor(R$Predictor,
                        levels = c("(Intercept)",
                                   object$tuning[[1]]$model$beta@Dimnames[[1]]))
  R$fold <- factor(R$fold, levels = unique(R$fold))

  R$Rep <- sapply(strsplit(as.character(R$fold), split = "\\."), "[", 2)

  if (type == "all") return(R)
  # Aggregate by Repetition:
  R <- aggregate(R$Coef,
                 by = list(Predictor = R$Predictor,
                           Rep = R$Rep),
                 mean)
  names(R)[3] <- "Coef"
  return(switch(type,
                rep = R,
                mean = stats::setNames(
                  aggregate(R$Coef,
                            by = list(Predictor = R$Predictor),
                            mean),
                  c("Predictor", "Coef")),
                median = stats::setNames(
                  aggregate(R$Coef,
                            by = list(Predictor = R$Predictor),
                            stats::median),
                  c("Predictor", "Coef")),
                stop(paste("type:",
                           type,
                           " - should be one of: all, rep, mean, median"))
  ))
}


coefficients.dCVnet <- function(object, ...) coef.dCVnet(object, ...)


#' @export
print.dCVnet <- function(x, ...) {

  parsed <- x$input$parsed

  stab <- table(parsed$y)

  nalphas <- length(x$input$alphas)
  nlambdas <- length(x$input$lambdas[[1]])

  typemeas <- attr(x$final$tuning$inner_results, "type.measure")

  outerfolds <- x$folds
  innerfolds <- x$tuning[[1]]$tuning$inner_folds

  outerkey <- data.frame(do.call(rbind,
                                 strsplit(names(outerfolds),
                                          split = ".", fixed = T)))
  colnames(outerkey) <- c("Folds", "Reps")

  outer_k <- length(unique(outerkey$Folds))
  outer_nrep <- length(unique(outerkey$Reps))

  inner_k <- max(innerfolds[[1]])
  inner_nrep <- length(innerfolds)

  # What options were specified?:
  opts.empirical <- x$input$call[["opt.empirical_cutoff"]]
  if ( is.null(opts.empirical) ) opts.empirical <- FALSE

  opts.ystratify <- x$input$call[["opt.ystratify"]]
  if ( is.null(opts.ystratify) ) opts.ystratify <- TRUE

  opts.uniquefolds <- x$input$call[["opt.uniquefolds"]]
  if ( is.null(opts.uniquefolds) ) opts.uniquefolds <- FALSE

  cat("A dCVnet object, from the dCVnet Package\n\n")
  print(x$input$call)
  cat("\n")
  cat(paste0("Runtime: ", as.numeric(x$input$runtime.mins), " mins\n\n"))

  cat("Input\n-----\n")
  cat("Features:\n")
  cat(paste0("\t", ncol(parsed$x_mat), "\n"))

  cat("Observations:\n")
  cat(paste0("\t", nrow(parsed$x_mat), " subjects\n"))
  cat(paste0("\t", stab[1], " of outcome: ", names(stab)[1], "\n"))
  cat(paste0("\t", stab[2], " of outcome: ", names(stab)[2], "\n"))

  cat("Hyperparameter Tuning:\n")
  cat(paste("\tOptimising: ", typemeas, "\n"))
  cat(paste0("\t", nalphas, " alpha values:\t"))
  cat(x$input$alphas)
  cat(paste0("\n\t", nlambdas, " lambda values\n"))

  cat("Cross-validation:\n")
  cat(paste("\tInner:\n"))
  cat(paste("\t\tk =\t", inner_k, "\n"))
  cat(paste("\t\tnrep =\t", inner_nrep, "\n"))
  cat(paste("\tOuter:\n"))
  cat(paste("\t\tk =\t", outer_k, "\n"))
  cat(paste("\t\tnrep =\t", outer_nrep, "\n"))

  cat("Options:\n")
  cat(paste("\tUse Empirical Thresholding for LogReg Classification: ",
            opts.empirical, "\n"))

  cat(paste("\tStratify k-fold sampling by outcome: ",
            opts.ystratify, "\n"))
  cat(paste("\tCheck random folds are unique: ",
            opts.uniquefolds, "\n\n"))

  invisible(x)
}



#' summary.dCVnet
#'
#' a summary of key options and results for a
#'     \code{\link{dCVnet}} object.
#'
#' Prints the following sections: \itemize{
#' \item{Input descriptives}
#' \item{Outer Loop classification performance}
#' \item{Inner Loop model stability}
#' \item{Inner Loop model performance}
#' \item{Production model performance (fit to *all* data)}
#' }
#'     The outer loop tests generalisation of the model.
#'
#' @param object a a \code{\link{dCVnet}} object.
#' @param ... NULL
#'
#' @return a dataframe of tuning results for each outer fold/rep.
#'     \[Note: this may change in future\]
#'
#' @export
summary.dCVnet <- function(object, ...) {
  .titlecat <- function(title) {
    n <- nchar(title)
    divider <- paste0(rep("-", length.out = n), collapse = "")
    cat(paste0(divider, "\n"))
    cat(paste0(title, "\n"))
    cat(paste0(divider, "\n"))
  }
  # Start with print.
  cat("Summary of ")
  print(object)

  # Outerloop CV results:
  outcv <- report_classperformance_summary(object)

  min_vars <- c("Accuracy", "Sensitivity",
                "Specificity", "Balanced Accuracy",
                "AUROC")

  outernreps <- length(unique(object$performance$label))

  if ( outernreps == 1 ) {
    min_outcv <- outcv[outcv$Measure %in% min_vars, "Rep1", drop = F]
    colnames(min_outcv) <- "Value"
  } else {
    min_outcv <- outcv[outcv$Measure %in% min_vars,
                       c("mean", "sd", "min", "max")]
  }

  min_outcv <- round(min_outcv, 3)
  row.names(min_outcv) <- min_vars

  # What do the 'best-fitting' results of the inner loops look like:
  R <- lapply(object$tuning, function(x) {
    summary.multialpha.repeated.cv.glmnet(x$tuning, print = F)
  })

  R <- data.frame(do.call(rbind, R))
  R <- R[R$best, names(R)[!names(R) %in% "best"]]
  rownames(R) <- names(object$tuning)

  R$Rep <- sapply(strsplit(rownames(R), split = ".", fixed = T), "[", 2)

  # summarise lambdas:
  lambda_summary <- c(mean = mean(R$lambda),
                      sd = sd(R$lambda),
                      min = min(R$lambda),
                      max = max(R$lambda))

  # Summarise inner loop performances (i.e. best cvm)  by rep:
  cvm_repmean <- aggregate(R$cvm, by = list(Rep = R$Rep), FUN = mean)

  # summarise 'final' model performance:
  fmp <- summary(object$final$performance)
  min_fmp <- fmp$Value[fmp$Measure %in% min_vars]
  min_fmp <- round(min_fmp, 3)
  names(min_fmp) <- min_vars

  # fInal model hyper parameters:
  fmp_hp <- summary(object$final$tuning, print = F)
  fmp_hp <- fmp_hp[fmp_hp$best, ]

  fmp_hp_str <- c(alpha = fmp_hp$alpha,
                  lambda = formatC(fmp_hp$lambda),
                  cvm = formatC(fmp_hp$cvm))

  # Write this out:
  cat("\n")
  .titlecat("Outer Loop CV Performance")
  print(min_outcv)
  cat("\n")
  .titlecat("Inner Loop Model Stability")
  cat(paste0("Tuned alphas (table):"))
  print(table(R$alpha))
  cat("Tuned lambdas (descriptives):\n")
  print(lambda_summary, digits = 3)
  cat("\n")

  .titlecat("Inner Loop Model Performance")
  cat("[Note: inner-loop performance is not independently cross-validated]\n\n")
  cat(paste0("metric (cvm): ",
             attr(object$final$tuning$inner_best, "type.measure"), "\n\n"))
  cat("cvm summary (per-rep, averaged over folds):\n")
  print(summary(cvm_repmean$x))
  cat("cvm summary (over folds):\n")
  print(summary(R$cvm))
  cat("\n")

  .titlecat("'Production' Model")
  cat("Production Performance (not cross-validated):\n")
  print(min_fmp)
  cat("Production Hyperparameter Tuning:\n")
  print(fmp_hp_str, quote = F)

  invisible(R)
}

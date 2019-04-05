
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
#' @param opt.empirical_cutoff Boolean.
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
  y,
  data,
  f = "~.", # nolint

  nrep_outer = 2,
  k_outer = 10,
  nrep_inner = 5,
  k_inner = 10,
  alphalist = c(0.2, 0.5, 0.8),
  nlambda = 100,
  type.measure = "deviance",
  family = "binomial",
  positive = 1,
  offset = NULL,

  opt.lambda.type = c("minimum", "se", "percentage"),
  opt.lambda.type.value = 1.0,
  opt.empirical_cutoff = FALSE,
  opt.uniquefolds = FALSE,
  opt.ystratify = TRUE,

  ...) {

  # Parse input -------------------------------------------------------------

  thecall <- match.call()
  callenv <- c(as.list(environment()), list(...))

  time_start <- force(Sys.time()) # for logging.
  parsed <- parse_dCVnet_input(f = f,
                               y = y,
                               data = as.data.frame(data),
                               family = family,
                               positive = positive)
  x <- parsed$x_mat
  y <- parsed$y

  if ( type.measure == "auc" ) {
    # magic number for innerloop is:
    #   Ncases * outer train proportion * inner test proportion
    auc_magic <- (nrow(x) * (1 - (1 / k_outer)) * (1 / k_inner))
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

  # are we working with empirical cutoffs?
  cutoff <- 0.5
  if ( opt.empirical_cutoff ) {
    cutoff <- (as.numeric(table(y)[1]) / sum(as.numeric(table(y))))
  }
  if ( family == "binomial" ) cat(paste("Cutoff: ", cutoff, "\n"))

  # Print some general run info to screen.
  startup_message(k_inner = k_inner, nrep_inner = nrep_inner,
                  k_outer = k_outer, nrep_outer = nrep_outer,
                  nalpha = nalpha, nlambda = nlambda,
                  parsed = parsed,
                  family = family,
                  time.start = time_start)

  # Create outer folds ------------------------------------------------------
  # Note: this is default stratified by y, we obtain unstratified sampling
  #         by giving caret::createMultiFolds a single-level factor/char.
  ystrat <- y
  if ( identical(opt.ystratify, FALSE) | family %in% c("cox", "mgaussian") ) {
    ystrat <- rep("x", length(y))
  }
  outfolds <- caret::createMultiFolds(y = ystrat,
                                      k = k_outer,
                                      times = nrep_outer)
  rm(ystrat)
  # Warn if repeated folds were obtained.
  if ( opt.uniquefolds ) checkForDuplicateCVFolds(outfolds)

  imax <- length(outfolds)
  names(outfolds) <- paste0("Out", names(outfolds)) # give names


  # Create Lambda paths -----------------------------------------------------
  .Lambdafinder <- function(x, y, family, nlambda, alphalist, ...) {
    # inspired by the source of cv.glmnet, we get the list of lambdas once,
    #   by running a glmnet at the required alpha for all data and
    #   extracting the list of lambdas used.
    # Previously a wider/more representative range of lambdas was obtained by
    #   by formulaic calculation of max lambda for random samplings of the
    #   data (per alpha).
    # However, 1) this is not done by Hastie et al who uses the run a
    #             model method in cv.glmnet.
    #          2) it was time consuming
    #          3) the calculation of maxlambda from the data is opaque
    #             and performed by the fortran code.
    #             there do not appear to be published formula for doing these
    #             calculations for mgaussian/cox family models.
    #          4) using fortran avoids rounding errors & excess precision.
    lapply(alphalist,
           function(aa) {
             m <- glmnet(x = x, y = y,
                         alpha = aa,
                         family = family,
                         nlambda = nlambda,
                         ...)
             m$lambda
           } )
  }

  lambdas <- .Lambdafinder(x = x, y = y,
                           nlambda = nlambda,
                           family = family,
                           alphalist = alphalist,
                           ...)


  # Outer Loop Start --------------------------------------------------------
  outers <- parallel::mclapply(
    seq(along = outfolds),
    mc.cores = getOption("mc.cores", 1L),
    function(i) {

      if ( getOption("mc.cores", 1L) > 1 ) {
        cat(paste0("Outerloop:", i, " of ", imax, "\n"))
      } else {
        cat(paste0("\nOuterloop:", i, " of ", imax, "\n"))
      }

      # split & preprocess data ---------------------------------------------

      of <- outfolds[[i]]

      sel_train <- seq(nrow(x)) %in% of
      sel_test <- !sel_train

      trainx <- subset(x, sel_train)
      trainy <- subset(y, sel_train)

      testx <- subset(x, sel_test)
      testy <- subset(y, sel_test)

      # scaling to mean=0, sd=1 (calculated on train data, applied to test data)
      PPx <- caret::preProcess(trainx, method = c("center", "scale"))
      trainx <- predict(PPx, trainx)
      testx <- predict(PPx, testx)


      # inner loop train data -------------------------------------------
      #   - receives no (outer) test data
      inners <- multialpha.repeated.cv.glmnet(
        nrep = nrep_inner,
        k = k_inner,
        alphalist = alphalist,
        lambdas = lambdas,
        y = trainy, x = trainx,
        type.measure = type.measure,
        family = family,
        standardize = FALSE,
        opt.lambda.type = opt.lambda.type,
        opt.lambda.type.value = opt.lambda.type.value,
        opt.ystratify = opt.ystratify,
        opt.uniquefolds = opt.uniquefolds,
        ...)

      # hyperparameter selection --------------------------------------------
      #   - based on best CV/out-of-sample performance in inner loop
      fit_alpha <- as.numeric(inners$inner_best$alpha)
      fit_lambda <- inners$inner_best$lambda

      # tuned outer-train model ---------------------------------------------
      #   - tuned wrt best inner alpha
      #   - full lambda path is fit with subsequent selection
      #       (see https://web.stanford.edu/~hastie/glmnet/glmnet_beta.html:
      #       "[single lambda fit] is not recommended and is not the spirit
      #        of the package."
      model <- glmnet::glmnet(x = trainx,
                              y = trainy,
                              family = family,
                              alpha = fit_alpha,
                              standardize = FALSE,
                              lambda = lambdas[[which(alphalist == fit_alpha)]],
                              ...)

      newx_performance <- tidy_predict.glmnet(
        mod = model,
        newx = testx,
        family = family,
        newy = testy,
        binomial_thresh = cutoff,
        offset = offset,
        label = strsplit(names(outfolds)[[i]],
                         split = ".", fixed = TRUE)[[1]][2],
        s = fit_lambda)

      return(list(tuning = inners,
                  model = model,
                  performance = structure(newx_performance,
                                          class = c("classperformance",
                                                    "data.frame"))))
    }
  )
  names(outers) <- names(outfolds)
  # Outer Loop End --------------------------------------------------------


  # Performance -------------------------------------------------------------

  # Gather performance from the outers into a single dataframe.
  performance <- do.call(rbind, lapply(outers, "[[", "performance"))
  rownames(performance) <- NULL
  performance <- structure(performance,
                           class = c("classperformance",
                                     "data.frame"))
  # and remove it in the source:
  for (ii in seq_along(outers)) {
    outers[[ii]]$performance <- NULL
  }


  # Production (Final) Model ------------------------------------------------
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
    family = family,
    type.measure = type.measure,
    opt.lambda.type = opt.lambda.type,
    opt.lambda.type.value = opt.lambda.type.value,
    opt.ystratify = opt.ystratify,
    opt.uniquefolds = opt.uniquefolds
  )

  final_model <- glmnet(
    x = xs,
    y = y,
    family = family,
    alpha = final_tuning$inner_best$alpha,
    lambda = lambdas[[which(alphalist == final_tuning$inner_best$alpha)]])

  final_performance <- tidy_predict.glmnet(
    mod = final_model,
    newx = xs,
    family = family,
    newy = y,
    binomial_thresh = cutoff,
    offset = offset,
    label = "Final",
    s = final_tuning$inner_best$lambda)

  final_performance <- structure(final_performance,
                                 class = c("classperformance",
                                           "data.frame"))

  final <- list(tuning = final_tuning,
                performance = final_performance,
                model = final_model)

  time_stop <- Sys.time()
  run_time <- difftime(time_stop, time_start, units = "hours")


  # Return object -----------------------------------------------------------
  obj <- structure(list(tuning = outers,
                        performance = performance,
                        folds = outfolds,
                        final = final,
                        input = list(callenv = callenv,
                                     runtime = run_time,
                                     lambdas = lambdas)),
                   class = c("dCVnet", "list"))

  cat("\n\n")
  print(thecall)
  cat("\n")
  cat(paste0("Finished.\n"))
  cat(paste0(time_stop, "\n"))
  cat(paste0("Runtime: ", sprintf("%.2f", run_time), " hours\n"))
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
  R <- lapply(seq_along(object$tuning), function(ii) {
    tt <- object$tuning[[ii]]
    coefs <- as.matrix(predict(tt$model,
                               type = "coef",
                               s = tt$tuning$inner_best$lambda))
    RR <- data.frame(Predictor = rownames(coefs),
                     Coef = c(coefs),
                     fold = names(object$folds)[[ii]],
                     stringsAsFactors = FALSE)
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


#' @export
print.dCVnet <- function(x, ...) {

  callenv <- x$input$callenv

  parsed <- parse_dCVnet_input(f = callenv$f,
                               y = callenv$y,
                               data = callenv$data,
                               family = callenv$family,
                               positive = callenv$positive)

  nalphas <- length(callenv$alphalist)
  nlambdas <- length(x$input$lambdas[[1]])

  typemeas <- attr(x$final$tuning$inner_results, "type.measure")

  outerfolds <- x$folds
  innerfolds <- x$tuning[[1]]$tuning$inner_folds

  outerkey <- data.frame(do.call(rbind,
                                 strsplit(names(outerfolds),
                                          split = ".", fixed = TRUE)))
  colnames(outerkey) <- c("Folds", "Reps")

  outer_k <- length(unique(outerkey$Folds))
  outer_nrep <- length(unique(outerkey$Reps))

  inner_k <- max(innerfolds[[1]])
  inner_nrep <- length(innerfolds)

  # What options were specified?:
  opts.empirical <- callenv[["opt.empirical_cutoff"]]
  opts.ystratify <- callenv[["opt.ystratify"]]
  opts.uniquefolds <- callenv[["opt.uniquefolds"]]
  opts.lambdatype <- callenv[["opt.lambda.type"]]
  if ( length(opts.lambdatype) > 1 ) opts.lambdatype <- opts.lambdatype[1]
  opts.lambdatypevalue <- callenv[["opt.lambda.type.value"]]

  cat("A dCVnet object, from the dCVnet Package\n\n")
  print(callenv[[1]]) # the call.
  cat("\n")
  cat(paste0("Runtime: ", sprintf("%.2f", x$input$runtime), " hours\n\n"))

  cat("Input\n-----\n")
  cat("Features:\n")
  cat(paste0("\t", ncol(parsed$x_mat), "\n"))

  cat("Observations:\n")
  cat(paste0("\t", nrow(parsed$x_mat), " subjects\n"))
  cat("Outcome:\n")
  if ( callenv$family %in% c("binomial", "multinomial")) {
    stab <- table(parsed$y)
    sapply(seq_along(stab), function(i) {
      cat(paste0("\t", stab[i], " of outcome: ", names(stab)[i], "\n"))
    })
  } else {
    print(summary(parsed$y)); cat("\n")
  }

  cat("Hyperparameter Tuning:\n")
  cat(paste("\tOptimising: ", typemeas, "\n"))
  cat(paste0("\t", nalphas, " alpha values:\t"))
  cat(callenv$alphalist)
  cat(paste0("\n\t", nlambdas, " lambda values\n"))

  cat("Cross-validation:\n")
  cat(paste("\tInner:\n"))
  cat(paste("\t\tk =\t", inner_k, "\n"))
  cat(paste("\t\tnrep =\t", inner_nrep, "\n"))
  cat(paste("\tOuter:\n"))
  cat(paste("\t\tk =\t", outer_k, "\n"))
  cat(paste("\t\tnrep =\t", outer_nrep, "\n"))

  cat("Options:\n")
  cat(paste0("\tOptimal Lambda selection (and parameter): \"",
            opts.lambdatype, "\" (", opts.lambdatypevalue, ")\n"
            ))
  cat(paste("\tUse Empirical Thresholding for LogReg Classification: ",
            opts.empirical, "\n"))
  cat(paste("\tStratify k-fold sampling by outcome: ",
            opts.ystratify, "\n"))
  cat(paste("\tCheck random folds are unique: ",
            opts.uniquefolds, "\n\n"))

  invisible(x)
}


#' selected_hyperparameters
#'
#' What hyperparameters were chosen for the outer loop of a dCVnet object.
#'
#' @param object a dCVnet object
#' @param what desired returns:
#'     \itemize{
#'     \item{\code{"summary"} - alpha/lambda summary over all outer folds/reps.}
#'     \item{\code{"data"} - return underlying data.}
#'     \item{\code{"both"} - both of the above in a list.}
#'     }
#' @param ... " "
#'
#' @name selected_hyperparameters
#'
#' @export
selected_hyperparameters <- function(object,
                                     what = c("summary", "data", "both")) {
  what <- match.arg(what)

  # what were the final hyperparams:
  FF <- as.data.frame(object$final$tuning$inner_best)
  FF.summary <- FF[, c("alpha", "lambda")]

  # What do the 'best-fitting' results of the inner loops look like:
  R <- lapply(object$tuning, function(x) {
    summary(x$tuning, print = FALSE)
  })
  R <- data.frame(do.call(rbind, R))
  R <- R[R$best, names(R)[!names(R) %in% "best"]]
  rownames(R) <- names(object$tuning)

  R$Rep <- sapply(strsplit(rownames(R), split = ".", fixed = TRUE), "[", 2)

  if ( what == "data") return(list(CVfolds = R, FinalModel = FF))

  alphas <- setNames(sort(unique(R$alpha)),
                     paste0("Alpha:", sort(unique(R$alpha))))

  # alpha, lambda and joint lambda|alpha descriptives:
  A <- table(R$alpha)
  L <- summary(R$lambda)
  J <- lapply(alphas, function(a) {
    data.frame(`n times selected` = sum(R$alpha == a),
               `lambda mean` = mean(R$lambda[R$alpha == a]),
               `lambda sd` = sd(R$lambda[R$alpha == a]),
               `lambda min` = min(R$lambda[R$alpha == a]),
               `lambda max` = max(R$lambda[R$alpha == a]))
  } )
  J <- do.call(rbind, J)

  if ( what == "summary") {
    return(list(alpha = A, lambda = L, joint = J, FinalModel = FF.summary))
  } else {
    return(list(data = R,
                summary = list(alpha = A, lambda = L, joint = J,
                               FinalModel = FF.summary)))
  }
}


#' coefficients_summary
#'
#' Describe the coefficients in the outer loop/final production model
#'     of a dCVnet object.
#'
#' @param object a dCVnet object
#' @param ... " "
#'
#' @name coefficients_summary
#'
#' @export
coefficients_summary <- function(object, ...) {

  Medians <- setNames(coef(object, type = "median"),
                      c("Predictor", "Median Coef"))

  Range <- coef(object, type = "all")
  Range.preds <- setNames(unique(Range$Predictor), unique(Range$Predictor))
  Range <- lapply(Range.preds, function(i) {
    kk <- Range$Coef[Range$Predictor == i]
    return(data.frame(min = min(kk),
                      max = max(kk),
                      propnz = sum(kk != 0) / length(kk)))
  } )
  names(Range) <- names(Range.preds)
  Range <- do.call(rbind, Range)

  FinalModel <- coef(object$final$model,
                     s = object$final$tuning$inner_best$lambda)
  FinalModel <- setNames(as.data.frame(as.matrix(FinalModel)), "FinalModel")

  return(data.frame(FinalModel,
                    OuterMedian = Medians[, 2],
                    Range))
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
    min_outcv <- outcv[outcv$Measure %in% min_vars, "Rep1", drop = FALSE]
    colnames(min_outcv) <- "Value"
  } else {
    min_outcv <- outcv[outcv$Measure %in% min_vars,
                       c("mean", "sd", "min", "max")]
  }

  min_outcv <- round(min_outcv, 3)
  row.names(min_outcv) <- min_vars

  # What do the 'best-fitting' results of the inner loops look like:
  R <- lapply(object$tuning, function(x) {
    summary.multialpha.repeated.cv.glmnet(x$tuning, print = FALSE)
  })

  R <- data.frame(do.call(rbind, R))
  R <- R[R$best, names(R)[!names(R) %in% "best"]]
  rownames(R) <- names(object$tuning)

  R$Rep <- sapply(strsplit(rownames(R), split = ".", fixed = TRUE), "[", 2)

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
  fmp_hp <- summary(object$final$tuning, print = FALSE)
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
  print(fmp_hp_str, quote = FALSE)

  invisible(R)
}

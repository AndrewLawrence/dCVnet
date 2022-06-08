
# Outer Loop functions ----------------------------------------------------

#' Fit a doubly cross-validated elastic-net regularised logistic regression
#'
#' Double (nested) repeated k-fold cross-validation used to:
#'     \itemize{
#'     \item{Produce unbiased estimates of out-of-sample
#'         classification performance (outer CV).}
#'     \item{Select optimal hyperparameters for the elasticnet (inner CV).}
#'     }
#'     Elasticnet hyperparameters are
#'     \bold{lambda} (the total regularisation penalty)
#'     and \bold{alpha} (the balance of L1 and L2 regularisation types).
#'
#' @inheritParams parse_dCVnet_input
#' @inheritSection parse_dCVnet_input Factor Outcomes
#' @inheritSection parse_dCVnet_input Notes
#' @inheritParams multialpha.repeated.cv.glmnet
#'
#' @param k_inner an integer, the k in the inner k-fold CV.
#' @param k_outer an integer, the k in the outer k-fold CV.
#' @param nrep_inner an integer, the number of repetitions (k-fold inner CV)
#' @param nrep_outer an integer, the number of repetitions (k-fold outer CV)
#' @param alphalist
#'     a numeric vector of values in (0,1]. This sets the search space for
#'     optimising hyperparameter alpha.
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
#'     Otherwise classify at 50% probability.
#' @param ... Arguments to pass through to cv.glmnet
#'     (may break things).
#' @return a dCVnet object.
#' @examples
#' \dontrun{
#'
#' # Iris example: Setosa vs. Virginica
#' #
#' # This example is fast to run, but not very informative because it is a
#' #  simple problem without overfitting and the predictors work 'perfectly'.
#' # `help(iris)` for more infomation on the data.
#'
#' # Make a two-class problem from the iris dataset:
#' siris <- droplevels(subset(iris, iris$Species != "versicolor"))
#' # scale the iris predictors:
#' siris[,1:4] <- scale(siris[,1:4])
#'
#' set.seed(1) # for reproducibility
#' model <- dCVnet(y = siris$Species,
#'                      f = ~ Sepal.Length + Sepal.Width +
#'                            Petal.Length + Petal.Width,
#'                      data = siris,
#'                      alphalist = c(0.2, 0.5, 1.0),
#'                      opt.lambda.type = "1se")
#'
#' #Note: in most circumstances non-default (larger) values of
#' #      nrep_inner and nrep_outer will be required.
#'
#' # Input summary:
#' dCVnet::parseddata_summary(model)
#'
#' # Model summary:
#' summary(model)
#'
#' # Detailed cross-validated model performance summary:
#' summary(performance(model))
#'
#' # hyperparameter tuning plot:
#' plot(model)
#' # as above, but zoomed in:
#' plot(model)$plot + ggplot2::coord_cartesian(ylim = c(0,0.03), xlim=c(-4, -2))
#'
#' # Performance ROC plot:
#' plot(model, type = "roc")
#'
#' # predictor importance (better with more outer reps):
#' dCVnet::coefficients_summary(model)
#' #    show variability over both folds and reps:
#' dCVnet::plot_outerloop_coefs(model, "all")
#'
#' # selected hyperparameters:
#' dCVnet::selected_hyperparameters(model, what = "data")
#'
#' # Reference logistic regressions (unregularised & univariate):
#' ref_model <- dCVnet::reflogreg(model)
#'
#' dCVnet::report_reference_performance_summary(ref_model)
#'
#'
#' }
#' @importFrom stats aggregate as.formula coef glm model.frame model.matrix
#' @importFrom stats predict sd terms var
#' @importFrom data.table rbindlist
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
  offset = NULL,

  opt.lambda.type = c("min", "1se"),
  opt.empirical_cutoff = FALSE,
  opt.uniquefolds = FALSE,
  opt.ystratify = TRUE,

  ...) {

  # Parse input -------------------------------------------------------------

  # call for printing:
  thecall <- match.call()

  # call for posterity (save objects from environment):
  callenv <- c(as.list(environment()), list(...))

  time_start <- force(Sys.time()) # for logging.
  parsed <- parse_dCVnet_input(f = f,
                               y = y,
                               data = data,
                               family = family,
                               passNA = FALSE)
  x <- parsed$x_mat
  y <- parsed$y

  # stop if AUC requested and magic number not met.
  if ( type.measure == "auc" ) {
    # magic number for glmnet (in the inner loop) is:
    #   Ncases * outer train proportion * inner test proportion
    auc_magic <- (nrow(x) * (1 - (1 / k_outer)) * (1 / k_inner))
    if (auc_magic < 11) {
      stop(
        paste("AUC is not possible due to small sample size!\n
              Estimated inner foldsize =", auc_magic,
              "\nuse different type.measure (e.g. 'class')."))
    }
  }

  # Check the alpha values:
  alphalist <- parse_alphalist(alphalist)
  # how many alphas did we feed in:
  nalpha <- length(alphalist)

  # are we reporting performance at 0.5 or the empirical cutoff?
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
  # Note: this is by default stratified by y, we obtain unstratified sampling
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

  # Initialise the call -----------------------------------------------------

  # create a call to multialpha.repeated.cv.glmnet which will be
  #   updated with new data and re-evaluated within each run of the outer loop.
  cl.marcvglm.argset <- unique(
    c(methods::formalArgs(multialpha.repeated.cv.glmnet),
      methods::formalArgs(repeated.cv.glmnet),
      methods::formalArgs(glmnet::glmnet),
      methods::formalArgs(glmnet::cv.glmnet))
  )
  cl.marcvglm <- callenv[names(callenv) %in% cl.marcvglm.argset]
  # force no standardisation:
  cl.marcvglm$standardize <- FALSE
  # force keep multialpha models:
  cl.marcvglm$opt.keep_models <- "best"
  # translate some dCVnet args to multialpha.repeated.cv.glmnet:
  cl.marcvglm$k <- k_inner
  cl.marcvglm$nrep <- nrep_inner

  # add the list of alpha values:
  cl.marcvglm$alphalist <- alphalist


  # Create Lambda paths -----------------------------------------------------

  # Lambda path notes:
  # Inspired by the source of cv.glmnet, we get the list of lambdas once,
  #   rather than per-CV fold. This is done by running a glmnet at each alpha
  #   for all data and extracting the list of lambdas used.
  # Previously a wider/more representative range of lambdas was obtained by
  #   by applying a formula for calculation of max lambda to random samplings
  #   of the data (per alpha).
  #
  # However, 1) The formula involved matrix multiplication and
  #             The multiple random samplings was time consuming.
  #          2) There is no published formula for these calculations for
  #             mgaussian/cox family models, so the previous method could not
  #             be extended.
  #          3) Using a call to fortran ensures the lambda sequence is formatted
  #             in the way glmnet expects annd avoids inconsistency in rounding
  #             and excess precision.
  #          4) Hastie et al use the whole-data fortran-call method in cv.glmnet
  #
  # Note: creating lambdas before standardising the data is OK because
  #   the internal glmnet code which uses internal standardisation and
  #   standardisation is idempotent. This would need to be revisited if
  #   other transformations are implemented - e.g. yeo-johnsson.

  # subset the above call such that it containts all valid arguments for glmnet:
  cl.lf <- cl.marcvglm[
    names(cl.marcvglm) %in% methods::formalArgs(glmnet::glmnet)
  ]
  # modify to force glmnet internal standardisation:
  cl.lf$standardize <- TRUE
  # add x and y
  cl.lf$x <- x
  cl.lf$y <- y

  # Loop over alpha values and extract lambda list for each alpha:
  lambdas <- lapply(alphalist,
                    function(aa) {
                      cl.lf$alpha <- aa
                      m <- do.call(glmnet::glmnet, cl.lf)
                      return(m$lambda)
                    })


  # add to the main call:
  cl.marcvglm$lambdas <- lambdas


  # Run Outer Loop (Start) --------------------------------------------------

  # Notes: x and y are added to the call within the outer loop.
  #        inner folds are generated by multialpha.repeated.cv.glmnet

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

      # scaling to mean=0, sd=1
      #   (scaling calculated on completecases train data, applied to test data)
      PPx <- caret::preProcess(trainx, method = c("center", "scale"))
      trainx <- predict(PPx, trainx)
      testx <- predict(PPx, testx)

      # inner loop train data -------------------------------------------
      #   - receives no (outer) test data

      cl.marcvglm$x <- trainx
      cl.marcvglm$y <- trainy

      inners <- do.call("multialpha.repeated.cv.glmnet",
                        cl.marcvglm)

      # hyperparameter selection --------------------------------------------
      #   - based on best CV/out-of-sample performance in inner loop
      fit_lambda <- inners$best$lambda


      # tuned outer-train model ---------------------------------------------
      #   - tuned wrt best inner alpha
      #   - full lambda path is fit with subsequent selection
      #       (see https://web.stanford.edu/~hastie/glmnet/glmnet_beta.html:
      #       "[single lambda fit] is not recommended and is not the spirit
      #        of the package."
      newx_performance <- tidy_predict.glmnet(
        mod = inners$models[[inners$bestmodel]],
        newx = testx,
        family = family,
        newy = testy,
        binomial_thresh = cutoff,
        offset = offset,
        label = strsplit(names(outfolds)[[i]],
                         split = ".", fixed = TRUE)[[1]][2],
        s = fit_lambda)

      return(list(tuning = drop_models.multialpha.repeated.cv.glmnet(inners),
                  model = inners$models[[inners$bestmodel]],
                  performance = structure(newx_performance,
                                          family = family,
                                          class = c("performance",
                                                    "data.frame"))))
    }
  )
  names(outers) <- names(outfolds)
  # Run Outer Loop (End) ----------------------------------------------------


  # Outer Loop Performance --------------------------------------------------

  # Gather performance from the outers into a single dataframe.
  performance <- as.data.frame(
    data.table::rbindlist(
      lapply(outers, "[[", "performance")
    ),
    stringsAsFactors = FALSE
  )
  rownames(performance) <- NULL
  performance <- structure(performance,
                           family = family,
                           class = c("performance",
                                     "data.frame"))
  # and remove it in the source:
  for (ii in seq_along(outers)) {
    outers[[ii]]$performance <- NULL
  }


  # Production (Final) Model ------------------------------------------------
  prod_PPx <- caret::preProcess(x, method = c("center", "scale"))
  xs <- predict(prod_PPx, x)

  # Run the inner loop for the production model:
  cat("\n\nProduction Model\n")

  cl.marcvglm$y <- y
  cl.marcvglm$x <- xs
  prod_tuning <- do.call("multialpha.repeated.cv.glmnet", cl.marcvglm)

  prod_performance <- tidy_predict.glmnet(
    mod = prod_tuning$models[[prod_tuning$bestmodel]],
    newx = xs,
    family = family,
    newy = y,
    binomial_thresh = cutoff,
    offset = offset,
    label = "Production",
    s = prod_tuning$best$lambda)

  prod_performance <- structure(prod_performance,
                                family = family,
                                 class = c("performance",
                                           "data.frame"))

  prod <- list(
    tuning = drop_models.multialpha.repeated.cv.glmnet(prod_tuning),
    performance = prod_performance,
    model = prod_tuning,
    preprocess = prod_PPx) # include the preprocessing for predict method.

  time_stop <- Sys.time()
  run_time <- difftime(time_stop, time_start, units = "hours")

  # Return object -----------------------------------------------------------
  obj <- structure(list(tuning = outers,
                        performance = performance,
                        folds = outfolds,
                        prod = prod,
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


# Extract coefficients from a dCVnet object:
#' coef.dCVnet
#'
#' Coefficients from a dCVnet object. If type is "production" this gives the
#' coefficients from the production model fit to all data. Otherwise
#' coefficients are returned from the outerloop of the cross-validation.
#'
#' @param object a dCVnet object
#' @param type how to return coefficients.
#'     \itemize{
#'     \item{\code{"production"} - the "production" model coefficients}
#'     \item{\code{"all"} - return separate coefficients
#'         for each fold & rep of CV}
#'     \item{\code{"mean"} - return mean over \code{"all"}.}
#'     \item{\code{"median"} - return median over \code{"all"}.}
#'     \item{\code{"byrep"} - return separate coefficients for each rep of CV
#'         (take mean average over folds).}
#'     \item{\code{byrep_mean} - mean over \code{"byrep"}}
#'     \item{\code{byrep_mean} - median over \code{"byrep"}}
#'     }
#' @param ... " "
#'
#' @return a data.frame of coefficient values (see type argument) containing
#'    columns: Predictor and Coef (as well as Rep and fold if applicable)
#'
#' @name coef.dCVnet
#'
#' @seealso coefficients_summary
#'
#' @export
coef.dCVnet <- function(object, type = "all", ...) {
  # Type can be:
  #   production - a.k.a. final model coefficients
  #   all - coefficients for each rep/fold.
  #   rep - mean average per rep.
  #   mean - mean of per-rep means.
  #   median = median of per-rep means.

  type <- match.arg(type, choices = c("all",
                                      "production",
                                      "mean",
                                      "median",
                                      "byrep",
                                      "byrep_mean",
                                      "byrep_median"))

  # first handle production:
  if ( type == "production" ) {
    return(tidy_coef.multialpha.repeated.cv.glmnet(object$prod$model))
  }

  # next handle outerloop coefficients:

  # works on the output of the outerloop.
  #   Given the 'best' alpha/lambda from the inner loop,
  #     what were the coefficients in a model selected with the best alpha.
  R <- lapply(seq_along(object$tuning), function(ii) {
    tt <- object$tuning[[ii]]
    coefs <- tidy_coef.glmnet(
      tt$model,
      s = tt$tuning$best$lambda
    )
    coefs$fold <- names(object$folds)[[ii]]
    return(coefs)
  })
  R <- as.data.frame(data.table::rbindlist(R), stringsAsFactors = FALSE)

  R$Predictor <- factor(R$Predictor, levels = unique(as.character(R$Predictor)))
  R$fold <- factor(R$fold, levels = unique(as.character(R$fold)))

  R$Rep <- vapply(X = strsplit(as.character(R$fold), split = "\\."),
                  FUN = "[",
                  FUN.VALUE = c(""),
                  2)

  if (type == "all") return(R)
  if (type == "mean") {
    R <- aggregate(R$Coef, by = list(Predictor = R$Predictor), mean)
    colnames(R)[2] <- "Coef"
    return(R)
  }
  if (type == "median") {
    R <- aggregate(R$Coef, by = list(Predictor = R$Predictor), stats::median)
    colnames(R)[2] <- "Coef"
    return(R)
  }
  # Remaining cases are byrep:

  # Aggregate over Repetitions:
  R <- aggregate(R$Coef,
                 by = list(Predictor = R$Predictor,
                           Rep = R$Rep),
                 mean)
  names(R)[3] <- "Coef"
  R <- R[, c(1, 3, 2)] # reorder columns
  return(switch(type,
                byrep = R,
                byrep_mean = stats::setNames(
                  aggregate(R$Coef,
                            by = list(Predictor = R$Predictor),
                            mean),
                  c("Predictor", "Coef")),
                byrep_median = stats::setNames(
                  aggregate(R$Coef,
                            by = list(Predictor = R$Predictor),
                            stats::median),
                  c("Predictor", "Coef")),
                stop(paste("type:",
                           type,
                           " - is non-sensical"))
  ))
}

#' @export
print.dCVnet <- function(x, ...) {

  callenv <- x$input$callenv

  parsed <- parse_dCVnet_input(f = callenv$f,
                               y = callenv$y,
                               data = callenv$data,
                               family = callenv$family)

  nalphas <- length(callenv$alphalist)
  nlambdas <- range(vapply(x$input$lambdas, length, c(1)))

  outerfolds <- x$folds
  innerfolds <- x$tuning[[1]]$tuning$folds

  outerkey <- data.frame(do.call(rbind,
                                 strsplit(names(outerfolds),
                                          split = ".", fixed = TRUE)),
                         stringsAsFactors = FALSE)
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
  print(describe_outcome(parsed$y, family = family(x)))
  cat("\n")

  cat("Hyperparameter Tuning:\n")
  cat(paste("\tOptimising: ", callenv$type.measure, "\n"))
  cat(paste0("\t", nalphas, " alpha values:\t"))
  cat(callenv$alphalist)
  cat(paste0("\n\trequested ", callenv$nlambda, " lambdas, ",
             " fitted between ", nlambdas[1],
             " and ", nlambdas[2],
             " (range over alphas)\n"))

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
                                     what = c("summary", "data", "both"),
                                     ...) {
  what <- match.arg(what)

  # what were the production hyperparams:
  FF <- as.data.frame(object$prod$tuning$best, stringsAsFactors = FALSE)
  FF.summary <- FF[, c("alpha", "lambda")]

  # What do the 'best-fitting' results of the inner loops look like:
  R <- lapply(object$tuning, function(x) {
    summary(x$tuning, print = FALSE)
  })
  R <- data.frame(do.call(rbind, R), stringsAsFactors = FALSE)
  R <- R[R$best, names(R)[!names(R) %in% "best"]]
  rownames(R) <- names(object$tuning)

  R$Rep <- vapply(X = strsplit(rownames(R), split = ".", fixed = TRUE),
                  FUN = "[",
                  FUN.VALUE = c(""),
                  2)

  if ( what == "data") return(list(CVfolds = R, ProductionModel = FF))

  alphas <- setNames(sort(unique(R$alpha)),
                     paste0("Alpha:", sort(unique(R$alpha))))

  # alpha, lambda and joint lambda|alpha descriptives:
  A <- table(R$alpha)
  L <- summary(R$lambda)
  J <- lapply(alphas, function(a) {
    data.frame(alpha = a,
               `n times selected` = sum(R$alpha == a),
               `lambda mean` = mean(R$lambda[R$alpha == a]),
               `lambda sd` = sd(R$lambda[R$alpha == a]),
               `lambda min` = min(R$lambda[R$alpha == a]),
               `lambda max` = max(R$lambda[R$alpha == a]),
               stringsAsFactors = FALSE)
  } )
  J <- as.data.frame(data.table::rbindlist(J), stringsAsFactors = FALSE)

  if ( what == "summary" ) {
    return(list(alpha = A, lambda = L, joint = J, ProductionModel = FF.summary))
  } else {
    return(list(data = list(CVfolds = R,
                            ProductionModel = FF),
                summary = list(alpha = A, lambda = L, joint = J,
                               ProductionModel = FF.summary)))
  }
}


#' coefficients_summary
#'
#' Extract the coefficients in the production model (fit to complete data)
#'     and provide descriptives for the coefficients in the outerloop
#'     cross-validation
#'
#' @param object a dCVnet object
#' @param ... " "
#'
#' @name coefficients_summary
#'
#' @return A data.frame with a row for each term/predictor in the model.
#'     Columns present coefficients of the Production Model
#'     (\code{ProductionModel}) and descriptives of coefficients from the
#'     outer-loop cross-validation. \code{OuterMedian}: the median over folds and
#'     repetitions of cross-validation. min, max = the smallest, largest values
#'     held by that coefficient over the cross-validation. propnz = the
#'     proportion of cross-validation folds for which the variable was non-zero
#'     (Note: this will be close to 1 for low alpha = ridge penalty)
#'
#' @seealso coef.dCVnet
#'
#' @export
coefficients_summary <- function(object, ...) {

  Medians <- setNames(coef(object, type = "median"),
                      c("Predictor", "Median Coef"))

  Range <- coef(object, type = "all")
  Range.preds <- setNames(unique(Range$Predictor),
                          unique(Range$Predictor))
  Range <- lapply(Range.preds, function(i) {
    kk <- Range$Coef[Range$Predictor == i]
    return(data.frame(min = min(kk),
                      max = max(kk),
                      propnz = sum(kk != 0) / length(kk),
                      stringsAsFactors = FALSE))
  } )
  names(Range) <- names(Range.preds)
  Range <- as.data.frame(data.table::rbindlist(Range),
                         stringsAsFactors = FALSE)

  ProductionModel <- coef(object, type = "production")
  colnames(ProductionModel)[2] <- "ProductionModel"

  return(data.frame(ProductionModel,
                    OuterMedian = Medians[, 2],
                    Range,
                    stringsAsFactors = FALSE))
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
  min_outcv <- report_performance_summary(object, short = TRUE)

  # What do the 'best-fitting' results of the inner loops look like:
  R <- lapply(object$tuning, function(x) {
    summary.multialpha.repeated.cv.glmnet(x$tuning, print = FALSE)
  })

  R <- data.frame(do.call(rbind, R), stringsAsFactors = FALSE)
  R <- R[R$best, names(R)[!names(R) %in% "best"]]
  rownames(R) <- names(object$tuning)

  R$Rep <- vapply(X = strsplit(rownames(R), split = ".", fixed = TRUE),
                  FUN = "[",
                  FUN.VALUE = c(""),
                  2)

  # summarise lambdas:
  lambda_summary <- c(mean = mean(R$lambda),
                      sd = sd(R$lambda),
                      min = min(R$lambda),
                      max = max(R$lambda))

  # Summarise inner loop performances (i.e. best cvm)  by rep:
  cvm_repmean <- aggregate(R$cvm, by = list(Rep = R$Rep), FUN = mean)

  # summarise 'production' model performance
  min_pmp <- summary(object$prod$performance, short = TRUE, label = "None")

  # final model hyper parameters:
  pmp_hp <- summary(object$prod$tuning, print = FALSE)
  pmp_hp <- pmp_hp[pmp_hp$best, ]

  pmp_hp_str <- c(alpha = pmp_hp$alpha,
                  lambda = formatC(pmp_hp$lambda),
                  cvm = formatC(pmp_hp$cvm))

  # Write this out:
  cat("\n")
  .titlecat("Outer Loop CV Performance")
  print(min_outcv, digits = 3)
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
             attr(object$prod$tuning$best, "type.measure"), "\n\n"))
  cat("cvm summary (per-rep, averaged over folds):\n")
  print(summary(cvm_repmean$x))
  cat("cvm summary (over folds):\n")
  print(summary(R$cvm))
  cat("\n")

  .titlecat("'Production' Model")
  cat("Production Performance (not cross-validated):\n")
  print(min_pmp, digits = 3)
  cat("Production Hyperparameter Tuning:\n")
  print(pmp_hp_str, quote = FALSE)

  invisible(R)
}

#' predict.dCVnet
#'
#' predict method for a \code{\link{dCVnet}} object.
#'     predictions come from the "production" model fitted to the full data.
#'     As a result they do not reflect cross-validated/out-of-sample
#'     performance.
#'
#' @param object a a \code{\link{dCVnet}} object.
#' @param ... passed to \code{\link{predict.multialpha.repeated.cv.glmnet}},
#'     then \code{\link[glmnet]{predict.glmnet}}
#' @inheritParams predict.multialpha.repeated.cv.glmnet
#' @inheritParams glmnet::predict.glmnet
#'
#' @return predictions of the requested type.
#'
#' @importFrom utils getFromNamespace
#'
#' @export
predict.dCVnet <- function(object,
                           newx,
                           ...) {
  # apply the preprocessing from the production model:

  preProcessPredict <- utils::getFromNamespace(x = "predict.preProcess",
                                               ns = "caret")
  newx <- as.matrix(newx) # coerce to matrix (this fixed a random bug...)
  newx <- preProcessPredict(object$prod$preprocessing, newx)
  # run the prediction:
  predict.multialpha.repeated.cv.glmnet(object$prod$model,
                                        newx = newx,
                                        ...)
}

#' family.dCVnet
#'
#' Extracts the model family from a dCVnet object as a string
#'
#' @param object a a \code{\link{dCVnet}} object.
#' @param ... further arguments passed to methods
#'
#' @return character vector (of length 1)
#'
#' @export
family.dCVnet <- function(object, ...) {
  object$input$callenv$family
}

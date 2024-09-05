
# dCVnet function ----------------------------------------------------

#' Double cross-validated elastic-net regularised regression
#'
#' Fits and cross-validates an elastic-net regularised regression model
#' using independent cross-validation to select optimal alpha and lambda
#' hyperparameters for the regularisation.
#'
#' The time-consuming double (i.e. nested) cross-validation (CV) is used because
#' single cross-validation - which both tunes hyperparameters and estimates
#' out-of-sample classification performance - will be optimistically biased.
#'
#' Cross-validation for both the inner and outer loop is *repeated k-fold*.
#'
#' Both alpha and lambda hyperparameters of the elastic-net can be tuned:
#' \itemize{
#'     \item{\bold{lambda} - the total regularisation penalty}
#'     \item{\bold{alpha} - the balance of L1(LASSO) and L2 (Ridge)
#'         regularisation types}
#' }
#'
#' @section Coefficients:
#'
#' Currently all coefficients reported / used by dCVnet are
#' *semi-standardised* for x, but not y. In other words, the predictor matrix x
#' is mean-centered and scaled by the standard deviation prior to calculations.
#'
#' The predict method for dCVnet stores these means/SDs, and will apply the same
#' standardisation to x on predicting new values.
#'
#' As a result the reported coefficients for dCVnet can be interpreted as
#' (semi-)standardised effect sizes. A coefficient of 0.5 is the effect of a 1SD
#' difference in that x matrix element. Note this is the same for all elements
#' of x, even factors, and so if coefficients for binary variables are
#' interpreted as an effect size this will include the impact of prevalence.
#'
#' When running cross-validation, standardisation is based on the means
#' and standard deviations of the training dataset, not the held-out test data.
#' This prevents leakage from train to test data.
#'
#' This approach can be contrasted with glmnet which internally standardises,
#' and then back-transforms, coefficients to the original scale. In contrast
#' dCVnet model coefficients are *always* presented as (semi-)standardised.
#'
#' Coefficients in the original scale can be recovered using the
#' standard deviations and means employed in the standardisation
#' (see: \url{https://stats.stackexchange.com/a/75025})
#' These means and standard deviations are retained for the production model
#' in the `preprocess` slots of the dCVnet object: `my_model$prod$preprocess`.
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
#' @param opt.use_imputation Boolean.
#'     Run imputation on missing predictors?
#' @param opt.imputation_method Which imputation method?: \itemize{
#'     \item \code{"mean"} - mean imputation (unconditional)
#'     \item \code{"knn"} - k-nearest neighbours imputation
#'         (uses \code{\link[caret]{preProcess}}).
#'     \item \code{"missForestPredict"} - use the
#'     missForestPredict package to impute missing values.
#'     }
#' @param opt.imputation_usey Boolean.
#'     Should conditional imputation methods use y in the imputation model?
#'     Note: no effect if \code{opt.use_imputation} is \code{FALSE}, or if
#'     \code{opt.imputation_method} is \code{"mean"}.
#' @param ... Arguments to pass through to cv.glmnet
#'     (may break things).
#' @return a dCVnet object containing:
#'  \itemize{
#'  \item{`input`: call arguments and input data}
#'  \item{`prod`: the production model and preprocessing information
#'      used in making new predictions.}
#'  \item{`folds`: outer-loop CV fold membership}
#'  \item{`performance`: Cross-validated \code{\link{performance}}}
#'  \item{`tuning`: Outer loop CV tuning information}
#'  }
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
#' ref_model <- dCVnet::refunreg(model)
#'
#' dCVnet::report_reference_performance_summary(ref_model)
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
  opt.random_seed = NULL,

  opt.use_imputation = FALSE,
  opt.imputation_usey = FALSE,
  opt.imputation_method = c("mean", "knn", "missForestPredict"),

  ...
) {

  # ~ Parse input -------------------------------------------------------------

  # capture the function call (for printing):
  thecall <- match.call()

  # check imputation arguments:
  opt.imputation_method <- match.arg(opt.imputation_method)

  if (!opt.use_imputation) {
    opt.imputation_method <- "mean" # will not be applied.
    opt.imputation_usey <- FALSE
  }

  if ( opt.imputation_usey && opt.imputation_method == "mean" ) {
    stopmsg <- paste(
      "Invalid option combination.",
      "Set opt.imputation_usey to FALSE if using mean imputation",
      sep = "\n"
    )
    stop(stopmsg)
  }

  # If the imputation method is missForestPredict then check we have the package
  #   and fail if not:
  if (opt.imputation_method == "missForestPredict" &&
        ! requireNamespace("missForestPredict", quietly = TRUE) ) {
    stop("Please install the missForestPredict package to use this option.")
  }

  # argument matching for model family:
  family <- match.arg(family, choices = supported_model_families())

  # save arguments for posterity
  #   (i.e. package objects from the calling environment):
  callenv <- c(as.list(environment()), list(...))

  time_start <- force(Sys.time()) # for logging.

  # Parse basic model input:
  parsed <- parse_dCVnet_input(f = f,
                               y = y,
                               data = data,
                               family = family,
                               passNA = opt.use_imputation,
                               offset = offset)
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
              "\nuse different type.measure (e.g. 'class').")
      )
    }
  }

  # Check the alpha values:
  alphalist <- parse_alphalist(alphalist)
  nalpha <- length(alphalist)

  # for classification problems:
  #   are we reporting performance at 0.5 or the empirical cutoff?
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

  # ~ Preprocessing functions -------------------------------------------------

  # Data is preprocessed repeatedly during cross-validation,
  #   the preprocessing differs depending on imputation method
  #   use a utility function to get appropriate preprocessing functions
  #   depending on opt.imputation_method.

  pp_fn <- preproc_imp_functions(opt.imputation_method = opt.imputation_method)

  # ~ Create outer folds ------------------------------------------------------
  # Note: this is by default stratified by y, we obtain unstratified sampling
  #         by giving caret::createMultiFolds a single-level factor/char.
  ystrat <- y
  if ( identical(opt.ystratify, FALSE) || family %in% c("cox", "mgaussian") ) {
    ystrat <- rep("x", length(y))
  }

  set.seed(opt.random_seed)
  # But don't use the same fixed seed for all the inner loop calls:
  opt.random_seed <- NULL

  outfolds <- caret::createMultiFolds(y = ystrat,
                                      k = k_outer,
                                      times = nrep_outer)
  rm(ystrat)
  # Warn if repeated folds were obtained.
  if ( opt.uniquefolds ) checkForDuplicateCVFolds(outfolds)

  imax <- length(outfolds)
  names(outfolds) <- paste0("Out", names(outfolds)) # give names

  # ~ Initialise the call -----------------------------------------------------

  # create a call to multialpha.repeated.cv.glmnet which will be
  #   used later as the basis of the outer loop.
  #
  #   For each outer fold this call will be updated with new data and
  #   evaluated.
  cl.marcvglm.argset <- unique(
    c(methods::formalArgs(multialpha.repeated.cv.glmnet),
      methods::formalArgs(repeated.cv.glmnet),
      methods::formalArgs(glmnet::glmnet),
      methods::formalArgs(glmnet::cv.glmnet))
  )
  # Initialise the call with all relevant arguments provided to dCVnet:
  cl.marcvglm <- callenv[names(callenv) %in% cl.marcvglm.argset]

  # Then make some key adjustments
  # Ensure no standardisation in glmnet
  #   (because standardisation will be handled independently by preProcess):
  cl.marcvglm$standardize <- FALSE
  # Ensure we keep multialpha models:
  if ( is.null(cl.marcvglm$opt.keep_models) ||
        cl.marcvglm$opt.keep_models == "none" ) {
    cl.marcvglm$opt.keep_models <- "best"
  }
  # Translate some dCVnet args to multialpha.repeated.cv.glmnet format:
  cl.marcvglm$k <- k_inner
  cl.marcvglm$nrep <- nrep_inner

  # Add the processed list of alpha values:
  cl.marcvglm$alphalist <- alphalist


  # ~ Production (Final) Model ------------------------------------------------
  #   start at the end. This is the model we are cross-validating.
  #   note that no output from this model are used in the cross-validation.
  #   excepting - the lambda list.
  #
  #   first, create the preprocessing object:
  prod_PPx <- structure(
    pp_fn$fit(x),
    fit = pp_fn$fit,
    apply = pp_fn$apply
  )
  xs <- attr(prod_PPx, "apply")(prod_PPx, x)

  # ~ Create Lambda paths -----------------------------------------------------

  # Lambda paths:
  # dCVnet works with a single, fixed list of lambdas applied to all folds/reps
  #   of the CV, rather than independent lambda lists.
  #
  # Lambda lists are obtained by running a glmnet at each alpha and extracting
  #   the list of lambdas used.
  #
  # In an earlier version of glmnet, aiming to be more robust, a single list
  #   was obtained by directly calculation of max-lambda for random samplings of
  #   the data. This used a formula for max lambda.
  #
  # However, 1) The formula involved matrix multiplication and
  #             The multiple random samplings was time consuming.
  #          2) For mgaussian/cox family models there is no published formula
  #             for these calculations , so the previous method could not
  #             be generalised.
  #          3) Using a call to glmnet ensures the lambda sequence is formatted
  #             exactly as glmnet expects, avoiding rounding inconsistency and
  #             excess precision.
  #          4) Hastie et al use the whole-data fortran-call method in cv.glmnet
  #
  # Note: creating lambdas before standardising the data is OK because
  #   the glmnet uses internal standardisation and standardisation is idempotent
  #   This would need to be revisited if other transformations are implemented
  #   - e.g. yeo-johnsson.

  # Subset the above call to retain all valid arguments for glmnet:
  cl.lf <- cl.marcvglm[
    names(cl.marcvglm) %in% methods::formalArgs(glmnet::glmnet)
  ]
  # add x and y
  cl.lf$x <- xs
  cl.lf$y <- y

  # Loop over alpha values and extract lambda list for each alpha:
  lambdas <- lapply(alphalist,
                    function(aa) {
                      cl.lf$alpha <- aa
                      m <- do.call(glmnet::glmnet, cl.lf)
                      return(m$lambda)
                    })


  # add lambdas to the main call:
  cl.marcvglm$lambdas <- lambdas

  # ~ Production (Final) Model (cont.) ----------------------------------------

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
    s = prod_tuning$best$lambda
  )

  prod_performance <- structure(prod_performance,
                                family = family,
                                class = c("performance",
                                          "data.frame"))

  prod <- list(
    tuning = drop_models.multialpha.repeated.cv.glmnet(prod_tuning),
    performance = prod_performance,
    model = prod_tuning,
    preprocess = prod_PPx # include the preprocessing for predict method.
  )



  # ~ Run Outer Loop (Start) --------------------------------------------------

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

      # ~ split & preprocess data ---------------------------------------------

      of <- outfolds[[i]]

      sel_train <- seq_len(nrow(x)) %in% of
      sel_test <- !sel_train

      trainx <- subset(x, sel_train)
      trainy <- subset(y, sel_train)

      testx <- subset(x, sel_test)
      testy <- subset(y, sel_test)

      # scaling to mean=0, sd=1
      #   (scaling calculated on completecases train data, applied to test data)
      PPx <- structure(
        pp_fn$fit(trainx),
        fit = pp_fn$fit,
        apply = pp_fn$apply
      )
      trainx <- attr(PPx, "apply")(PPx, trainx)
      testx  <- attr(PPx, "apply")(PPx, testx)

      # ~ inner loop train data -------------------------------------------
      #   - receives no (outer) test data

      cl.marcvglm$x <- trainx
      cl.marcvglm$y <- trainy

      inners <- do.call("multialpha.repeated.cv.glmnet",
                        cl.marcvglm)

      # ~ hyperparameter selection --------------------------------------------
      #   - based on best CV/out-of-sample performance in inner loop
      fit_lambda <- inners$best$lambda


      # ~ tuned outer-train model ---------------------------------------------
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
        s = fit_lambda
      )

      return(list(tuning = drop_models.multialpha.repeated.cv.glmnet(inners),
                  model = inners$models[[inners$bestmodel]],
                  performance = structure(newx_performance,
                                          family = family,
                                          class = c("performance",
                                                    "data.frame"))))
    }
  )
  names(outers) <- names(outfolds)
  # ~ Run Outer Loop (End) ----------------------------------------------------


  # ~ Outer Loop Performance --------------------------------------------------

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


  time_stop <- Sys.time()
  run_time <- difftime(time_stop, time_start, units = "hours")

  # ~ Return object -----------------------------------------------------------
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


# dCVnet methods ----------------------------------------------------------


# Extract coefficients from a dCVnet object:
#' coef.dCVnet
#'
#' Coefficients from a dCVnet object. If type is "production" this gives the
#' coefficients from the production model fit to all data. Other options for
#' type return coefficients from the outerloop of the cross-validation which
#' are summarised/presented in different ways.
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

  # first: "production"
  if ( type == "production" ) {
    return(tidy_coef.multialpha.repeated.cv.glmnet(object$prod$model))
  }

  # next: outerloop coefficients

  # works on the output of the outerloop.
  #   Given the 'best' alpha/lambda from the inner loop,
  #     what were the coefficients in a model selected with the best alpha?
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
  return(switch(
    type,
    byrep = R,
    byrep_mean = stats::setNames(aggregate(R$Coef,
                                           by = list(Predictor = R$Predictor),
                                           mean),
                                 c("Predictor", "Coef")),
    byrep_median = stats::setNames(
      aggregate(R$Coef,
                by = list(Predictor = R$Predictor),
                stats::median),
      c("Predictor", "Coef")
    ),
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
                               family = callenv$family,
                               passNA = callenv$opt.use_imputation,
                               offset = callenv$offset)

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
#'     outer-loop cross-validation. \code{OuterMedian}: the median over folds
#'     and repetitions of cross-validation.
#'     \code{min}, \code{max}: the smallest, largest values
#'     held by that coefficient over the cross-validation. \code{propnz}: the
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
  preProcessPredict <- attr(object$prod$preprocess,
                            "apply")
  newx <- as.matrix(newx) # coerce to matrix (this fixed a random bug...)
  newx <- preProcessPredict(object$prod$preprocess, newx)
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

#' variable_importance
#'
#' Variable importance for dCVnet/glmnet models does not require permutation
#'     methods, because coefficients are directly interpretable.
#'
#'     This VI function follows caret's example (see \code{\link[caret]{varImp}}
#'     function) and simply returns the absolute values of the coefficients.
#'
#'     As variable importance is inferential this is done for the tuned
#'     "production" model rather than the cross-validated outer-loop.
#'
#' @param x a dCVnet object
#' @param scale Boolean. Should the return values be scaled so the most
#' important value is 1?
#' @param percentage Boolean. Should the return values be scaled so the most
#'     important value is 100?
#' @return a data.frame of variable names "Predictor"
#'     and variable importance "varImp".
#' @seealso \code{\link[caret]{varImp}}
#' @export
variable_importance <- function(x,
                                scale = FALSE,
                                percentage = FALSE) {
  vi <- tidy_coef.multialpha.repeated.cv.glmnet(x$prod$model)

  names(vi)[names(vi) == "Coef"] <- "varImp"
  vi$varImp <- abs(vi$varImp)

  if ( scale && !percentage ) {
    vi$varImp <- vi$varImp / max(vi$varImp)
    names(vi)[names(vi) == "varImp"] <- "varImp_scaled"
  }
  if ( percentage ) {
    vi$varImp <- 100 * vi$varImp / max(vi$varImp)
    names(vi)[names(vi) == "varImp"] <- "varImp_pc"
  }
  vi
}

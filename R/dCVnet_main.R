
# Single Alpha Functions: ------------------------------------------------------
#   These functions wrap cv.glmnet to run over multiple specified folds and
#   pass back results along with selected best performing lambdas.
#   as suggested by the name, only one alpha is present.

repeated.cv.glmnet <- function(folds,
                               lambdas,
                               alpha,
                               x, y,
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
    x = results[,!names(results) %in% c("lambda", "rep",
                                        "lambda.min", "lambda.1se")],
    by = list(lambda = results$lambda),
    FUN = mean)

  # What optimising function should we use?:
  optfun <- ifelse("auc" %in% measure_name, max, min)

  # Which is the optimum lambda?
  av$lambda.min <- av$cvm == optfun(av$cvm, na.rm = T)

  # Deal with ties (multiple lambdas with the same minimum).
  nties <- sum(av$lambda.min)
  if ( nties > 1 ) {
    # warning(paste("Multiple lambdas:", ties))
    tiedlambdas <- av$lambda[av$lambda.min]
    # For odd length vectors the median exists in the data.
    #   otherwise return the larger lambda of the middle two.
    if ( length(tiedlambdas) %% 2 == 1 ) {
      av$lambda.min <- av$lambda == median(tiedlambdas)
    } else {
      larger_middle_lambda <- (length(tiedlambdas)/2) + 1
      av$lambda.min <- av$lambda == sort(tiedlambdas)[larger_middle_lambda]
    }
  }

  av[order(av$lambda, decreasing = T),]
  attr(av, 'type.measure') <- measure_name
  attr(av, 'class') <- c("repeated.cv.glmnet", "data.frame")

  # models and results could be useful, but we do not usually want these.
  if ( !debug ) {
    return(av)
  } else {
    return(list(averaged = av,
                results = results,
                models = models))
  }

}


summary.repeated.cv.glmnet <- function(obj) {
  alpha <- unique(obj$alpha)
  lambda <- range(obj$lambda)
  #best <- obj$cvm[obj$lambda.min]
  best_up <- obj$cvup[obj$lambda.min]
  best_lo <- obj$cvlo[obj$lambda.min]

  close <- sapply(unique(obj$lambda), function(x) {
    probe <- obj$cvm[obj$lambda == x]
    return(probe < best_up & probe > best_lo)
  })

  cat(paste("Tuning metric (cvm):", attr(obj, 'type.measure') , "\n"))
  cat(paste0("\talpha:\t\t", prettyNum(alpha), "\n"))
  cat(paste0("\tnlambda:\t", length(unique(obj$lambda)),"\n"))
  cat(paste0("\tlambda range:\t", paste(prettyNum(lambda), collapse = ", "), "\n\n"))
  cat("Best fit:\n")
  print.data.frame(obj[obj$lambda.min, ], row.names = F)
  cat(paste0("\n",sum(close), "/", length(close),
             " (", prettyNum(round(mean(close),2)), "%) ",
             "lambdas with CV fit within 1 S.E. of best fitting\n"))

  invisible(obj[obj$lambda.min, ])
}


# Multiple Alpha Functions: --------------------------------------------------

multialpha.repeated.cv.glmnet <- function(alphalist,
                                          lambdas,
                                          y,
                                          k,
                                          nrep,
                                          ...) {
  # Fold generation:
  folds <- lapply(1:nrep, function(i) {
    caret::createFolds(y = y, k = k, list = FALSE, returnTrain = FALSE)
  })

  malist <- lapply(1:length(alphalist),
                   function(i) {
                     cat(paste("Inner Alpha", i, "of",
                               length(alphalist), Sys.time(),"\n"))
                     a <- alphalist[[i]]

                     repeated <- repeated.cv.glmnet(alpha = a,
                                                    lambdas = lambdas[[i]],
                                                    y = y,
                                                    folds = folds,
                                                    ...)
                     repeated$alpha <- as.character(a)

                     return(repeated)
                   })
  tmeas <- attr(malist[[1]], "type.measure")
  malist <- do.call(rbind, malist)
  attr(malist, "type.measure") <- tmeas

  bestfun <- ifelse(tmeas == "auc", max, min)
  bestcandidates <- malist[malist$lambda.min,]
  best <- bestcandidates[bestcandidates$cvm == bestfun(bestcandidates$cvm),]

  # ties are broken by smaller cvsd followed by
  #   sparser solutions.
  if ( nrow(best) > 1 ) {
    best <- best[order(best$cvsd, best$nzero),]
    best <- best[1,]
  }

  # add column to the multialpha list.
  malist$best <- (malist$alpha == best$alpha) & (malist$lambda == best$lambda)

  return(structure(list(inner_results = malist,
                        inner_best = best,
                        inner_folds = folds),
                   class = "multialpha.repeated.cv.glmnet"))
}

summary.multialpha.repeated.cv.glmnet <- function(obj, print = T) {
  .get_nsimilar <- function(marc) {
    kk <- lapply(unique(marc$alpha), function(A){

      tdf <- marc[marc$alpha == A,]
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

  best <- obj$inner_best
  obj <- obj$inner_results

  alpha <- unique(obj$alpha)
  ## lambda <- range(obj$lambda)

  R <- obj[obj$lambda.min, c("alpha", "lambda", "nzero", "cvm", "cvsd", "cvup", "cvlo")]
  R <- merge(R, .get_nsimilar(obj), by = "alpha")

  selected <- (R$alpha == best$alpha) & (R$lambda == best$lambda)

  if ( print ) {
    cat("Summary of a multiple-alpha repeated k-fold cross-validated glmnet\n\n")
    cat(paste0("Alphas: ", paste(alpha, collapse = ", "), "\n\n"))

    cat("Best fitting alpha:\n")
    print(R[selected,], row.names = F)

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

#' Fit a doubly cross-validated elastic-net regularised (generalised) linear model
#'
#' repeated k-fold cross-validation used to:
#'     \itemize{
#'     \item{Produce unbiased estimates of out-of-sample classification performance (outer CV).}
#'     \item{Select optimal hyperparameters for the elasticnet (inner CV).}}
#'     Elasticnet hyperparameters are
#'     \bold{lambda} (the total regularisation penalty)
#'     and \bold{alpha} (the balance of L1 and L2 regularisation types).
#'
#' @param f a two sided formula. LHS must be binary and formula must include an intercept.
#' @param data data.frame containing all terms in f.
#' @param k_inner an integer, the k in the inner k-fold CV.
#' @param k_outer an integer, the k in the outer k-fold CV.
#' @param nrep_inner an integer, the number of repetitions (k-fold outer CV)
#' @param nrep_outer an integer, the number of repetitions (k-fold outer CV)
#' @param tuning_searchsize_lambda an integer, number of gradations between
#'     lambda.min and lambda.max to search.
#'     See \code{glmnet} argument \code{nlambda}.
#' @param alphalist a numeric vector of values in [0,1].
#'     Values of alpha to evaluate in inner cross-validation.
#' @param type.measure passed to \code{cv.glmnet}.
#'     The loss to use for hyperparameter selection in the inner cross-validation.
#'     Options: \code{"deviance"}, \code{"class"}, \code{"mse"}, \code{"mae"}
#' @param option.selection Method for hyperparameter selection in the inner cross-validation.
#'     One of \code{"default"} (lambda at best loss), \code{"lambda.3pc"} (add 3% to the default lambda),
#'     \code{"lambda.1se"} (add 1SE to the default lambda).
#'     Selection between alpha values is not affected by \code{option.selection}.
#' @param positive What level of the outcome is a 'positive' result
#'     (in the sense of a diagnostic test)
#' @param option.empirical_cutoff Boolian.
#'     Use the empirical proportion of cases as the cutoff for outer CV classification
#'     (affects outer CV performance only). Otherwise classify at 50\% probability.
#' @param debug Boolian.
#'     In 'debug' mode inner loop models are retained. This can produce very large return sizes.
#' @param speedup Boolian.
#'     Requests the faster \code{"modified.Newton"} method for logistic \code{glmnet}.
#'     See: \code{glmnet}'s \code{type.logistic} argument.
#' @return a dCVnet object.
#' @examples
#' \dontrun{
#' iris_class <- dCVreg(f = Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
#'                      data = subset(iris, iris$Species != c("versicolor")),
#'                      alphalist = 0.5)}
#' @importFrom stats aggregate as.formula coef glm model.frame model.matrix
#' @importFrom stats predict sd terms var
#' @export
dCVnet <- function(
  f, data,
  nrep_outer = 2,
  k_outer = 10,
  nrep_inner = 5,
  k_inner = 10,
  alphalist,
  type.measure,
  positive = 1,
  option.empirical_cutoff,
  nlambda,
  ...) {

  thecall <- match.call()

  time.start <- Sys.time() # for logging.

  parsed <- parse_input(f = f, data = data, positive = positive)
  x <- parsed$x_mat
  y <- parsed$y

  min_lambda_ratio <- ifelse(ncol(parsed$x_mat) > nrow(parsed$x_mat),
                             0.01, 0.0001)

  if ( type.measure == "auc" ) {
    # magic number for innerloop is:
    #   Ncases * outer train proportion * inner test proportion
    auc_magic <- (nrow(parsed$x_mat) * (1 - (1 / k_outer)) * (1 / k_inner))
    if (auc_magic < 11) {
      stop(
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
  startup_message(k_inner, nrep_inner,
                  k_outer, nrep_outer,
                  nalpha, nlambda,
                  parsed, time.start)
  # Main work starts here:

  # Step 1: make repeated outer folds in x & y according to nrep_outer, k_outer.
  outfolds <- caret::createMultiFolds(y = y,
                                      k = k_outer,
                                      times = nrep_outer)
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
                              length.out = nlambda)) })

  # are we working with empirical cutoffs?
  cutoff <- 0.5
  if ( option.empirical_cutoff ) {
    cutoff <- (as.numeric(table(y)[1]) / sum(as.numeric(table(y))))
  }
  cat(paste("Cutoff: ", cutoff, "\n"))

  # Main outer loop - for each outer fold run an inner cv loop
  outers <- parallel::mclapply(
    seq(along = outfolds),
    mc.cores = getOption("mc.cores", 1L),

    function(i) {
      cat(paste0("\n\nOuterloop:", i, " of ", imax, "\n"))
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
                                              ...)

      # extract the best alpha/lambda based on out of sample performance:
      fit_alpha <- as.numeric(inners$inner_best$alpha)
      fit_lambda <- inners$inner_best$lambda

      # fit a model to ALL the train data using the selected alpha/lambda:
      model <- glmnet(x = trainx,
                      y = trainy,
                      family = "binomial",
                      alpha = fit_alpha,
                      lambda = lambdas[[which(alphalist == fit_alpha)]])

      # how well did it do?:
      newx_prediction <- predict.lognet(model,
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
                                levels = c(F,T),
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
  performance <- do.call(rbind, lapply(outers, '[[', "performance"))
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
    type.measure = type.measure
  )

  final_model <- glmnet(
    x = xs,
    y = y,
    family = "binomial",
    alpha = final_tuning$inner_best$alpha,
    lambda = lambdas[[which(alphalist == final_tuning$inner_best$alpha)]])

  final_prediction <- predict.lognet(final_model,
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
                            levels = c(F,T),
                            labels = final_model$classnames))

  rownames(final_performance) <- rownames(final_prediction)

  final_performance <- structure(final_performance,
                                 class = c("classperformance",
                                           "data.frame"))

  final <- list(tuning = final_tuning,
                performance = final_performance,
                model = final_model)

  # reference models:
  # reference <- reference_models.dCVnet(parsed)
  #   these are large and fast to calculate.
  #   so don't store this.

  # The final object:
  obj <- structure(list(tuning = outers,
                        performance = performance,
                        folds = outfolds,
                        final = final,
                        #reference = reference,
                        input = list(parsed = parsed,
                                     lambdas = lambdas,
                                     alphas = alphalist)),
                   class = c("dCVnet","list"))


  time.stop <- Sys.time()

  cat("\n\n")
  print(thecall)
  cat("\n")
  cat(paste0("Finished.\n"))
  cat(paste0(time.stop, "\n"))
  cat(paste0((time.stop - time.start) / 60, " mins\n"))
  cat("------------")

  return(obj)
}

coefficients.dCVnet <- function(...) { coef.dCVnet(...) }

# Extract logreg coefficients from the outerloop best models:
coef.dCVnet <- function(obj, type = "all") {
  # Type can be:
  #   all - coefficients for each rep/fold.
  #   rep - mean average per rep.
  #   mean - mean of per-rep means.
  #   median = median of per-rep means.


  # # First handle reference cases:
  #   #   refglm = coefficients from the reference glm
  #   #   refuni = coefficients from the reference glmlist
  # if ( type == "refglm" ) { return(coef(obj$reference$glm)) }
  #
  # if ( type == "refuni" ) {
  #   CC <- vapply(obj$reference$univariate,
  #                function(x) { coef(x)[2] },
  #                numeric(1))
  #   names(CC) <- names(obj$reference$univariate)
  #   return(CC) }

  # next handle dCVnet cases:

  # works on the output of the outerloop.
  #   Given the 'best' alpha/lambda from the inner loop,
  #     what were the coefficients in a model selected with the best alpha.
  R <- lapply(1:length(obj$tuning), function(ii) {
    tt <- obj$tuning[[ii]]
    coefs <- as.matrix(glmnet::predict.lognet(tt$model,
                                              type = "coef",
                                              s = tt$tuning$inner_best$lambda))
    RR <- data.frame(Predictor = rownames(coefs),
                     Coef = c(coefs),
                     fold = names(obj$folds)[[ii]],
                     stringsAsFactors = F)
    return(RR)
  })
  R <- do.call(rbind, R)

  R$Predictor <- factor(R$Predictor,
                        levels = c("(Intercept)",
                                   obj$tuning[[1]]$model$beta@Dimnames[[1]]))
  R$fold <- factor(R$fold, levels = unique(R$fold))

  R$Rep <- sapply(strsplit(as.character(R$fold), split = "\\."), '[', 2)

  if (type == "all") { return(R) }
  # Aggregate by Repetition:
  R <- aggregate(R$Coef,
                 by = list(Predictor = R$Predictor,
                           Rep = R$Rep),
                 mean)
  names(R)[3] <- "Coef"
  return(switch(type,
                rep = R,
                mean = setNames(
                  aggregate(R$Coef,
                            by = list(Predictor = R$Predictor),
                            mean),
                  c("Predictor", "Coef")),
                median = setNames(
                  aggregate(R$Coef,
                            by = list(Predictor = R$Predictor),
                            median),
                  c("Predictor", "Coef")),
                stop(paste("type:", type ,
                           " - should be one of: all, rep, mean, median"))
  ))
}


summary.dCVnet <- function(obj) {

  # What do the 'best-fitting' results of the inner loops look like:
  R <- lapply(obj$tuning, function(x) {
    summary.multialpha.repeated.cv.glmnet(x$tuning, print = F)
  })

  R <- data.frame(do.call(rbind, aa))
  R <- R[R$best, names(R)[!names(R) %in% "best"]]
  rownames(R) <- names(obj$tuning)
  return(R)

}


# Code:
# dCVnet_functions.R
#
# Author          : Andrew J. Lawrence
# Version Number  : 0.6.1
# Date            : 26/07/2018

# ---------------------------------------------------------------------------- #
# Parallel Processing -------------------------------------------------------
# ---------------------------------------------------------------------------- #

# Set R option mc.cores to the number of cores desired for parallel processing
#   (if running on an linux/posix environment)
# to set:     options(mc.cores=3)
# to check:   getOption("mc.cores")
# to check with default 2L if not set:
#             getOption("mc.cores", 2L)

# ---------------------------------------------------------------------------- #
# Utility Functions -------------------------------------------------------
# ---------------------------------------------------------------------------- #

# Utility function to check packages, install those missing and then attach.
#   this function has no return and is employed for side effects.
# .initialise.dCVreg <- function() {
#   # Checks system packages and install any missing ones
#   wants <- c("caret",        # cross-validation functions; confusion matrices
#              "ggplot2",      # plotting
#              "glmnet",       # fit the glm with elasticnet regularisation
#              "mgcv",         # GAM smoother for summary of GoF over repeated CV
#              "ModelMetrics", # efficient AUC
#              "parallel",     # mclapply for parallel on linux.
#              "psych",        # used for psych::describe
#              "ROCR",         # used for ROC plots
#              "openxlsx"      # write results to logfile.
#   )
#   has <- wants %in% rownames(installed.packages())
#   if (any(!has)) install.packages(wants[!has])
#
#   # Attach the required packages.
#   suppressPackageStartupMessages(
#     lapply(wants, function(x) { library(x, character.only = TRUE) })
#   )
#   invisible()
# }


# Given a formula and dataframe:
#     - filter non-complete cases
#     - separate y variable
#     - make x_mat  a model.matrix input for cv.glmnet
#                   i.e. process interactions & dummy coding
#                   also strip intercept.
#     - make f0     a flattened formula (with intercept)
#                   based on the regressands in x_mat
parse_input <- function(f, data) {
  df <- model.frame(f, data = data,
                    drop.unused.levels = T) # does not include interaction.
  if (identical(nrow(df), 0L)) { stop("Error: no complete cases found.") }
  y <- df[,1] # extract y variable.
  x_mat <- model.matrix(f, data = data)[,-1] # drop intercept.
  f0.string <- paste(f[[2]],
                     "~",
                     paste(colnames(x_mat),
                           collapse = "+"))
  f0 <- as.formula(f0.string)
  return(list(y = y,
              x_mat = x_mat,
              f0 = f0))
}

# Utility function: calculate a list of lambdas given the data + alpha.
get_lambdalist <- function(x, y, alpha, nlambdas, min_lambda_ratio) {
  # Internal utility function. Given a dataset and alpha level
  #                     Return the maximum lambda value.
  .get_maxlambda <- function(x,y, alpha) {
    # https://stats.stackexchange.com/questions/166630/glmnet-compute-maximal-lambda-value
    # y must be 0/1 coded.
    # x must be a matrix.
    # If a zero alpha is fed in maxlambda is infinity - because of maths!
    n <- length(y) # number of subjects.

    # Alternative standardisation for the data:
    sdn <- function(kk) { sqrt(sum((kk-mean(kk))^2)/length(kk)) }
    x_sds <- apply(x, 2, sdn)

    # Scale x and y:
    x <- scale(x, scale = x_sds)
    y <- y - mean(y)*(1 - mean(y))

    # calculate the maximum lambda.
    max(abs(t(y) %*% x )) / (alpha * n)
  }
  maxl <- .get_maxlambda(x,y,alpha)
  lambdas <- exp(seq(from = log(maxl),
                 to = log(maxl * min_lambda_ratio),
                 length.out = nlambdas))
  return(lambdas)
}


# utility function to provide coefficients for the best fitting
#     repeated.cv.alpha.glmnet object
describe_best_rcag <- function(obj, option.selection = 'default') {
  # a result2 object (i.e. from repeated.cv.alpha.glmnet)
  mods <- obj$models
  alpharef <- mods$Rep1$alphalist

  if ( option.selection == "default" ) {
    BEST <- obj$best.param
  } else if ( option.selection == "lambda.3pc" ) {
    BEST <- obj$best.3pc
  } else if ( option.selection == "lambda.1se" ) {
    BEST <- obj$best.1se
  } else {
    stop("Problem with option.selection")
  }

  best.param.alpha.index <- which(alpharef == BEST$alpha)
  best.param.mod <- mods$Rep1$alphamodels[[best.param.alpha.index]]

  COEF <- coef(best.param.mod,
               s = BEST$lambda)

  return(list(coef = COEF,
              alpha = BEST$alpha,
              lambda = BEST$lambda))
}


# Utility function to extract cv.glmnet k-fold cross-validation information:
cv.glmnet.modelsummary <- function(mod,      # a cv.glmnet model
                                   alpha=NA, # optional: label with alpha
                                   rep=NA) { # optional: label with rep
  return(data.frame(lambda = mod$lambda,
                    cvm = mod$cvm,
                    cvsd = mod$cvsd,
                    alpha = alpha,
                    rep = rep))
}


# Utility function to return the best hyperparameters
#                         from a pooled cv.glmnet run:
get_best_hyperparams <- function(y, y_se, x, alpha, maximise = FALSE,
                                 type = "default") {
  # 3 types of 'best' hyperparameters:
  #     default - the conventional : lambda at the optimal model fitness.
  #     lambda.1se : largest lambda within 1se of fit of the optimal model fit
  #     lambda.3pc : default lambda + 3pc.

  # Internal function:
  #   return 1se given a standardised dataframe of results: y, y_se, x
  #     (note alpha is for label.)
  .get_1se <- function(df, alpha, maximise = FALSE, tol = 1e-6) {
    # df has y, x and y_se -- ALPHA MUST BE CONSTANT.

    # arrange the df vars:
    df <- df[, c("y","y_se","x")]

    # select the functions to be used depending on minimising or maximising:
    # optimise (min / max):
    optf <- match.fun(
      ifelse(identical(maximise,TRUE), "which.max", "which.min"))
    # create the 1se boundary value:
    boundf <- match.fun(ifelse(identical(maximise,TRUE), "-", "+"))
    # compare the values with the boundary:
    compf <- match.fun(ifelse(identical(maximise,TRUE), ">=", "<="))

    # Which is the 'greedy' best result:
    optrow <- optf(df$y)
    if ( is.na(optrow) ) { "Error: no optimum..." }

    # What are the contents of the optimum row:
    best <- df[optrow,]
    # what is the largest lambda assessed?
    max_x <- max(df$x)

    # Check if the optimum and maximum lambda are similar
    #   if so we are at the edge of parameter space and
    #   1se should just return the 'best':
    if ( abs(best$x - max_x) < tol ) {
      R <- best[,c("y","x")]
      names(R) <- c("fit","lambda")
      R$alpha <- alpha
      return(R[,c("lambda","alpha","fit")])
    }

    # filter to lambdas at least as penalised as the best:
    sefilter1 <- (df$x >= best$x)
    #   with model fits within 1se of the best fit
    sefilter2 <- compf(df$y, boundf(best$y,best$y_se))

    sefilter <- sefilter1 & sefilter2

    if ( !any(sefilter) ) {
      # this code runs where the sefilter is empty.
      #   an empty sefilter occurs when there are
      #   no models penalised worse than the best fitting model.
      R <- best[,c("y","x")]
      names(R) <- c("fit","lambda")
      R$alpha <- alpha
      # warning("ropey filters...") debug.
      return(R[,c("lambda","alpha","fit")])
    }

    # apply the filter:
    df <- df[sefilter, ]

    # return the largest lambda in the remaining dataset:
    R <- df[which.max(df$x),c("y","x")]
    names(R) <- c("fit","lambda")
    R$alpha <- alpha
    return(R[,c("lambda","alpha","fit")])
  }

  # get_best_hyperparameters - main:
  type <- match.arg(type, c("default", "lambda.1se", "lambda.3pc"))
  bestfun <- match.fun(
    ifelse(identical(maximise,TRUE), "which.max", "which.min"))

  df <- data.frame(y = y, y_se = y_se, x = x, alpha = alpha)

  # Two cases are simple as we just want the minimum or
  #     a linear xfm of the minumum:
  if ( type != "lambda.1se" ) {
    # do default and 3% here:
    R <- df[,c("x","alpha","y")]
    names(R) <- c("lambda", "alpha", "fit")
    # can just add 3% to lambda values and take the best fit:
    if ( type == "lambda.3pc" ) { R$lambda <- R$lambda*1.03 }
    return(R[bestfun(R$fit),])
  }

  # The last case is more complex:
  #   lambda.1se is the largest value of lambda (for a given alpha) which
  #   gives a fit within 1se of the optimum lambda.
  #   so we must iterate through all alphas, find the lambda.1se for each and
  #   then return the best fitting lambda.1se.
  R <- do.call(rbind,
               lapply(unique(df$alpha), function(tt) {
                 rdf <- df[df$alpha == tt, ]
                 .get_1se(rdf,
                          alpha = tt,
                          maximise = maximise)
               }))
  # return the best fitting 1se between the alphas:
  return(R[bestfun(R$fit),])
}


# utility function to return a list of univariate glm coefficients.
get_uni_fits <- function(f, df) {
  # f is the model formula, df is the data frame.
  # jiggle about the data to get a scaled model frame & design matrix:
  yname <- as.character(f[[2]]) # what is the outcome name.

  df <- model.frame(f, df)
  Y <- df[,yname]

  M <- model.matrix(f, df)
  M[,-1] <- scale(M[,-1]) # standardise the design matrix (not intercept)

  Mnames <- colnames(M)

  Ms <- lapply(Mnames, function(X) {
    xdf <- data.frame(y = Y, x = M[,X])
    return(glm( y ~ x, data = xdf, family = "binomial" ))
  } )
  result <- sapply(Ms, function(X) coef(X)[2])
  result[[1]] <- coef(Ms[[1]])[1]
  return(as.numeric(result))
}


# returns a flat string of ns and percents for levels of a factor.
factor_summary <- function(VAR, name) {
  if ( ! "factor" %in% class(VAR) ) { stop("Error: VAR must be a factor.") }
  TAB <- table(VAR) # frequencies.
  TOT <- sum(TAB)
  PERC <- round(100 * as.numeric(TAB / TOT), 1)
  paste0(
    name, "(n=", TOT, "): ",
    paste0(names(TAB), " = ", as.numeric(TAB),
           "(", PERC, "%)", collapse = "; "),
    collapse = ""
  )
}


# Calculate range, return as single string for display.
range_textsummary <- function(x, digits = 4, ...) {
  if ( !is.numeric(x) ) { stop("Error: input not numeric") }
  R <- prettyNum(range(x, na.rm = T), digits = digits, ...)
  return(paste0("[", R[[1]], ", ", R[[2]], "]"))
}
# Make a vectorised version:
range_textsummary_v <- base::Vectorize(range_textsummary)


# Utility function to produce list of prediction/label tables
#   one table for each repetition summarising data over the
#   k test folds in that k-fold repetition.
collate_results_over_reps <- function(obj) {
  # Input can either be a 'fit' object, or a 'outers' object.
  #   This code only needs the outers object.
  #   For convenience it checks if the input list has the name 'outers'
  #     then we treat it as a fit object and extract the outers.
  #       otherwise assume it is a 'outers' object.
  if ( "outers" %in% names(obj) ) { obj <- obj$outers }

  # Inner Utility fxn to get the indices assocated with each Rep.
  # returns: list of dataframes w/ preds, probs & labs
  .extract_reps <- function(obj) {
    # Input: a dCVreg object
    # Return: list of indices - RES
    #           RES[[n]] are the indices associated with the nth rep.
    # names(example.results$outers)[RES[[1]]] then gives
    #     the FoldIDs of a given rep.
    N <- data.frame(FoldID = names(obj), stringsAsFactors = F)
    N$RepID <- sapply(strsplit(N$FoldID, split = ".", fixed = T), '[',2)
    Reps <- unique(N$RepID)
    N$Index <- 1:nrow(N)
    RES <- lapply(Reps, function(x) { N$Index[N$RepID == x] } )
    names(RES) <- Reps
    return(RES)
  }

  # Which outers belong to which rep:
  RepInds <- .extract_reps(obj)
  # gather all test results per rep:
  collate <- lapply(RepInds, function(kk) {
    RES <- do.call(rbind, lapply(obj[kk],
                                 function(sobjs) {
                                   data.frame(A = sobjs$best_preds,
                                              B = sobjs$best_probs,
                                              C = sobjs$true_labs)
                                 }))
    rownames(RES) <- NULL
    colnames(RES) <- c("preds", "probs", "labs")
    return(RES)
  })
  return(collate)
}

# Utility function to produce summary tables
#   for a list of caret::ConfusionMatrices
summarise_cms <- function(cmlist, digits = 3) {
  # Input is a list of confusion matrices and a precision value
  c.t <- sapply(cmlist, function(x) as.data.frame(x$table)[,3])
  c.t.s <- apply(c.t, 1, sum, na.rm = T)

  # a classification table showing outer loop predicted and actual labels:
  classification.table <- data.frame(
    as.data.frame(cmlist[[1]]$table)[,1:2],
    Proportion = round(c.t.s / sum(c.t.s), digits),
    OuterRuns = rep("-"),
    c.t)

  # extract the confusion matrix statistics:
  c.r <- sapply(cmlist, function(x) { c(x$overall, x$byClass)})
  c.r.mean <- round(apply(c.r, 1, mean, na.rm = T), digits)
  c.r.sd <- round(apply(c.r, 1, sd, na.rm = T), digits)
  c.r.range <- range_textsummary_v(data.frame(t(c.r)), digits)

  # form the classification results showing the statistics
  classification.results <- data.frame(Mean = c.r.mean,
                                       SD = c.r.sd,
                                       Range = c.r.range,
                                       OuterRuns = rep("-"),
                                       c.r)
  return(list(
    classification.table = classification.table,
    classification.results = classification.results
  ))
}

# Utility function: produce a dataframe of model coefficients
#   employed for each fold/rep of the outer loop.
#   coefficients are extracted at the best lambda/alpha
#   chosen by the inner loop.
per_fold_model_coefficients <- function(obj) {
  if ( "outers" %in% names(obj) ) { obj <- obj$outers }
  # for above see collate_results_over_reps.
  R <- data.frame(sapply(obj, function(x) {
    as.matrix(coef(x$best_mod, s = x$best_params$lambda))
    } ))
  rownames(R) <- rownames(coef(obj[[1]]$best_mod))
  colnames(R) <- names(obj)
  return(R)
}

# Utility: return the fold membership frame for a dCVnet object (outers).
outerfoldmembership.dCVnet <- function(obj) {
  if ( "outers" %in% names(obj) ) { obj <- obj$outers } # work with outers.
  R <- lapply(obj, function(x) x$outerfolds)
  Nmax <- max(sapply(R, max))

  R <- do.call(data.frame, lapply(R, function(x) {
    factor(1:Nmax %in% x,
           levels = c(TRUE, FALSE),
           labels = c("Train", "Test")) }))
  return(R)
}

# ---------------------------------------------------------------------------- #
# AUC function for object -----------------------------------------------
# ---------------------------------------------------------------------------- #

# extract test auc values and final model auc (train+test).
auc.dCVreg <- function(x) {
  # extract vector of AUCs for the outer loops:
  .get_outer_test_aucs <- function(obj) {
    RES <- sapply(obj$outers, function(x) {
      ModelMetrics::auc(actual = x$true_labs,
                        predicted = x$best_probs)
    } )
    return(RES)
  }
  outer_test <- as.numeric(.get_outer_test_aucs(x))
  outer_test <- pmax(outer_test, 1 - outer_test)

  # reformat to include summary stats:
  outer_test <- data.frame(
    Mean = round(mean(outer_test, na.rm = T), digits = 4),
    SD = round(sd(outer_test, na.rm = T), digits = 4),
    Range = range_textsummary(outer_test, 4),
    OuterRuns = "-",
    t(round(outer_test, digits = 4)), stringsAsFactors = F)
  colnames(outer_test)[5:ncol(outer_test)] <- names(x$outers)
  rownames(outer_test) <- "AUC"
  # extract auc for final model:
  .get_final_auc <- function(obj) {
    ModelMetrics::auc(actual = obj$log$data[,as.character(obj$log$f[[2]])],
                      predicted = obj$final.model$best_probs)
  }
  final_all <- .get_final_auc(x)
  final_all <- pmax(final_all, 1 - final_all)

  # Daniel Stahl Suggestion:
  #   AUC for the pooled predictions from all test-folds in a repetition.
  .get_rep_AUC <- function(obj) {
    # function to calculate per-repetition AUC
    #   where AUC is calculated given the predictions from the k-folds.
    Collated <- collate_results_over_reps(obj)
    sapply(Collated,
           function(X) {
             ModelMetrics::auc(
               actual = X$labs,
               predicted = X$probs) })
  }
  rep_AUC <- .get_rep_AUC(x) # vector of AUCS:
  rep_AUC <- pmax(rep_AUC, 1 - rep_AUC)

  rep_AUC <- data.frame(
    Mean = round(mean(rep_AUC, na.rm = T), digits = 4),
    SD = round(sd(rep_AUC, na.rm = T), digits = 4),
    Range = range_textsummary(rep_AUC, digits = 4),
    OuterRuns = "-",
    t(rep_AUC),
    row.names = "AUC")

  # Return object:
  return(list(outer_folds = outer_test,
              outer_reps = rep_AUC,
              final_all = final_all))
}


# ---------------------------------------------------------------------------- #
# Log file function for object --------
# ---------------------------------------------------------------------------- #

# Generate a logfile given a model object - no return, used for side effect
#  of writing to file.
log.dCVreg <- function(X, fileroot) {
  # X == a dCVreg results object
  # fileroot == a path / root to write to:
  #           e.g. 'path/to/example' gives: 'path/to/example_log.xlsx'
  logfile <- paste0(fileroot, "_log.xlsx")

  nit_inner <- X$log$nrep_inner*X$log$k_inner*X$log$tuning_searchsize_alpha
  nit_outer <- X$log$nrep_outer*X$log$k_outer
  nit_total <- nit_inner * nit_outer

  # First: write the run notes:
  out.Notes <- data.frame(ModelNotes = c(
    paste(X$log$cmd[[1]],
          "glmnet-regularised logistic regression with estimation of generalisation error by repeated k-fold cross-validation - includes fully nested selection of optimal hyperparameters by repeated k-fold crossvalidation."),
    "",
    paste("Start:", X$log$start.time),
    paste("Finish:", X$log$stop.time),
    paste("Runtime(mins):",
          round(as.numeric(difftime(X$log$stop.time,
                                    X$log$start.time,
                                    units = "mins")),1)),
    "",
    paste("Data:", X$log$cmd$data),
    #paste("Output (Results):", example.outfile),
    #paste("Save file:", example.outsave),
    #paste("Tuning Grid Plots:", example.tuningplots),
    "",
    paste("Model Formula:", deparse(X$log$f, width.cutoff = 20L)),
    "",
    paste("Tuning metric: ", X$log$type.measure),
    paste("Hyperparameter Selection Method:", X$log$option.selection),
    paste("Use Empirical Cutoff for Classification?:", X$log$option.empirical_cutoff),
    "",
    paste("Outer # Folds:", X$log$k_outer),
    paste("Outer # Repetitions:", X$log$nrep_outer),
    paste("Inner # Folds:", X$log$k_inner),
    paste("Inner # Repetitions:", X$log$nrep_inner),
    paste("Lambda Searchsize:", X$log$tuning_searchsize_lambda),
    paste("Alpha Searchsize:", X$log$tuning_searchsize_alpha),
    paste("Alpha Values:", paste0(X$log$alphalist, collapse = ", ")),
    "",
    paste("# Runs per Inner Loop:", nit_inner),
    paste("# Runs per Outer Loop:", nit_outer),
    paste("# Runs Total:", nit_total)
  ),
  stringsAsFactors = F)

  # Descriptives for input
  # 1. Factors (at least outcome must be a factor)
  var.type <- sapply(X$log$data, function(xx) "factor" %in% class(xx))

  factor.vars <- names(which(var.type))
  cont.vars <- names(which(!var.type))

  # make a datasheet for factors:
  out.xfactors <- do.call(rbind,
                          lapply(factor.vars,
                                 function(VV) {
                                   factor_summary(X$log$data[[VV]], VV)
                                 } ))
  colnames(out.xfactors) <- "Factor Summary"

  out.cont <- psych::describe(X$log$data[,cont.vars])

  out.key <- data.frame(Sheets = c(
    "Notes - Analysis parameters timestamps etc",
    "Descriptives1 - Predictor Descriptives (factors)",
    "Descriptives2 - Predictor Descriptives (continuous)",
    "ClassStats - Outer loop results (per repetition)",
    "ClassTables - Confusion matrices for the outer loops (per repetition)",
    "FoldClassStats - Outer loop results (per fold)",
    "FoldClassTables - Confusion matrices for the outer loops (per fold)",
    "TuningParams - optimal tuning parameters selected each outer loop.",
    "ModelCoefficients - coefficient vector for each fold of the outer loop",
    "FoldMembership - Which subjects were in which fold.",
    "FinalModCoefficients - regularised model coefficients of a 'production' model alongside unregularised logistic regression coefficients for reference. All coefficients are standardised regression betas.",
    "FinalModTuningParams - the optimum tuning parameters selected for the 'production' model.",
    "FinalModClassStats - final model classification statistics  (Caution: Optimistic!)",
    "FinalModClassTable - Confusion matrices for the final model (Caution: Optimistic!)",
    "SheetKey - a key for the output contents."
  ))

  # Gather output into a named list:

  # Do we include row names?
  log.opts <- list(
    Notes = F,
    Descriptives1 = F,
    Descriptives2 = T,
    ClassStats = T,
    ClassTables = F,
    FoldClassStats = T,
    FoldClassTables = F,
    TuningParams = T,
    ModelCoefficients = T,
    FoldMembership = T,
    FinalModCoefficients = T,
    FinalModTuningParams = T,
    FinalModClassStats = T,
    FinalModClassTable = F,
    SheetKey = F
  )

  log <- list(
    Notes = out.Notes,
    Descriptives1 = out.xfactors,
    Descriptives2 = out.cont,
    ClassStats = X$rep.class.results,
    ClassTables = X$rep.class.table,
    FoldClassStats = X$fold.class.results,
    FoldClassTables = X$fold.class.table,
    TuningParams = X$tuning.results,
    ModelCoefficients = X$fold.coefs,
    FoldMembership = outerfoldmembership.dCVnet(X),
    FinalModCoefficients = X$final.model$coef,
    FinalModTuningParams = X$final.model$tuning,
    FinalModClassStats = X$final.model$classification.results,
    FinalModClassTable = X$final.model$classification.table,
    SheetKey = out.key
  )


  # Write all this out:
  wb <- openxlsx::createWorkbook()

  Map(function(x,RN,Name) {
    openxlsx::addWorksheet(wb, sheetName = Name)
    openxlsx::writeData(wb, sheet = Name, x = x, rowNames = RN)
    # Hack: RN is logical and coerced to 0/1 in sum for column indices:
    openxlsx::setColWidths(wb, sheet = Name, widths = "auto", cols = 1:(ncol(x) + RN))
    #print(Name)
    return(NULL)
  }, x = log, RN = log.opts, Name = names(log))

  saveWorkbook(wb, file = logfile, overwrite = T)

  invisible()
}


# ---------------------------------------------------------------------------- #
# Plot Hyperparams function for object --------
# ---------------------------------------------------------------------------- #

# Make a multipage pdf showing inner loop performance
#   and selected hyperparameters per each outer loop.
plot.hp.dCVreg <- function(X, fileroot) {
  # X == a dCVreg results object
  # fileroot == a path / root to write to:
  #           e.g. 'path/to/example' gives: 'path/to/example_log.xlsx'

  plotfile <- paste0(fileroot, "_plots.pdf")

  # extract the pooled model fitness (averaged over repetitions for each alpha):
  tsearch <- lapply(X$outers, function(x) {
    kk <- x$inners$pooledresults
    kk$alpha <- as.factor(kk$alpha)
    return(kk)
  })
  names(tsearch) <- names(X$outers)

  # Get global limits for the plots:
  tsearch_lambdarange <- t(sapply(tsearch, function(x) range(x$lambda)))
  tsearch_lambdarange <- c(min(tsearch_lambdarange[,1]),
                           max(tsearch_lambdarange[,2]))
  tsearch_cvmrange <- t(sapply(tsearch, function(x) {
    range(c(x$cvm + x$cvsd,x$cvm - x$cvsd))
  }))
  tsearch_cvmrange <- c(min(tsearch_cvmrange[,1]),max(tsearch_cvmrange[,2]))

  tbests <- lapply(X$outers, function(x) {
    tbest <- x$best_params
    tbest$lab_cvm <- tsearch_cvmrange[2] # put lab at max.
    tbest$lab_lambda <- tsearch_lambdarange[2]
    tbest$label <- paste0("Chosen Param.\nfit: ",
                          round(tbest$fit,3), "\nalpha: ",
                          tbest$alpha, "\nlambda: ",
                          round(tbest$lambda,3))
    return(tbest)
  })
  names(tbests) <- names(X$outers)

  pdf(file = plotfile, width = 7, height = 7, onefile = T)
  for (i in seq(along = tsearch)) {
    cat( paste0("Plotting: ", i, " of ", length(tsearch), "\n") )
    best <- tbests[[i]]
    names(best)[names(best) == "fit"] <- "cvm"
    p <- ggplot2::ggplot(tsearch[[i]],
                         ggplot2::aes(y = cvm,
                                      ymin = cvm - cvsd,
                                      ymax = cvm + cvsd,
                                      x = lambda,
                                      colour = alpha,
                                      group = alpha)) +
      #ggplot2::stat_smooth(alpha = 1, geom = "line", method = "loess",
      #                     se = F, span = 0.1) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = best, ggplot2::aes(y = cvm, x = lambda),
                          colour = "black", shape = 24, inherit.aes = F) +
      ggplot2::geom_text(data = best,
                         ggplot2::aes(y = lab_cvm,
                                      x = lab_lambda,
                                      label = label),
                         colour = "black", inherit.aes = F,
                         vjust = "inward", hjust = "inward") +
      ggplot2::scale_x_log10(limits = tsearch_lambdarange) +
      ggplot2::scale_y_continuous(limits = tsearch_cvmrange) +
      ggplot2::xlab("lambda (log 10)") +
      ggplot2::ylab(paste0("Mean +/- 1SD model ", X$log$type.measure)) +
      ggplot2::theme_light()
    print(p)
  }
  dev.off()
}


# ---------------------------------------------------------------------------- #
# Plot ROC curves function for object --------
# ---------------------------------------------------------------------------- #

# make a multipage pdf showing ROC curves for:
#   1) final model all data (train + test)
#   2) outer loops test data.
plot.roc.dCVreg <- function(X, fileroot, dumpraw = F) {
  # Utility subfunction convert a ROC performance object to a dataframe.
  .performance_to_data_frame <- function(perf, names) {
    ns <- sapply(perf@y.values, length)
    runs <- rep(names, ns)
    # convert to dataframe for ggplot:
    outer.df <- data.frame(y = do.call(c, perf@y.values),
                           x = do.call(c, perf@x.values),
                           alpha = do.call(c, perf@alpha.values),
                           run = runs)
    return(outer.df)
  }
  Collated <- collate_results_over_reps(X)
  # multipred object:
  outer.pred <- ROCR::prediction(
    predictions = lapply(Collated, '[[', "probs"),
    labels = lapply(Collated, '[[', "labs"))
  outer.perf <- performance(outer.pred, "tpr", "fpr")
  # single for final performance:
  final.pred <- ROCR::prediction(
    predictions = X$final.model$best_probs,
    labels = X$log$data[,as.character(X$log$f[[2]])])
  final.perf <- performance(final.pred, "tpr", "fpr")

  outer.df <- .performance_to_data_frame(outer.perf, names(Collated))
  final.df <- .performance_to_data_frame(final.perf, "FinalModel")

  # Plot the outer results:
  p <- ggplot2::ggplot(outer.df, aes(y = y,
                                     x = x,
                                     group = run,
                                     colour = run)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "black", size = 1.1) +
    ggplot2::geom_line() +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    ggplot2::guides(colour = F) +
    ggplot2::theme_classic()
  ggplot2::ggsave(filename = paste0(fileroot, "_ROC_outers_test.pdf"))

  # Plot the final results:
  p <- ggplot2::ggplot(final.df, aes(y = y,
                                     x = x,
                                     group = run,
                                     colour = run)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "black", size = 1.1) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    ggplot2::theme_classic()
  ggplot2::ggsave(filename = paste0(fileroot, "_ROC_final_traintest.pdf"))

  # return the raw data:
  return(list(outer = outer.df, final = final.df))
}


# ---------------------------------------------------------------------------- #
# Subfunctions of dCVreg -----------------------------------------
# ---------------------------------------------------------------------------- #
# 3 nested subfunctions which do the main body of work behind dCVreg

# Lowest layer runs cv.glmnet over a number of alphas (not repeated)
cv.alpha.glmnet <- function(alphalist,
                            foldids,
                            y,
                            x,
                            nfolds=10,
                            rep=NA, # optional label for output
                            type.measure,
                            nlambda,
                            min.lambda.ratio,
                            ...) {
  # function is a wrapper to cv.glmnet
  #  now *have* to specify previously optional nfolds and foldids
  #  dots (...) pass rest to cv.glmnet
  names(alphalist) <- paste0("Alpha",seq(along = alphalist))

  # Create a per-alpha list of lambda values to use
  #       - this will be the same over repetitions.
  lambdalist <- lapply(alphalist, function(aa) {
    # for each alpha calculate the maximum lambda
    # return a list of lambda geometrically progressing to maxlambda * minlambdaratio.
    get_lambdalist(x = x,
                   y = as.numeric(y)-1,
                   alpha = aa,
                   nlambdas = nlambda,
                   min_lambda_ratio = min.lambda.ratio)
  } )

  # Basic data blob
  mods <- lapply(seq(along = alphalist),
                 function(kkk) {
                   A <- alphalist[[kkk]]
                   L <- lambdalist[[kkk]]
                   glmnet::cv.glmnet(y = y, x = x, alpha = A,
                                     type.measure = type.measure,
                                     nfolds = nfolds,
                                     foldid = foldids,
                                     lambda = L,
                                     ...)
                 })
  names(mods) <- names(alphalist)

  # Apply an extraction function to return all crossvalidated fits over
  #     all the alphas
  result_frame <- do.call(rbind,
                          mapply(cv.glmnet.modelsummary,
                                 mod = mods,
                                 alpha = alphalist,
                                 rep = rep(rep, length(alphalist)),
                                 SIMPLIFY = F, USE.NAMES = F))
  # returns a copy of alphalist and a list of associated results
  return(list(alphamodels = mods,
              alphalist = alphalist,
              result = result_frame))
}


# Middle function repeats application of cv.alpha.glmnet and pools
#   results to give a repeated k-fold crossvlalidation.
repeated.cv.alpha.glmnet <- function(nreps,
                                     nfolds,
                                     y,
                                     type.measure,
                                     ...) {
  # generate list random folds for each repetition / fold:
  fixfolds <- lapply(1:nreps, function(i) {
    caret::createFolds(y = y, k = nfolds, list = FALSE, returnTrain = FALSE)
  })

  # Run cv.alpha.glmnet for each fold (rep/cv).
  blob <- lapply(seq(along = fixfolds), function(i) {
    f <- fixfolds[[i]]
    cv.alpha.glmnet(y = y, rep = i,
                    type.measure = type.measure, foldids = f, ...)
  })
  names(blob) <- paste0("Rep",1:nreps)

  # extract the full table:
  pooledresults <- do.call(rbind, lapply(blob, '[[', "result"))

  # save the original pooledresults if debug:
  rawpooledresults <- "NA: Debug is off"
  if ( isTRUE(debug) ) { rawpooledresults <- pooledresults }

  # Mean average over reps:
  pooledresults <- aggregate(
    x = list(cvm = pooledresults$cvm,
             cvsd = pooledresults$cvsd),
    by = list(alpha = pooledresults$alpha,
              lambda = pooledresults$lambda),
    FUN = mean)

  # Are we minimising or maximising:
  optfun <- "which.min"
  if ( type.measure %in% c("auc") ) {
    optfun <- "which.max"
  }

  # what are the best hyperparameters
  #       (generate all 3 methods and choose in the outer loop.):
  best.param <- get_best_hyperparams(
    y = pooledresults$cvm,
    y_se = pooledresults$cvsd,
    x = pooledresults$lambda,
    alpha = pooledresults$alpha,
    maximise = ifelse(type.measure %in% c("auc"),
                      TRUE, FALSE),
    type = "default")

  best.1se <- get_best_hyperparams(
    y = pooledresults$cvm,
    y_se = pooledresults$cvsd,
    x = pooledresults$lambda,
    alpha = pooledresults$alpha,
    maximise = ifelse(type.measure %in% c("auc"),
                      TRUE, FALSE),
    type = "lambda.1se")

  best.3pc <- get_best_hyperparams(
    y = pooledresults$cvm,
    y_se = pooledresults$cvsd,
    x = pooledresults$lambda,
    alpha = pooledresults$alpha,
    maximise = ifelse(type.measure %in% c("auc"),
                      TRUE, FALSE),
    type = "lambda.3pc")

  return(list(models = blob,
              pooledresults = pooledresults,
              rawpooledresults = rawpooledresults,
              best.param = best.param,
              best.3pc = best.3pc,
              best.1se = best.1se,
              alphalist = blob[[1]]$alphalist, # copy of alphalist.
              innerfolds = fixfolds))
}


# This is the manager of the outer loop.
#   it runs the outer loop and evaluates the results.
nested.repeated.cv.alpha.glmnet <- function(
  nrep_outer = 5,
  k_outer = 10,
  nrep_inner = 5,
  k_inner = 10,
  y, x,
  type.measure,
  option.selection = option.selection,
  option.empirical_cutoff = option.empirical_cutoff,
  report.precision = 3L,
  debug = FALSE,
  ...) {

  # Step 1: make repeated outer folds in x & y according to nrep_outer, k_outer.
  outfolds <- caret::createMultiFolds(y = y, k = k_outer, times = nrep_outer)
  imax <- length(outfolds)
  names(outfolds) <- paste0("Out", names(outfolds)) # give names

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
      cat(paste0("Outerloop:",i," of ", imax, "\n"))
      of = outfolds[[i]]
      # first generate the train and test data from the fold:
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

      inners <- repeated.cv.alpha.glmnet(nreps = nrep_inner,
                                         nfolds = k_inner,
                                         y = trainy, x = trainx,
                                         type.measure = type.measure,
                                         ...)
      # What method to pick the best hyperparam.
      if ( option.selection == "default" ) {
        BEST <- inners$best.param
      } else if ( option.selection == "lambda.3pc" ) {
        BEST <- inners$best.3pc
      } else if ( option.selection == "lambda.1se" ) {
        BEST <- inners$best.1se
      } else {
        stop("Problem with option.selection")
      }

      # Extract the model with the inner loop optimal alpha/lambda
      #   per each outer loop:
      best_alpha_ref <- inners$models$Rep1$alphalist
      best_alpha_index <- which(best_alpha_ref == BEST$alpha)
      best_lambda <- BEST$lambda
      best_mod <- inners$models$Rep1$alphamodels[[best_alpha_index]]

      # Extract predictions from the test data:
      #   continuous:
      best_probs <- glmnet::predict.cv.glmnet(best_mod,
                                              newx = testx,
                                              s = best_lambda,
                                              type = "response")
      #   categorical:
      best_preds <- factor(ifelse(best_probs > cutoff,
                                  levels(y)[2],
                                  levels(y)[1]),
                           levels = levels(y))

      # strip inner models
      #  to reduce the size (anon function used for side effect.):
      if (!isTRUE(debug)) { inners$models <- NULL }

      return(list(inners = inners,
                  best_preds = best_preds,
                  best_probs = best_probs,
                  true_labs = testy,
                  best_params = BEST,
                  best_mod = best_mod,
                  outerfolds = of))
    })
  names(outers) <- names(outfolds) # add back names


  # gather the tuning results (best selected) for each outer loop.
  tuning.results <- t(sapply(outers, function(x) { x$best_params }))

  # gather confusion matrices of best predictions on the outer folds/reps:
  fold_results <- lapply(outers,
                         function(x) {
                           caret::confusionMatrix(
                             data = x$best_preds,
                             reference = x$true_labs) })

  # and collated over folds, so 1 set per repetition.
  rep_results <- lapply(collate_results_over_reps(outers),
                        function(x) {
                          caret::confusionMatrix(
                            data = x$preds,
                            reference = x$labs) } )

  # Apply a function to summarise confusion matrices:
  fold_results <- summarise_cms(fold_results, report.precision)
  rep_results <- summarise_cms(rep_results, report.precision)

  # Also extract...
  fold_cs <- per_fold_model_coefficients(outers)
  # ...and summarise the model coefficients:
  fold_coefs <- data.frame(
    Mean = round(apply(fold_cs, 1, mean), digits = report.precision),
    SD = round(apply(fold_cs, 1, sd), digits = report.precision),
    PerFold = rep("-"),
    fold_cs,
    row.names = rownames(fold_cs)
  )

  # Finally produce a 'grand' model over all the data (train and test):
  #   this is the 'production model' useful for coefficients.
  #   it should not be evaluated in terms of classification accuracy.
  final.PPx <- caret::preProcess(x)   # preprocessed data (centered & scaled)
  xs <- predict(final.PPx, x)

  # Run the inner loop for the final model:
  final.model <- repeated.cv.alpha.glmnet(
    nreps = nrep_inner,
    nfolds = k_inner,
    y = y, x = xs,
    type.measure = type.measure,
    ...
  )
  # what were the coefs for the best model?
  final.model.coef <- describe_best_rcag(final.model,
                                         option.selection = option.selection)
  # and the hyperparameters:
  final.model.tuning <- data.frame(alpha = final.model.coef$alpha,
                                   lambda = final.model.coef$lambda)
  rownames(final.model.tuning) <- "final model"

  # extract best bits of final model:
  final.best_alpha_ref <- final.model$models$Rep1$alphalist
  final.best_alpha_index <- which(final.best_alpha_ref ==
                                    final.model.tuning$alpha)
  final.best_lambda <- final.model.tuning$lambda
  final.best_mod <- final.model$models$Rep1$alphamodels[[final.best_alpha_index]]

  # class probabilities
  final.best_probs <- glmnet::predict.cv.glmnet(final.best_mod,
                                                newx = xs,
                                                s = final.best_lambda,
                                                type = "response")
  # class predictions
  final.best_preds <- factor(ifelse(final.best_probs > cutoff,
                                    levels(y)[2], levels(y)[1]),
                             levels = levels(y))

  final.cm <- caret::confusionMatrix(data = final.best_preds,
                                     reference = y)
  # Final Classification Table:
  final.cm.table <- as.data.frame(final.cm$table)
  final.cm.table$Prop <- final.cm.table$Freq / sum(final.cm.table$Freq,
                                                   na.rm = T)

  # Final Classification Stats:
  final.cm.stats <- data.frame(c(final.cm$overall, final.cm$byClass))
  names(final.cm.stats) <- "FinalModel"

  final.model <- list(outers = final.model,
                      PPx = final.PPx,
                      tuning = final.model.tuning,
                      coef = final.model.coef[[1]],
                      best_mod = final.best_mod,
                      best_probs = final.best_probs,
                      best_preds = final.best_preds,
                      classification.results = final.cm.stats,
                      classification.table = final.cm.table)

  return(list(outers = outers,
              tuning.results = tuning.results,
              rep.class.results = rep_results$classification.results,
              rep.class.table = rep_results$classification.table,
              fold.class.results = fold_results$classification.results,
              fold.class.table = fold_results$classification.table,
              fold.coefs = fold_coefs,
              final.model = final.model))
}


# ---------------------------------------------------------------------------- #
# Main function --------
# ---------------------------------------------------------------------------- #
#   runs a regularised logistic regression with
#     nested crossvalidation of elasticnet parameters alpha/lambda.
dCVreg <- function(f,
                    data,
                    nrep_outer = 5,
                    k_outer = 10,
                    nrep_inner = 5,
                    k_inner = 10,
                    type.measure = "deviance",
                    option.selection = 'default',
                    option.empirical_cutoff = FALSE,
                    tuning_searchsize_lambda = 1000,
                    alphalist = seq(from = 0, to = 1.0, by = 0.1),
                    debug = FALSE,
                    speedup = FALSE
                   ) {
  # FUNCTION NOTES.
  #
  # Required arguments:
  #   f - an r formula
  #     Note: this must include an intercept as per R default behaviour.
  #           if you specify y ~ 0 + x1 + x2 or y ~ x1 + x2 - 1
  #           then results will be unpredictable.
  #           This is true even if a variable (e.g. x1) is a hard-coded
  #           intercept variable. The appropriate attribute of terms(f)
  #           must be set correctly.
  #
  #   data - an r dataframe.
  #     Note1: only complete cases for terms of f in data will be analysed.
  #     Note2: after removing incomplete cases unused factor levels of both
  #            outcome and predictors will be dropped.
  #
  # Optional arguments:
  #   k (outer + inner)     - the number of folds for k-fold cv
  #   nrep (outer + inner)  - number of times to repeat k-fold cv
  #   type.measure          - which metric to use to select hyperparameters
  #                           (inner loop).
  #   option.selection      - how to choose the best hyperparameters
  #                           (default, lambda.1se, lambda.3pc)
  #
  #   tuning_searchsize_lambda
  #     - determines the number of values of lambda to consider
  #       (per alpha value).
  #       lambda is the amount of regularisation of coefficients
  #       and lambda values are selected automatically between
  #       (high to low): lambda.max to (0.0001 * lambda.max)
  #       [or 0.01 if more vars than cases].
  #
  #       high to low with the largest lambda being
  #       the smallest lambda to produce a model
  #         - This script tunes 2 parameters: alpha & lambda
  #               lambda = the amount of regularisation penalty which
  #                         scales between zero (no penalty) and
  #                         the smallest lambda which results in all
  #                         coefficients being zero (set to 1).
  #               alpha  = the elastic net mixing paramter between the
  #                         L2 (ridge) and L1 (lasso) penalties.
  #                         alpha = 0 is ridge, alpha = 1 is lasso.
  #       tuning_searchsize_lambda determines the number of parameter values
  #         to consider per parameter.
  #
  #   tuning_metric
  #               - on what basis do we choose the best tuning parameters.
  #   debug - if true, we will keep all the inner results.

  ## Internal functions to wrap cv.glmnet and format results:
  ##  1. cv.alpha.glmnet
  ##        run glmnet over a list of alphas
  ##        needs a list of alphas.
  ##        passes up a performance dataframe.
  ##  2. repeated.cv.alpha.glmnet
  ##        repeatedly call 1.
  ##        to run glmnet over a list of alphas
  ##        generating fixedfolds / repetitions as required.
  ##        merges performance dataframes and returns
  ##        best parameters.
  ##  3. nested.repeated.cv.alpha.glmnet
  ##        repeatedly call 2.
  ##        manage the outer loop folds/reps.

  # FUNCTION BODY:

  # Record call for posterity:
  thecall <- match.call()

  # Check input 1: formula is formula and has intercept.
  f <- as.formula(f) # attempt to coerce in case of character string.
  if ( !identical(attr(terms(f),"intercept"), 1L) ) {
    stop("Error: formula must have an intercept. See terms(f)")
  }

  # Check input 2: data is a data.frame:
  if ( !"data.frame" %in% class(data) ) {
    stop("Error: data must be data.frame.")
  }

  parsed <- parse_input(f, data)

  min_lambda_ratio <- ifelse(ncol(parsed$x_mat) > nrow(parsed$x_mat), 0.01, 0.0001)

  # extract outcome and xvars from the input:
  #   NOTE: as a result only complete cases are included...
  df <- model.frame(f, data = data,
                    drop.unused.levels = T) # does not include interaction.

  # check the y-variable:
  if ( !"factor" %in% class(parsed$y) ) {
    # allows "ordered" as well as "factor".
    stop("Error: outcome is not a factor.")
  }
  if ( !identical(nlevels(parsed$y), 2L) ) {
    # *Must* be binary.
    stop("Error: outcome is not binary.")
  }

  # Confirm than all predictors have non-zero variance after dropping
  #   incomplete cases
  if ( any(sapply(as.data.frame(parsed$x_mat), var) == 0.0)) {
    stop("One or more predictors has zero variance")
  }

  # Check type.measure argument:
  #   We need at least 10 subjects per inner test fold to use test fold
  #   AUC as the cost function for hyperparameter selection.
  #   If AUC is requested and req isn't met: print a warning and use deviance.
  if ( type.measure == "auc" ) {
    # magic number for innerloop is:
    #   Ncases * outer train proportion * inner test proportion
    auc_magic <- (nrow(parsed$x_mat) * (1 - (1/k_outer)) * (1 / k_inner))
    if (auc_magic < 11) {
      warning(
        paste("Warning: using model deviance - not AUC - as sample size small!
               Estimated inner foldsize =", auc_magic))
      type.measure <- "deviance"
    }
  }

  # how many alphas did we feed in:
  tuning_searchsize_alpha = length(alphalist)
  # Check alpha values are OK:
  if ( any(is.na(alphalist)) ) { stop("Error: missing value(s) in alphalist") }
  if ( min(alphalist) < 0.0 ) { stop("Error: alphas must be positive") }
  if ( max(alphalist) > 1L ) { stop("Error: alphas must be <= 1") }

  # substitute zero-length alphas:

  # N.B. There appears to currently be a bug in the glmnet package (2.0.13)
  #   it appears to occur in the conjunction of the following 3 conditions:
  #         * use of fixed folds
  #         * pure ridge regression (alpha = 0.0)
  #         * large number of different lambda values
  # Number 1 is non-negotiable for our use case (repeated cv).
  # Number 3 is very useful to give stable results for hyperparameter selection.
  #
  # So we will change number 2 as follows:
  #   Replace any 'pure' ridge with a strongly ridge-like elastic-net
  #   (either alpha = 0.01, or half the smallest nonzero alpha
  #     if this is lower than 0.01.)

  if ( any(alphalist == 0.0) ) {
    replacement_alpha <- min(c(0.01, 0.5 * min(alphalist[alphalist > 0.0])))
    alphalist[alphalist == 0.0] <- replacement_alpha
    warning(paste0("BUGFIX: alpha=0.0 replaced with: ", replacement_alpha))
  }
  # Note: often low values of alpha (<0.1) are very slow to fit.

  # Print some general run info to screen.
  nit_inner <- k_inner*nrep_inner*tuning_searchsize_alpha
  nit_outer <- k_outer*nrep_outer
  nit_total <- nit_inner * nit_outer

  cat(paste0("Analysing n=",nrow(df)," subjects\n"))
  cat(paste0(table(parsed$y)[1], " of outcome: ", levels(parsed$y)[1], "\n"))
  cat(paste0(table(parsed$y)[2], " of outcome: ", levels(parsed$y)[2], "\n"))
  cat(paste0("running ", nit_inner , " inner runs\n"))
  cat(paste0("running ", nit_outer, " outer runs\n"))
  cat(paste0("total: ", nit_total, " runs\n"))

  # Log start time for model fitting:
  time.start <- Sys.time()

  # Outerloop:
  fit <- nested.repeated.cv.alpha.glmnet(
    nrep_outer = nrep_outer,
    k_outer = k_outer,
    nrep_inner = nrep_inner,
    k_inner = k_inner,
    y = parsed$y,
    x = parsed$x_mat,
    nlambda = tuning_searchsize_lambda,
    family = "binomial",
    type.measure = type.measure,
    option.selection = option.selection,
    option.empirical_cutoff = option.empirical_cutoff,
    alphalist = alphalist,
    min.lambda.ratio = min_lambda_ratio,
    type.logistic = ifelse(speedup, "modified.Newton", "Newton"),
    debug = debug)

  # ---
  # Add in a reference glm to help interpretation of final model coefficients.
  reference_data <- data.frame(y = parsed$y,
                               as.data.frame(
                                 predict(fit$final.model$PPx, parsed$x_mat)))
  # copy in outcome names from LHS of formula.
  names(reference_data)[1] <- as.character(f[[2]])

  #   Wrap the reference glm in suppressWarnings as they
  #     almost always fail to converge for ML problems.
  suppressWarnings(
    reference_glm <- glm(formula = parsed$f0, # use the 'flattened' formula.
                         family = "binomial",
                         data = reference_data)
  )

  # copy glm results into final model fit:
  fit$final.model$coef <-
    data.frame(FinalModel = as.numeric(fit$final.model$coef),
               UnregularisedReference = coef(reference_glm))
  rownames(fit$final.model$coef) <- names(coef(reference_glm))

  # also add in reference univarate fits:
  fit$final.model$coef$UnivariateReference <- get_uni_fits(f, data)
  # ---

  # All model fitting completed:
  time.stop <- Sys.time()

  # Make a log object.
  # Times:
  fit$log$start.time <- time.start
  fit$log$stop.time <- time.stop

  # Arguments:
  fit$log$f <- f
  fit$log$data <- df # after removal of non-complete cases.
  fit$log$nrep_outer <- nrep_outer
  fit$log$k_outer <- k_outer
  fit$log$nrep_inner <- nrep_inner
  fit$log$k_inner <- k_inner
  fit$log$type.measure <- type.measure
  fit$log$option.selection <- option.selection
  fit$log$option.empirical_cutoff <- option.empirical_cutoff
  fit$log$tuning_searchsize_lambda <- tuning_searchsize_lambda
  fit$log$tuning_searchsize_alpha <- tuning_searchsize_alpha
  fit$log$alphalist <- alphalist
  fit$log$speedup <- speedup

  # The call:
  fit$log$cmd <- thecall

  # add AUC data to classification results:
  auc.results <- auc.dCVreg(fit)
  fit$rep.class.results <- rbind(fit$rep.class.results,
                                 auc.results$outer_reps)

  fit$fold.class.results <- rbind(fit$fold.class.results,
                                 auc.results$outer_folds)

  fit$final.model$AUC <- auc.results$final_all

  return(fit)
}


# Todo list ---------------------------------------------------------------

# Version 0.6.1
#   - averaging of the ROC curves
#   - develop a plot of per-subject test-group membership vs. performance.
#   - separate documentation from script.
#   - tests/evaluation of data where p>>n.
#   - write up and submit bug to glmnet package authors.
#       It turns out this has been done:
#         https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R
#   - add report of performance within the hyperparameter optimisation loops,
#       arguably this allows diagnosis of overfitting.
#   - currently the preprocessing (scaling) is performed at the level of the
#       outer loop before entering the inner loop.
#       Thus there is some theoretical leakage between the test/train data
#       for the inner loop.
#       This may detriment hyperparameter selection.
#       However, any fix will result in a additional processing overhead as
#       quadratically more calls to caret::preProcess will be required.
#       Investigate this and potentially change the implmentation.
#   - make a dCVnet S3/S4 class (which?), add print, predict, summary methods
#       and convert existing log, plot functions to methods.
#   - rename helper functions from helper() to _helper().
#   - convert to package.

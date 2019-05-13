
# Utility Functions -------------------------------------------------------

#' parse_dCVnet_input
#'
#' Collate a formula and dataset into a standardised object ready for dCVnet
#'
#' The outcome is coerced to a binary factor which is releveled (if necessary)
#' such that the test=positive level comes first.
#' @name parse_dCVnet_input
#'
#' @param data a data.frame containing all terms in f.
#' @param y the outcome (can be numeric vector,
#'      a factor (for binomial / multinomial) or a matrix for cox/mgaussian)
#' @param f a one sided formula.
#'     The RHS must refer to columns in \code{data} and may include
#'     interactions, transformations or expansions (like \code{\link{poly}}, or
#'     \code{\link{log}}).
#'     The formula *must* include an intercept.
#' @param family the model family (see \code{\link{glmnet}})
#' @param offset optional model offset (see \code{\link{glmnet}})
#' @param positive What level of the outcome is a 'positive' result
#'     (in the sense of a diagnostic test).
#'     Can be a numeric indicating which of \code{levels(y)} to use
#'     (i.e. 1 | 2). Alternatively a character specifying the exact level
#'     (e.g. \code{"patient"}).
#' @param yname an optional label for the outcome / y variable.
#' @param passNA should NA values be excluded (FALSE) or passed through (TRUE)?
#'
#' @return \itemize{
#'     \item{ \code{y} - outcome factor
#'         with level ordered according to \code{positive}}
#'     \item{ \code{x_mat} - predictor matrix
#'         including expansions, interaction terms specified in \code{f}}
#'     }
#'
#' @export
parse_dCVnet_input <- function(data,
                               y,
                               family,
                               f = "~.", # nolint
                               positive = 1,
                               offset = NULL,
                               yname = "y",
                               passNA = FALSE) {
  # Check input:
  f <- as.formula(f)
  data <- as.data.frame(data)
  fterms <- stats::terms(f, data = data)

  if ( !identical(attr(fterms, "intercept"), 1L) ) {
    stop("Error: formula must have an intercept. See stats::terms") # nolint
  }
  if ( !identical(attr(fterms, "response"), 0L) ) {
    stop("Error: use a RHS formula to specify in data")
  }

  data <- as.data.frame(data)
  vars <- attr(fterms, "term.labels")
  data <- data[, vars, drop = FALSE]

  # remove missing based on data (unless imputing):
  if ( !passNA && any(!stats::complete.cases(data)) ) {
    cat(paste0("Removing ", sum(!stats::complete.cases(data)),
               " of ", NROW(data),
               " subjects due to missing data.\n"))
    y <- subset(y, stats::complete.cases(data))
    data <- subset(data, stats::complete.cases(data))

  }
  # always remove missing based on y:
  if ( any(!stats::complete.cases(y)) ) {
    cat(paste0("Removing ", sum(!stats::complete.cases(y)),
               " of ", NROW(y),
               " subjects due to missing y.\n"))
    data <- subset(data, stats::complete.cases(y))
    y <- subset(y, stats::complete.cases(y))
  }

  # special treatment if the outcome should be factor:
  if ( family %in% c("binomial", "multinomial")) {
    y <- as.factor(y) # extract y variable & coerce to factor.

    # relevel s.t. 'positive' is first:
    lvl <- levels(y)
    if ( is.numeric(positive) ) {
      positive <- lvl[positive]
    } else {
      positive <- as.character(positive)
    }
    if ( lvl[1] != positive ) y <- relevel(y, positive)
  }

  # Make a model matrix of RHS variables
  #   i.e. parse dummy coding / interaction terms & drop intercept:
  xf <- model.frame(formula = f, data = data, na.action = "na.pass")
  x_mat <- model.matrix(f, data = xf)[, -1]

  # return the outcome, predictor matrix and flattened formula.
  return(list(y = y,
              x_mat = x_mat,
              yname = yname,
              family = family))
}


#' parseddata_summary
#'
#' Simple descriptives for dCVnet's input (parsed dataset).
#'
#' This acts on a model matrix, so factor variables will be indicator coded.
#'     nnz should be informative for such variables.
#'
#' @name parseddata_summary
#' @param object either a \code{\link{dCVnet}} object, or a parsed dataset
#'      output from \code{\link{parse_dCVnet_input}}
#'
#' @return a list containing: \itemize{
#'     \item{ \code{OutcomeData} - frequencies, proportions and level names
#'         for the outcome factor}
#'     \item{ \code{PredictorData} - descriptives for the predictor matrix}
#'     }
#'
#' @export
parseddata_summary <- function(object) {
  # we want to operate either directly on a 'parsed' input,
  #   or extract one from a dCVnet object.
  if ( "dCVnet" %in% class(object) ) {
    object <- parse_dCVnet_input(f = object$input$callenv$f,
                                 y = object$input$callenv$y,
                                 data = object$input$callenv$data,
                                 family = object$input$callenv$family,
                                 positive = object$input$callenv$positive)
  }
  # First describe the target:
  if ( object$family %in% c("binomial", "multinomial") ) {
    ytab <- table(object$y)
    yptab <- round(prop.table(ytab), 3) * 100
    stry <- paste0(names(ytab),
                   ": ",
                   sprintf(ytab, fmt = "%i"),
                   " (", yptab, "%)")
  } else {
    if ( object$family == "cox") {
      stry <- aggregate(list(Time = object$y[, 1]),
                        by = list(Outcome = object$y[, 2]),
                        summary)
    } else {
      stry <- summary(object$y)
    }
  }

  # Next the predictor matrix:
  xdes <- lapply(as.data.frame(object$x_mat),
                 function(x) {
                   data.frame(mean = mean(x, na.rm = TRUE),
                              sd = sd(x, na.rm = TRUE),
                              min = min(x, na.rm = TRUE),
                              max = max(x, na.rm = TRUE),
                              skew = e1071::skewness(x, na.rm = TRUE, type = 2),
                              kurt = e1071::kurtosis(x, na.rm = TRUE, type = 2),
                              nnz = sum(x != 0))
                 })
  xdes <- do.call(rbind, xdes)
  return(list(OutcomeData = stry,
              PredictorData = xdes))
}

#' parse_alphalist
#'
#' Check and process a numeric of alpha values.
#'
#' \itemize{
#' \item{Zero alphas (i.e. pure L2 / ridge penalty) do not currently work
#'     due to a bug in glmnet. Zeros are replaced with a small non-zero value.
#'     For a 'pure' L2/ridge regression use a different package.}
#' \item{it is often the case that small alphas (<0.1) are slow to fit.}
#' }
#'
#' @name parse_alphalist
#' @param alphalist a numeric of alpha values.
#'
#' @return non-duplicated, non-missing alpha values between (0, 1] -
#'     i.e. excluding zero.
#'
#' @export
parse_alphalist <- function(alphalist) {
  # no duplicates:
  alphalist <- unique(alphalist)

  # check no missing values & within [0,1]
  if ( any(is.na(alphalist)) ) stop("Error: missing value(s) in alphalist")
  if ( min(alphalist) < 0.0 ) stop("Error: alphas must be positive")
  if ( max(alphalist) > 1L ) stop("Error: alphas must be <= 1")

  # substitute zero-length alphas due to bug.:
  if ( any(alphalist == 0.0) ) {
    replacement_alpha <- min(c(0.01, 0.5 * min(alphalist[alphalist > 0.0])))
    alphalist[alphalist == 0.0] <- replacement_alpha
    warning(paste0("BUGFIX: alpha=0.0 replaced with: ", replacement_alpha))
  }
  return(alphalist)
}

checkForDuplicateCVFolds <- function(folds) {
  cat("Checking for duplicate folds...\n")
  R <- length(folds) == length(unique(folds))
  if ( !R ) {
    ndup <- sum(duplicated(folds))
    warning(paste("Outer Loop CrossValidation - there were",
                  ndup,
                  "duplicated folds.\n",
                  "Consider reducing number of folds and/or",
                  "turning off stratification"))
    #TODO: could implement a while loop to run until a non-dup set of the
    #       desired size is recovered.
  }
  invisible(R)
}


cvlambdafinder <- function(lambda, cvm, cvsd,
                           minimise = TRUE,
                           type = c("minimum", "se", "percentage"),
                           type.value = 1) {
  # Adapted from glmnet.
  #   'minimum' will just return lambda.min.
  #   'se' will return lambda + type.value*SE.
  #   'percentage' will return lambda * type.value.
  #     e.g. 'percentage' & type.value = 1.03 gives lambda + 3%.
  #     e.g. 'se' & type.value = 1.0 gives the 'standard' lambda+1se.
  #
  type <- match.arg(type)
  if ( !minimise ) cvm <- cvm * -1

  cvmin <- min(cvm, na.rm = TRUE)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin], na.rm = TRUE)

  if ( type == "minimum" ) {
    return(list(lambda.min = lambda.min))
  } else {
    if ( type == "se" ) {
      idmin <- match(lambda.min, lambda)
      semin <- (cvm + type.value * cvsd)[idmin]
      idmin <- cvm <= semin
      lambda.1se <- max(lambda[idmin], na.rm = TRUE)
      return(list(lambda.min = lambda.1se))
    } else {
      # type must be percentage:
      return(list(lambda.min = lambda.min * type.value))
    }
  }
}

#' lambda_rangefinder
#'
#' What range of lambda penalty amounts should we consider for a dataset?
#'     The maximum for a given dataset is determined by a simple formula,
#'     but this can vary from subset to subset. This code bootstraps
#'     the calculation of max-lambda for subsets of the input data.
#'
#' @name lambda_rangefinder
#' @param y binary outcome (from \code{\link{parse_dCVnet_input}})
#' @param x predictor matrix (from \code{\link{parse_dCVnet_input}})
#' @param alphalist alphas to consider, e.g. \code{c(0.2, 0.5, 0.8)}
#' @param prop proportion of the input to bootstrap. This should be the
#'     fraction of observations remaining in the inner loop CV.
#'     for example: if outer and inner loops are both 5-fold,
#'     prop should be 0.8 * 0.8 = 0.64
#' @param niter the number of times to bootstrap.
#'
#' @return \code{maxlambda} - the largest bootstrapped lambda value
#'     for each alpha requested.
#'
#' @export
lambda_rangefinder <- function(y, x,
                               alphalist,
                               prop,
                               niter = 10000) {
  # Bruteforce selection of an overall max lambda per alpha
  #   based on random sampling.
  # prop should be the inner loop sample size for fitting.
  #   e.g. if 5-fold / 5-fold then prop = (1 * 0.8 * 0.8) = 0.64
  #       if 10-fold / 5 fold then prop = (1 * 0.9 * 0.8) = 0.81
  y <- as.numeric(y) - 1
  # Utility function:
  .get_maxlambda <- function(x, y, alphalist) {
    # https://stats.stackexchange.com/questions/166630
    # y must be 0/1 coded.
    # x must be a matrix.
    # If a zero alpha is fed in maxlambda is infinity - because of maths!
    n <- length(y) # number of subjects.

    # Alternative standardisation for the data:
    .sdn <- function(kk) {
      sqrt( sum( (kk - mean(kk)) ^ 2 ) / length(kk) )
    }
    x_sds <- apply(x, 2, .sdn)

    # Scale x and y:
    x <- scale(x, scale = x_sds)
    y <- y - mean(y) * (1 - mean(y))

    # calculate the maximum lambda.
    max(abs(t(y) %*% x )) / (alphalist * n)
  }

  subsize <- round(length(y) * prop)

  result <- parallel::mclapply(1:niter,
                               mc.cores = getOption("mc.cores", 1L),
                               function(i) {
                                 subsamp <- sample(seq_along(y), size = subsize)
                                 .get_maxlambda(x = x[subsamp, ],
                                                y = y[subsamp],
                                                alphalist = alphalist)
                               })
  result <- do.call(rbind, result)

  if ( length(alphalist) == 1 ) {
    return(max(result))
  } else {
    result <- setNames(data.frame(result), as.character(alphalist))
    return(apply(result, 2, max))
  }
}

#' lambda_seq_list
#'
#' produce descending lambda sequences given a maximum, minimum fraction and
#'     the number of steps
#'
#' @name lambda_seq_list
#'
#' @param maxlambdas max lambdas vector,
#'     output by \code{\link{lambda_rangefinder}}
#' @param nlambda how many steps to take
#' @param lambda_min_ratio what fraction of the maximum to take as the minimum.
#'     no default is set, glmnet recommends
#'     \code{ifelse(nobs<nvars,0.01,0.0001)}.
#'
#' @export
lambda_seq_list <- function(maxlambdas, nlambda, lambda_min_ratio) {
  # Utility function, given one or more max_lambdas
  #     (possibily for a list of alphas)
  # as well as a number of outputs and a minimum ratio
  #   returns the lambda sequence.
  lapply(maxlambdas, function(ml) {
    exp(seq(from = log(ml),
            to =   log(ml * lambda_min_ratio),
            length.out = nlambda))
  })
}

# Print some run info to screen.
startup_message <- function(k_inner, nrep_inner,
                            k_outer, nrep_outer,
                            nalpha, nlambda,
                            parsed, time.start,
                            family) {
  nit_inner <- k_inner * nrep_inner * nalpha
  nit_outer <- k_outer * nrep_outer
  nit_total <- nit_inner * nit_outer

  cat(paste0("-- dCVnet --\n"))
  cat(paste0(time.start, "\n"))
  cat(paste0("Model family: ", family, "\n"))

  cat("Features:\n")
  cat(paste0("\t", ncol(parsed$x_mat), "\n"))

  cat("Observations:\n")
  cat(paste0("\t", nrow(parsed$x_mat), " subjects\n"))

  cat("Outcome:\n")
  if ( family %in% c("binomial", "multinomial")) {
    stab <- table(parsed$y)
    vapply(X = seq_along(stab),
           FUN = function(i) {
             cat(paste0("\t", stab[i], " of outcome: ", names(stab)[i], "\n"))
           },
           FUN.VALUE = c(""))
  } else {
    print(summary(parsed$y)); cat("\n")
  }

  cat("Tuning:\n")
  cat(paste0("\t", nalpha, " alpha values\n"))
  cat(paste0("\t", nlambda, " lambda values\n"))

  cat("Models:\n")
  cat(paste0("\t", nit_inner, " inner models\n"))
  cat(paste0("\t", nit_outer, " outer models\n"))
  cat(paste0("\t", nit_total, " models in total\n"))

  cat("------------\n\n")
}


#' cv.glmnet.modelsummary
#'
#' return a dataframe of cv.glmnet results.
#'
#' @param mod a cv.glmnet object
#' @param alpha an optional label
#' @param rep an optional label
#'
#' @name cv.glmnet.modelsummary
#'
#' @export
cv.glmnet.modelsummary <- function(mod,
                                   alpha=NA,
                                   rep=NA) {
  return(data.frame(lambda = mod$lambda,
                    cvm = mod$cvm,
                    cvsd = mod$cvsd,
                    cvup = mod$cvup,
                    cvlo = mod$cvlo,
                    nzero = mod$nzero,
                    lambda.min = mod$lambda == mod$lambda.min,
                    lambda.1se = mod$lambda == mod$lambda.1se,
                    alpha = alpha,
                    rep = rep,
                    stringsAsFactors = FALSE))
}

#' tidy_predict.glmnet
#'
#' return a dataframe of glmnet predictions associated with outcomes (when these
#'     are provided.)
#'
#' @param mod a fitted glmnet object (alpha is determined by the object)
#' @param newx new values of x for which predictions are desired.
#' @param s the value of penalty parameter lambda at which predictions are
#'     required.
#' @param family the glmnet model family
#'     (this determines the format of the return)
#' @param label an optional label (value is added in column "label")
#' @param newy outcome associated with newx. If provided these will be included
#'     in the output (useful for subsequent performance assessment).
#' @param newoffset if an offset is used in the fit, then one must be supplied
#'     for making predictions.
#' @param binomial_thresh this allows non-default thresholds to
#'     be used for classification. This is only relevant for binary
#'     classification. E.g. for an imbalanced binary outcome
#'     with 70:30 allocation, setting the decision threshold to 0.7
#'     gives a better balance of sensitivity and specificity
#'     without requiring threshold tuning (as in AUC optimal threshold).
#' @param ... passed to \code{\link{predict.glmnet}}
#'
#' @return a \code{\link{data.frame}} containing column(s) for:
#'     \itemize{
#'     \item{'prediction' - result of
#'         \code{predict.glmnet(.., type = "response")}.
#'         The interpretation depends on the model family.
#'         }
#'     \item{'rowid' - the rownames of newx}
#'     \item{'label' - a column of a single label used when merging predictions}
#'     }
#'     Optionally the data.frame will contain:
#'     \itemize{
#'     \item{'classification' - the predicted class (for nominal outcomes).}
#'     \item{'reference' - the response being predicted (if newy specified).}
#'     }
#'
#' @name tidy_predict.glmnet
#'
#' @export
tidy_predict.glmnet <- function(mod,
                                newx,
                                s,
                                family,
                                newy = NULL,
                                newoffset = NULL,
                                label = "",
                                binomial_thresh = 0.5,
                                ...) {
  # always specify a value for lambda (s)
  # if rownames were used in fitting the model they will be carried through
  p <- predict(object = mod,
               newx = newx,
               type = "response",
               s = s,
               exact = FALSE,
               newoffset = newoffset,
               ...)
  if ( class(p) == "array" ) {
    p <- p[,,1] # nolint # extra dims not needed.
  }

  p <- as.data.frame(p)

  class(p) <- c("dcvntidyperf", "data.frame")
  attr(p, "family") <- family

  if ( !is.null(newy)) {
    a <- as.data.frame(newy)
  }

  # different rules for column names of p:
  if ( family %in% c("mgaussian", "multinomial") ) {
    colnames(p) <- paste0("prediction", colnames(p))
  } else {
    colnames(p) <- "prediction"
  }

  p$rowid <- rownames(p)

  if ( family %in% c("gaussian", "poisson") ) {
    if ( !is.null(newy) ) p$reference <- a[[1]]
    p$label <- label
    return(p)
  }

  if ( family == "cox" ) {
    if ( !is.null(newy) ) p$reference.Time <- a[, 1]
    if ( !is.null(newy) ) p$reference.Status <- a[, 2]
    p$label <- label
    return(p)
  }

  if ( family == "binomial" ) {
    lvl <- mod$classnames

    p$classification <- factor(p$prediction > binomial_thresh,
                               levels = c(FALSE, TRUE),
                               labels = lvl)

    if ( !is.null(newy) ) {
      p$reference <- factor(newy, levels = lvl)
    }

    p$label <- label
    return(p)
  }

  if ( family == "mgaussian" ) {
    if ( !is.null(newy) ) {
      colnames(a) <- paste0("reference", colnames(a))
      p <- data.frame(p, a)
    }
    p$label <- label
    return(p)
  }

  if ( family == "multinomial" ) {
    p$classification <- c(predict(object = mod,
                                  newx = newx,
                                  type = "class",
                                  s = s,
                                  exact = FALSE,
                                  newoffset = newoffset))
    if ( !is.null(newy) ) p$reference <- a[[1]]
    p$label <- label
    return(p)
  }
  stop(paste0("family error: ", family))
}


#' tidy_coef.glmnet
#'
#' return a dataframe of glmnet predictions associated with outcomes (when these
#'     are provided.)
#'
#' @inheritParams tidy_predict.glmnet
#'
#' @name tidy_coef.glmnet
#' @export
tidy_coef.glmnet <- function(mod,
                             newx,
                             newy,
                             s,
                             family,
                             label = "",
                             newoffset = NULL,
                             binomial_thresh = 0.5) {
  # always specify a value for lambda (s)
  # if rownames were used in fitting the model they will be carried through
  p <- predict(object = mod,
               newx = newx,
               type = "coefficients",
               s = s,
               exact = FALSE,
               newoffset = newoffset)

  if ( is.null(dim(p)) ) {
    nm <- names(p)
    p <- mapply(function(x, n) {
      x <- as.matrix(x)
      rownames(x) <- paste0(rownames(x), "_", n)
      return(x)
    },
    x = p, n = nm, SIMPLIFY = FALSE)
    p <- do.call(rbind, p)
  } else {
    p <- as.matrix(p)
  }
  colnames(p) <- label
  return(p)
}



#' tidy_confusionmatrix
#'
#' return contents of a \code{\link[caret]{confusionMatrix}} as a
#'    'tidy' one column data.frame.
#'
#' @param mat output from \code{confusionMatrix.default}.
#'
#' @name tidy_confusionmatrix
#' @return a one column data.frame
#' @export
tidy_confusionmatrix <- function(mat) {
  tab <- as.data.frame(mat$table)
  tablabs <- paste("Predicted", tab[, 1], "Actual", tab[, 2])
  tab <- data.frame(Measure = tablabs, Value = tab[, 3])

  overall <- data.frame(Measure = names(mat$overall),
                        Value = mat$overall,
                        stringsAsFactors = FALSE)

  byclass <- data.frame(Measure = names(mat$byClass),
                        Value = mat$byClass,
                        stringsAsFactors = FALSE)

  tab <- rbind(tab, overall)
  tab <- rbind(tab, byclass)
  rownames(tab) <- NULL
  return(tab)
}

#' predict_cat.glm
#'
#' Category predictions from a binomial glm object.
#'
#' Note: this is not implemented as an S3 generic despite the name.
#'     It has no dispatch.
#'
#' @param glm a binomial family glm object
#' @param threshold the prediction threshold
#' @export
predict_cat.glm <- function(glm, threshold = 0.5) {
  if ( !("glm" %in% class(glm)) ) stop("input must be of class 'glm'")
  if ( !("binomial" %in% glm$family$family) ) stop("input glm must be binomial")

  # extract the outcome variable.
  outcome <- model.frame(glm)[[1]]

  # Return categorical predictions from a glm given a threshold (default = 0.5):
  lvl <- levels(outcome)

  if (is.null(lvl)) stop("outcome is missing levels")

  R <- as.numeric(stats::fitted(glm) > threshold) + 1
  R <- lvl[R]
  return(factor(R, levels = lvl))
}

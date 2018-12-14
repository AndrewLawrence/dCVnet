
# Utility Functions -------------------------------------------------------

#' parse_dCVnet_input
#'
#' Collate a formula and dataset into a standardised object ready for dCVnet
#'
#' The outcome is coerced to a binary factor which is releveled (if necessary)
#' such that the test=positive level comes first.
#' @name parse_dCVnet_input
#' @param f a two sided formula.
#'     The LHS must be a single binary variable.
#'     The RHS must refer to columns in \code{data}, but may include
#'     interactions, transformations or expansions (like \code{\link{poly}}, or
#'     \code{\link{log}}).
#'     The formula *must* include an intercept.
#' @param data a data.frame containing all terms in f.
#' @param positive What level of the outcome is a 'positive' result
#'     (in the sense of a diagnostic test).
#'     Can be a numeric indicating which of \code{levels(y)} to use
#'     (i.e. 1 | 2). Alternatively a character specifying the exact level
#'     (e.g. \code{"patient"}).
#'
#' @return \itemize{
#'     \item{ \code{y} - outcome factor
#'         with level ordered according to \code{positive}}
#'     \item{ \code{x_mat} - predictor matrix
#'         including expansions, interaction terms specified in \code{f}}
#'     }
#'
#' @export
parse_dCVnet_input <- function(f, data, positive = 1) {
  # Ordering levels of y for classification:
  # We will follow a convention that the condition we are aiming to
  #   detect/diagnose ('i.e. the disease') should be the first
  #   level of the factor.
  # If this is not the case for the input then 'positive'
  #   can be used to specify.

  # Check input:
  f <- as.formula(f)
  if ( !identical(attr(terms(f, data = data), "intercept"), 1L) ) {
    stop("Error: formula must have an intercept. See terms(f)")
  }

  if ( !"data.frame" %in% class(data) ) {
    stop("Error: data must be data.frame.")
  }

  df <- model.frame(f, data = data,
                    drop.unused.levels = T) # does not include interaction.
  if (identical(nrow(df), 0L)) {
    stop("Error: no complete cases found.")
  }
  # TODO:: missing data imputation?
  y <- as.factor(df[, 1]) # extract y variable & coerce to factor.

  if ( length(levels(y)) != 2 ) stop("LHS (i.e. outcome) must have 2 levels")

  # Make a model matrix of RHS variables
  #   i.e. parse dummy coding / interaction terms & drop intercept:
  x_mat <- model.matrix(f, data = data)[, -1]

  # Recode levels.
  lvl <- levels(y)
  if ( is.numeric(positive) ) {
    positive <- lvl[positive]
  } else {
    positive <- as.character(positive)
  }
  if ( lvl[1] != positive ) y <- relevel(y, positive)

  # return the outcome, predictor matrix and flattened formula.
  return(list(y = y,
              x_mat = x_mat))
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


#' lambda_rangefinder
#'
#' What range of lambda penalty amounts should we consider for a dataset?
#'     The maximum for a given dataset is detemined by a simple formula,
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
    # https://stats.stackexchange.com/questions/166630/glmnet-compute-maximal-lambda-value
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
                                 subsamp <- sample(1:length(y), size = subsize)
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


startup_message <- function(k_inner, nrep_inner,
                            k_outer, nrep_outer,
                            nalpha, nlambda,
                            parsed, time.start) {
  # Print some general run info to screen.
  nit_inner <- k_inner * nrep_inner * nalpha
  nit_outer <- k_outer * nrep_outer
  nit_total <- nit_inner * nit_outer

  stab <- table(parsed$y)

  cat(paste0("-- dCVnet --\n"))
  cat(paste0(time.start, "\n"))
  cat("Features:\n")
  cat(paste0("\t", ncol(parsed$x_mat), "\n"))
  cat("Observations:\n")
  cat(paste0("\t", nrow(parsed$x_mat), " subjects\n"))
  cat(paste0("\t", stab[1], " of outcome: ", names(stab)[1], "\n"))
  cat(paste0("\t", stab[2], " of outcome: ", names(stab)[2], "\n"))
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
                    stringsAsFactors = F))
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
                        stringsAsFactors = F)

  byclass <- data.frame(Measure = names(mat$byClass),
                        Value = mat$byClass,
                        stringsAsFactors = F)

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
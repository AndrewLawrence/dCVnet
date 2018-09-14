
# Utility Functions -------------------------------------------------------


parse_input <- function(f, data, positive = 1) {
  # Ordering levels of y for classification:
  # We should follow convention that the condition we are aiming to
  #   detect/diagnose ('i.e. the disease') should be the first
  #   level of the factor.
  # If this is not the case then 'positive' can be used to specify.
  #   if positive is numeric it will be used as the position of the
  #     level of y to use as positive.
  #   if positive is character then this level will be used as positive.

  # Check input 1: formula is formula and has intercept.
  f <- as.formula(f) # attempt to coerce in case of character string.
  if ( !identical(attr(terms(f, data = df), "intercept"), 1L) ) {
    stop("Error: formula must have an intercept. See terms(f)")
  }

  # check input 2:
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

  if ( length(levels(y)) != 2 ) { stop("LHS (y) must have 2 levels") }

  # Make a model matrix of RHS variables
  #   i.e. parse dummy coding / interaction terms & drop intercept:
  x_mat <- model.matrix(f, data = data)[, -1]

  # produce a flattened formula based on x_mat:
  f0.string <- paste(f[[2]],
                     "~",
                     paste(colnames(x_mat),
                           collapse = "+"))
  f0 <- as.formula(f0.string)

  # Recode levels.
  lvl <- levels(y)
  if ( is.numeric(positive) ) {
    positive = lvl[positive]
  } else {
    positive = as.character(positive)
  }
  if ( lvl[1] != positive ) { y <- relevel(y, positive) }

  # return the outcome, predictor matrix and flattened formula.
  return(list(y = y,
              x_mat = x_mat,
              f0 = f0))
}


parse_alphalist <- function(alphalist) {
  # Check alpha values are OK:
  if ( any(is.na(alphalist)) ) { stop("Error: missing value(s) in alphalist") }
  if ( min(alphalist) < 0.0 ) { stop("Error: alphas must be positive") }
  if ( max(alphalist) > 1L ) { stop("Error: alphas must be <= 1") }

  # substitute zero-length alphas due to bug.:
  if ( any(alphalist == 0.0) ) {
    replacement_alpha <- min(c(0.01, 0.5 * min(alphalist[alphalist > 0.0])))
    alphalist[alphalist == 0.0] <- replacement_alpha
    warning(paste0("BUGFIX: alpha=0.0 replaced with: ", replacement_alpha))
  }
  # Note: it is often the case that small alphas (<0.1)
  #       are also slower to fit.
  return(alphalist)
}


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

  subsize <- round(length(y)*prop)

  result <- sapply(1:niter, function(i) {
    subsamp <- sample(1:length(y), size = subsize)
    .get_maxlambda(x = x[subsamp,], y = y[subsamp], alphalist = alphalist)
  })


  if ( length(alphalist) == 1 ) {
    return(max(result))
  } else {
    result <- setNames(data.frame(t(result)), as.character(alphalist))
    return(apply(result, 2, max))
  }
}


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


cv.glmnet.modelsummary <- function(mod,      # a cv.glmnet model
                                   alpha=NA, # optional: add a label with alpha
                                   rep=NA) { # optional: add a label with rep
  return(data.frame(lambda = mod$lambda,
                    cvm = mod$cvm,
                    cvsd = mod$cvsd,
                    cvup = mod$cvup,
                    cvlo = mod$cvlo,
                    nzero = mod$nzero,
                    lambda.min = mod$lambda == mod$lambda.min,
                    lambda.1se = mod$lambda == mod$lambda.1se,
                    alpha = alpha,
                    rep = rep))
}


tidy_confusionmatrix <- function(mat) {
  tab <- as.data.frame(mat$table)
  tablabs <- paste("Predicted", tab[,1], "Actual", tab[,2])
  tab <- data.frame(Measure = tablabs, Value = tab[,3])

  overall <- data.frame(Measure = names(mat$overall),
                        Value = mat$overall)

  byClass <- data.frame(Measure = names(mat$byClass),
                        Value = mat$byClass)

  tab <- rbind(tab, overall)
  tab <- rbind(tab, byClass)
  rownames(tab) <- NULL
  return(tab)
}


predict_cat.glm <- function(glm, threshold = 0.5) {
  # Return categorical predictions from a glm given a threshold (default = 0.5):
  lvl <- levels(glm$data[,1])
  R <- as.numeric(fitted(glm) > threshold) + 1
  R <- lvl[R]
  return(factor(R, levels = lvl))
}


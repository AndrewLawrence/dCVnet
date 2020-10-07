
# Utility Functions -------------------------------------------------------

check_categorical_outcome <- function(cat) {

  .is_alphabetical <- function(factor) {
    lvl <- levels(factor)
    identical(lvl, sort(lvl))
  }

  if ( is.numeric(cat) ) return(cat)
  if ( is.factor(cat) ) {
    if ( ! .is_alphabetical(cat) ) {
      stop("The order of factor levels is not alphabetical.
            glmnet ignores levels and treats factor levels in y
            alphabetically.")
    }
    return(cat)
  }
  if ( is.character(cat) ) return(factor(cat, levels = sort(unique(cat))))
  # at this point we have already returned for:
  #     is.numeric, is.factor and is.character.
  # failure mode: co-erce factor
  return(as.factor(cat))
}


#' parse_dCVnet_input
#'
#' Collate an outcome (y) predictor matrix (x) into a standardised object ready
#'    for dCVnet functions. Optionally x can be a dataframe and a
#'    one-sided formula (f) can be provided to allow interactions,
#'    transformations and expansions using R formula notation.
#'
#' @section Factor Outcomes:
#'    For categorical families (binomial, multinomial) input can be:
#'    \itemize{
#'        \item{numeric (integer): c(0,1,2)}
#'        \item{factor: factor(1:3, labels = c("A", "B", "C")))}
#'        \item{character: c("A", "B", "C")}
#'        \item{other}
#'    }
#'    These are treated differently.
#'
#'    Numeric data is used as provided.
#'    Character data will be coerced to a factor:
#'        \code{factor(x, levels = sort(unique(x)))}.
#'    Factor data will be used as provided, but *must* have levels in
#'    alphabetical order.
#'
#'    In all cases *the reference category must be ordered first*,
#'    this means for the binomial family the 'positive' category is second.
#'
#'    Why alphabetical? Previously bugs arose due to different handling
#'    of factor levels between functions called by dCVnet. These appear to be
#'    resolved in the latest versions of the packages, but this restriction will
#'    stay until I can verify.
#'
#' @section Notes:
#'    Sparse matrices are not supported by dCVnet.
#'
#' @name parse_dCVnet_input
#'
#' @param data a data.frame containing variables needed for the formula (f).
#' @param y the outcome (can be numeric vector,
#'      a factor (for binomial / multinomial) or a matrix for cox/mgaussian)
#'      For factors see Factor Outcomes section below.
#' @param f a one sided formula.
#'     The RHS must refer to columns in \code{data} and may include
#'     interactions, transformations or expansions (like \code{\link{poly}}, or
#'     \code{\link{log}}).
#'     The formula *must* include an intercept.
#' @param family the model family (see \code{\link{glmnet}})
#' @param offset optional model offset (see \code{\link{glmnet}})
#' @param yname an optional label for the outcome / y variable.
#' @param passNA should NA values in data be excluded (FALSE)
#'                   or passed through (TRUE)?
#'
#' @return a list containing
#'     \itemize{
#'     \item{ \code{y} - outcome}
#'     \item{ \code{x_mat} - predictor matrix
#'         including expansions, interaction terms specified in \code{f}}
#'     \item{ \code{yname} - a variable name for the y-variable }
#'     \item{ \code{family} - the model family }
#'     }
#'
#' @export
parse_dCVnet_input <- function(data,
                               y,
                               family,
                               f = "~.", # nolint
                               offset = NULL,
                               yname = "y",
                               passNA = FALSE) {
  f <- as.formula(f)
  data <- as.data.frame(data)
  fterms <- stats::terms(f, data = data)
  vars <- all.vars(fterms)

  # check formula doesn't exclude the intercept:
  if ( !identical(attr(fterms, "intercept"), 1L) ) {
    stop("Error: formula must have an intercept. See stats::terms") # nolint
  }
  # check formula is has no response variable:
  if ( !identical(attr(fterms, "response"), 0L) ) {
    stop("Error: use a RHS formula to specify in data")
  }

  data <- data[, vars, drop = FALSE]

  # Paired removal of incomplete data from x/y
  #     missing y is never passed,
  #     missing x (data) is optionally passed according to passNA
  ycomplete <- stats::complete.cases(y)
  dcomplete <- stats::complete.cases(data)
  # only worry about missing data if y is not missing:
  dmiss <- any( (!dcomplete) & ycomplete )

  # special behaviour for missing in data not y if passNA:
  if ( passNA && dmiss ) {
    dcomplete <- rep(TRUE, length(dcomplete))
    dmiss <- FALSE
    warning("Passing NAs - missing data (non-outcome)")
  }

  # joint removal:
  complete <- ycomplete & dcomplete
  miss <- any(!complete)

  if ( miss ) {
    y <- subset(y, complete)
    data <- subset(data, complete)
    warning(paste0("Removing ", sum(!complete),
               " of ", length(complete),
               " subjects due to missing data.\n"))
  }

  # coerce y into factor and check, because...
  #   - glmnet (for its own infernal reasons) does not honour factor ordering,
  #     the factor levels are converted to character and sorted alphabetically
  #   - previously in dCVnet there was an argument ('positive') to change factor
  #     ordering from the call to dCVnet. This was not honoured by glmnet.
  #   - if y is numeric it is assumed the user knows what they are doing.
  #   NOTE: now not certain if this was my bug rather than glmnet.
  if ( family %in% c("binomial", "multinomial") ) {
    y <- check_categorical_outcome(y)
  }

  # Make a model matrix of RHS variables
  #   i.e. parse dummy coding / interaction terms & drop intercept:
  mf <- stats::model.frame(formula = f,
                           data = data,
                           na.action = ifelse(passNA,
                                              yes = stats::na.pass,
                                              no = stats::na.omit))
  x_mat <- model.matrix(f, data = mf)[, -1]

  # force matrices to vectors:
  if ( family %in% c("gaussian", "poisson") && inherits(y, "matrix") ) {
    y <- as.vector(y)
  }

  # check nothing was mangled:
  if ( NROW(y) != NROW(x_mat) ) stop("Error x & y lengths do not match.")

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
  if ( inherits(object, "dCVnet") ) {
    object <- parse_dCVnet_input(f = object$input$callenv$f,
                                 y = object$input$callenv$y,
                                 data = object$input$callenv$data,
                                 family = object$input$callenv$family)
  }
  # First describe the target:
  if ( object$family %in% c("binomial", "multinomial") ) {
    ytab <- table(object$y)
    yptab <- round(prop.table(ytab), 3) * 100
    stry <- paste0(names(ytab),
                   ": ",
                   sprintf(ytab, fmt = "%i"),
                   " (", yptab, "%)")
  } else if ( object$family == "cox") {
    stry <- aggregate(list(Time = object$y[, 1]),
                      by = list(Outcome = object$y[, 2]),
                      summary)
  } else {
    # should be gaussian (1d mat / vector),
    #           poisson (1d mat / vector) or
    #           mgaussian (data.frame)
    stry <- summary(object$y)
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
                              nnz = sum(x != 0),
                              stringsAsFactors = FALSE)
                 })
  xdes <- as.data.frame(data.table::rbindlist(xdes), stringsAsFactors = FALSE)
  rownames(xdes) <- colnames(object$x_mat)
  return(list(OutcomeData = stry,
              PredictorData = xdes))
}

#' parse_alphalist
#'
#' Check / standardise a numeric vector of alpha values for a multi-alpha
#'     glmnet model.
#'
#' \itemize{
#' \item{Very small alphas are typically slow to fit.}
#' \item{At last time of testing alpha=0.0 (i.e. a pure L2 / ridge penalty)
#'     did not work due to a bug related to using fixed folds in glmnet.
#'     As a workaround Zeros are replaced with a small non-zero value.
#'     In effect different software is required to get results for a
#'     'pure' L2/ridge regression model.}
#' }
#'
#' @name parse_alphalist
#' @param alphalist a numeric of alpha values.
#' @param stripNA missing values either throw an error or are discarded.
#' @param dedupe duplicate alpha values either throw an error or are discarded.
#'
#' @return a named numeric vector of, de-duplicated,
#'     alpha values between (0, 1] - i.e. excluding zero.
#'
#' @export
parse_alphalist <- function(alphalist,
                            stripNA = FALSE,
                            dedupe = FALSE) {
  # coerce numeric:
  alphalist <- base::unname(as.numeric(alphalist))
  # Missing values:
  if ( stripNA ) {
    alphalist <- alphalist[is.na(alphalist)]
  } else {
    if ( any(is.na(alphalist)) ) stop("Error: missing value(s) in alphalist")
  }
  # duplicates:
  if ( dedupe ) {
    alphalist <- unique(alphalist)
  } else {
    if ( length(alphalist) != length(unique(alphalist)) ) {
      stop("Error: alphalist has duplicates")
    }
  }

  # check within [0,1]
  if ( min(alphalist) < 0.0 ) stop("Error: alphas must be positive")
  if ( max(alphalist) > 1L ) stop("Error: alphas must be <= 1")

  # substitute zero-length alphas due to bug.:
  if ( any(alphalist == 0.0) ) {
    replacement_alpha <- min(c(0.01, 0.5 * min(alphalist[alphalist > 0.0])))
    alphalist[alphalist == 0.0] <- replacement_alpha
    warning(paste0("BUGFIX: alpha=0.0 replaced with: ", replacement_alpha))
  }
  # add sensible names:
  nm <- prettyNum(alphalist)
  # fallback mode - longer names:
  if ( max(table(nm)) != 1 ) nm <- as.character(nm)
  return(stats::setNames(alphalist, nm))
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
  result <- as.data.frame(data.table::rbindlist(result),
                          stringsAsFactors = FALSE)

  if ( length(alphalist) == 1 ) {
    return(max(result))
  } else {
    result <- setNames(data.frame(result, stringsAsFactors = FALSE),
                       as.character(alphalist))
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
    for ( i in seq_along(stab) ) {
      cat(paste0("\t", stab[i], " of outcome: ", names(stab)[i], "\n"))
    }
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
  return(data.frame(s = names(mod$nzero),
                    lambda = mod$lambda,
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


#' extract_glmnet_alpha
#'
#' glmnet model objects do not explicitly include the alpha value which was
#' used when the model was fit. This utility function tries to extract the
#' value of alpha used to fit a glmnet object - this is not always
#' possible (e.g. if alpha was set by a named variable no-longer in the
#' environment)
#'
#' Note: dCVnet will quietly extend the glmnet models it fits to include
#'    additional information. This information includes the alpha value
#'    (see \code{\link{set_glmnet_alpha}}). This dCVnet specified
#'    alpha value will be extracted with priority over alpha values implied
#'    in the model without checking.
#'    The dCVnet extension of the (cv.)glmnet S3 class is very basic.
#'
#' @param mod a \code{\link[glmnet]{glmnet}} or \code{\link[glmnet]{cv.glmnet}}
#'     model
#' @export
#' @seealso \code{\link{set_glmnet_alpha}}
extract_glmnet_alpha <- function(mod) {
  if ( "cv.glmnet" %in% class(mod) ) mod <- mod$glmnet.fit
  # If alpha has been set by dCVnet then use this value:
  if ( ! is.null(mod$alpha) ) return(mod$alpha)
  # Otherwise we must use information provided by glmnet.
  # If there is no alpha in the call, then return the fxn default:
  if ( is.null(mod$call$alpha) ) return( base::formals(glmnet::glmnet)$alpha )
  alpha <- mod$call$alpha
  # If alpha evaluates to a valid numeric then return this:
  if ( is.numeric(alpha) && identical(length(alpha), 1L) ) return(alpha)
  # otherwise look for the named object in the environment and if found
  #   warn and return what it points to:
  if ( is.name(alpha) ) {
    warning("alpha in model call is a named object,
              an object of this name was found in the environment,
              but it is not guaranteed to be the correct alpha.")
    try( return( get(as.character(alpha)) ) )
  }
  # as a fallback, stop and print name of missing object:
  stop(
    paste0("alpha in model call is a named object
              which could not be found in the environment.
             object name was: ", alpha))
  return(NULL)
}




#' set_glmnet_alpha
#'
#' glmnet model objects do not explicitly include the alpha value which was
#' used when the model was fit. This utility function extends the glmnet S3
#' class by adding a slot for 'alpha'. This can be manually specified (setalpha)
#' or extracted from the model object (\code{\link{extract_glmnet_alpha}}).
#'
#' Note: This extension of the (cv.)glmnet S3 class is very basic.
#'
#' @inheritParams extract_glmnet_alpha
#' @param setalpha Specify an alpha value
#'
#' @export
#' @seealso \code{\link{extract_glmnet_alpha}}, \code{\link[glmnet]{glmnet}}
set_glmnet_alpha <- function(mod, setalpha = NULL) {
  if ( "cv.glmnet" %in% class(mod) ) {
    if ( missing(setalpha) ) setalpha <- extract_glmnet_alpha(mod$glmnet.fit)
    mod$glmnet.fit$alpha <- setalpha
  } else {
    if ( missing(setalpha) ) setalpha <- extract_glmnet_alpha(mod)
    mod$alpha <- setalpha
  }
  return(mod)
}


#' amalgamate_cv.glmnet
#'
#' Gathers results from a list of \code{\link[glmnet]{cv.glmnet}} objects
#'     and returns a merged, averaged object.
#'
#' The arithmetic mean k-fold cross-validated loss (i.e. type.measure) is taken
#' over the models (with the sd averaged via variance).
#' The cv SE upper and lower limits (used in lambda.1se calculation) are then
#' calculated from on the averaged data and finally the cv optimal lambda.1se
#' and lambda.min values calculated for the averaged performance.
#'
#' Consistent with cv.glmnet, the model coefficients within folds are not
#' made available, averaged or otherwise investigable, but a whole data model
#' is returned in the \code{glmnet.fit} slot.
#'
#' The cvglmlist must contain cv.glmnet models suitable for averaging together.
#'     This typically means all models having the same:
#'     \itemize{
#'     \item{family}
#'     \item{x and y data}
#'     \item{alpha value}
#'     \item{lambda sequence}
#'     \item{type.measure}
#'     \item{number of k-fold CV folds}
#'     \item{other cv.glmnet options}
#'     }
#'     in order for the amalgamated results to "make sense".
#'     Essentially the models in the list should only differ on the random
#'     allocation of folds to cases (usually specified in foldid).
#'
#' Some limited checks are implemented to ensure alpha, lambda and type.measure
#'     are identical. There is an option to turn these checks off, but this is
#'     not recommended.
#'
#' This function presently does not honour the "keep" argument of cv.glmnet and
#'     all additional arrays/vectors are silently dropped.
#'
#' @param cvglmlist a list of cv.glmnet models
#' @param checks which input checks to run
#' @inherit glmnet::cv.glmnet return
#' @examples
#' \dontrun{
#' data("CoxExample", package = "glmnet") # x and y
#' # folds for unstratified 10x-repeated 5-fold cv:
#' foldlist <- replicate(10,
#' sample(1:5, size = NROW(x), replace = TRUE),
#' simplify = FALSE)
#' names(foldlist) <- paste0("Rep", 1:10) # label the replications.
#' lambdaseq <- glmnet::cv.glmnet(x=x, y=y, family = "cox")$lambda
#' # create a list of models:
#' modellist <- lapply(foldlist, function(ff) {
#' glmnet::cv.glmnet(x = x, y = y, family = "cox", foldid = ff,
#'     lambda = lambdaseq) } )
#'
#' # use amalgamate to average results:
#' mod <- amalgamate_cv.glmnet(modellist)
#'
#' # compare rep-rep performance variability with the average performance:
#' # rep1:
#' glmnet::plot.cv.glmnet(modellist[[1]], main = "rep1")
#' # rep2:
#' glmnet::plot.cv.glmnet(modellist[[2]], main = "rep2")
#' # etc...
#' # mean:
#' glmnet::plot.cv.glmnet(mod, main = "averaged")
#' }
#' @seealso \code{\link[glmnet]{cv.glmnet}}
#' @export
amalgamate_cv.glmnet <- function(cvglmlist,
                                 checks = list(alpha = TRUE,
                                               lambda = FALSE,
                                               type.measure = TRUE)) {
  # apply with base::Reduce to check all element of a list are identical:
  .reducing_identical <- function(x, y) if (identical(x, y)) x else FALSE
  # Sometimes the glmnet path fit stops before all lambdas are tested.
  #  we will filter down to lambdas found in every repetition.
  lambda <- base::Reduce(.reducing_identical, lapply(cvglmlist, "[[", "lambda"))
  # we have common lambdas unless lambda == FALSE, handle the non-common case:
  if ( identical(lambda, FALSE) ) {
    # get commmon lambdas:
    lambda <- base::Reduce(base::intersect, lapply(cvglmlist, "[[", "lambda"))
    if ( checks$lambda ) {
      warning("lambda lists are not identical")
      cat("lambda count per repetition:")
      vapply(cvglmlist, function(x) length(x$lambda), FUN.VALUE = c(1))
    }
    # filter results to common set:
    common <- lapply(cvglmlist, function(x) x[["lambda"]] %in% lambda)
    cvglmlist <- mapply(function(mod, sel) {
      mod$lambda <- mod$lambda[sel]
      mod$cvm <- mod$cvm[sel]
      mod$cvsd <- mod$cvsd[sel]
      mod$cvup <- mod$cvup[sel]
      mod$cvlo <- mod$cvlo[sel]
      mod$nzero <- mod$nzero[sel]
      return(mod) },
      mod = cvglmlist,
      sel = common,
      SIMPLIFY = FALSE)
  }

  # check alphas
  alphas <- vapply(cvglmlist,
                   FUN = extract_glmnet_alpha,
                   FUN.VALUE = numeric(1))
  if (checks$alpha && !identical(length(unique(alphas)), 1L)) {
    stop("alphas must be identical")
  }
  # check names (i.e. type.measure)
  cvname <- base::Reduce(.reducing_identical, lapply(cvglmlist, "[[", "name"))
  if ( checks$type.measure && identical(cvname, FALSE) ) {
    stop("all models must use the same name/type.measure")
  }

  # merge the list of results:
  cvm <- rowMeans(as.data.frame(base::Map("[[", cvglmlist, "cvm"),
                                stringsAsFactors = FALSE))
  cvsd <- as.data.frame(base::Map("[[", cvglmlist, "cvsd"),
                        stringsAsFactors = FALSE)
  # for sd take variance, average and return to sd.
  cvsd <- sqrt(rowMeans(cvsd ^ 2))
  # Note: cv.glmnet (2.0.16) gets nzero from the reference glmnet (prior to
  #       cross-validation), it doesn't measure nzero independently for each
  #       fold.
  #       This makes the below averaging unneccesary when dealing with
  #       a list of calls to cv.glmnet with different folds.
  nz <- rowMeans(as.data.frame(base::Map("[[", cvglmlist, "nzero"),
                               stringsAsFactors = FALSE))
  # format nzero:
  nm <- paste0("s", seq.int(from = 0, length.out = length(lambda)))
  names(nz) <- nm
  # lookup lambda.min/lambda.1se
  lamin <- if (cvname == "AUC") {
    glmnet_getmin(lambda, -cvm, cvsd)
  } else {
    glmnet_getmin(lambda, cvm, cvsd)
  }
  return(structure(c(list(lambda = lambda,
                          cvm = cvm,
                          cvsd = cvsd,
                          cvup = cvm + cvsd,
                          cvlo = cvm - cvsd,
                          nzero = nz,
                          name = cvname,
                          glmnet.fit = cvglmlist[[1]]$glmnet.fit,
                          call = cvglmlist[[1]]$call),
                     as.list(lamin)),
                   nrep = length(cvglmlist),
                   class = "cv.glmnet"))
}


#' tidy_predict.glmnet
#'
#' return a dataframe of glmnet predictions associated with outcomes (when these
#'     are provided). Standardises return over different model families.
#'
#' @param mod a fitted glmnet object (alpha is determined by the object)
#' @param newx new values of x for which predictions are desired.
#' @param s the value of penalty parameter lambda at which predictions are
#'     required.
#' @param family the glmnet model family
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
  # Check mod's class is glmnet:
  if ( !("glmnet" %in% class(mod)) ) {
    msg <- paste0(
      "Input must be a glmnet model! ",
      "(i.e. not cv.glmnet, multialpha.repeated.cv.glmnet)",
      collapse = " "
    )
    stop(msg)
  }

  # always specify a value for lambda (s)
  # if rownames were used in fitting the model they will be carried through
  # see formals(glmnet::predict.glmnet)
  p <- predict(object = mod,
               newx = newx,
               type = "response",
               s = s,
               exact = FALSE,
               newoffset = newoffset,
               ...)
  # remove excess dims, convert to data.frame
  p <- as.data.frame(drop(p), stringAsFactors = FALSE)

  if ( !is.null(newy)) {
    a <- as.data.frame(newy, stringsAsFactors = FALSE)
  }

  cn <- "prediction"
  # different rules for column names of p:
  if ( family %in% c("mgaussian", "multinomial") ) {
    cn <- paste0("prediction", colnames(p))
  }

  colnames(p) <- cn

  p$rowid <- rownames(p)

  if ( family %in% c("gaussian", "poisson") ) {
    if ( !is.null(newy) ) p$reference <- a[[1]]
    p$label <- label
  }

  if ( family == "cox" ) {
    if ( !is.null(newy) ) p$reference.Time <- a[, 1]
    if ( !is.null(newy) ) p$reference.Status <- a[, 2]
    p$label <- label
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
  }

  if ( family == "mgaussian" ) {
    if ( !is.null(newy) ) {
      colnames(a) <- paste0("reference", colnames(a))
      p <- data.frame(p, a, stringsAsFactors = FALSE)
    }
    p$label <- label
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
  }

  class(p) <- c("dcvntidyperf", "data.frame")
  attr(p, "family") <- family
  return(p)
}

#' tidy_predict.multialpha.repeated.cv.glmnet
#'
#' return a dataframe of glmnet predictions associated with outcomes (when these
#'     are provided). Standardises return over different model families.
#'
#' @param alpha specify an alpha, or leave NULL to use the optimal alpha
#'     identified by \code{\link{multialpha.repeated.cv.glmnet}}.
#' @param s specify a lambda, or leave NULL to use the optimal lambda
#'     identified by \code{\link{multialpha.repeated.cv.glmnet}}.
#'
#' @inheritParams tidy_predict.glmnet
#'
#' @export
tidy_predict.multialpha.repeated.cv.glmnet <- function(mod,
                                                       newx,
                                                       s = NULL,
                                                       alpha = NULL,
                                                       newy = NULL,
                                                       newoffset = NULL,
                                                       label = "",
                                                       binomial_thresh = 0.5,
                                                       ...) {
  # multialpha objects have attributes:
  #   family, type.measure, type.lambda, opt.keep_models
  if ( !is_multialpha.repeated.cv.glmnet(mod) ) {
    stop("Input must be a multialpha.repeated.cv.glmnet model")
  }

  if ( attr(mod, "opt.keep_models") == "none" ) {
    stop(paste0("The object ", deparse(substitute(object)),
                " does not include the fitted models required for predict.\n",
                "Rerun multialpha.repeated.cv.glmnet with ",
                "opt.keep_models = best",
                "or opt.keep_models = all",
                "to make predictions"))
  }
  if ( !is.null(alpha) && attr(mod, "opt.keep_models") == "best" ) {
    stop(paste0("The object ", deparse(substitute(object)),
                " does not include all fitted models required for predict.\n",
                "Rerun multialpha.repeated.cv.glmnet with ",
                "opt.keep_models = all",
                "in order to predict at the non-best alpha."))
  }

  if ( is.null(s) ) s <- mod$best$lambda
  if ( is.null(alpha) ) {
    # use best alpha:
    alpha <- mod$best$alpha
    gmod <- mod$models[[mod$bestmodel]]
  } else {
    # look for the specified alpha.
    sel <- alpha %in% mod$alphas
    if ( !any(sel) ) stop("Error requested alpha not found in model.")
    gmod <- mod$models[[which(sel)]]
  }

  tidy_predict.glmnet(gmod,
                      newx = newx,
                      s = s,
                      family = attr(mod, "family"),
                      newy = newy,
                      newoffset = newoffset,
                      label = label,
                      binomial_thresh = binomial_thresh,
                      ...)
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
    p <- as.data.frame(data.table::rbindlist(p), stringsAsFactors = FALSE)
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
  tab <- as.data.frame(mat$table, stringsAsFactors = FALSE)
  tablabs <- paste("Predicted", tab[, 1], "Actual", tab[, 2])
  tab <- data.frame(Measure = tablabs,
                    Value = tab[, 3],
                    stringsAsFactors = FALSE)

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
  if ( !inherits(glm, "glm") ) stop("input must be of class 'glm'")
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


#' cv_performance_glm
#'
#' Cross-validated estimates of model performance by
#' repeated k-fold cross-validation.
#'
#' This function is nothing revolutionary. The idea is to
#' extend \code{\link{boot}{cv.glm}} with an interface that better matches
#' the other functions in this package.
#'
#' The additions are:
#' \itemize{
#' \item{Repeated k-fold rather than single k-fold}
#' \item{Option to provide the fold membership}
#' \item{Default use of stratified sampling by outcome class}
#' \item{Performance assessed with \code{\link{summary.performance}}}
#' }
#'
#' @param y outcome vector (numeric or factor)
#' @param data predictors in a data.frame
#' @param f a formula to apply to x
#' @param return_summary bool. return summarised performance (default), or
#'     \code{\link{performance}} objects for further analysis (set to FALSE)
#' @param ... other arguments
#' @inheritParams multialpha.repeated.cv.glmnet
#' @inheritParams repeated.cv.glmnet
#' @return A list containing the following:
#'     \itemize{
#'     \item{glm.performance - summary(performance(x))
#'         for the uncrossvalidated model}
#'     \item{cv.performance - report_performance_summary(cv.fits)
#'         for the crossvalidated model}
#'     \item{folds - the folds used in cross-validation}
#'     \item{call - the function call}
#'     }
#' @seealso \code{\link[boot]{cv.glm}}, \code{\link{performance}}
#' @export
cv_performance_glm <- function(y,
                               data,
                               f = "~.", # nolint
                               folds = NULL,
                               k = 10,
                               nrep = 2,
                               family = "binomial",
                               opt.ystratify = TRUE,
                               opt.uniquefolds = FALSE,
                               return_summary = TRUE,
                               ...) {
  cl <- as.list(match.call())[-1]

  parsed <- parse_dCVnet_input(data = data, y = y, f = f, family = family)

  x <- data.frame(parsed$x)
  y <- parsed$y

  # observed model:
  m0 <- glm(y ~ .,
            data = data.frame(y = y, x, stringsAsFactors = FALSE),
            family = family)

  # prediction performance:
  p0 <- summary(performance(m0), label = "observed")

  # fold generation
  if ( missing(folds) ) {
    strat_y <- rep(".", times = NROW(x))
    if ( opt.ystratify ) strat_y <- y

    folds <- lapply(seq.int(nrep),
                    function(i) {
                      caret::createFolds(y = strat_y,
                                         k = k,
                                         list = FALSE)
                    })
  } else {
    nrep <- length(folds)
    k <- max(folds[[1]])
  }

  if ( identical(opt.uniquefolds, TRUE) ) checkForDuplicateCVFolds(folds)

  # cross-validation loop:
  pp <- lapply(seq_along(folds), function(i) {
    cat(paste0("rep ", i, " of ", nrep, "\n"))
    rep <- folds[[i]]
    ppp <- lapply(1:max(rep), function(j) {
      xtrain <- x[rep != j, ]
      ytrain <- y[rep != j]
      xtest <- x[rep == j, ]
      ytest <- y[rep == j]

      m <- glm(y ~ .,
               data = data.frame(y = ytrain, xtrain, stringsAsFactors = FALSE),
               family = family)
      p <- performance(m,
                       newdata = data.frame(y = ytest,
                                            xtest,
                                            stringsAsFactors = FALSE))
      return(p)
    } )
    # merge the folds:
    ppp <- structure(as.data.frame(data.table::rbindlist(ppp),
                                   stringsAsFactors = FALSE),
                     class = c("performance", "data.frame"))
    ppp$label <- paste("rep", as.character(i))
    return(ppp)
  } )

  # merge the repetitions:

  pp <- structure(as.data.frame(data.table::rbindlist(pp),
                                stringsAsFactors = FALSE),
                  class = c("performance", "data.frame"))


  ppp <- report_performance_summary(pp)

  if ( return_summary ) {
    return(list(
      glm.performance = p0,
      cv.performance = ppp,
      folds = folds,
      call = cl))
  }

  return(list(
    glm.performance = performance(m0),
    cv.performance = pp,
    folds = folds,
    call = cl
  ))
}

# function was removed from glmnet in 3.0 update.
glmnet_getmin <- function(lambda, cvm, cvsd) {
  cvmin <- min(cvm, na.rm = TRUE)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin], na.rm = TRUE)
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin], na.rm = TRUE)
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

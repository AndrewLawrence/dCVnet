
# Reference Models --------------------------------------------------------


#' refunreg
#'
#' Unregularised prediction models for comparison
#'
#' Calculate uncrossvalidated performance for reference unregularised
#'     models to help interpret dCVnet performance.
#'     given \eqn{n} observations of \eqn{p} predictors, the models calculated
#'     are:
#'     \itemize{
#'     \item{a regression using all variables (if \eqn{n > 5 * p}),
#'         otherwise one using the first \eqn{round(n / 5)}
#'         principal components}
#'     \item{a series of models, one for each column in the
#'         design matrix - i.e. a mass-univariate approach}
#'     }
#'
#' Dev Note: the mass-univariate component has a class ('glmlist')
#'     used in some summary functions, but this is not currently
#'     fully implemented with methods etc.
#'
#' @name refunreg
#'
#' @param object an object to calculate reference logistic regressions for
#' @param ... arguments to pass on
#'
#' @return a list containing: \itemize{
#' \item{\code{glm} - the multiple-predictor model (possibly PCA-reduced)}
#' \item{\code{univariate} - an optional \code{glmlist} of models,
#'                 one for each predictor}
#' }
#'
#' @examples
#' \dontrun{
#' data(QuickStartExample, package = "glmnet")
#' m <- dCVnet(y = QuickStartExample$y,
#'             data = QuickStartExample$x,
#'             family = "gaussian")
#' report_reference_performance_summary(refunreg(m))
#' }
#' @export
refunreg <- function(object, ...) {
  UseMethod("refunreg", object)
}

#' refunreg.default
#'
#' @describeIn refunreg refunreg for \code{\link{dCVnet}} object
#' @export
refunreg.default <- function(object, ...) {
  stop("This function is only implemented for dCVnet class")
}


#' refunreg.dCVnet
#'
#' @describeIn refunreg refunreg for \code{\link{dCVnet}} object
#'
#' @param univariate boolean: also calculate per-variable
#'     models?
#' @param doPCA first run PCA on the features can be "auto" or a boolean.
#'                  \itemize{
#'                  \item{\code{"auto"} determines based on ratio of
#'                             observations to predictors}
#'                  \item{\code{TRUE|FALSE} forces pca/no-pca.}
#'                  }
#' @param ncomp specify how many components for pca (integer).
#'                  \code{"auto"} (use n / 10)
#'
#' @export
refunreg.dCVnet <- function(object,
                            univariate = TRUE,
                            doPCA = "auto",
                            ncomp = "auto",
                            ...) {
  cl <- c(as.list(environment()), list(...))

  f <- family(object)

  switch(
    f,
    binomial = do.call("refunreg_glm", args = cl),
    gaussian = do.call("refunreg_glm", args = cl),
    poisson = do.call("refunreg_glm", args = cl),
    stop(paste0(
      "family ",
      f,
      " is currently unsupported for refunreg"
    ))
  )
}

.estimate_ncomp <- function(y,
                            xmat,
                            familystring,
                            cases_per_predictor = 10) {
  if ( familystring %in% c("binomial", "multinomial")) {
    n <- min(table(y))
  } else {
    n <- NROW(y)
  }
  k <- NCOL(xmat)

  return(list(k = k, estncomp = base::floor(n / cases_per_predictor)))
}

# Use 1 function for all model families handled by glm:
refunreg_glm <- function(object,
                         univariate = TRUE,
                         doPCA = "auto",
                         ncomp = "auto",
                         ...) {

  # extract parsed input
  parsed <- parse_dCVnet_input(f = object$input$callenv$f,
                               y = object$input$callenv$y,
                               data = object$input$callenv$data,
                               family = object$input$callenv$family,
                               passNA = object$input$callenv$opt.use_imputation,
                               offset = object$input$callenv$offset)
  opt.imputation_method <- object$input$callenv$opt.imputation_method

  pp_fn <- preproc_imp_functions(opt.imputation_method = opt.imputation_method)

  prod_PPx <- structure(
    pp_fn$fit(parsed$x_mat),
    fit = pp_fn$fit,
    apply = pp_fn$apply
  )
  parsed$x_mat <- attr(prod_PPx, "apply")(prod_PPx, parsed$x_mat)

  # estimate reasonable dimensionality:
  estncomp <- .estimate_ncomp(y = parsed$y,
                              xmat = parsed$x_mat,
                              familystring = parsed$family)

  # Settings:
  doPCA <- ((doPCA == "auto") && (estncomp$k > estncomp$estncomp)) ||
    identical(doPCA, TRUE)

  if ( identical(ncomp, "auto") ) ncomp <- estncomp$estncomp

  # Is PCA required?:
  if ( identical(doPCA, TRUE) ) {
    #   reduce the number of variables via a (svd)PCA on the data.
    x_pca <- prcomp(parsed$x_mat,
                    rank. = ncomp,
                    retx = TRUE)
    pglm.data <- data.frame(y = parsed$y,
                            as.data.frame.matrix(x_pca$x),
                            stringsAsFactors = FALSE)

    pglm <- glm(y ~ .,
                data = pglm.data,
                family = parsed$family)
  } else {
    pglm.data <- data.frame(y = parsed$y,
                            parsed$x_mat,
                            stringsAsFactors = FALSE)
    pglm <- glm(y ~ .,
                data = pglm.data,
                family = parsed$family)
  }

  # univariate models:
  if ( !univariate ) {
    punivariate <- NA
  } else {
    punivariate <- lapply(seq.int(estncomp$k), function(i) {
      f <- as.formula(paste0("y ~ ", colnames(parsed$x_mat)[i]))
      glm(f, data = data.frame(y = parsed$y,
                               parsed$x_mat,
                               stringsAsFactors = FALSE),
          family = parsed$family)
    } )
    names(punivariate) <- colnames(parsed$x_mat)
    punivariate <- structure(punivariate,
                             class = c("glmlist",
                                       class(punivariate)))
  }
  # Merge:
  rlr <- structure(
    list(
      glm = pglm,
      univariate = punivariate,
      object = object,
      preproc = prod_PPx,
      options = list(
        doPCA = doPCA,
        ncomp = ncomp,
        univariate = univariate
      )
    ),
    class = c("refunreg", "list")
  )
  return(rlr)
}

# Simple print function for \code{\link{refunreg}} objects.
#' @export
print.refunreg <- function(x, ...) {
  cat(paste("\nReference models fit to", deparse(substitute(x$object)), "\n\n"))
  cat(paste("\tncases:", length(x$glm$y), "\n"))
  cat(paste("\tnpreds:", length(x$glm$coefficients) - 1, "\n\n"))
  cat(paste("options:\n"))
  print(x$options)
  cat("\n")
}


#' report_reference_performance_summary
#'
#' Returns a report of performance summary information for a set of dCVnet
#'     reference models.
#'
#' @param refobj a set of reference models provided by
#'     \code{\link{refunreg.dCVnet}}
#'
#' @name report_reference_performance_summary
#' @inherit report_performance_summary return
#' @export
report_reference_performance_summary <- function(refobj) {

  glm <- summary(performance(refobj$glm), "GLM")

  glm <- glm[, -3]
  names(glm)[2] <- "Reference GLM"

  if ( refobj$options$univariate ) {
    univ <- report_performance_summary(refobj$univariate)

    names(univ)[2:5] <- paste("UnivPred", names(univ)[2:5])

    return(data.frame(glm, univ[, -1], stringsAsFactors = FALSE))
  } else {
    return(glm)
  }
}


#' coef_refunreg
#'
#' Returns model coefficients (betas) for glm/univariate reference models.
#'
#' @param refobj a set of reference models provided by
#'     \code{\link{refunreg}}
#'
#' @return a data.frame with columns for:
#'     \itemize{
#'     \item{\code{glm} - unless PCA was performed}
#'     \item{\code{univariate} - unless univariate was not requested}
#'     }
#' @name coef_refunreg
#' @export
coef_refunreg <- function(refobj) {
  a <- !refobj$options$doPCA
  b <- refobj$options$univariate

  if ( a ) {
    G <- coef(refobj$glm)
  }
  if ( b ) {
    U <- c(NA,
           vapply(X = refobj$univariate,
                  FUN = function(x) coef(x)[2],
                  FUN.VALUE = c(1.0)))
    names(U) <- c("(Intercept)", names(refobj$univariate))
  }

  if ( a && b ) return(data.frame(glm = G,
                                  univariate = U,
                                  stringsAsFactors = FALSE))
  if ( a ) return(data.frame(glm = G, stringsAsFactors = FALSE))
  if ( b ) return(data.frame(univariate = U, stringsAsFactors = FALSE))
  return(NA)
}

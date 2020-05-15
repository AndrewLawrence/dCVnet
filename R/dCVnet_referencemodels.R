
# Reference Models --------------------------------------------------------


#' reflogreg
#'
#' Unregularised logistic regression models for comparison
#'
#' Calculate reference logistic regressions to help interpret performance.
#'     models for \eqn{n} observations of \eqn{p} predictors calculated are:
#'     \itemize{
#'     \item{a logistic regression using all variables (if \eqn{n > 5 * p}),
#'         otherwise a logistic regression on the first \eqn{round(n / 5)}
#'         principal components}
#'     \item{a series of logistic regressions, one for each column in the
#'         design matrix - a mass univariate approach}
#'     }
#'
#' The univariate component has a class ('glmlist') used in some summary
#'     functions. This is not currently correctly implemented.
#'
#' @name reflogreg
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
#' @export
reflogreg <- function(object, ...) {
  UseMethod("reflogreg")
}

#' reflogreg.default
#'
#' @describeIn reflogreg reflogreg for \code{\link{dCVnet}} object
#' @export
reflogreg.default <- function(object, ...) {
  stop("This function is only implemented for dCVnet class")
}


#' reflogreg.dCVnet
#'
#' @describeIn reflogreg reflogreg for \code{\link{dCVnet}} object
#'
#' @param univariate calculate per-variable logistic-regression models.
#' @param doPCA first run PCA on the features can be "auto" or a boolean.
#'                  \itemize{
#'                  \item{\code{"auto"} determines based on ratio of
#'                             observations to predictors}
#'                  \item{\code{TRUE|FALSE} forces pca/no-pca.}
#'                  }
#' @param ncomp specify how many components for pca (integer).
#'                  \code{"auto"}
#'
#' @export
reflogreg.dCVnet <- function(object,
                             univariate = TRUE,
                             doPCA = "auto",
                             ncomp = "auto",
                             ...) {
  # extract parsed input
  parsed <- parse_dCVnet_input(f = object$input$callenv$f,
                               y = object$input$callenv$y,
                               data = object$input$callenv$data,
                               family = object$input$callenv$family)

  n <- min(table(parsed$y))
  # above we assume effective N for a logistic regression is n minority cases.
  k <- ncol(parsed$x_mat)
  estncomp <- round(n / 5) # rule of thumb: 5 cases per predictor.

  # Settings:
  doPCA <- ( ( (doPCA == "auto") && (k > estncomp) ) || identical(doPCA, TRUE) )
  if ( identical(ncomp, "auto") ) ncomp <- estncomp

  # Is PCA required?:
  if ( identical(doPCA, TRUE) ) {
    #   reduce the number of variables via a (svd)PCA on the data.
    x_pca <- prcomp(parsed$x_mat,
                    rank. = ncomp,
                    retx = TRUE)
    pglm.data <- data.frame(y = parsed$y,
                            as.data.frame.matrix(x_pca$x),
                            stringsAsFactors = FALSE)

    pglm <- glm(y ~ ., data = pglm.data, family = "binomial")
  } else {
    pglm.data <- data.frame(y = parsed$y, parsed$x_mat)
    pglm <- glm(y ~ ., data = pglm.data, family = "binomial")
  }

  # univariate models:
  if ( !univariate ) {
    punivariate <- NA
  } else {
    punivariate <- lapply(1:k, function(i) {
      f <- as.formula(paste0("y ~ ", colnames(parsed$x_mat)[i]))
      glm(f, data = data.frame(y = parsed$y,
                               parsed$x_mat,
                               stringsAsFactors = FALSE),
          family = "binomial")
    } )
    names(punivariate) <- colnames(parsed$x_mat)
    punivariate <- structure(punivariate,
                             class = c("glmlist",
                                       class(punivariate)))
  }
  # Merge:
  rlr <- structure(list(glm = pglm,
                        univariate = punivariate,
                        object = substitute(object),
                        options = list(doPCA = doPCA,
                                       ncomp = ncomp,
                                       univariate = univariate)
  ),
  class = c("reflogreg", "list"))
  return(rlr)
}

# Simple print function for \code{\link{reflogreg}} objects.
#' @export
print.reflogreg <- function(x, ...) {
  cat(paste("\nReference models fit to", x$object, "\n\n"))
  cat(paste("\tncases:", length(x$glm$y), "\n"))
  cat(paste("\tnpreds:", length(x$glm$coefficients) - 1, "\n\n"))
  cat(paste("options:\n"))
  print(x$options)
}


#' report_reference_performance_summary
#'
#' Returns a report of performance summary information for a set of dCVnet
#'     reference models.
#'
#' @param refobj a set of reference models provided by
#'     \code{\link{reflogreg.dCVnet}}
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


#' coef_reflogreg
#'
#' Returns model coefficients (betas) for glm/univariate reference models.
#'
#' @param refobj a set of reference models provided by
#'     \code{\link{reflogreg}}
#'
#' @return a data.frame with columns for:
#'     \itemize{
#'     \item{\code{glm} - unless PCA was performed}
#'     \item{\code{univariate} - unless univariate was not requested}
#'     }
#' @name coef_reflogreg
#' @export
coef_reflogreg <- function(refobj) {
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

  if ( a & b ) return(data.frame(glm = G, univariate = U, stringsAsFactors = FALSE))
  if ( a ) return(data.frame(glm = G, stringsAsFactors = FALSE))
  if ( b ) return(data.frame(univariate = U, stringsAsFactors = FALSE))
  return(NA)
}

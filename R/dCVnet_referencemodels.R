
# Reference Models --------------------------------------------------------


#' dCVnet_refmodels
#'
#' calculate some reference models to help interpret dCVnet performance.
#'     models for \eqn{n} observations of \eqn{p} predictors calculated are:
#'     \itemize{
#'     \item{a logistic regression using all variables (if \eqn{n > 5 * p}),
#'         otherwise a logistic regression on the first \eqn{round(n / 5)}
#'         principal compoments}
#'     \item{a series of logistic regressions, one for each column in the
#'         design matrix - a mass univariate approach}
#'     }
#'
#' The univariate component has a class ('glmlist') used in some summary
#'     functions. This is not currently correctly implemented.
#'
#' @name dCVnet_refmodels
#'
#' @param object a dCVnet object
#'
#' @return a list containing: \itemize{
#' \item{\code{glm} - the multiple-predictor model (possibly PCA-reduced)}
#' \item{\code{univariate} - a \code{glmlist} of models, one for each predictor}
#' }
#'
#' @export
dCVnet_refmodels <- function(object) {
  parsed <- object$input$parsed
  n <- min(table(parsed$y))
  # above we assume effective N for a logistic regression is n minority cases.
  k <- ncol(parsed$x_mat)
  ideal_nvars <- round(n / 5) # ideally we want at least 5 cases per predictor.

  if ( k > ideal_nvars ) {
    # if we have more variables than subjects then first
    #   reduce the number of variables via a (svd)PCA on the data.
    x_pca <- prcomp(parsed$x_mat,
                    rank. = ideal_nvars,
                    retx = T)
    pglm.data <- data.frame(y = parsed$y, as.data.frame.matrix(x_pca$x))

    pglm <- glm(y ~ ., data = pglm.data, family = "binomial")
  } else {
    pglm.data <- data.frame(y = parsed$y, parsed$x_mat)
    pglm <- glm(y ~ ., data = pglm.data, family = "binomial")
  }
  # next univariate prediction:
  punivariate <- lapply(1:k, function(i) {
    f <- as.formula(paste0("y ~ ", colnames(parsed$x_mat)[i]))
    glm(f, data = data.frame(y = parsed$y, parsed$x_mat), family = "binomial")
  } )
  names(punivariate) <- colnames(parsed$x_mat)
  return(list(glm = pglm,
              univariate = structure(punivariate,
                                     class = c("glmlist",
                                               class(punivariate)))))
}


#' report_reference_classperformance_summary
#'
#' Returns a report of classperformance summary information for a set of dCVnet
#'     reference models.
#'
#' @param refobj a set of reference models provided by
#'     \code{\link{dCVnet_refmodels}}
#'
#' @name report_reference_classperformance_summary
#' @inherit report_classperformance_summary return
#' @export
report_reference_classperformance_summary <- function(refobj) {

  glm <- summary(classperformance(refobj$glm), "GLM")

  glm <- glm[, -3]
  names(glm)[2] <- "Reference GLM"

  univ <- report_classperformance_summary(refobj$univariate)

  names(univ)[2:5] <- paste("UnivPred", names(univ)[2:5])

  return(data.frame(glm, univ[, -1]))
}

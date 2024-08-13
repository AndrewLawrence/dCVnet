# documented NULL function for ROxygen and meta-info.


#' dCVnet: Doubly Cross-Validated Elastic-Net Regularised (Generalised)
#'     Linear Models
#'
#' dCVnet fits and cross-validates regression with elastic-net regularisation.
#' Optimal hyperparameters (lambda & alpha) are selected by double (nested)
#' cross-validation with model performance evaluated in an independent outer
#' cross-validation loop.
#'
#' The \code{lambda} and \code{alpha} hyperparameters of the elastic-net
#' allow models to freely
#' range from effectively unregularised, to heavily-regularised (via lambda).
#' With a balance of two types of regularisation: L2 (ridge) and L1 (LASSO)
#' (alpha). This regularisation produces dimensionality reduction and
#' variable selection in the predictors.
#'
#' The values of the lambda and alpha hyperparameters (and so the amount
#' and type of regularisation) are selected on the basis of cross-validated
#' model performance - i.e. a data-driven approach.
#'
#' Using a fully nested cross-validation (instead of reporting the CV
#' performance of the selected hyperparameters) keeps the cross-validated
#' performance estimates 'honest' by reducing optimism bias related to
#' hyperparameter tuning.
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats aggregate as.formula coef glm
#' @importFrom stats model.frame model.matrix family
#' @importFrom stats predict fitted predict.glm
#' @importFrom caret confusionMatrix
#' @importFrom graphics plot
#' @importFrom stats median prcomp relevel setNames
#' @importFrom methods formalArgs
#'
#' @name dCVnet-package
"_PACKAGE"


# model families are being added gradually to this package.
#   supported model families and methods for checking are set here
supported_model_families <- function() {
  c("binomial", "multinomial", "gaussian", "poisson", "mgaussian", "cox")
}

check_model_family <- function(family) {
  if ( ! family %in% supported_model_families() ) {
    cat(paste("Supported families are",
              paste0(supported_model_families(), collapse = " "),
              "\n"))
    stop(paste("family not supported by dCVnet:", family))
  }
}

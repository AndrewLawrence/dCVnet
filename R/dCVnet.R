# documented NULL function for ROxygen.


#' dCVnet: Doubly Cross-Validated Elastic-Net Regularised (Generalised) Linear Models
#'
#' The dCVnet package provides functions to conduct and evaluate classification
#' performance using binary logistic regression with elastic-net regulatisation.
#'
#' The alpha and lambda hyperparameters of the elastic-net permit principled
#' data-driven variable selection and dimensionality reduction in the
#' predictors. While returning a model more readily interpretable than
#' 'black-box' techniques. The nested cross-validation keeps both hyperparameter
#' selection and performance estimates 'honest' by reducing optimism bias.
#'
#' dCVnet adds nested cross-validation functionality to \code{glmnet::cv.lognet}.
#'
#' @docType package
#' @name dCVnet-package
NULL



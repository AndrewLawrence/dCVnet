# documented NULL function for ROxygen.


#' dCVnet: Doubly Cross-Validated Elastic-Net Regularised (Generalised)
#'     Linear Models
#'
#' dCVnet runs binary logistic regression with elastic-net regularisation.
#' Optimal hyperparameters are selected by double (nested) cross-validation with
#' model performance evaluated in an independent outer loop.
#'
#' The \code{lambda} and \code{alpha} hyperparameters of the elastic-net
#' allow models to
#' range from effectively unregularised, to heavily-regularised. Further a
#' mixture of two types of regularisation: L2 (ridge) and L1 (LASSO) are
#' possible. This regularisation produces dimensionality reduction and
#' variable selection in the predictors.
#'
#' The values of the hyperparameters and thereby the amount and type of
#' regularisation is selected on the basis of what
#' hyperparameter combinations produce performance which generalises the best
#' in the inner cross-validation.
#'
#' The fully nested cross-validation keeps the final performance estimates
#' 'honest' by reducing optimism bias.
#'
#' @importFrom glmnet glmnet cv.glmnet predict.lognet
#' @importFrom stats aggregate as.formula coef glm model.frame model.matrix
#' @importFrom stats predict fitted predict.glm
#' @importFrom caret confusionMatrix
#' @importFrom graphics plot
#' @importFrom stats median prcomp relevel setNames
#'
#' @docType package
#' @name dCVnet-package
NULL

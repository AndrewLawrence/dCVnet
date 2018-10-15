# documented NULL function for ROxygen.


#' dCVnet: Doubly Cross-Validated Elastic-Net Regularised (Generalised) Linear Models
#'
#' dCVnet runs binary logistic regression with elastic-net regulatisation.
#' Optimal hyperparameters are selected by double (nested) cross-validation with
#' model performance evaluated in an independent outer loop.
#'
#' The lambda and alpha hyperparameters of the elastic-net permit models ranging
#' from effectively unregularised to heavily-regularised models with a mixture
#' of two types of regularisation. The amount and type of regularisation will be
#' selected on the basis of what generalises best in the inner cross-validation.
#' L2 (ridge) and L1 (LASSO) effectively produce dimensionality reduction and
#' variable selection in the predictors.
#'
#' The fully nested cross-validation keeps the final performance estimates
#' 'honest' by reducing optimism bias.
#'
#' dCVnet is a wrapper around \code{glmnet::cv.lognet}.
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

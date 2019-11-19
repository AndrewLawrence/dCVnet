#' Simulated binary logistic regression data
#'
#' blrsim is a single realisation of data under the following scheme:
#' x contains 20 random, correlated gaussian predictors (see \code{\link[MASS]{mvrnorm}});
#' the binary outcome (y) is set to 1 if the sum of x1:x10 (+ error) is positive, and zero otherwise.
#'
#' The data is split into train and test sets.
#' For reference, the classification accuracy of a model using the generative rule
#' in the complete data is provided.
#'
#' @format list
#' \describe{
#' \item{x}{20 Gaussian predictors. Training data, 50 observations}
#' \item{y}{Class membership for x. Training data, 50 observations}
#' \item{x.test}{20 Gaussian predictors. Test data, 1000 observations}
#' \item{y.test}{Class membership for x. Test data, 1000 observations}
#' \item{bayes.rate}{Classification accuracy in the complete sample given the
#' rule used to generate the data (see: raw-data/blrsim.R)}
#' }
"blrsim"

#' Simulated binary logistic regression data
#'
#' blrsim is a single realisation of data under the following scheme:
#' x contains 40 random, correlated gaussian predictors
#' (see \code{\link[MASS]{mvrnorm}});
#' the binary outcome (y) is set to 1 if the sum of x1:x10 (+ error)
#' is positive, and zero otherwise. In effect variables 11-40 are not-relevant
#' to the outcome.
#'
#' The data is split into train and test sets.
#' For reference, the classification accuracy of a model using the generative rule
#' in the complete data is provided.
#'
#' @format list
#' \describe{
#' \item{x}{40 Gaussian predictors. Training data, 50 observations}
#' \item{y}{Class membership for x. Training data, 50 observations}
#' \item{x.test}{20 Gaussian predictors. Test data, 1000 observations}
#' \item{y.test}{Class membership for x. Test data, 1000 observations}
#' \item{bayes.rate}{Classification accuracy in the complete sample given the
#' rule used to generate the data (see: raw-data/blrsim.R).
#'
#' Note: calling this the bayes rate is a bit of a liberty
#' as bayes rate is an error/loss while classification accuracy
#' is the complement of such a error/loss.}
#' }
"blrsim"
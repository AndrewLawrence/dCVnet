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
#' For reference, the classification accuracy of a model
#' (using the generative rule) in the complete data is provided.
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

#' Prostate cancer data from Stamey et al. (1989)
#'
#' This dataset is used as an example in Hastie, Tibshirani & Freedman's
#' Elements of Statistical Learning. It was included in the ElemStatLearn
#' package which (at time of writing) is orphaned and no-longer available on
#' CRAN.
#'
#' There are 8 predictors (columns 1:8), one outcome (column 9) and a
#' marker for test/train data used in the textbook examples (column 10).
#'
#' Observations are from 97 men who underwent prostatectomy. The original
#' paper investigates the post-surgical characteristics that predict
#' pre-surgical prostate-specific antigen (PSA) score (variable: lpsa).
#'
#' Variables prefixed with 'l' have been log transformed.
#'
#' The following descriptions have been adapted from Ryan Tibshirani's
#' \href{https://www.stat.cmu.edu/~ryantibs/statcomp-F16/lectures/exploratory_data_slides.html#(3)}{lecture notes on EDA} #nolint
#'
#' \describe{
#' \item{lpsa:}{log PSA score}
#' \item{lcavol:}{log cancer volume}
#' \item{lweight:}{log prostate weight}
#' \item{age:}{age of patient}
#' \item{lbph:}{log of the amount of benign prostatic hyperplasia}
#' \item{svi:}{seminal vesicle invasion}
#' \item{lcp:}{log of capsular penetration}
#' \item{gleason:}{Gleason score}
#' \item{pgg45:}{percent of Gleason scores 4 or 5}
#' }
#'
#' The dataset is provided in the original units, a scaled version can be
#' obtained with
#' \code{sprostate <- data.frame(scale(prostate[,-10]), train = prostate[,10])}.
#'
#' Observations are ordered by outcome.
"prostate"
# This is a link for the paper:
# \href{https://www.ncbi.nlm.nih.gov/pubmed/2468795}{Stamey et al (1989)}

#' Used for vignette: dCVnet-limitations
"modsel_perf"

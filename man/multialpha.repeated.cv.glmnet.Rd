% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dCVnet_innerloop.R
\name{multialpha.repeated.cv.glmnet}
\alias{multialpha.repeated.cv.glmnet}
\title{multialpha.repeated.cv.glmnet}
\usage{
multialpha.repeated.cv.glmnet(
  x,
  y,
  alphalist = round(seq(0.2, 1, len = 6)^exp(1), 2),
  lambdas = NULL,
  k = 10,
  nrep = 5,
  opt.lambda.type = c("min", "1se"),
  opt.ystratify = TRUE,
  opt.uniquefolds = FALSE,
  opt.random_seed = NULL,
  family,
  opt.keep_models = c("best", "none", "all"),
  ...
)
}
\arguments{
\item{x}{input matrix, of dimension nobs x nvars; each row is an observation
vector. Can be in sparse matrix format (inherit from class
\code{"sparseMatrix"} as in package \code{Matrix}).
Requirement: \code{nvars >1}; in other words, \code{x} should have 2 or more columns.}

\item{y}{response variable. Quantitative for \code{family="gaussian"}, or
\code{family="poisson"} (non-negative counts). For \code{family="binomial"}
should be either a factor with two levels, or a two-column matrix of counts
or proportions (the second column is treated as the target class; for a
factor, the last level in alphabetical order is the target class). For
\code{family="multinomial"}, can be a \code{nc>=2} level factor, or a matrix
with \code{nc} columns of counts or proportions. For either
\code{"binomial"} or \code{"multinomial"}, if \code{y} is presented as a
vector, it will be coerced into a factor. For \code{family="cox"}, preferably
a \code{Surv} object from the survival package: see Details section for
more information. For \code{family="mgaussian"}, \code{y} is a matrix
of quantitative responses.}

\item{alphalist}{a vector of alpha values to search.}

\item{lambdas}{a list of lambda sequence lists
(corresponding to alphas given in alphalist)}

\item{k}{the number of folds for k-fold cross-validation.}

\item{nrep}{the number of repetitions}

\item{opt.lambda.type}{Method for selecting optimum lambda. One of
\itemize{
\item{\code{"min"} - returns the lambda with best
CV score.}
\item{\code{"1se"} - returns the +1 se lambda}
}}

\item{opt.ystratify}{Boolean.
Outer and inner sampling is stratified by outcome.
This is implemented with \code{\link[caret]{createFolds}}}

\item{opt.uniquefolds}{Boolean.
In most circumstances folds will be unique. This requests
that random folds are checked for uniqueness in inner and outer loops.
Currently it warns if non-unique values are found.}

\item{opt.random_seed}{Interpreted as integer. This is used to control
the generation of random folds.}

\item{family}{Either a character string representing
one of the built-in families, or else a \code{glm()} family object. For more
information, see Details section below or the documentation for response
type (above).}

\item{opt.keep_models}{The models take up memory. What should we return?
\itemize{
\item{ best - model with the alpha value selected as optimal. }
\item{ none - no models, just cv results. }
\item{ all - list of models at all alphas. }
}}

\item{...}{arguments passed to \code{\link[glmnet]{cv.glmnet}}}
}
\value{
an object of class \code{\link{multialpha.repeated.cv.glmnet}}.
Containing: \itemize{
\item{results - merged \code{\link{repeated.cv.glmnet}} with
additional columns indicating \emph{alpha} and logical for \emph{best} overall}
\item{best - best selected row from results}
\item{folds - record of folds used}
\item{models - models requested by opt.keep_models.}
\item{bestmodel - index of the best model such that
\code{models[[bestmodel]]} returns the model selected as optimal.}
}
}
\description{
Runs \code{\link{repeated.cv.glmnet}} for a list of alpha values and
returns averaged results, selects the 'best' alpha.
One key difference between (repeated.)cv.glmnet and this function is
that a single 'best' lambda/alpha combination is identified
based on opt.lambda.type.
\emph{This is intended to be a dCVnet internal function}
}
\seealso{
\code{\link{repeated.cv.glmnet}}
}

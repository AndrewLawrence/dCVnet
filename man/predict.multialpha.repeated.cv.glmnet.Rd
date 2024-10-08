% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dCVnet_innerloop.R
\name{predict.multialpha.repeated.cv.glmnet}
\alias{predict.multialpha.repeated.cv.glmnet}
\title{predict.multialpha.repeated.cv.glmnet}
\usage{
\method{predict}{multialpha.repeated.cv.glmnet}(object, newx, alpha = NULL, s = NULL, ...)
}
\arguments{
\item{object}{a a \code{\link{multialpha.repeated.cv.glmnet}} object.}

\item{newx}{matrix of new values for x at which predictions are required.
Note: no impact when type is "coefficients", "class" or "nonzero".}

\item{alpha}{the penalty type alpha at which prediction is required.
Leave NULL to use the cv-optimised value.}

\item{s}{The penalty amount lambda at which prediction is required.
Leave NULL to use the cv-optimised value.}

\item{...}{passed to \code{\link[glmnet]{predict.glmnet}}}
}
\description{
obtain predictions from a
\code{\link{multialpha.repeated.cv.glmnet}} object.
Uses the 'best' alpha & lambda hyperparameters determined by the
internal cross-validation results. For lambda this will be lambda.min
or lambda.1se (determined at model runtime) unless requested otherwise.
}
\seealso{
\code{\link[glmnet]{predict.cv.glmnet}},
\code{\link[glmnet]{predict.glmnet}}
}

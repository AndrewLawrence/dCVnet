% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dCVnet_main.R
\name{coefficients_summary}
\alias{coefficients_summary}
\title{coefficients_summary}
\usage{
coefficients_summary(object, ...)
}
\arguments{
\item{object}{a dCVnet object}

\item{...}{" "}
}
\value{
A data.frame with a row for each term/predictor in the model.
Columns present coefficients of the Production Model
(\code{ProductionModel}) and descriptives of coefficients from the
outer-loop cross-validation. \code{OuterMedian}: the median over folds
and repetitions of cross-validation.
\code{min}, \code{max}: the smallest, largest values
held by that coefficient over the cross-validation. \code{propnz}: the
proportion of cross-validation folds for which the variable was non-zero
(Note: this will be close to 1 for low alpha = ridge penalty)
}
\description{
Extract the coefficients in the production model (fit to complete data)
and provide descriptives for the coefficients in the outerloop
cross-validation
}
\seealso{
coef.dCVnet
}

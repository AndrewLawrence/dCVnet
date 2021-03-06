% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dCVnet_utilities.R
\name{parse_alphalist}
\alias{parse_alphalist}
\title{parse_alphalist}
\usage{
parse_alphalist(alphalist, stripNA = FALSE, dedupe = FALSE)
}
\arguments{
\item{alphalist}{a numeric of alpha values.}

\item{stripNA}{missing values either throw an error or are discarded.}

\item{dedupe}{duplicate alpha values either throw an error or are discarded.}
}
\value{
a named numeric vector of, de-duplicated,
alpha values between (0, 1] - i.e. excluding zero.
}
\description{
Check / standardise a numeric vector of alpha values for a multi-alpha
glmnet model.
}
\details{
\itemize{
\item{Very small alphas are typically slow to fit.}
\item{At last time of testing alpha=0.0 (i.e. a pure L2 / ridge penalty)
did not work due to a bug related to using fixed folds in glmnet.
As a workaround Zeros are replaced with a small non-zero value.
In effect different software is required to get results for a
'pure' L2/ridge regression model.}
}
}

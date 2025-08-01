% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dCVnet_performance.R
\name{report_performance_summary}
\alias{report_performance_summary}
\title{report_performance_summary}
\usage{
report_performance_summary(
  dCVnet_object,
  short = FALSE,
  somersD = FALSE,
  pvprevalence = "observed"
)
}
\arguments{
\item{dCVnet_object}{result from a call to \code{\link{dCVnet}}}

\item{short}{(bool) return a core set of performance measures.}

\item{somersD}{(bool) calculate the computationally expensive Somers' D for
certain families (gaussian, poisson)}

\item{pvprevalence}{allows calculation of PPV/NPV at different prevalences.
set to "observed" to use the prevalence of the data.
For binomial data use a single value, for multinomial use a
named vector of prevalences with names as per the levels of y.
Note: does not affect the presented prevalence value in the table.}
}
\value{
a data.frame of summarised and raw performance statistics.
}
\description{
extracts performance from a dCVnet object
calculates classification statistics and
provides a summary of
\code{\link[base]{mean}},
\code{\link[stats]{sd}},
\code{\link[base:Extremes]{min}}, and \code{\link[base:Extremes]{max}}
}

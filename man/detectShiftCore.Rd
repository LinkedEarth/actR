% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanShift.R
\name{detectShiftCore}
\alias{detectShiftCore}
\title{Detect a shift in a single time-value pair of paleodata}
\usage{
detectShiftCore(
  time,
  vals,
  minimum.segment.length = 1,
  cpt.fun = changepoint::cpt.mean,
  gaussianize = TRUE,
  calc.deltas = TRUE,
  ...
)
}
\arguments{
\item{time}{a vector of time data}

\item{vals}{a vector paleodata}

\item{minimum.segment.length}{the minimum allowed length of a detected segment (in time units)}

\item{cpt.fun}{which function from the changepoint package to use, changepoint::cpt.mean, changepoint::cpt.var or changepoint::cpt.meanvar}

\item{gaussianize}{Force vals to gaussian distribution before analysis. Default (TRUE). Most (all?) methods in the changepoint package assume gaussian distributions, so this is strongly recommended.}

\item{calc.deltas}{Calculate the difference in means between change point sections. Default FALSE}

\item{...}{options to pass to cpt.fun . See changepoint function documentation for details.}
}
\value{
A tibble of output data and metadata
}
\description{
Detect a shift in a single time-value pair of paleodata
}

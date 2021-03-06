% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_functions.R
\name{pexp_truncated}
\alias{pexp_truncated}
\title{Distribution function for truncated exponential distribution}
\usage{
pexp_truncated(x, rate, upper)
}
\arguments{
\item{x}{Vector of values.}

\item{rate}{Rate parameter.}

\item{upper}{Upper (right-hand) limit of the distribution.}
}
\description{
For each value in the vector \code{x}, returns the probability
for the interval (0, x]. Any values of x less than zero will have
probability zero. Any values of x greater than \code{upper}
will have probability 1.0.
}

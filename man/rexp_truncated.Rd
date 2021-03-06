% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_functions.R
\name{rexp_truncated}
\alias{rexp_truncated}
\title{Draw random values from a truncated exponential distribution}
\usage{
rexp_truncated(n, rate, upper)
}
\arguments{
\item{n}{Number of values.}

\item{rate}{Vector of one or more values for the rate parameter.}

\item{upper}{Upper (right-hand) limit of the distribution.}
}
\value{
A vector of random values.
}
\description{
Draws \code{n} random values from an exponential distribution with
rate \code{rate}. Unlike the unbounded exponential distribution,
the expected mean value will not be \code{1 / rate}. For example,
with \code{rate = 1/20, upper = 10} the expected mean is 4.59
rather than 20. The mean of a truncated exponential distribution is
given by:
\code{(1 / rate) - upper / (exp(rate * upper) - 1)}.
The function \code{mean_exp_truncated} does this calculation.
}

#' Draw random values from a truncated exponential distribution
#'
#' Draws \code{n} random values from an exponential distribution with
#' rate \code{rate}. Unlike the unbounded exponential distribution,
#' the expected mean value will not be \code{1 / rate}. For example,
#' with \code{rate = 1/20, upper = 10} the expected mean is 4.59
#' rather than 20. The mean of a truncated exponential distribution is
#' given by:
#' \code{(1 / rate) - upper / (exp(rate * upper) - 1)}.
#' The function \code{mean_exp_truncated} does this calculation.
#'
#' @param n Number of values.
#' @param rate Vector of one or more values for the rate parameter.
#' @param upper Upper (right-hand) limit of the distribution.
#'
#' @return A vector of random values.
#'
#' @export
#'
rexp_truncated <- function(n, rate, upper) {
  r <- runif(n) * (1 - exp(-upper * rate))
  -log(1 - r) / rate
}

#' Density function for truncated exponential distribution
#'
#' Any values \code{<=0} or \code{>upper} will have density zero.
#'
#' @param x Vector of values.
#' @param rate Rate parameter.
#' @param upper Upper (right-hand) limit of the distribution.
#'
#' @export
#'
dexp_truncated <- function(x, rate, upper) {
  d <- x
  ii <- x > 0 & x <= upper

  d[ii] <- rate * exp(-rate * x[ii]) / (1 - exp(-rate * upper))
  d[!ii] <- 0

  d
}


#' Distribution function for truncated exponential distribution
#'
#' For each value in the vector \code{x}, returns the probability
#' for the interval (0, x]. Any values of x less than zero will have
#' probability zero. Any values of x greater than \code{upper}
#' will have probability 1.0.
#'
#' @param x Vector of values.
#' @param rate Rate parameter.
#' @param upper Upper (right-hand) limit of the distribution.
#'
#' @export
#'
pexp_truncated <- function(x, rate, upper) {
  p <- x
  ii <- x > 0 & x <= upper

  p[ii] <- 1 - (exp(rate * (upper - x[ii])) - 1) / (exp(rate * upper) - 1)
  p[x <= 0] <- 0
  p[x > upper] <- 1

  p
}


#' Mean value of a truncated exponential distribution
#'
#' @param rate Vector of one or more values for the rate parameter.
#' @param upper Upper (right-hand) limit of the distribution.
#'
#' @export
#'
mean_exp_truncated <- function(rate, upper) {
  (1 / rate) - upper / (exp(rate * upper) - 1)
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risk_functions.R
\name{summarize_block_risk}
\alias{summarize_block_risk}
\title{Summarize scan line risk values for landscape blocks}
\usage{
summarize_block_risk(line.risk, blocks, quantiles = c(0.25, 0.75),
  strict.crs = FALSE, quiet = FALSE)
}
\arguments{
\item{line.risk}{Either a \code{risk} data frame as returned by
\code{calculate_risk} or a character string giving the path and filename
of a SQLite database of risk data.}

\item{blocks}{Polygons (or multi-polygons) to profile. These can be provided
as an \code{sf} data frame or as the path to a file of vector data
(e.g. ESRI shapefile).}

\item{quantiles}{A vector of quantiles (probabilities with values in [0,1])
to calculate in addition to the mean. If \code{NULL} or an empty vector
no quantiles are calculated.}

\item{strict.crs}{If TRUE, the blocks and risk objects must have
exactly the same coordinate reference system defined. If FALSE,
the function will make relaxed assumptions about the reference
systems as described for function \code{lines_in_blocks}.}

\item{quiet}{If TRUE, suppress messages and warnings about coordinate
reference systems of objects.}
}
\value{
An \code{sf} object with data from the input \code{risk} object
  plus columns giving the number of scan lines intersecting each block;
  number of locations to which the intersected scan lines belong;
  the mean value of scan line probabilities, and (optionally) specified
  quantiles.
}
\description{
This function takes a set of risk values for scan lines,
as returned by \code{calculate_risk}, and summarizes the probability values
of lines intersecting with landscape blocks. Blocks are polygons representing
management units (e.g. areas for prescribed burning). For each block the
function calculates the mean probability of all intersecting lines and,
(optionally) specified quantiles.
}
\examples{
\dontrun{
risk.lines <- calculate_risk(scanlines, tsf, forest, 100)
blocks <- st_read("c:/foo/bar/data/burning_blocks.shp")
risk.blocks <- summarize_block_risk(risk.lines, blocks, quantiles = c(0.1, 0.9))
}

}
\seealso{
\code{\link{lines_in_blocks}}
}

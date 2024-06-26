% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiationtools.R
\name{mean_slope}
\alias{mean_slope}
\title{Calculates the mean slope to the horizon}
\usage{
mean_slope(dtm, steps = 36, reso = 100)
}
\arguments{
\item{dtm}{a SpatRaster object, two-dimensional array or matrix of elevations (m).}

\item{steps}{an optional integer. The mean slope is calculated from the horizon angle in specified directions. Steps defines the total number of directions used. If the default 36 is specified, then the horizon angle is calculated at 10º intervals.}

\item{reso}{a single numeric value representing the spatial resolution of \code{dtm} (m).}
}
\value{
a SpatRaster object or a two-dimensional array of the mean slope angle to the horizon in all directions (º).
}
\description{
\code{mean_slope} is used to calculates the mean slope to the horizon in all directions.
}
\details{
If \code{dtm} is a SpatRaster object, a SpatRaster object is returned.
The projection system associated with \code{dtm} must be such that
units of x, y and z are identical. Use \code{\link[terra:project]{terra::project()}} to convert the
projection to a Universal Transverse Mercator type projection system.
}
\examples{
library(terra)
ms <- mean_slope(rast(dtm100m), reso = 100)
plot(ms, main = "Mean slope to horizon")
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiationtools.R
\name{skyviewveg}
\alias{skyviewveg}
\title{Calculates a sky view correction factor underneath vegetation}
\usage{
skyviewveg(dtm, l, x, steps = 36, reso = 1)
}
\arguments{
\item{dtm}{a SpatRaster object, two-dimensional array or matrix of elevations (m).}

\item{l}{a SpatRaster object, two-dimensional array or matrix of leaf area index values as returned by \code{\link[=lai]{lai()}}.}

\item{x}{a SpatRaster object, two-dimensional array of numeric values representing the ratio of vertical to horizontal projections of leaf foliage as returned by \code{\link[=leaf_geometry]{leaf_geometry()}}.}

\item{steps}{an optional integer. The sky view is calculated from the horizon angle in specified directions. Steps defines the total number of directions used. If the default 36 is specified, then the horizon angle is calculated at 10º intervals.}

\item{reso}{a single numeric value representing the spatial resolution of \code{dtm} (m).}
}
\value{
a SpatRaster object or a two-dimensional array of numeric values representing the proportion of isotropic radiation received by a surface partially obscured by topography relative to the full hemisphere underneath vegetation.
}
\description{
\code{skyviewveg} is used to calculate a coefficient to correct for the
proportion of sky obscured by topography when calculating net shortwave or
longwave radiation above the canopy.
}
\details{
If \code{dtm} is a SpatRaster object, a SpatRaster object is returned.
The projection system associated with \code{dtm} must be such that
units of x, y and z are identical. Use \code{\link[terra:project]{terra::project()}} to convert the
projection to a Universal Transverse Mercator type projection system.
If a proportion of the sky of partially obscured, then the isotropic
radiation flux received by a surface underneath canopy can be determined by
integrating the single direction radiation transmission over the proportion
of sky in view. This function returns a computationally efficient
approximation of the integrated transmission over the proportion of sky in
view expressed as a proportion of the integrated transmission over the
entire hemisphere.
}
\examples{
library(terra)
l <- lai(aerial_image[,, 3], aerial_image[,, 4])
x <- leaf_geometry(rast(veg_hgt))
sv <- skyviewveg(rast(dtm1m), l, x)
plot(sv, main = "Sky view factor")
}
\seealso{
The function \code{\link[=skyviewtopo]{skyviewtopo()}} calculates a sky view correction factor above vegetation.
}

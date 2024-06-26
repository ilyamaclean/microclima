% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiationtools.R
\name{canopy}
\alias{canopy}
\title{Calculates canopy cover}
\usage{
canopy(l, ref = 0)
}
\arguments{
\item{l}{a SpatRaster object, two-dimensional array or matrix of leaf area index values as returned by \code{\link[=lai]{lai()}}.}

\item{ref}{leaf reflectance.}
}
\value{
a SpatRaster object or a two-dimensional array of numeric values representing fractional canopy cover estimated as the proportion of isotropic radiation transmitted through the canopy.
}
\description{
\code{canopy} is used to calculate fractional canopy cover.
}
\details{
Canopy cover calculated by this function is defined as 1 - the proportion of isotropic
radiation transmitted through the canopy. This is effectively te same as canopy
cover if \code{ref} = 0, but if \code{ref} is greater than zero, scattered light is also considered.
If \code{l} is a raster object, a raster object is returned.
}
\examples{
library(terra)
l <- lai(aerial_image[,,3], aerial_image[,,4])
l <- if_raster(l, dtm1m) # convert to raster
fr <- mask(canopy(l),rast(dtm1m))
plot(fr, main = "Fractional canopy cover")
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solartools.R
\name{latlongfromraster}
\alias{latlongfromraster}
\title{Derives latitude and longitude of centre of SpatRaster object}
\usage{
latlongfromraster(r)
}
\arguments{
\item{r}{a SpatRaster object with the coordinate reference system defined by \code{\link[=crs]{crs()}}}
}
\value{
a data.frame with the latitude and longitude of the centre of the raster
}
\description{
\code{latlongfromraster} is used to calculate the latitude and longitude of
the centre of a raster object.
}
\examples{
latlongfromraster(rast(dtm1m))
latlongfromraster(rast(dtm100m))
}

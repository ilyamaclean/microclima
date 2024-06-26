% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiationtools.R
\name{shortwavetopo}
\alias{shortwavetopo}
\title{Downscales shortwave radiation accounting for topographic effects}
\usage{
shortwavetopo(
  dni,
  dif,
  julian,
  localtime,
  lat = NA,
  long = NA,
  dtm = array(0, dim = c(1, 1)),
  slope = NA,
  aspect = NA,
  svf = 1,
  alb = 0.23,
  albr = 0.23,
  ha = 0,
  reso = 100,
  merid = round(long/15, 0) * 15,
  dst = 0,
  shadow = TRUE,
  component = "sw",
  difani = TRUE
)
}
\arguments{
\item{dni}{a single numeric value, SpatRaster object, two-dimensional array or matrix of coarse-resolution direct radiation perpendicular to the solar beam (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}}).}

\item{dif}{a single numeric value, SpatRaster object, two-dimensional array or matrix of diffuse radiation horizontal ot the surface (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}}).}

\item{julian}{a single integer representing the Julian as returned by \code{\link[=julday]{julday()}}.}

\item{localtime}{a single numeric value representing local time (decimal hour, 24 hour clock).}

\item{lat}{a single numeric value representing the mean latitude of the location for which downscaled radiation is required (decimal degrees, -ve south of equator).}

\item{long}{a single numeric value representing the mean longitude of the location for which downscaled radiation is required (decimal degrees, -ve west of Greenwich meridian).}

\item{dtm}{an optional SpatRaster object, two-dimensional array or matrix of elevations (m), orientated as if derived using \code{\link[=is_raster]{is_raster()}}. I.e. \verb{[1, 1]} is the NW corner.}

\item{slope}{a single value, SpatRaster object, two-dimensional array or matrix of slopes (º). If an array or matrix, then orientated as if derived using \code{\link[=is_raster]{is_raster()}}. I.e. \verb{[1, 1]} is the NW corner.}

\item{aspect}{a single value, SpatRaster object, two-dimensional array or matrix of aspects (º). If an array or matrix, then orientated as if derived using \code{\link[=is_raster]{is_raster()}}. I.e. \verb{[1, 1]} is the NW corner.}

\item{svf}{an optional single value, SpatRaster object, two-dimensional array or matrix of values representing the proportion of isotropic radiation received by a partially obscured surface relative to the full hemisphere as returned by \code{\link[=skyviewtopo]{skyviewtopo()}}.}

\item{alb}{an optional single value, SpatRaster object, two-dimensional array or matrix of surface albedo(s) (range 0 - 1) derived using \code{\link[=albedo]{albedo()}} or \code{\link[=albedo_adjust]{albedo_adjust()}}.}

\item{albr}{an optional single value, SpatRaster object, two-dimensional array or matrix of values of albedo(s) of adjacent surfaces (range 0 - 1) as returned by \code{\link[=albedo_reflected]{albedo_reflected()}}.}

\item{ha}{an optional SpatRaster object, two-dimensional array or matrix of values representing the mean slope to the horizon (decimal degrees) of surrounding surfaces from which radiation is reflected for each cell of \code{dtm} as returned by \code{\link[=mean_slope]{mean_slope()}}.}

\item{reso}{a single numeric value representing the spatial resolution of \code{dtm} (m).}

\item{merid}{an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is \code{round(long / 15, 0) * 15}}

\item{dst}{an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if \code{merid} = 0).}

\item{shadow}{an optional logical value indicating whether topographic shading should be considered (False = No, True = Yes).}

\item{component}{an optional character string of the component of radiation to be returned. One of "sw" (net shortwave radiation, i.e. accounting for albedo), "sw2" (total incoming shortwave radiation), "dir" (direct), "dif" (diffuse), "iso" (isotropic diffuse), "ani" (anistopic diffuse), "ref" (reflected).}

\item{difani}{an optional logical indicating whether to treat a proportion of the diffuse radiation as anistropic (see details).}
}
\value{
If component is "sw", a SpatRaster object or two-dimensional array of numeric values representing net shortwave radiation (MJ m^-2 hr^-1).

If component is "sw2", a SpatRaster object or two-dimensional array of numeric values representing total incoming shortwave radiation (MJ m^-2 hr^-1).

If component is "dir", a SpatRaster object or two-dimensional array of numeric values representing direct shortwave radiation (MJ m^-2 hr^-1).

If component is "dif", a SpatRaster object or two-dimensional array of numeric values representing diffuse shortwave radiation (MJ m^-2 hr^-1).

If component is "iso", a SpatRaster object or two-dimensional array of numeric values representing isotropic diffuse shortwave radiation (MJ m^-2 hr^-1).

If component is unspecified, then the default "sw" is returned.
}
\description{
\code{shortwavetopo} is used to downscale components of the flux density of
shortwave radiation received at the surface of the Earth using a
high-resolution digital elevation dataset, ignoring canopy effects.
}
\details{
If \code{slope} is unspecified, and \code{dtm} is a SpatRaster, \code{slope} and \code{aspect} are calculated from
the raster. If \code{slope} is unspecified, and \code{dtm} is not a SpatRaster, the slope and aspect
are set to zero. If \code{lat} is unspecified, and \code{dtm} is a SpatRaster with a coordinate reference
system defined, \code{lat} and \code{long} are calculated from the SpatRaster. If \code{lat} is unspecified,
and \code{dtm} is not a SpatRaster, or a SpatRaster without a coordinate reference system defined, an
error is returned. If \code{dtm} is specified, then the projection system used must be such that
units of x, y and z are identical. Use \code{\link[terra:project]{terra::project()}} to convert the projection to a
Universal Transverse Mercator type projection system. If \code{dtm} is a SpatRaster object, a raster
object is returned. If \code{dtm} is a SpatRaster object, a SpatRaster object is returned. If \code{dni} or
\code{dif} are SpatRaster objects, two-dimensional arrays or matrices, then it is assumed that they
have been derived from coarse-resolution data by interpolation, and have the same extent as \code{dtm}.
If no value for \code{ha} is provided, the mean slope to the horizon is assumed to be 0. If no value
for \code{svv} is provided, then the entire hemisphere is assumed to be in view. If no value for \code{svf}
is provided, then the entire hemisphere is assumed to be in view. If values of \code{alb} and \code{albr}
are not specified, then a default value of 0.23, typical of well-watered grass is
assumed. If single values of \code{alb} and \code{albr} are given, then the entire
area is assumed to have the same albedo. If \code{dtm} is specified, then the
projection system used must be such that the units of x, y and z are
identical. Use \code{\link[terra:project]{terra::project()}} to convert the projection to a Universal
Transverse Mercator type projection system. If no value for \code{dtm} is
provided, radiation is downscaled by deriving values on the inclined
surfaces specified in \code{slope} and \code{aspect} and topographic shadowing is
ignored. If single values are provided for \code{slope} and \code{aspect} single
values of components of shortwave radiation for an inclined surface are
returned. Only single values of \code{lat} and \code{long} are taken as inputs. Under partially
cloudy conditions, a proportion of diffuse radiation is typically anistropic
(directional). If \code{difani} is TRUE (the default), then the assumption is made that
hourly direct radiation transmission can define the portions of the diffuse
radiation to be treated as anistropic and isotropic.
If \code{difani} is FALSE, all diffuse radiation is treated as isotropic.
If \code{dtm} covers a large extent, the \code{dtm} is best divided into blocks and
seperate calculations performed on each block. Since horizon angles,
topographic shading and sky view correction factors may be influenced by
locations beyond the extent of \code{dtm}, it is best to ensure \code{dtm} covers a
larger extent than that for which radiation values are needed, and to
ensure sub-divided blocks overlap in extent. Calculations are faster if values
for all inputs are provided.
}
\examples{
library(terra)
# =================================
# Extract data for 2010-05-24 11:00
# =================================
dni <-dnirad[,,3444]
dif <-difrad[,,3444]
# ===========================
# Resample to 100m resolution
# ===========================
e<-ext(-5.40, -5.00, 49.90, 50.15)
dnir <- rast(dni)
difr <- rast(dif)
ext(dnir) <- e
ext(difr) <- e
crs(dnir) <- '+init=epsg:4326'
crs(difr) <- '+init=epsg:4326'
dnir <- project(dnir, crs(rast(dtm100m)))
difr <- project(difr, crs(rast(dtm100m)))
dni <- resample(dnir, rast(dtm100m))
dif <- resample(difr, rast(dtm100m))
sv <- skyviewtopo(rast(dtm100m))
jd <- julday(2010, 5, 24)
ha <- mean_slope(rast(dtm100m))
# ================================================================
# Calculate and plot net shortwave radiation for 2010-05-24 11:00
# ================================================================
netshort100m <- shortwavetopo(dni, dif, jd, 11, dtm = rast(dtm100m),
                              svf = sv, ha = ha)
plot(mask(netshort100m, rast(dtm100m)),
     main = "Net shortwave radiation")
}
\seealso{
Function \code{\link[=shortwaveveg]{shortwaveveg()}} returns net shortwave radiation below a canopy.
}

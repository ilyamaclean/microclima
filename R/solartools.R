#' Flexible conversion to raster object
#'
#' @description
#' `if.raster` is used to permit flexibility in the use of rasters, matrices or arrays in many functions.
#'
#' @param x an R object
#' @param r an R object
#'
#' @return if `r` is a raster, `x` is converted to a raster with the same attributes as `r`, otherwise returns `x`
#' @export
#' @import raster
#'
#' @examples
#' r <- is.raster(dtm100m)
#' r1 <- if.raster(r, dtm100m)
#' r2 <- if.raster(r, r)
#' class(r1) # is a RasterLayer
#' class(r2) # is a matrix
if.raster <- function(x, r) {
  if (class(r) == "RasterLayer")
    x <- raster(x, template = r)
  x
}
#' Checks whether object is a raster and returns a matrix if yes

#' @description
#' `is.raster` is used to permit flexibility in the use of rasters, matrices or arrays in many functions.
#'
#' @param r an R object
#'
#' @return if `r` is a raster, returns a matrix containing all values of `r`, otherwise returns `r`
#' @export
#' @import raster
#'
#' @examples
#' r <- is.raster(dtm100m)
#' class(dtm100m) # is a RasterLayer
#' class(r) # is a matrix
#' plot(r) # not a raster
#' plot(raster(r)) # converts to raster
is.raster <- function(r) {
  if (class(r) == "RasterLayer")
    r <- getValues(r, format = "matrix")
  r
}
#' Derives latitude and longitude of centre of raster object
#'
#' @description `latlongfromraster` is used to calculate the latitude and longitude of
#' the centre of a raster object.
#'
#' @param r a raster object with the coordinate reference system defined by [crs()]
#'
#' @return a data.frame with the latitude and longitude of the centre of the raster
#' @export
#' @import raster rgdal
#' @importFrom sp coordinates
#'
#' @examples
#' latlongfromraster(dtm1m)
#' latlongfromraster(dtm100m)
latlongfromraster <- function(r) {
  e <- extent(r)
  xy <- data.frame(x = (e@xmin + e@xmax) / 2,
                   y = (e@ymin + e@ymax) / 2)
  coordinates(xy) = ~x+y
  proj4string(xy) = crs(r)
  ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
  ll <- data.frame(lat = ll$y, long = ll$x)
  ll
}
#' Calculates the astronomical Julian day
#'
#' @description `julian` is used to calculate the astronomical Julian day (days since since January 1, 4713 BCE at noon UTC) from a given year, month and day.
#'
#' @param year year (AD).
#' @param month month in numeric form (1-12).
#' @param day days of the month (1-31).
#' @param hour hours (decimal, 0-23).
#' @param min minutes (decimal, 0-59).
#' @param sec seconds (decimal, 0-59).
#' @param tz an optional numeric value specifying the time zones expressed as hours different from GMT (-ve to west).
#'
#' @return Julian Day. I.e. the number of days since January 1, 4713 BCE at noon UTC.
#' @export
#'
#' @examples
#' jd1 <- julday(2010, 1, 31)
#' jd2 <- julday(2010, 1, 31, 11, 0, 0)
#' jd1 - jd2
julday <- function(year, month, day, hour = 12, min = 0, sec = 0, tz = 0) {
  day_decimal <- day + (hour - tz + (min + sec / 60) / 60) / 24
  monthadj <- month + (month < 3) * 12
  yearadj <- year + (month < 3) * -1
  julian_day <- trunc(365.25 * (yearadj + 4716)) + trunc(30.6001 *
                (monthadj + 1)) + day_decimal - 1524.5
  B <- (2 - trunc(yearadj / 100) + trunc(trunc(yearadj / 100) / 4))
  julian_day <- julian_day + (julian_day > 2299160) * B
  julian_day
}
#' Calculates the solar time
#'
#' @description `solartime` is used to calculate the solar time. I.e. the time that would be measured by a sundial.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param long longitude of the location for which the solar time is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian Julian day expressed as an integer as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for UK time).
#' @param dst an optional numeric value representing the local summer time adjustment (hours, e.g. +1 for BST).
#'
#' @return the solar time. I.e. the times that would be measured by a sundial (hours).
#' @export
#'
#' @details
#' ‘solartime’ accounts for two factors: firstly, east or west component of the analemma, namely
#' the angular offset of the Sun from its mean position on the celestial sphere as viewed from
#' Earth due the eccentricity of the Earth's orbit and the obliquity due to tilt of the Earth's
#' rotational axis. These two factors have different wavelengths, amplitudes and phases,
#' that vary over geological timescales. The equations used here are those derived by Milne.
#'
#' @examples
#' jd <- julday (2010, 6, 21) # Julian day
#' solartime(12, -5, jd) # solartime at noon on 21 June 2010, 5ºW
solartime <- function(localtime, long, julian, merid = 0, dst = 0) {
  m <- 6.24004077 + 0.01720197 * (julian -  2451545)
  eot <- -7.659 * sin(m) + 9.863 * sin (2 * m + 3.5932)
  st <- localtime + (4 * (long + merid) + eot) / 60 - dst
  st
}
#' Calculates the solar azimuth
#'
#' @description `solazi` is used to calculate the solar azimuth at any given location from the local time.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the solar azimuth is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar azimuth is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian Julian day expressed as an integer as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for UK time).
#' @param dst an optional numeric value representing the local summer time adjustment (hours, e.g. +1 for BST).
#'
#' @return a numeric value representing the solar azimuth (decimal degrees).
#' @export
#'
#' @examples
#' # solar azimuth at noon on 21 June 2010, Porthleven, Cornwall, UK
#' jd <- julday (2010, 6, 21) # Julian day
#' solazi(12, 50.08, -5.31, jd)
solazi <- function(localtime, lat, long, julian, merid = 0, dst = 0) {
  stime <- solartime(localtime, long, julian, merid, dst)
  tt <- 0.261799 * (stime - 12)
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((julian - 171) / 365.25))
  Sinh <- sin(declin) * sin(lat * pi / 180) + cos(declin) *
          cos(lat * pi / 180) * cos(tt)
  hh <- (atan(Sinh / sqrt(1 - Sinh * Sinh)))
  sinazi <- cos(declin) * sin(tt) / cos(hh)
  cosazi <- (sin(lat * pi / 180) * cos(declin) * cos(tt) -
            cos(pi * lat / 180) * sin(declin)) / sqrt((cos(declin) *
            sin(tt)) ^ 2 + (sin(pi * lat / 180) * cos(declin) * cos(tt) -
            cos(pi * lat / 180) * sin(declin)) ^ 2)
  sqt <- 1 - sinazi * sinazi
  sqt[sqt < 0] <- 0
  solz <- 180 + (180 * atan(sinazi / sqrt(sqt))) / pi
  if (cosazi < 0) {
    if (sinazi < 0) solz <- 180 - solz else solz <- 540 - solz
  }
  solz
}
#' Calculates the solar altitude
#'
#' @description `solalt` is used to calculate the solar altitude at any given location from the local time.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the solar altitude is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar altitude is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian Julian day expressed as an integer as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for UK time).
#' @param dst an optional numeric value representing the local summer time adjustment (hours, e.g. +1 for BST).
#'
#' @return a numeric value representing the solar altitude (º).
#' @export
#'
#' @examples
#' # solar altitude at noon on 21 June 2010, Porthleven, Cornwall
#' jd <- julday (2010, 6, 21) # Julian day
#' solalt(12, 50.08, -5.31, jd)
solalt <- function(localtime, lat, long, julian, merid = 0, dst = 0) {
  stime <- solartime(localtime, long, julian, merid, dst)
  tt <- 0.261799 * (stime - 12)
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((julian - 171) / 365.25))
  sinh <- sin(declin) * sin(lat * pi / 180) + cos(declin) *
          cos(lat * pi / 180) * cos(tt)
  sa <- (180 * atan(sinh / sqrt(1 - sinh * sinh))) / pi
  sa
}
#' Calculates the tangent of the horizon angle
#'
#' @description `horizonangle` is used to calculate the tangent of the angle to the horizon in a specified direction.
#'
#' @param dtm a raster object, two-dimensional array or matrix of elevations (m). If not a raster, orientated as if derived from a raster using [is.raster()]. I.e. `[1, 1]` is the NW corner.
#' @param azimuth a numeric value representing the direction of the horizon as, for example, returned by [solazi()] (º from north).
#' @param res a single numeric value representing the spatial resolution of `dtm` (m).
#'
#' @return a raster object or two-dimensional array of numeric values representing the tangent of the angle to the horizon in a specified direction.
#' @import raster
#' @export
#'
#' @details
#' To enable calculation of horizon angles near the edge of `dtm` a 100 pixel buffer is of
#' zeros is placed around it. NAs in `dtm` are converted to zeros.
#' The projection system used must be such that units of x, y and z are identical. Use
#' [projectRaster()] to convert the projection to a Universal Transverse Mercator type
#' projection system. If `dtm` is a raster object, a raster object is returned.

#' @examples
#' ha <- horizonangle(dtm1m, 0)
#' plot(ha, main = "Tangent of angle to horizon")
horizonangle <- function(dtm, azimuth, res = 1) {
  r <- dtm
  dtm <- is.raster(r)
  dtm[is.na(dtm)] <- 0
  dtm <- (dtm * 5) / res
  azimuth <- azimuth - 90
  azi <- azimuth * (pi / 180)
  horizon <- array(0, dim(dtm))
  dtm3 <- array(0, dim(dtm) + 200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x + 100), 101:(y + 100)] <- dtm
  for (step in 1:10) {
    horizon[1:x, 1:y] <- pmax(horizon[1:x, 1:y], (dtm3[(101 + sin(azi) *
                         step ^ 2):(x + 100 + sin(azi) * step ^ 2),
                         (101 + cos(azi) * step ^ 2):(y + 100 + cos(azi) *
                         step ^ 2)] - dtm3[101:(x + 100), 101:(y + 100)]) /
                         (5 * step ^ 2))
  }
  horizon <- if.raster(horizon, r)
  horizon
}
#' Calculates the solar index
#'
#' @description `solarindex` is used to calculate the proportion of direct beam radiation incident on an inclined surface at a specified time and location.
#'
#' @param slope a single value, raster object, two-dimensional array or matrix of slopes (º). If an array or matrix, then orientated as if derived using [is.raster()]. I.e. `[1, 1]` is the NW corner.
#' @param aspect a single value, raster object, two-dimensional array or matrix of aspects (º). If an array or matrix, then orientated as if derived using [is.raster()]. I.e. `[1, 1]` is the NW corner.
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock).
#' @param lat a single numeric value representing the mean latitude of the location for which the solar index is required (decimal degrees, -ve south of the equator).
#' @param long a single numeric value representing the mean longitude of the location for which the solar index is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian a single integer representing the Julian day as returned by [julday()].
#' @param dtm an optional raster object, two-dimensional array or matrix of elevations (m). If not a raster, orientated as if derived from a raster using [is.raster()]. I.e. `[1, 1]` is the NW corner.
#' @param res a single numeric value representing the spatial resolution of `dtm` (m).
#' @param merid an optional single numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for UK time).
#' @param dst an optional single numeric value representing the local summer time adjustment (hours, e.g. +1 for BST).
#' @param shadow an optional logical value indicating whether topographic shading should be considered (TRUE = Yes, FALSE = No).
#'
#' @return If shadow is `TRUE`, a raster object or a two-dimensional array of numeric values representing the proportion of direct beam radiation incident on an inclined surface, accounting for topographic shading.
#' @return If shadow is `FALSE`, a raster object or a two-dimensional array of numeric values representing the proportion of direct beam radiation incident on an inclined surface, not accounting for topographic shading.
#' @return If no `dtm` is provided, a vector, array or single numeric value of the proportion of direct beam radiation incident on the inclined surfaces specified by `slope` and `aspect`, and topographic shading is ignored.
#' @import raster
#'
#' @details
#' If `slope` is unspecified, and `dtm` is a raster, `slope` and `aspect` are calculated from
#' the raster. If `slope` is unspecified, and `dtm` is not a raster, the slope and aspect
#' are set to zero. If `lat` is unspecified, and `dtm` is a raster with a coordinate reference
#' system defined, `lat` and `long` are calculated from the raster. If `lat` is unspecified,
#' and `dtm` is not a raster, or a raster without a coordinate reference system defined, an
#' error is returned. If `dtm` is specified, then the projection system used must be such that
#' units of x, y and z are identical. Use [projectRaster()] to convert the projection to a
#' Universal Transverse Mercator type projection system. If `dtm` is a raster object, a raster
#' object is returned.
#'
#' @export
#'
#' @seealso the raster package function [terrain()] can be used to derive slopes and aspects from `dtm` (see example).
#'
#' @examples
#' jd <- julday (2010, 6, 21) # Julian day
#' # slope, aspect, lat & long calculated from raster
#' si1 <- solarindex(localtime = 8, julian = jd, dtm = dtm1m)
#' si2 <- solarindex(localtime = 8, julian = jd, dtm = dtm1m, shadow = FALSE)
#' par(mfrow = c(2, 1))
#' plot(si1, main = "Solar index with topographic shadowing")
#' plot(si2, main = "Solar index without topographic shadowing")
#' ll <- latlongfromraster(dtm1m)
#' solarindex(0, 0, 8, lat = ll$lat, long = ll$long, jd)
solarindex <- function(slope = NA, aspect, localtime, lat = NA, long, julian,
                      dtm = array(0, dim = c(1, 1)), res = 1, merid = 0,
                      dst = 0, shadow = TRUE) {
  r <- dtm
  if (class(slope) == "logical" & class(r) == "RasterLayer") {
    slope <- terrain(r, opt = "slope", unit = "degrees")
    aspect <- terrain(r, opt = "aspect", unit = "degrees")
  }
  if (class(slope) == "logical" & class(r) != "RasterLayer") {
    slope <- 0
    aspect <- 0
  }
  if (class(lat) == "logical" & class(crs(r)) == "CRS") {
    lat <- latlongfromraster(r)$lat
    long <- latlongfromraster(r)$long
  }
  if (class(lat) == "logical" & class(crs(r)) != "CRS")
    stop("Latitude not defined and cannot be determined from raster")
  dtm <- is.raster(dtm)
  slope <- is.raster(slope)
  aspect <- is.raster(aspect)
  saltitude <- solalt(localtime, lat, long, julian, merid, dst)
  alt <- saltitude * (pi / 180)
  zen <- pi / 2 - alt
  sazimuth <- solazi(localtime, lat, long, julian, merid, dst)
  azi <- sazimuth * (pi / 180)
  sl <- slope * (pi / 180)
  asp <- aspect * (pi / 180)
  shadowmask <- array(1, dim(dtm))
  if (shadow) shadowmask[horizonangle(dtm, sazimuth, res) > tan(alt)] <- 0
  index <- array(0, dim(dtm))
  index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
  index[index < 0] <- 0
  index <- index * shadowmask
  index <- if.raster(index, r)
  index
}
#' Calculates the airmass coefficient
#'
#' @description `airmasscoef` is used to calculate, for a given location and time, the direct optical path length of a solar beam through the atmosphere of the Earth, expressed as a ratio relative to the path length vertically upwards.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the airmass coefficient is required (decimal degrees, -ve south of equator).
#' @param long longitude of the location for which the airmass coefficient is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian a numeric value representing the Julian day as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for UK time).
#' @param dst an optional numeric value representing the local summer time adjustment (hours, e.g. +1 for BST).
#'
#' @return the airmass coefficient, i.e. the direct optical path length of a solar beam through the Earth’s atmosphere, expressed as a ratio relative to the path length vertically upwards for a given location and time.
#' @export
#'
#' @examples
#' # airmass coefficient at noon on 21 June 2010, Porthleven, Cornwall
#' jd <- julday (2010, 6, 21) # Julian day
#' airmasscoef(12, 50.08, -5.31, jd)
airmasscoef <- function(localtime, lat, long, julian, merid = 0, dst = 0) {
  sa <- solalt(localtime, lat, long, julian, merid, dst)
  z <- 90 - sa
  thickness <- 1 / (cos(z * pi / 180) + 0.50572 * (96.07995 - z) ^ (-1.6364))
  thickness[sa < -1.757] <- NA
  thickness
}

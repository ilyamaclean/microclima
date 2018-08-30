#' Converts between different measures of humidity
#'
#' @description `humidityconvert` is used to convert between different measures of humidity, namely relative, absolute or specific. Vapour pressure is also returned.
#'
#' @param h humidity value(s). Units as follows: specific humidity (\ifelse{html}{\out{kg kg<sup>-1</sup>}}{\eqn{kg kg^{-1}}}), absolute humidity (\ifelse{html}{\out{kg m<sup>-3</sup> }}{\eqn{kg m^{-3}}}), relative humidity (\%), vapour pressure (kPa).
#' @param intype a character string description of the humidity type of `h`. One of "relative", "absolute" or "specific".
#' @param tc A numeric value specifying the temperature (ºC).
#' @param p An optional numeric value specifying the atmospheric pressure (Pa).
#'
#' @details This function converts between vapour pressure and specific,
#' relative and absolute humidity, based on sea-level pressure and
#' temperature. It returns a list of relative, absolute and specific
#' humidity and vapour pressure. If `intype` is unspecified, then `h` is
#' assumed to be relative humidity. If `p` is unspecified, pressure assumed
#' to be 101300, a typical value for sea-level pressure. If one or more of the
#' relative humidity values exceeds 100\% a warning is given.
#'
#' @return a list of numeric humidity values with the following components:
#' @return `relative` relative humidity (\%).
#' @return `absolute`  absolute humidity (\ifelse{html}{\out{kg m<sup>-3</sup> }}{\eqn{kg m^{-3}}}).
#' @return `specific` specific humidity (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}}).
#' @return `vapour_pressure` vapour pressure (kPa).
#' @export
#'
#' @examples
#' humidityconvert(90, 'relative', 20)
#' humidityconvert(0.01555486, 'absolute', 20)
#' humidityconvert(0.01292172, 'specific', 20)
humidityconvert <- function(h, intype = "relative", tc = 20, p = 101300) {
  tk <- tc + 273.15
  pk <- p / 1000
  if (intype != "specific" & intype != "relative" & intype != "absolute") {
    warning ("No valid input humidity specified. Humidity assumed to be
             relative")
    intype <- "relative"
  }
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ws <- 0.622 * e0 / pk
  if (intype == "specific") {
    e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
    hr <- (h / ws) * 100
  }
  if (intype == "absolute") {
    ea <- (tk * h) / 2.16679
    hr <- (ea / e0) * 100
  }
  if (intype == "relative") hr <- h
  if (max(hr, na.rm = T) > 100) warning(paste("Some relative humidity values > 100%",
                            max(hr, na.rm = T)))
  ea <- e0 * (hr / 100)
  hs <- (hr / 100) * ws
  ha <- 2.16679 * (ea / tk)
  return(list(relative = hr, absolute = ha, specific = hs,
              vapour_pressure = ea))
}
#' Calculates land to sea ratio in upwind direction
#'
#' @description `invls` is used to calculate an inverse distance\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}} weighted ratio of land to sea in a specified upwind direction.
#'
#' @param landsea A raster object with NAs (representing sea) or any non-NA value (representing land). The object should have a larger extent than that for which land-sea ratio values are needed, as the calculation requires land / sea coverage to be assessed upwind outside the target area.
#' @param e an extent object indicating the region for which land-sea ratios are required.
#' @param direction a single numeric value specifying the direction (decimal degrees) from which the wind is blowing.
#' @param maxdist The maximum distance (in the same units as `landsea`) over which land-sea ratios should be calculated.  If the maximum distance is greater than the extent of `landsea`, then the area outside `landsea` is ignored.
#'
#' @details This function calculates a coefficient of the ratio of land to
#' sea pixels in a specified upwind direction, across all elements of a
#' raster object, weighted using an inverse distance squared function,
#' such that nearby pixels have a greater influence on the coefficient.
#' It returns a raster object representing values ranging between zero
#' (all upwind pixels sea) to one (all upwind pixels land).
#'
#' @return a raster object with distance-weighted proportions of upwind land pixels
#' @import raster rgdal
#' @importFrom sp coordinates
#' @export
#'
#' @examples
#' ls1 <- invls(dtm100m, extent(dtm1m), 180)
#' ls2 <- invls(dtm100m, extent(dtm1m), 270)
#' par(mfrow=c(2,1))
#' plot(ls1, main = "Land to sea weighting, southerly wind")
#' plot(ls2, main = "Land to sea weighting, westerly wind")
invls <- function(landsea, e, direction, maxdist = 2e+05) {
  res <- xres(landsea)
  slr <- landsea * 0 + 1
  s <- c((8:1000) / 8) ^ 2 * res
  s <- s[s <= maxdist]
  lss <- crop(slr, e)
  lsm <- getValues(lss, format = "matrix")
  lsw <- array(NA, dim = dim(lsm))
  for (yy in 1:dim(lsm)[1]) {
    for (xx in 1:dim(lsm)[2]) {
      if (is.na(lsm[yy, xx]) == F) {
        x <- xx * res + e@xmin - res / 2
        y <- e@ymax + res / 2 - yy * res
        xdist <- round(s * sin(direction * pi / 180), 0)
        ydist <- round(s * cos(direction * pi / 180), 0)
        xy <- data.frame(x = x + xdist, y = y + ydist)
        coordinates(xy) <- ~x + y
        lsc <- extract(slr, xy)
        lsw[yy, xx] <- sum(lsc, na.rm = T) / length(lsc)
      }
    }
  }
  lsr <- raster(lsw, template = lss)
  lsr
}
#' Calculates coastal effects using thin-plate spline
#'
#' @description `coastalTps` uses thin-plate spline interpolation to estimate
#' the effects of coastal buffering of land-temperatures by the sea.
#'
#' @param dT a coarse-resolution raster of sea - land temperatures (ºC).
#' @param lsw a fine-resolution raster of coastal exposure upwind, as produced by [invls()].
#' @param lsa a fine-resolution raster of mean coastal exposure in all directions.
#' @return a fine-resolution raster of sea - land temperature differences (ºC).
#' @export
#' @importFrom rgcvpack fitTps
#' @importFrom rgcvpack predict.Tps
#'
#' @examples
#' # =========================================
#' # Calculate land-sea temperature difference
#' # =========================================
#' temp <- tas[,,1]
#' sst <- 10.665
#' dT <- if_raster(sst - temp, dtm1km)
#' # ============================
#' # Obtain coastal exposure data
#' # ============================
#' lsw <- landsearatios[,,7]  # upwind
#' lsa <- apply(landsearatios, c(1, 2), mean) # mean, all directions
#' lsw <- if_raster(lsw, dtm100m)
#' lsa <- if_raster(lsa, dtm100m)
#' # ==========================================================
#' # Calculate coastal effects using thin-plate spline and plot
#' # ==========================================================
#' dTf <- coastalTps(dT, lsw, lsa)
#' par(mfrow = c(2, 1))
#' plot(sst - dT, main = expression(paste("Temperature ",(~degree~C))))
#' plot(sst - dTf, main = expression(paste("Temperature ",(~degree~C))))
coastalTps <- function(dT, lsw, lsa) {
  lswc <- resample(lsw, dT)
  lsac <- resample(lsa, dT)
  xy <- data.frame(xyFromCell(lswc, 1:ncell(dT)))
  z1 <- extract(lswc, xy)
  z2 <- extract(lsac, xy)
  v <- extract(dT, xy)
  xyz <- cbind(xy, z1, z2)
  sel <- which(is.na(v) == F)
  v <- v[is.na(v) == F]
  xyz <- xyz[sel, ]
  tps <- fitTps(xyz, v, m = 3)
  xy <- data.frame(xyFromCell(lsw, 1:ncell(lsw)))
  z1 <- extract(lsw, xy)
  z2 <- extract(lsa, xy)
  xyz <- cbind(xy, z1, z2)
  sel <- which(is.na(z1) == FALSE)
  xyz <- xyz[sel, ]
  xy$z <- NA
  xy$z[sel] <- predict.Tps(tps, xyz)
  r <- rasterFromXYZ(xy)
  r
}
#' Calculates the moist adiabatic lapse rate
#'
#' @description `lapserate` is used to calculate changes in temperature with height.
#'
#' @param tc a single numeric value, raster object, two-dimensional array or matrix of temperature (ºC).
#' @param h a single numeric value, raster object, two-dimensional array or matrix of specific humidity (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}}).
#' @param p an optional single numeric value, raster object, two-dimensional array or matrix of atmospheric pressure (Pa).
#'
#' @return the lapse rate (\ifelse{html}{\out{º m<sup>-1</sup> }}{\eqn{ º m^{-1}}}).
#' @export
#' @import raster
#'
#' @details if tc is a raster, a raster object is returned. This function calculates the
#' theoretical lapse rate. Environmental lapse rates can vary due to winds.
#'
#' @examples
#' lapserate(20, 0) * 1000 # dry lapse rate per km
#' h <- humidityconvert(100, intype = "relative", 20)
#' lapserate(20, h$specific) * 1000 # lapse rate per km when relative humidity is 100%
lapserate <- function(tc, h, p = 101300) {
  r <- tc
  tc <- is_raster(tc)
  h <- is_raster(h)
  p <- is_raster(p)
  pk <- p / 1000
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ws <- 0.622 * e0 / pk
  rh <- (h / ws) * 100
  rh[rh > 100] <- 100
  ea <- e0 * (rh / 100)
  rv <- 0.622 * ea / (pk - ea)
  lr <- 9.8076 * (1 + (2501000 * rv) / (287 * (tc + 273.15))) / (1003.5 +
        (0.622 * 2501000 ^ 2 * rv) / (287 * (tc + 273.15) ^ 2))
  lr <- lr * -1
  if_raster(lr, r)
}
#' Applies height correction to wind speed measurements
#'
#' @description `windheight` is used to to apply a height correction to wind speed measured at a specified height above ground level to obtain estimates of wind speed at a desired height above the ground.
#'
#' @param ui numeric value(s) of measured wind speed (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{ m s^{-1}}}) at height `zi` (m).
#' @param zi a numeric value idicating the height (m) above the ground at which `ui` was measured.
#' @param zo a numeric value indicating the height (m) above ground level at which wind speed is desired.
#'
#' @details Thus function assumes a logarithmic height profile to convert
#' wind speeds. It performs innacurately when `uo` is lower than 0.2
#' and a warning is given. If `uo` is below ~0.08 then the logairthmic height
#' profile cannot be used, and `uo` is converted to 0.1 and a warning given.
#'
#' @seealso The function [windcoef()] calculates a topographic or vegetation sheltering effect.
#'
#' @return numeric values(s) of wind speed at height specified by `zo` (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{m s^{-1}}}).
#' @export
#'
#' @examples
#' windheight(3, 10, 1) # good
#' windheight(3, 10, 0.15) # performs poorly. Warning given
#' windheight(3, 10, 0.05) # cannot calculate. ui converted and warning given
#'
windheight <- function(ui, zi, zo) {
  if (zo < 0.2 & zo > (5.42 / 67.8)) {
    warning("Wind-height profile function performs poorly when wind
            height is below 20 cm")
  }
  if (zo <= (5.42 / 67.8)) {
    warning(paste("wind-height profile cannot be calculated for height ",
            zo * 100, " cm"))
    print("Height converted to 10 cm")
    zo <- 0.1
  }
  uo <- ui * log(67.8 * zo - 5.42) / log(67.8 * zi - 5.42)
  uo
}
#' Calculates wind shelter coefficient
#'
#' @description `windcoef` is used to apply a topographic shelter coefficient to wind data.
#'
#' @param dsm raster object, two-dimensional array or matrix of elevations (m) derived either from a digital terrain or digital surface model, and orientated as if derived using [is_raster()]. I.e. `[1, 1]` is the NW corner.
#' @param direction a single numeric value specifying the direction from which the wind is blowing (º).
#' @param hgt a single numeric value specifying the height (m) at which wind speed is derived or measured. The wind speeds returned are also for this height, and account for the fact topography affords less shelter to wind at greater heights.
#' @param res a single numeric value specifying the the resolution (m) of `dsm`.
#'
#' @details If dsm is a raster object, then a raster object is returned.
#' If elevations are derived from a digital terrain model, then
#' the sheltering effects of vegetation are ignored. If derived from a
#' digital surface model, then the sheltering effects of vegetation are
#' accounted for. If `res` is unspecified `dtm` is assumed to have a resolution of one m.
#'
#' @seealso The function [windheight()] converts measured wind heights to a standard reference height.
#'
#' @return a raster object, or two-dimensional array of shelter coefficients. E.g. a shelter coefficient of 0.5 indicates that wind speed at that location is 0.5 times that measured at an unobscured location.
#' @import raster
#' @export
#'
#' @examples
#' dsm <- dtm1m + veg_hgt
#' wc <- windcoef(dsm, 0)
#' plot(mask(wc, dtm1m), main ="Northerly wind shelter coefficient")
windcoef <- function(dsm, direction, hgt = 1, res = 1) {
  r <- dsm
  dsm <- is_raster(dsm)
  dsm[is.na(dsm)] <- 0
  dtm <- dsm / res
  hgt <- hgt / res
  direction <- direction - 90
  azi <- direction * (pi / 180)
  horizon <- array(0, dim(dtm))
  dtm3 <- array(0, dim(dtm) + 200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x + 100), 101:(y + 100)] <- dtm
  for (step in 1:10) {
    horizon[1:x, 1:y] <- pmax(horizon[1:x, 1:y], (dtm3[(101 + sin(azi) *
                         step ^ 2):(x + 100 + sin(azi) * step ^ 2),
                         (101 + cos(azi) * step ^ 2 ):(y + 100 + cos(azi) *
                         step ^ 2)] - dtm3[101:(x + 100), 101 : (y + 100)]) /
                         step ^ 2)
    horizon[1:x, 1:y] <- ifelse(horizon[1:x, 1:y] < (hgt / step ^ 2), 0,
                                horizon[1:x, 1:y])
  }
  index <- 1 - atan(0.17 * 100 * horizon) / 1.65
  index <- if_raster(index, r)
  index
}
#' Calculates sunrise, sunset and daylength
#'
#' @description `suntimes` is used to calculate, sunrise, sunset and daylength on any given data at a specified location.
#'
#' @param julian the Julian day as returned by [julday()].
#' @param lat latitude of the location for which `suntime` is required (decimal degrees, -ve south of equator).
#' @param long longitude of the location for which `suntime` is required (decimal degrees, -ve west of Greenwich meridian).
#' @param tz an optional numeric value specifying the time zones expressed as hours different from GMT (-ve to west).
#' @param dst an optional numeric value representing the local summer time adjustment (hours, e.g. +1 for BST).
#'
#' @return a data.frame with three components:
#' @return `sunrise` a vector of sunrise (hours in 24 hour clock).
#' @return `sunrise` a vector of sunset (hours in 24 hour clock).
#' @return `sunrise` a vector of daylengths (hours).
#' @export
#'
#' @details if the sun is above or below the horizon for 24 hours on any
#' days, a warning is given, and the sunrise and sunset returned are the
#' appoximate times at which that the sun is closest to the horizon.
#'
#' @examples
#' jd <- julday(2018, 1, 16)
#' suntimes(jd, 50.17, -5.12) # Cornwall, UK
#' suntimes(jd, 78.22, 15.64, 1) # Longyearbyen, Svalbad
#' suntimes(-54.94, -67.61, -3) # Puerto Williams, Cabo de Hornos, Chile
suntimes <- function(julian, lat, long, tz = 0, dst = 0) {
  tz <- tz + dst
  lw <- long * -1
  n <- julian - 2451545 - 0.0009 - (lw / 360)
  n <- floor (n) + 0.5
  sn <- 2451545 + 0.0009 + (lw / 360) + n
  msa <- (357.5291 + 0.98560028 * (sn - 2451545))%%360
  eoc <- 1.9148 * sin(msa * pi / 180) + 0.02 * sin(2 * msa * pi / 180) +
    0.0003 * sin(3 * msa * pi / 180)
  ecl <- (msa + 102.9372 + eoc + 180)%%360
  st <- sn + (0.0053 * sin(msa * pi / 180)) - (0.0069 *
                                               sin(2 * ecl * pi / 180))
  d <- asin(sin(ecl * pi / 180) * sin(23.45 * pi / 180))
  coshas <- (sin(-0.83 * pi / 180) - sin(lat * pi / 180) * sin(d)) /
    (cos(lat * pi / 180) * cos(d))
  if (max(coshas) > 1) {
    coshas <- ifelse(coshas > 1, 1, coshas)
    warning("sun above horizon for 24 hours on some days")
  }
  if (min(coshas) < -1) {
    coshas <- ifelse(coshas < -1, -1, coshas)
    warning("sun below horizon for 24 hours on some days")
  }
  has <- acos(coshas)
  jset <- 2451545 + 0.0009 + (((has * 180 / pi + lw) / 360) + n +
                                0.0053 * sin(msa * pi / 180)) - 0.0069 *
                                sin(2 * ecl * pi / 180)
  jrise <- st - (jset - st)
  hset <- jset%%1 * 24 + tz
  hrise <- jrise%%1 * 24 + tz
  dl <- (jset - jrise) * 24
  sunvars <- data.frame(sunrise = hrise, sunset = hset, daylight = dl)
  sunvars
}
#' derived fraction of solar day
.solarday <- function(julian, localtime, lat, long, tz = 0, dst = 0) {
  src <- suntimes(julian, lat, long, tz, dst)$sunrise
  ssc <- suntimes(julian, lat, long, tz, dst)$sunset
  ssp <- suntimes(julian - 1, lat, long, tz, dst)$sunset
  srn <- suntimes(julian + 1, lat, long, tz, dst)$sunrise
  st <- ifelse(localtime >= src & localtime <= ssc,
               ((localtime - src) / (ssc - src)) * 12 + 6, localtime)
  st <- ifelse(localtime > ssc, ((localtime - ssc) / (srn + 24 - ssc)) *
                 12 + 18, st)
  st <- ifelse(localtime < src, ((localtime + 24 - ssp) / (src + 24 - ssp)) *
                 12 - 6, st)
  st <- (st / 24)%%1
  st
}
#' derived proportion of maximum radiation
.propmaxrad <- function(dif, dct, am) {
  theta <- 1.1 * (0.7 ^ (am ^ 0.678))
  mxr <- 4.87 * theta
  pr <- (dif + dct) / mxr
  pr <- ifelse(pr > 1, 1, pr)
  pr
}
#' calculates emsisivity from humidity etc
.emissivity <- function(h, tc, n, p = 100346.13) {
  pk <- p / 1000
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ws <- 0.622 * e0 / pk
  rh <- (h / ws) * 100
  rh[rh > 100] <- 100
  ea <- e0 * rh / 100
  emcs <- 0.23 + 0.433 * (ea / (tc + 273.15)) ^ (1 / 8)
  em <- emcs * (1 - n ^ 2) + 0.976 * n ^ 2
  em
}


#' Derives hourly temperatures from daily data
#'
#' @description `hourlytemp` is used to derive hourly temperatures from daily maxima and minima.
#'
#' @param julian vector of julian days expressed as integers for every day for which `mintemp` and `maxtemp` are provided, as returned by function [julday()].
#' @param em an optional vector of hourly emissivities. If not provided, calculated from `h`, `n` and `p`.
#' @param h a vector of hourly specific humidities (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}}). Ignored if `em` provided.
#' @param n a vector of hourly fractional cloud cover values (range 0 - 1). Ignored if `em` provided.
#' @param p an optional vector of hourly atmospheric pressure values (Pa). Ignored if `em` provided.
#' @param dni a vector of hourly direct radiation values normal to the solar beam (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}}).
#' @param dif a vector of hourly diffuse radiation values (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}}).
#' @param mintemp a vector of daily minimum temperatures (ºC).
#' @param maxtemp a bector of daily maximum temperatures (ºC).
#' @param lat a single numeric value representing the latitude of the location for which hourly temperatures are required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the longitude of the location for which hourly temperatures are required (decimal degrees, -ve west of Greenwich meridian).
#' @param merid an optional numeric value representing the longitude of the local time zone meridian (º) (0 for UK time).
#' @param tz an optional single  numeric value or vector of values specifying the time zones expressed as hours different from GMT for each day (-ve to west).
#' @param dst a single numeric value or vector of values representing the local summer time adjustment for each day (hours, e.g. +1 for BST).
#'
#' @return a vector of hourly temperatures (ºC).
#' @export
#'
#' @details A warning is returned if any of following conditions are not
#' met. (1)  `h`, `n`, `dif` and `dct` differ in length.
#' (2) `julian`, `mintemp` and `maxtemp` differ in length. (3) E.g. `h` / 24 is not
#' an integer. (4) the length of e.g. `h` is not equal to the length of
#' e.g. `mintemp` x 24.
#'
#' @examples
#' jd <- julday(2010, 5 , c(1,2))
#' ht <- hourlytemp(jd, NA, microvars$humidity[1:48], microvars$cloudcover[1:48],
#'                  microvars$pressure[1:48], dni = microvars$dni[1:48],
#'                  dif = microvars$dif[1:48], c(7.5, 7.2), c(14.6, 15.2),
#'                  49.968, -5.216)
#' plot(ht ~ c(0:47),type="l", xlab = "Hour", ylab = "Temperature",
#'      main = paste("tmins:", 7.2, 7.5, "tmaxs:", 14.6, 15.2))
hourlytemp <- function(julian, em = NA, h, n, p = 100346.13, dni, dif,
                       mintemp, maxtemp, lat, long, merid=0, tz=0, dst=0) {
  l1 <- unlist(lapply(list(dni, dif), length))
  l2 <- unlist(lapply(list(julian, mintemp, maxtemp), length))
  if (length(unique(l1)) != 1) {
    warning("Number of hourly values not identical")
  }
  if (length(unique(l2)) != 1) {
    warning("Number of daily values not identical")
  }
  if (length(dni)%%24 != 0) warning ("Hourly values must be multiple of 24")
  if (l1[1] != l2[1] * 24) {
    warning("Number of hourly values not 24 times number of daily values")
  }
  o <- order(rep(c(1:length(julian)), 24))
  jd <- rep(julian, 24)[o]
  localtime <- rep(c(0:23), length(mintemp))
  solfrac <- .solarday(julian, localtime, lat, long, tz, dst)
  am <- airmasscoef(localtime, lat, long, julian, merid, dst)
  dct <- dni * siflat(localtime, lat, long, jd, merid, dst)
  pr <- .propmaxrad(dif, dct, am)
  tfrac <- 110.42201 * sin(solfrac) - 38.64802 * cos(solfrac) - 79.82963 *
    sin(2 * solfrac) + 65.39122 * cos(2 * solfrac) + 15.54387 *
    sin(3 * solfrac) - 26.30047 * cos(3 * solfrac)
  mns <- rep(mintemp, 24)[o]
  mxs <- rep(maxtemp, 24)[o]
  tc <- (mxs - mns) * tfrac + mns
  if (is.na(em[1]))
    em <- .emissivity(h, tc, n, p)
  day <- which(is.na(pr) == F)
  ngt <- which(is.na(pr))
  tfrac[day] <- -0.5154639 + 1.1060481 * tfrac[day] + 6.1147121 * pr[day] -
    7.5826432 * tfrac[day] * pr[day]
  tfrac[ngt] <- -4.614759 + 5.909548 * tfrac[ngt] + 4.308512 * em[ngt] -
    5.020017 * tfrac[ngt] * em[ngt]
  tfrac <- 1 / (1 + exp(-1 * tfrac))
  td <- array(tfrac, dim = c(24, length(tfrac) / 24))
  tmns <- apply(td, 2, min)
  tmxs <- apply(td, 2, max)
  for (i in 1:dim(td)[2]) {
    td[, i] <- (td[, i] - tmns[i]) / (tmxs[i] - tmns[i])
  }
  tfrac <- as.vector(td)
  tc <- (mxs - mns) * tfrac + mns
  tc
}
#' Applies a spline function to an array of values
#'
#' @description `arrayspline` is used to derive e.g. hourly climate data from e.g. daily values.
#'
#' @param a a three-dimensional array (row, column, time)
#' @param tme an object of class POSIXct of times for `a`. I.e. `length(tme) = dim(a)[3]`
#' @param nfact indicates the time interval for which outputs are required. E.g to derive hourly from daily data `nfact = 24`, or derive six-hourly from daily data `nfact = 4`
#' @param out an optional character vector indicating the time for which the output is required. Format must be as for `tme`.
#'
#' @return
#' If out is unspecified, a three dimensional array of size `c(a[1], a[2], (length(tme) - 1) * nfact + 1)` is returned.
#' If out is specified, a three dimensional array of size `c(a[1], a[2], length(out))` is returned.
#'
#' @importFrom stats splinefun
#' @export
#'
#' @details
#' `arrayspline` uses the Forsythe, Malcolm and Moler method of splining, as specified by
#' `"fmm"` in [spline()]. If `a[i, j, ]` is a vector of NAs, then tthe corresponding vector in
#' the output array is also a vector of NAs. It is assumed that all spatial data within `a` have
#' equivelent times. I.e. the time of both `a[1, 1, 1]` and `a[1000, 1000, 1]` is identical and
#' equal to `tme[1]`.
#'
#' @examples
#' tme <- as.POSIXct(c(0:364) * 24 * 3600, origin="2010-01-01", tz = "GMT")
#' h <- arrayspline(huss, tme, out = "2010-05-01 11:00")
#' plot(raster(h, template = dtm1km),
#'      main = "Specific humidity 2010-05-01 11:00")
arrayspline <- function(a, tme, nfact = 24, out = NA) {
  n <- (length(tme) - 1) * nfact + 1
  sf <- function(y, tme, nfact) {
    n <- (length(tme) - 1) * nfact + 1
    if (is.na(mean(y)) == FALSE) {
      xy <- spline(as.numeric(tme), y, n = n)
      xy$y
    }
    else rep(NA, n)
  }
  ao <-  aperm(apply(a, c(1, 2), sf, tme, nfact), c(2, 3, 1))
  if (is.na(out[1]) == FALSE) {
    out <- as.POSIXct(out, tz = format(tme[1], "%Z"))
    tout <- spline(as.numeric(tme), c(1:length(tme)), n)$x
    tout <- as.POSIXct(tout, origin = "1970-01-01 00:00",
                       tz = format(tme[1], "%Z"))
    s <- 0
    for (i in 1:length(out)) {
      s[i] <- which(tout == out[i])
    }
    ao <- ao[,,s]
  }
  ao
}
#' Fits micro- or mesoclimate model
#'
#' @description
#' `fitmicro` is used to fit a micro- or mesoclimate model using field temperature readings, and estimates of reference temperature, net radiation
#' and wind speed at the locations of those readings.
#'
#' @param microfitdata a data.frame with at least the following columns (see, for example, microfit data):
#' \describe{
#'   \item{temperature}{microclimate temperature readings}
#'   \item{reftemp}{Reference (e.g. coarse-scale or weather station) temperatures}
#'   \item{wind}{Wind speeds}
#'   \item{netrad}{Net radiation values}
#' }
#' @param alldata an optional logical value indicating whether to fit the model using all data (TRUE) or using a randomization procedure (FALSE). See details.
#' @param windthresh an optional single numeric value indicating the threshold wind speed above which an alternative linear relationship between net radiation the microclimate temperature anomoly is fitted. See details.
#' @param iter a single integer specifying the iterations to perform during randomization. Ignored if `alldata` = TRUE.
#'
#' @return a data,frame with the following columns:
#' \describe{
#'   \item{Estimate}{parameter estimates and `windthresh`}
#'   \item{Std.Dev}{Standard deviation of parameter estimates}
#'   \item{P}{Two-tailed p-value}
#' }
#' @export
#' @details
#' If modelling mesoclimate, it is assumed that altitudinal, coastal and cold-air
#' drainage effects have already been accounted for in the calculation of `reftemp`.
#' It is therefore assumed that the most important energy fluxes determining near-surface
#' temperature are those due to the radiation flux and convection that occurs at
#' the surface-atmosphere boundary. Heat fluxes into the soil and latent heat
#' exchange are considered to be small and proportional to the net radiation
#' flux, and the heat capacity of the vegetation is considered to be small so
#' that, compared to the time-scale of the model, surface temperature rapidly
#' reach equilibrium. In consequence, the difference between the near-ground
#' temperature and the ambient temperature is a linear function of `netrad`.
#' The gradient of this linear relationship is a measure of the thermal
#' coupling of the surface to the atmosphere. If this relationship is applied
#' to vegetation, assuming the canopy to act like a surface, while air density
#' and the specific heat of air at constant pressure are constant, the slope
#' varies as a function of a wind speed factor, such that different slope values
#' are assumed under high and low wind conditions.  Hence, `fitmicro` fits a
#' linear model of the form `lm((temperature - reftemp) ~ netrad * windfact)`
#' where windfact is given by `ifelse(wind > windthresh, 1, 0)`. If all data is
#' FALSE, random subsets of the data are selected and the analyses repeated
#' `iter` times to reduce the effects of of temporal autocorrelation. Parameter
#' estimates are derived as the median of all runs. If no value is provided for `windthresh`,
#' it is derived by iteratively trying out different values, and selecting that which
#' yields the best fit. The gradient of the relationship is also dependent on
#' vegetation structure, and in some circumstances it may therefore be advisable
#' to fit seperate models for each vegetation type.
#'
#' @examples
#' fitmicro(microfitdata)
#' fitmicro(mesofitdata, alldata = TRUE)
fitmicro <- function(microfitdata, alldata = FALSE, windthresh = NA,
                     iter = 999) {
  pvals <- function(x) {
    iter <- length(x)
    d <- floor(log(iter, 10))
    p <- length(x > 0) / iter
    p <- ifelse(p > 0.5, 1 - p, p)
    p <- round(p * 2, d)
    p
  }
  microfitdata$anom <- microfitdata$temperature - microfitdata$reftemp
  if (is.na(windthresh)) {
    wrange <- floor(max(microfitdata$wind) - min(microfitdata$wind))
    rsq <- 0
    for (i in 1:100 * wrange) {
      wind.fact <- ifelse(microfitdata$wind > (i / 100), 1, 0)
      m <- lm(microfitdata$anom ~ microfitdata$netrad * wind.fact)
      rsq[i] <- summary(m)$r.squared
    }
    wthresh <- which(rsq == max(rsq, na.rm = TRUE)) / 100
  }
  if (alldata) {
    wf <- ifelse(microfitdata$wind > wthresh, 1, 0)
    m1 <- summary(lm(microfitdata$anom ~ wf * microfitdata$netrad))
    mdns <- as.vector(m1$coef[, 1])
    sds <-  as.vector(m1$coef[, 2]) * sqrt(dim(microfitdata)[1])
    pvs <- as.vector(m1$coef[, 4])
  }
  else {
    ln <- round(dim(microfitdata)[1] / 100, 0) * 10
    pint <- 0
    pwf <-0
    pnr <- 0
    pnrwf <- 0
    for (i in 1:iter) {
      u <- round(runif(ln, 1, dim(microfitdata)[1]), 0)
      one.iter <- microfitdata[u, ]
      wf <- ifelse(one.iter$wind > wthresh, 1, 0)
      m1 <- lm(one.iter$anom ~ wf * one.iter$netrad)
      pint[i] <- m1$coef[1]
      pwf[i] <- m1$coef[2]
      pnr[i] <- m1$coef[3]
      pnrwf[i] <- m1$coef[4]
    }
    mdns <- c(median(pint, na.rm = T), median(pwf, na.rm = T),
              median(pnr, na.rm = T), median(pnrwf, na.rm = T))
    sds  <- c(sd(pint, na.rm = T), sd(pwf, na.rm = T), sd(pnr, na.rm = T),
              sd(pnrwf, na.rm = T))
    pvs <- c(pvals(pint), pvals(pwf), pvals(pnr), pvals(pnrwf))
  }
  params <- data.frame(Estimate = c(mdns, wthresh), Std.Dev = c(sds, ""),
                       P = c(pvs, ""))
  row.names(params) <- c("Intercept", "Wind factor", "Net radiation",
                         "Wind factor:Net radiation", "Wind threshold")
  params
}
#' Runs micro- or mesoclimate model
#'
#' @description `microrun` produces a high-resolution dataset of downscaled
#' temperatures for one time interval
#'
#' @param params a data.frame of parameter estimates as produced by [fitmicro()]
#' @param netrad a raster object, two-dimensional array or matrix of downscaled net radiation as produced by [shortwaveveg()] - [longwaveveg()] or [shortwavetopo()] - [longwavetopo()].
#' @param wind a raster object, two-dimensional array or matrix of downscaled wind speed, as produced by reference wind speed x the output of [windcoef()].
#'
#' @return a raster object, two-dimensional array or matrix of temperature anomolies from reference temperature, normally in ºC, but units depend on those used in [fitmicro()].
#' @import raster
#' @export
#'
#' @details
#' If `netrad` is a raster object, a raster object is returned.
#' If modelling mesoclimate, it is assumed that altitudinal, coastal and cold-air
#' drainage effects have already been accounted for in the calculation of reference
#' temperature (see example). It is assumed that the most important energy fluxes
#' determining near-surface temperature are those due to the radiation flux and convection
#' that occurs at the surface–atmosphere boundary. Heat fluxes into the soil and latent heat
#' exchange are considered to be small and proportional to the net radiation
#' flux, and the heat capacity of the vegetation is considered to be small so
#' that, compared to the time-scale of the model, surface temperature rapidly
#' reach equilibrium. In consequence, the difference between the near-ground
#' temperature and the ambient temperature is a linear function of `netrad`.
#' The gradient of this linear relationship is a measure of the thermal
#' coupling of the surface to the atmosphere. If this relationship is applied
#' to vegetation, assuming the canopy to act like a surface, while air density
#' and the specific heat of air at constant pressure are constant, the slope
#' varies as a function of a wind speed factor, such that different slope values
#' are assumed under high and low wind conditions.
#'
#' @examples
#' # =======================================================================
#' # Run microclimate model for 2010-05-24 11:00 (one of the warmest hours)
#' # =======================================================================
#' params <- fitmicro(microfitdata)
#' netrad <- netshort1m - netlong1m
#' tempanom <- runmicro(params, netrad, wind1m)
#' reftemp <- raster(temp100[,,564])
#' extent(reftemp) <- extent(dtm1m)
#' reftemp <- resample(reftemp, dtm1m)
#' temps <- tempanom + getValues(reftemp, format = "matrix")
#' plot(if_raster(temps, dtm1m), main =
#'      expression(paste("Temperature ",(~degree~C))))
#'
#' # ======================================================================
#' # Run mesoclimate model for 2010-05-01 11:00 from first principles
#' # ======================================================================
#'
#' # -------------------------
#' # Resample raster function
#' # -------------------------
#' resampleraster <- function(a, ro) {
#'   r <- raster(a)
#'   extent(r) <- c(-5.40, -5.00, 49.90, 50.15)
#'   crs(r) <- "+init=epsg:4326"
#'   r <- projectRaster(r, crs = "+init=epsg:27700")
#'   r <- resample(r, ro)
#'   as.matrix(r)
#' }
#'
#' # --------------------------
#' # Resample raster: 24 hours
#' # --------------------------
#' get24 <- function(a) {
#'   ao <- array(NA, dim = c(dim(dtm1km)[1:2], 24))
#'   for (i in 1:24) {
#'     ai <- a[,,2880 + i]
#'     ao[,,i] <- resampleraster(ai, dtm1km)
#'   }
#'   ao
#' }
#'
#' # ----------------------------
#' # Derive hourly temperatures
#' # ----------------------------
#' tmax <- tas[,,121] + dtr[,,121] / 2
#' tmin <- tas[,,121] - dtr[,,121] / 2
#' tme <- as.POSIXct(c(0:364) * 24 * 3600, origin="2010-01-01", tz = "GMT")
#' out <- as.POSIXct(c(0:23) * 3600, origin="2010-05-01", tz = "GMT")
#' h <- arrayspline(huss, tme, out = out)
#' p <- arrayspline(pres, tme, out = out)
#' n <- get24(cfc)
#' dni <- get24(dnirad)
#' dif <- get24(difrad)
#' jd <- julday(2010, 5 , 1)
#' tc <- h[,,1] * NA
#' lr <- h[,,1] * NA # Also calculates lapse rate
#' for (i in 1:19) {
#'   for (j in 1:22) {
#'     if (is.na(tmax[i,j]) == F)
#'     {
#'       ht <- hourlytemp(jd, h[i, j, ], n[i, j, ], p[i, j, ], dni[i, j, ],
#'                        dif[i, j, ], tmin[i,j], tmax[i,j], 50.02, -5.20)
#'       tc[i, j]<- ht[12]
#'       lr[i, j] <- lapserate(tc[i, j], h[i, j, 12], p[i, j, 12])
#'     }
#'   }
#' }
#' # ----------------------------
#' # Calculate coastal effects
#' # ----------------------------
#' sst <- 10.771
#' dT <- if_raster(sst - tc, dtm1km)
#' lsw <- if_raster(landsearatios[,,28], dtm100m)  # upwind
#' lsa <- if_raster(apply(landsearatios, c(1, 2), mean), dtm100m) # mean, all directions
#' dTf <- coastalTps(dT, lsw, lsa)
#'
#' # ------------------------------
#' # Calculate altitudinal effects
#' # ------------------------------
#' lrr <- if_raster(lr, dtm1km)
#' lrr <- resample (lrr, dtm100m)
#' tc <- sst - dTf + lrr * dtm100m
#'
#' # ------------------------------
#' # Downscale radiation
#' # ------------------------------
#' dni <- resampleraster(dnirad[,,2891], dtm100m)
#' dif <- resampleraster(difrad[,,2891], dtm100m)
#' n <- resampleraster(cfc[,,2891], dtm100m)
#' h <- resample(if_raster(h[,,12], dtm1km), dtm100m)
#' p <- resample(if_raster(p[,,12], dtm1km), dtm100m)
#' sv <- skyviewtopo(dtm100m)
#' netshort <- shortwavetopo(dni, dif, jd, 11, dtm = dtm100m, svf = sv)
#' netlong <- longwavetopo(h, tc, p, n, sv)
#' netrad <- netshort - netlong
#'
#' # ------------------
#' # Downscale wind
#' # ------------------
#' ws <- array(windheight(wind2010$wind10m, 10, 1), dim = c(1, 1, 8760))
#' wh <- arrayspline(ws, as.POSIXct(wind2010$obs_time), 6, "2010-05-01 11:00")
#' ws <- windcoef(dtm100m, 270, res = 100) * wh
#'
#' # ------------------
#' # Fit and run model
#' # ------------------
#' params <- fitmicro(mesofitdata)
#' anom <- runmicro(params, netrad, ws)
#' tc <- tc + anom
#' plot(mask(tc, dtm100m), main =
#'      expression(paste("Mesoclimate temperature ",(~degree~C))))
runmicro <- function(params, netrad, wind) {
  r <- netrad
  netrad <- is_raster(netrad)
  wind <- is_raster(wind)
  wf <- ifelse(wind > params[5, 1], 1, 0)
  temp <- params[1, 1] + params[2, 1] * wf + params[3, 1] * netrad +
    params[4, 1] * wf * netrad
  if_raster(temp, r)
}

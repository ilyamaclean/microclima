#' Gets raster of elevations
#'
#' @description Several web services provide access to raster elevations. This function
#' provides access to the Registry of Open Data on AWS to retrieve a raster object of
#' elevations.
#'
#' @param r optional raster object covering the location for which elevations are
#' required
#' @param lat the latitude (in decimal degrees) of the location for which elevations
#' are required. Not used if `r` is a raster.
#' @param long the longitude (in decimal degrees) of the location for which elevations
#' are required. Not used if `r` is a raster
#' @param resolution the resolution (in metres) of the required dem (see details).
#' @param zmin assumed sea-level height. Values below this are set to zmin.
#'
#' @return a raster object of elevations in metres.
#' @import raster
#' @import elevatr
#' @export
#'
#' @details Currently, `get_dem` utilises the Amazon Web Services (AWS)
#' (https://aws.amazon.com/public-datasets/terrain/) to retrieve dems.  The spatial
#' resolution of the data available vary by location
#' (see https://mapzen.com/documentation/terrain-tiles/data-sources/ for details).
#' Global coverage is ~30m resolution and if the resolution of data requested is
#' higher than that available, elevations are derived by interpolation. If `r` is
#' unspecified, then `lat` and `long` are used to retrieve a 200 x 200 cell raster
#' centred on `lat` and `long`. Bathymetry is also returned automatically form AWS,
#' so `zmin` allows the user to specify a sea-level height. Elevations below `zmin`
#' are set at `zmin`. For use with `microclima` the returned raster has a Mercator
#' Coordinate Reference System if the latitude specified by `lat` is between 80
#' degrees S and 84 degrees N. In polar regions WGS 84 / NSIDC Sea Ice Polar
#' Stereographic projections are used. X, y and Z units are all in metres. If `r` is
#' specified, the Coordinate Reference System of `r` is used.
#'
#' @examples
#' dem50m <- get_dem(dtm100m, resolution = 50)
#' plot(dem50m) # 50m raster, OSGB projection system
#'
#' dem5m <- get_dem(lat = 49.97, long = -5.22, resolution = 5)
#' plot(dem5m) # 5m raster, Mercator projection system
get_dem <- function(r = NA, lat, long, resolution = 30, zmin = 0) {
  if (resolution < 30) {
    warning("Higher resolution data only available for some locations. DEM
            may be derived by interpolation")
  }
  if (class(r) != "RasterLayer") {
    xy <- data.frame(x = long, y = lat)
    coordinates(xy) = ~x + y
    proj4string(xy) = "+init=epsg:4326"
    if (lat >= -80 & lat <= 84)
      xy <- as.data.frame(spTransform(xy, CRS("+init=epsg:3395")))
    if (lat > 84)
      xy <- as.data.frame(spTransform(xy, CRS("+init=epsg:3413")))
    if (lat < -80)
      xy <- as.data.frame(spTransform(xy, CRS("+init=epsg:3976")))
    e <- extent(c(xy$x - 100 * resolution, xy$x + 100 * resolution,
                  xy$y - 100 * resolution, xy$y + 100 * resolution))
    r <- raster(e)
    res(r) <- resolution
    if (lat >= -80 & lat <= 84)
      crs(r) <- "+init=epsg:3395"
    if (lat > 84)
      crs(r) <- "+init=epsg:3413"
    if (lat < -80)
      crs(r) <- "+init=epsg:3976"
  } else {
    lat <- latlongfromraster(r)$lat
    long <- latlongfromraster(r)$long
    res(r) <- resolution
  }
  z = ceiling(log((cos(lat * pi/180) * 2 * pi * 6378137) / (256 * resolution), 2))
  z <- ifelse(z > 14, 14, z)
  r2 <-get_elev_raster(r, z = z, src = "aws")
  r2<- resample(r2, r)
  m2 <- getValues(r2, format = "matrix")
  m2[m2 < zmin] <- zmin
  m2[is.na(m2)] <- zmin
  r2 <- if_raster(m2, r2)
  return(r2)
}
#' Obtains NCEP data required to run microclima
#'
#' @param tme a POSIXlt object covering the duration for which data are required.
#' Hourly data are returned irrespective of the time interval of `tme`.
#' @param lat the latitude of the location for which data are required
#' @param long the longitude of the location for which data are required
#' @param reanalysis2 Logical. Should data be obtained from the Reanalysis II dataset (default) or
#' from Reanalysis I (data prior to 1979).
#'
#' @return a dataframe with the following variables: (i) obs_time:times in YYY-MM-DD HH:MM format (UTC),
#' (ii) Tk: mean temperatures at 2m (K), (iii) Tkmin:  minimum temperatures at 2m (K) (iv) Tkmax:
#' maximum temperatures at 2m (K), (v) sh: specific humidity at 2m (Kg / Kg),
#' (vi) pr: surface pressure (Pa), (vii) wu u wind vector at 10m (m/s) (viii) wv: v wind vector
#' at 10m (m/s) (ix) dlw: downward longwave radiation flux (\ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^-2}}),
#' (x) ulw: upward longwave radiation flux (\ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^-2}}), (xi)
#' dsw: downward shortwave radiation flux (\ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^-2}}),
#' (xii) tcdc: Total cloud cover for the entire atmosphere.

#' @import RNCEP
#' @export
#'
#' @seealso [hourlyNCEP()]
#'
#' @details The NCEP-NCAR and NCEP–DOE Atmospheric Model Intercomparison Project (Kanamitso et al 2002)
#' provides estimates of a range of climate and radiation data for the period 1948 to present.
#' Data at six-hourly intervals are available at a spatial resolution of ~ 1.8 to 2.5º. The function downloads
#' the following datasets: mean, min and max air temperature at 2m,  specific humidity and 2m, surface pressure,
#' the u and v wind vectors at 10, the downward and upward longwave radiation fluxes at the surface
#' and the downward shortwave flux at the surface. Data for days either side of `tme` are downloaded to aid with
#' interpolation by [hourlyNCEP()]. Requires internet connection.
#'
#' @references Kanamitsu M, Ebisuzaki W, Woollen J, Yang SK, Hnilo JJ, Fiorino M, Potter GL (2002)
#'Ncep–doe amip-ii reanalysis (r-2). Bulletin of the American Meteorological Society 83: 1631-1644.
#'
#' @examples
#' tme <- as.POSIXlt(c(1:15) * 24 * 3600, origin = "2015-01-15", tz = 'UTC')
#' head(get_NCEP(50, -5, tme))
get_NCEP <- function(lat, long, tme, reanalysis2 = TRUE) {
  sorttimes <- function(tme) {
    tme2 <- as.POSIXct(tme)
    tme2 <- c((tme2 - 24 * 3600), tme2, (tme2 + 24 * 3600), (tme2 + 48 * 3600))
    tme2 <- as.POSIXlt(tme2)
    yrs <- unique(tme2$year + 1900)
    mths <- unique(tme2$mon + 1)
    tma <- 0
    dm <- c(31,28,31,30,31,30,31,31,30,31,30,31) * 4 - 1
    for (yr in min(yrs):max(yrs)) {
      dm[2] <- ifelse(yr%%4 == 0, 115, 111)
      for (mth in min(mths):max(mths)) {
        tmym <- as.POSIXct(c(0:dm[mth]) * 3600 * 6,
                           origin = paste0(yr,"-",mth,"-01 00:00"),
                           tz = 'UTC')
        tma <- c(tma, tmym)
      }
    }
    tma <- as.POSIXlt(tma[-1], origin = '1970-01-01', tz = 'UTC')
    sel <- which(tma >= min(tme2) & tma <= max(tme2))
    sel <- sel[-length(sel)]
    return(list(tme = tma, sel = sel))
  }
  ncepget1 <- function(climvar, tme2, ll, sel) {
    yrs <- unique(tme2$year + 1900)
    mths <- unique(tme2$mon + 1)
    v <- NCEP.gather(climvar, level = 'gaussian',
                     years.minmax = c(min(yrs),max(yrs)),
                     months.minmax = c(min(mths):max(mths)),
                     lat.southnorth = c(ll$y,ll$y), lon.westeast = c(ll$x,ll$x),
                     reanalysis2 = reanalysis2, return.units = FALSE, status.bar = FALSE)
    if (is.null(dim(v)) == F) {
      latdif <- ((ll$y - as.numeric(rownames(v)))^2)^0.5
      londif <- ((ll$x%%360 - as.numeric(colnames(v)))^2)^0.5
      v <- v[which.min(londif), which.min(latdif),]
    }
    v[sel]
  }
  ll <- data.frame(x = long, y = lat)
  tme2 <- sorttimes(tme)$tme
  sel <- sorttimes(tme)$sel
  # These variables are forecasts valid 6 hours after the reference time.
  Tk <- ncepget1('air.2m', tme2, ll, sel)
  Tkmin <- ncepget1('tmin.2m', tme2, ll, sel)
  Tkmax <- ncepget1('tmax.2m', tme2, ll, sel)
  sh <- ncepget1('shum.2m', tme2, ll, sel)
  pr <- ncepget1('pres.sfc', tme2, ll, sel)
  wu <- ncepget1('uwnd.10m', tme2, ll, sel)
  wv <- ncepget1('vwnd.10m', tme2, ll, sel)
  # These variables are 6 hour averages starting at the reference time.
  dlw <- ncepget1('dlwrf.sfc', tme2, ll, sel)
  ulw <- ncepget1('ulwrf.sfc', tme2, ll, sel)
  dsw <- ncepget1('dswrf.sfc', tme2, ll, sel)
  tcdc <- ncepget1('tcdc.eatm', tme2, ll, sel)
  dfout <- data.frame(obs_time = tme2[sel], Tk, Tkmin, Tkmax, sh, pr, wu, wv, dlw, ulw, dsw, tcdc)
  rownames(dfout) <- NULL
  return(dfout)
}

#' Calculates the solar index for a flat surface
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
#' # solar index at noon on 21 June 2010, Porthleven, Cornwall
#' jd <- julday (2010, 6, 21) # Julian day
#' siflat(12, 50.08, -5.31, jd)
#'
siflat <- function(localtime, lat, long, julian, merid = 0, dst = 0){
  saltitude <- solalt(localtime, lat, long, julian, merid, dst)
  alt <- saltitude * (pi/180)
  index <- cos(pi/2 - alt)
  index[index < 0] <- 0
  index
}
#' Calculates the diffuse fraction from incoming shortwave radiation (time series)
#' @export
.difprop2 <- function (rad, julian, localtime, lat, long, hourly = FALSE,
                       watts = TRUE, merid = 0, dst = 0)
{
  if (watts)
    rad <- rad * 0.0036
  lt <- c(localtime:(localtime + 5)) + 0.5
  jd <- rep(julian, 6)
  sel <- which(lt > 24)
  lt[sel] <- lt[sel] -24
  jd[sel] <- jd[sel] + 1
  sa <- solalt(lt, lat, long, jd, merid, dst)
  alt <- sa * (pi/180)
  k1 <- 0.83 - 0.56 * exp(-0.06 * sa)
  si <- cos(pi/2 - alt)
  k <- rad/(4.87 * si)
  k[is.na(k)] <- 0
  k <- ifelse(k > k1, k1, k)
  k[k < 0] <- 0
  rho <- k/k1
  sigma3a <- 0.021 + 0.397 * rho - 0.231 * rho^2 - 0.13 *
    exp(-1 * (((rho - 0.931)/0.134)^2)^0.834)
  sigma3b <- 0.12 + 0.65 * (rho - 1.04)
  sigma3 <- ifelse(rho <= 1.04, sigma3a, sigma3b)

  k2 <- 0.95 * k1
  d1 <- ifelse(sa > 1.4, 0.07 + 0.046 * (90 - sa)/(sa + 3),
               1)
  d1 <- mean(d1)
  K <- 0.5 * (1 + sin(pi * (k - 0.22)/(k1 - 0.22) - pi/2))
  d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 *
                           K^2))
  d3 <- (d2 * k2) * (1 - k)/(k * (1 - k2))
  alpha <- (1/sin(alt))^0.6
  kbmax <- 0.81^alpha
  kmax <- (kbmax + d2 * k2/(1 - k2))/(1 + d2 * k2/(1 - k2))
  dmax <- (d2 * k2) * (1 - kmax)/(kmax * (1 - k2))
  d4 <- 1 - kmax * (1 - dmax)/k
  d <- ifelse(k <= kmax, d3, d4)
  d <- ifelse(k <= k2, d2, d)
  d <- ifelse(k <= 0.22, 1, d)
  kX <- 0.56 - 0.32 * exp(-0.06 * sa)
  kL <- (k - 0.14)/(kX - 0.14)
  kR <- (k - kX)/0.71
  delta <- ifelse(k >= 0.14 & k < kX, -3 * kL^2 * (1 - kL) *
                    sigma3^1.3, 0)
  delta <- ifelse(k >= kX & k < (kX + 0.71), 3 * kR * (1 -
                                                         kR)^2 * sigma3^0.6, delta)
  d[sigma3 > 0.01] <- d[sigma3 > 0.01] + delta[sigma3 > 0.01]
  d[si <= 0] <- NA
  d[d > 1] <- 1
  mean(d, na.rm = T)
}
#' Interpolate NCEP data for running microclima to hourly
#'
#' @description `hourlyNCEP` optionally downloads the required NCEP climate and radiation forcing data
#' required for running microclima and interpolates 4x daily data to hourly.
#'
#' @param ncepdata an optional  data frame of climate variables as returned by [get_ncep()].
#' @param tme a POSIXlt object covering the duration for which data are required. Ignored if `ncepdata``
#' provided.
#' Hourly data are returned irrespective of the time interval of `tme`.Ignored if `ncepdata` provided.
#' @param lat the latitude of the location for which data are required. Ignored if `ncepdata` provided.
#' @param long the longitude of the location for which data are required. Ignored if `ncepdata` provided.
#' @param reanalysis2 Logical. Should data be obtained from the Reanalysis II dataset (default) or
#' from Reanalysis I (data prior to 1979). Ignored if `ncepdata` provided.
#'
#' @return a dataframe with the following variables: (i) obs_time: a POSIXlt oobject of times,
#' (ii) temperature: temperatures at 2m (ºC), (iii) humidity: specific humidity at 2m (Kg / Kg),
#' (iv) pressure: surface pressure (Pa), (v) windspeed: wind speed at 2m (metres per
#' second), (vi) winddir: wind direction (degrees from N), (vii) emissivity: emissivity of the
#' atmosphere (0 - 1), (viii) netlong: Net longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})
#' (ix) uplong: upward longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}}),
#' (x) downlong: downward longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}}),
#' (xi) rad_dni: Direct radiation normal to the solar beam (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}}),
#' (xii) rad_dif: Diffuse radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}}),
#' (xiii) szenith: the zenith angle (degrees), (ixx) cloudcover: cloud cover (%).

#' @export
#' @import zoo
#'
#' @seealso [get_NCEP()]
#' @details If `ncepdata` is not provided, then [get_ncep()] is called and data are downloaded from NCEP Atmospheric
#' Model Intercomparison Project (Kanamitso et al 2002). Six-hourly data are interpolated as follows.
#' Pressure, humidity and the u and v wind vectors are converted to hourly using spline interpolation. Wind speeed and diretcion and then
#' calculated and adjusted to give values at 1 m using [windheight()]. The diffuse radiation
#' proportion is calculated using [difprop()], and hourly obtained by adusting for the direction
#' of the solar beam and airmass thickness using functions [siflat()] and [airmasscoef()].
#' Hourly temperature are derived using [hourlytemp()] and longwave radiation by splining emissivity
#' values.
#' @examples
#' tme <- as.POSIXlt(c(0:30) * 24 * 3600, origin ="2015-01-15 00:00", tz = "UTC")
#' # NB takes a while to download data
#' hdata<- hourlyNCEP(NA, 50, -5, tme)
#' head(hdata)
#' plot(temperature ~ as.POSIXct(obs_time), data = hdata, type = "l", xlab = "Month")
#' plot(rad_dif ~ as.POSIXct(obs_time), data = hdata, type = "l", xlab = "Month")
hourlyNCEP <- function(ncepdata = NA, lat, long, tme, reanalysis2 = TRUE) {
  bound <- function(x, mn = 0, mx = 1) {
    x[x > mx] <- mx
    x[x < mn] <- mn
    x
  }
  if (class(ncepdata) != "logical") {
    tme <- as.POSIXlt(ncepdata$obs_time)
    tme <- tme[5:(length(tme) - 4)]
  }
  int <- as.numeric(tme[2]) - as.numeric(tme[1])
  lgth <- (length(tme) * int) / (24 * 3600)
  tme2 <- as.POSIXlt(c(0:(lgth - 1)) * 3600 * 24, origin = min(tme), tz = 'UTC')
  if (class(ncepdata) == "logical") {
    ncepdata <- get_NCEP(lat, long, tme2, reanalysis2)
  }
  # *** NB sort out times (take form tme6 rather than tme)
  tme6 <- as.POSIXlt(ncepdata$obs_time)
  n <- (length(tme6) - 1) * 6 + 1
  h_pr <- spline(tme6, ncepdata$pr, n = n)$y
  h_sh <- spline(tme6, ncepdata$sh, n = n)$y
  h_cc <- bound(spline(tme6, ncepdata$tcdc, n = n)$y, mx = 100)
  tmo  <- spline(tme6, ncepdata$sh, n = n)$x
  tmo <- as.POSIXlt(tmo, origin = "1970-01-01 00:00",
                    tz = "UTC")
  tmo2 <- as.POSIXlt(seq(0,(length(tme6) / 4), by = 1/24) * 3600 * 24, origin = min(tme6), tz = 'UTC')
  jd <- julday(tmo2$year + 1900,  tmo2$mon + 1, tmo2$mday)
  si <- siflat(tmo2$hour, lat, long, jd)
  am <- airmasscoef(tmo2$hour, lat, long, jd)
  si_m <- 0
  am_m <- 0
  for (i in 1:(length(tme6))) {
    st <- (i - 1) * 6 + 1
    ed <- st + 5
    si_m[i] <- mean(si[st:ed], na.rm = T)
    am_m[i] <- mean(am[st:ed], na.rm = T)
  }
  # Calculate diffuse proportion
  jd <- julday(tme6$year + 1900, tme6$mon + 1, tme6$mday)
  dp <- 0
  for (i in 1:length(jd)) {
    dp[i] <- .difprop2(ncepdata$dsw[i], jd[i], tme6$hour[i], lat, long)
  }
  am_m[am_m > 5] <- 5
  # Adjust direct and diffuse
  dp[ncepdata$dsw == 0] <- NA
  dnir <- (ncepdata$dsw * (1 - dp)) / si_m
  dnir[si_m == 0] <- NA
  difr <- (ncepdata$dsw * dp)
  rmx <- 4.87 / 0.0036 * (1 - dp)
  rmx2 <- 4.87 / 0.0036 * (dp)
  sdni <- dnir - rmx
  sdni[sdni < 0] <- 0
  sdni[si_m < 0.0348995] <- dnir[si_m < 0.0348995]
  dnir <- dnir - sdni
  difr <- difr + sdni
  # Calculate extinction coefficient
  edni <- dnir / ((4.87 / 0.0036) * (1 - dp))
  edif <- difr / ((4.87 / 0.0036) * dp)
  #Calculate optical depths
  odni <- bound((log(edni) / -am_m))
  odif <- bound((log(edif) / -am_m))
  dnir <- dnir / exp(-am_m * odni)
  difr <- difr / exp(-am_m * odif)
  dnir <- ifelse(dnir > rmx, rmx, dnir)
  difr <- ifelse(difr > rmx2, rmx2, difr)
  # Interpolate to hourly
  globr <- dnir * si_m + difr
  globr <- bound(globr, mx = 4.87 / 0.0036)
  nd <- length(globr)
  sel <-which(is.na(globr * dp * odni * odif) == F)
  globr[1] <- globr[min(sel)]
  dp[1] <- dp[min(sel)]
  odni[1] <- odni[min(sel)]
  odif[1] <- odif[min(sel)]
  globr[nd] <- globr[max(sel)]
  dp[nd] <- dp[max(sel)]
  odni[nd] <- odni[max(sel)]
  odif[nd] <- odif[max(sel)]
  globr <- na.approx(globr, na.rm = F)
  dp <- na.approx(dp, na.rm = F)
  odni <- na.approx(odni, na.rm = F)
  odif <- na.approx(odif, na.rm = F)
  h_gr <- bound(spline(tme6, globr, n = n)$y, mx = 4.87 / 0.0036)
  h_dp <- bound(spline(tme6, dp, n = n)$y)
  h_oi <- bound(spline(tme6, odni, n = n)$y)
  h_od <- bound(spline(tme6, odif, n = n)$y)
  tmorad <- as.POSIXlt(tmo + 3600 * 3)
  jd <- julday(tmorad$year + 1900, tmorad$mon + 1, tmorad$mday)
  szenith <- 90 - solalt(tmorad$hour, lat, long, jd)
  si <- siflat(tmorad$hour, lat, long, jd)
  am <- airmasscoef(tmorad$hour, lat, long, jd)
  afi <- exp(-am * h_oi)
  afd <- exp(-am * h_od)
  h_dni <- (1 - h_dp) * afi * h_gr
  h_dif <- h_dp * afd * h_gr
  h_dni[si == 0] <- 0
  h_dif[is.na(h_dif)] <- 0
  em <- ncepdata$dlw / ncepdata$ulw
  em_h <- spline(tme6, em, n = n)$y
  t6sel <- c(4:(length(tme6)-5))
  thsel <-c(22:(length(tmorad)-22))
  thsel2 <-c(19:(length(tmo) - 25))
  tcmn <- ncepdata$Tkmin[t6sel] - 273.15
  tcmx <- ncepdata$Tkmax[t6sel] - 273.15
  tcmn <-t(matrix(tcmn, nrow = 4))
  tcmx <-t(matrix(tcmx, nrow = 4))
  tmin <- apply(tcmn, 1, min)
  tmax <- apply(tcmx, 1, max)
  jd <- julday(tme2$year + 1900, tme2$mon + 1, tme2$mday)
  h_tc<-suppressWarnings(hourlytemp(julian = jd, em = em_h[thsel], dni = h_dni[thsel],
                   dif = h_dif[thsel], mintemp = tmin, maxtemp = tmax,
                   lat = lat, long = long))
  hlwu <- 2.043e-10 * (h_tc + 273.15)^4
  hlwd <- em_h[thsel] *  hlwu
  h_nlw <- hlwu - hlwd
  h_uw <- spline(tme6, ncepdata$wu, n = n)$y
  h_vw <- spline(tme6, ncepdata$wv, n = n)$y
  h_ws <- sqrt(h_uw^2 + h_vw^2)
  h_ws <- windheight(h_ws, 10, 1)
  h_wd <- atan2(h_uw, h_vw) * 180/pi + 180
  h_wd <- h_wd%%360
  hourlyout <- data.frame(obs_time = tmorad[thsel], temperature = h_tc,
                          humidity = h_sh[thsel2], pressure = h_pr[thsel2],
                          windspeed = h_ws[thsel2], winddir = h_wd[thsel2],
                          emissivity = em_h[thsel2], cloudcover = h_cc[thsel],
                          netlong = h_nlw, uplong = hlwu, downlong = hlwd,
                          rad_dni = h_dni[thsel] * 0.0036,
                          rad_dif = h_dif[thsel] * 0.0036,
                          szenith = szenith[thsel])
  return(hourlyout)
}
#' Obtain daily precipitation from NCEP
#'
#' @param tme a POSIXlt object covering the duration for which data are required.
#' Intervals should be daily.
#' Hourly data are returned irrespective of the time interval of `tme`.
#' @param lat the latitude of the location for which data are required
#' @param long the longitude of the location for which data are required
#' @param reanalysis2 Logical. Should data be obtained from the Reanalysis II dataset (default) or
#' from Reanalysis I (data prior to 1979).
#'
#' @return a vector of daily precipitations (mm / day) for the period covered by tme
#' @export
#'
#' @examples
#' tme <- as.POSIXlt(c(1:15) * 24 * 3600, origin = "2015-01-15", tz = 'UTC')
#' dailyprecipNCEP(50, -5, tme)
dailyprecipNCEP <- function(lat, long, tme, reanalysis2 = TRUE) {
  int <- as.numeric(tme[2]) - as.numeric(tme[1])
  lgth <- (length(tme) * int) / (24 * 3600)
  tme <- as.POSIXlt(c(0:(lgth - 1)) * 3600 * 24, origin = min(tme), tz = 'UTC')
  yrs <- unique(tme$year + 1900)
  mths <- unique(tme$mon + 1)
  ll <- data.frame(x = long, y = lat)
  pre <- NCEP.gather('prate.sfc', level = 'gaussian',
                     years.minmax = c(min(yrs),max(yrs)),
                     months.minmax = c(min(mths):max(mths)),
                     lat.southnorth = c(ll$y,ll$y), lon.westeast = c(ll$x,ll$x),
                     return.units = FALSE, status.bar = FALSE, reanalysis2 = reanalysis2)
  if (is.null(dim(pre)) == F) {
    latdif <- ((ll$y - as.numeric(rownames(pre)))^2)^0.5
    londif <- ((ll$x%%360 - as.numeric(colnames(pre)))^2)^0.5
    pre <- pre[which.min(londif), which.min(latdif),]
  }
  pre <- pre * 6 * 3600
  tma <- 0
  dm <- c(31,28,31,30,31,30,31,31,30,31,30,31) * 4 - 1
  for (yr in min(yrs):max(yrs)) {
    dm[2] <- ifelse(yr%%4 == 0, 115, 111)
    for (mth in min(mths):max(mths)) {
      tmym <- as.POSIXct(c(0:dm[mth]) * 3600 * 6,
                         origin = paste0(yr,"-",mth,"-01 00:00"),
                         tz = 'UTC')
      tma <- c(tma, tmym)
    }
  }
  tma <- as.POSIXlt(tma[-1], origin = '1970-01-01', tz = 'UTC')
  sel <- which(tma >= min(tme) & tma < (max(tme) + 24 * 3600))
  pre <- pre[sel]
  dpre <- t(matrix(pre, nrow = 4))
  dpre <- apply(dpre, 1, sum)
  return(dpre)
}
#' get radiation and wind for use with NicheMapR
#' @export
.pointradwind <- function(hourlydata, dem, lat, long, l, x, albr = 0.15, zmin = 0, slope = NA, aspect = NA) {
  m <- is_raster(dem)
  m[is.na(m)] <- zmin
  m[m < zmin] <- zmin
  dem <- if_raster(m, dem)
  xy <- data.frame(x = long, y = lat)
  coordinates(xy) = ~x + y
  proj4string(xy) = "+init=epsg:4326"
  xy <- as.data.frame(spTransform(xy, crs(dem)))
  reso <- res(dem)[1]
  wsc36 <- 0
  wsc36atground <- 0
  for (i in 0:35) {
    wscr <- windcoef(dem, i*10, hgt = 1, res = reso)
    wscr2 <- windcoef(dem, i*10, hgt = 0, res = reso)
    wsc36[i + 1] <- extract(wscr, xy)
    wsc36atground[i + 1] <- extract(wscr2, xy)
  }
  wshelt <- 0
  wsheltatground <- 0
  hourlydata$winddir <- hourlydata$winddir%%360
  for (i in 1:length(hourlydata$winddir)) {
    daz <- round(hourlydata$winddir[i] / 10, 0) + 1
    daz[daz > 36] <- daz[daz > 36] - 36
    wshelt[i] <- wsc36[daz]
    wsheltatground[i] <- wsc36atground[daz]
  }
  if (class(slope) == "logical") {
    slope <- terrain(dem, unit = 'degrees')
    slope <- extract(slope, xy)
  }
  if (class(aspect) == "logical") {
    aspect <- terrain(dem, opt = 'aspect', unit = 'degrees')
    aspect <- extract(aspect, xy)
  }
  svf <- skyviewveg(dem, array(l, dim = dim(dem)[1:2]),
                    array(x, dim = dim(dem)[1:2]), res = reso)
  fr <- canopy(array(l, dim = dim(dem)[1:2]), array(x, dim = dim(dem)[1:2]))
  svf <- extract(svf, xy)
  fr <- mean(fr)
  mslope <- mean_slope(dem, res = reso)
  mha <- extract(mslope, xy)
  ha36 <- 0
  for (i in 0:35) {
    har <- horizonangle(dem, i*10, reso)
    ha36[i + 1] <- atan(extract(har, xy)) * (180/pi)
  }
  tme <- as.POSIXlt(hourlydata$obs_time)
  ha <- 0
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  for (i in 1:length(tme)) {
    saz <- solazi(tme$hour[i], lat, long, jd[i])
    saz <- round(saz / 10, 0) + 1
    saz <- ifelse(saz > 36, 1, saz)
    ha[i] <- ha36[saz]
  }
  si <- siflat(tme$hour, lat, long, jd)
  sa <- solalt(tme$hour, lat, long, jd)
  si[sa < ha] <- 0
  dirr <- si * hourlydata$rad_dni
  a <- slope * (pi/180)
  k <- hourlydata$rad_dni / 4.87
  k <- ifelse(k > 1, 1, k)
  isor <- 0.5 * hourlydata$rad_dif * (1 + cos(a)) * (1 - k)
  cisr <- k * hourlydata$rad_dif * si
  sdi <- (slope + mha) * (pi/180)
  refr <- 0.5 * albr * (1 - cos(sdi)) * hourlydata$rad_dif
  fd <- dirr + cisr
  fdf <- isor + refr
  kk <- ((x^2 + 1/(tan(sa * (pi/180))^2))^0.5)/(x + 1.774 * (x + 1.182)^(-0.733))
  trd <- exp(-kk * l)
  trf <- (1 - fr)
  fgd <- fd * trd
  fged <- fdf * trf * svf
  swrad <- fgd + fged
  windsp <- windheight(hourlydata$windspeed, 10, 1)
  cfc <- ((1 - trd) * fd + fr * fdf) / (fd + fdf)
  cfc[is.na(cfc)] <- ((1 - trd[is.na(cfc)]) * 0.5 + fr * 0.5) / (0.5 + 0.5)
  hourlyrad <- data.frame(swrad = swrad, skyviewfact = svf, canopyfact = cfc,
                          whselt = wsheltatground, windspeed = wshelt *  windsp,
                          slope = slope, aspect = aspect)
  hourlyrad
}
#' Cold air drainage direct from emissivity
#' @export
.cadconditions2 <- function (em, wind, startjul, lat, long,
                             starttime = 0, hourint = 1, windthresh = 4.5, emthresh = 0.725,
                             tz = 0, dst = 0, con = TRUE)
{
  jd <- floor(c(1:length(em)) * hourint/24 - hourint/24 + startjul +
                starttime/24)

  st <- suntimes(jd, lat, long, tz, dst)
  hrs <- (c(1:length(em)) * hourint - hourint + starttime)%%24
  dn <- ifelse(hrs > (st$sunrise + 3) & hrs < st$sunset, 1,
               0)
  cad <- ifelse(dn < 1 & wind < windthresh & em < emthresh,
                1, 0)
  if (hourint == 1 & con) {
    cad[2] <- ifelse(cad[1] & cad[2] == 1, 1, 0)
    cad <- c(cad[1:2], cad[3:length(cad)] * cad[2:(length(cad) -
                                                     1)] * cad[1:(length(cad) - 2)])
  }
  cad
}
# Calculates elevation effects
.eleveffects <- function(hourlydata, dem, lat, long, windthresh = 4.5,
                         emthresh = 0.78) {
  load(demworld)
  xy <- data.frame(x = long, y = lat)
  elevncep <- extract(demworld, xy)
  coordinates(xy) = ~x + y
  proj4string(xy) = "+init=epsg:4326"
  xy <- as.data.frame(spTransform(xy, crs(dem)))
  elev <- extract(dem, xy)
  # elevation effect
  lr <- lapserate(hourlydata$temperature, hourlydata$humidity,
                  hourlydata$pressure)
  elevt <- lr * (elev - elevncep) + hourlydata$temperature
  tme <- as.POSIXlt(hourlydata$obs_time)
  jds <- julday(tme$year[1] + 1900, tme$mday[1] + 1, tme$mday[1])
  cad <- .cadconditions2(hourlydata$emissivity, hourlydata$windspeed,
                         jds, lat, long, starttime = tme$hour[1], hourint = 1,
                         windthresh = windthresh, emthresh = emthresh)
  pxls <- dim(dem)[1] * dim(dem)[2]
  if (pxls > 300 * 300) {
    basins <- basindelin_big(dem)
  } else {
    basins <- basindelin(dem)
  }
  basins <- basinmerge(dem, basins, 2)
  basins <- basinsort(dem, basins)
  fa <- flowacc(dem)
  pfa <- is_raster(fa) * 0
  td <- is_raster(fa) * 0
  bm <- is_raster(basins)
  dm <- is_raster(dem)
  for (b in 1:max(bm, na.rm = TRUE)) {
    sel <- which(bm == b)
    fao <- log(is_raster(fa)[sel])
    pfa[sel] <- fao/max(fao, na.rm = TRUE)
    ed <- max(dm[sel], na.rm = TRUE) - dm[sel]
    td[sel] <- ed
  }
  cdif <- pfa * td
  cdif <- if_raster(cdif, dem)
  cdif <- extract(cdif, xy)
  cadt <- cdif * lr * cad
  tout <- data.frame(tref = hourlydata$temperature,
                     elev = elev, elevncep = elevncep,
                     telev = lr * elev,
                     tcad = cadt)
  tout
}
#' Calculates coastal exposure automatically
.invls.auto <- function(r, steps = 8, use.raster = T, zmin = 0, plot.progress = T) {
  adjust.lsr <- function(lsr, rs) {
    m <- is_raster(lsr)
    m[m < 0] <- 0
    s <- c((8:10000) / 8) ^ 2 * 30
    sl <- which(s <= rs)
    mval <- length(sl) / 1453
    m2 <- m * (1 - mval) + mval
    m2
  }
  ll <- latlongfromraster(r)
  xs <- seq(0, by = 1.875, length.out = 192)
  ys <- seq(-88.542, by = 1.904129, length.out = 94)
  xc <-xs[which.min(abs(xs - ll$long%%360))]
  yc <-ys[which.min(abs(ys - ll$lat))]
  rll <- raster(matrix(0, nrow = 3, ncol = 3))
  extent(rll) <- extent(xc -  2.8125, xc +  2.8125,
                        yc - 2.856193, yc + 2.856193)
  crs(rll) <- "+init=epsg:4326"
  rmer <- projectRaster(rll, crs = "+init=epsg:3395")
  ress <- c(30, 90, 500, 1000, 10000)
  ress <- ress[ress >= mean(res(r))]
  ress <- rev(ress)
  # Create a list of dems
  cat("Downloading land sea data \n")
  dem.list <- list()
  rmet <- projectRaster(rll, crs = crs(r))
  dem <- get_dem(rmet, resolution = ress[1], zmin = zmin)
  dem.list[[1]] <- projectRaster(dem, crs = crs(r))
  dc <- ceiling(max(dim(r)) / 2)
  rres <- mean(res(r)[1:2])
  for (i in 2:length(ress)) {
    d <- 50
    tst <- rres / ress[i] * dc * 2
    if (ress[i] <= mean(res(r)[1:2]))
      d <- dc
    if (tst > 50) d <- tst
    e <- extent(r)
    xy <- c((e@xmin + e@xmax) / 2, (e@ymin + e@ymax) / 2)
    e2 <- extent(c(xy[1] - d * ress[i], xy[1] + d * ress[i],
                   xy[2] - d * ress[i], xy[2] + d * ress[i]))
    rr <- raster()
    extent(rr) <- e2
    crs(rr) <- crs(r)
    dem <- suppressWarnings(get_dem(rr, resolution = ress[i], zmin = zmin))
    dem.list[[i]] <- projectRaster(dem, crs = crs(r))
  }
  cat("Computing coastal exposure \n")
  lsa.array <- array(NA, dim = c(dim(r)[1:2], steps))
  for (dct in 0:(steps - 1)) {
    direction <- dct * (360 / steps)
    cat(paste("Direction:", direction, "degrees"), "\n")
    lsa.list <- list()
    for (i in 1:length(ress)) {
      dem <- dem.list[[i]]
      m <- is_raster(dem)
      m[m == zmin] <- NA
      dem <- if_raster(m, dem)
      lsa <- invls(dem, extent(r), direction)
      lsm <- is_raster(lsa)
      if (max(lsm, na.rm = T) > min(lsm, na.rm = T)) {
        lsa.list[[i]] <- resample(lsa, r)
      } else  {
        lsa.list[[i]] <- if_raster(array(mean(lsm, na.rm = T),
                                         dim = dim(r)[1:2]), r)
      }
    }
    # find min vals
    for (i in 1:length(ress))
      lsa.list[[i]] <- adjust.lsr(lsa.list[[i]], ress[i])
    lsa <- is_raster(lsa.list[[1]])
    for (i in 2:length(ress))
      lsa <- pmin(lsa, is_raster(lsa.list[[i]]))
    lsa <- if_raster(lsa, r)
    if (use.raster) lsa <- mask(lsa, r)
    if (plot.progress) plot(lsa, main = paste("Direction:",direction))
    lsa.array[,,dct + 1] <- is_raster(lsa)
  }
  # Compute land-sea ratio of ncep grid cell
  cat("Computing mean coastal exposure of ncep grid cell \n")
  cncep <-0
  eone <- extent(xc - 1.875 / 2, xc + 1.875 / 2,
                 yc - 1.904129 / 2, yc + 1.904129 / 2)
  rllo <- crop(rll, eone)
  rone <- projectRaster(rllo, crs = crs(r))
  dem <- dem.list[[1]]
  m <- is_raster(dem)
  m[m == zmin] <- NA
  dem <- if_raster(m, dem)
  for (dct in 0:(steps - 1)) {
    direction <- dct * (360 / steps)
    lsa <- invls(dem, extent(rone), direction)
    cncep[dct + 1] <- mean(is_raster(lsa), na.rm = T)
  }
  cat("Adjusting coastal exposure by ncep means \n")
  for (i in 1:(steps-1)) lsa.array[,,i] <- lsa.array[,,i] / cncep[i]
  return(lsa.array)
}
#' Downloads sea-surface temperature data
.get_sst <- function(lat, long, tme) {
  sel1 <- which(tme$year == min(tme$year))
  sel2 <- which(tme$year == max(tme$year))
  dms <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  amonth <- min(tme$mon[sel1])
  amonth <- ifelse(amonth == 0, 12, amonth)
  pmonth <- max(tme$mon[sel1]) + 2
  pmonth <- ifelse(pmonth > 12, pmonth -12, pmonth)
  dm <- dms[max(tme$mon[sel2]) + 1]
  omn <- paste0(min(tme$year) + 1900, "-", min(tme$mon[sel1]) + 1, "-01 00:00:00")
  omx <- paste0(max(tme$year) + 1900, "-", max(tme$mon[sel2]) + 1, "-", dm, " 23:59:59")
  tmn <- as.POSIXlt(0, origin = omn, tz = "GMT") - dms[amonth] / 2 * 24 * 3600
  tmx <- as.POSIXlt(1, origin = omx, tz = "GMT") + dms[pmonth] / 2 * 24 * 3600
  tme2 <-  as.POSIXlt(seq(tmn, tmx, by = 'hour'))
  yrs <- unique(tme2$year) + 1900
  tmes <- as.POSIXlt(0, origin = "1900-01-01 00:00", tz = "GMT")
  sstdf <- data.frame(obs_time = tmes, sst = -99)
  for (yr in 1:length(yrs)) {
    dms[2] <- 28
    if (yrs[yr]%%4 == 0) dms[2] <- 29
    if (yrs[yr]%%100 == 0 & yrs[yr]%%400 != 0) dms[2] <- 28
    sel <- which(tme2$year + 1900 == yrs[yr])
    mths <- unique(tme2$mon[sel]) + 1
    for (mt in 1:length(mths)) {
      nc <- ersst(yrs[yr], mths[mt])
      sst <- ncvar_get(nc, 'sst')
      sst <- t(sst)
      sst <- apply(sst, 2, rev)
      sst <- na.approx(sst, na.rm = F)
      sstr <- raster(sst)
      extent(sstr) <- extent(-1, 359, -89, 89)
      long2 <- long%%360
      long2 <- ifelse(long2 > 359, long2 - 360, long2)
      s <- extract(sstr, data.frame(long2, lat))
      orn <- paste0(yrs[yr], "-", mths[mt], "-01 00:00")
      tmem <- as.POSIXlt(0 + 3600 * 24 * dms[mths[mt]] / 2, origin = orn, tz = "GMT")
      sstdo <- data.frame(obs_time = tmem, sst = s)
      sstdf <- rbind(sstdf, sstdo)
    }
  }
  sstdf <- sstdf[-1,]
  ssth <- spline(sstdf$sst~as.POSIXct(sstdf$obs_time), n = length(tme2))$y
  tmeh <- spline(sstdf$sst~as.POSIXct(sstdf$obs_time), n = length(tme2))$x
  tmeh <- as.POSIXlt(tmeh, origin = "1970-01-01 00:00", tz = "GMT")
  sel <- which(tmeh >=min(tme) & tmeh <=max(tme))
  return(ssth[sel])
}
#' Calculates coastal effects from NCEP data
#'
#' @param landsea landsea A raster object with NAs (representing sea) or any non-NA value and a projection system defined.
#' @param ncephourly a dataframe of hourly climate data as returned by [hourlyNCEP()].
#' @param steps  an optional integer. Coastal effects are calculated in specified directions upwind. Steps defines the total number of directions used. If the default 8 is specified, coastal effects are calculated at 45º intervals.
#' @param use.raster an optional logical value indicating whether to mask the output values by `landsea`.
#' @param zmin optional assumed sea-level height. Values below this are set to zmin
#' @param plot.progress logical value indicating whether to produce plots to track progress.
#'
#' @details coastalNCEP downloads digital elevation and varying resolutions to calculate
#' coastal effects by applying a variant of [invls()] over a greater extent then that specified by landsea, but the resulting outputs are
#' cropped to the same extent as landsea. Land temperatyres as a function of coastal exposure
#' and sea-surface temperature data are then calculated. Sea surface temperature data are
#' automatically downloaded from the NOAA.
#'
#' @return a three-dimension array of hourly temperature values for each pixel of `landsea`
#' @export
#' @import raster
#' @import sp
#' @import zoo
#' @import rnoaa
#' @import ncdf4
#' @examples
#' # Download NCEP data
#' ll <- latlongfromraster(dtm100m)
#' tme <-as.POSIXlt(c(0:31) * 3600 * 24, origin = "2015-03-15", tz = "GMT")
#' ncephourly<-hourlyNCEP(NA, ll$lat, ll$long, tme)
#' aout <- coastalNCEP(dtm100m, ncephourly)
#' # Calculate mean temperature and convert to raster
#' mtemp <- if_raster(apply(aout, c(1, 2), mean), dtm100m)
#' plot(mtemp, main = "Mean temperature")
coastalNCEP <- function(landsea, ncephourly, steps = 8, use.raster = T, zmin = 0, plot.progress = T) {
  bound <- function(x, mn = 0, mx = 1) {
    x[x < mn] <- mn
    x[x > mx] <- mx
    x
  }
  lsr1 <- .invls.auto(landsea, steps, use.raster, zmin, plot.progress)
  lsr2 <- lsr1
  for (i in 1:steps) {
    mn <- i - 1
    mn <- ifelse(mn == 0, 8, mn)
    mx <- i + 1
    mx <- ifelse(mx >8, mx - 8, mx)
    lsr2[,,i] <- 0.25 * lsr1[,,mn]  + 0.5 * lsr1[,,i] + 0.25 * lsr1[,,mx]
  }
  lsrm <- apply(lsr1, c(1,2), mean)
  ll <- latlongfromraster(landsea)
  tme <- as.POSIXlt(ncephourly$obs_time)
  cat("Downloading sea-surface temperature data \n")
  sst <- .get_sst(ll$lat, ll$long, tme)
  wd <- round(ncephourly$winddir%%360 / (360 / 8), 0) + 1
  wd <- ifelse(wd > 8, wd - 8, wd)
  lws <- log(ncephourly$windspeed)
  b1 <- 11.003 * lws - 9.357
  b1[b1 < 0] <- 0
  b2 <- 0.6253 * lws - 3.5185
  dT <- sst - ncephourly$temperature
  p1 <- 1.12399 - 0.39985 * dT - 0.73361 * lws
  p2 <- -0.011142 + 0.01048 * dT + 0.043311 * lws
  d2 <- bound((1 - lsrm + 2) / 3)
  aout <- array(NA, dim = c(dim(landsea)[1:2], length(tme)))
  cat("Applying coastal effects \n")
  for (i in 1:length(tme)) {
    d1 <- bound((1 - lsr2[,,wd[i]] + 2) / 3)
    pdT <- dT[i] + p1[i] * d1^b1[i] + p2[i] * d2^b2[i]
    aout[,,i] <- (pdT - sst[i]) * (-1)
    if (plot.progress == T & i%%500 == 0) {
      plot(if_raster(aout[,,i], landsea), main = tme[i])
    }
  }
  return(aout)
}

#' Extract NCEP data and run microclima for use with NicheMapR
#'
#' @param lat the latitude (in decimal degrees) of the location for point data are required.
#' @param long the longitude (in decimal degrees) of the location for point data are required.
#' @param dstart start date as character of time period required in format DD/MM/YYYY
#' @param dfinish end date as character of time period required in format DD/MM/YYYY
#' @param l  a single numeric value of the leaf area index for location (see details).
#' @param x  a single numeric value representing the ratio of vertical to horizontal
#' projections of leaf foliage for the location (see details).
#' @param hourlydata an optional data frame of ourly weather and radiation data as returned
#' by function [hourlyNCEP()]. Function called if data not provided.
#' @param dailyprecip an optional data frame of ourly weather and radiation data as returned
#' by function [dailyprecipNCEP()]. Function called if data not provided.
#' @param dem optional raster object of elevations covering the location for which point
#' data are required. If not provided, one is downloaded from the registry of Open
#' Data on AWS.
#' @param albr albedo of ground surface from which radiation is reflected (see details)
#' @param resolution the resolution (in metres) of the dem used for downscaling (see details).
#' @param zmin assumed sea-level height. Values below this are set to zmin (see details).
#' @param slope an optional slope value for the location in decimal degrees. If not specified, then
#' the slope is calculated from the retrieved dem.
#' @param aspect an optional aspect value for the location in decimal degrees. If not specified, then
#' the slope is calculated from the retrieved dem.
#' @param windthresh an optional threshold wind value (m /s) above which cold air drainage is
#' assumed not to occur
#' @param emthresh an optional threshold emissivity value above which cold air drainage is
#' assumed not to occur
#' @param reanalysis2 Logical. Should data be obtained from the Reanalysis II dataset (default) or
#' from Reanalysis I (data prior to 1979).
#'
#' @return a list with the following objects:
#'
#' (1) hourlydata: a data frame of hourly climate data from NCEP,
#' as returned by [hourlyNCEP()].
#'
#' (2) hourlyradwind: A data frame with seven columns:
#' (i) `swrad`: total downward shortwave radiation as a funcion of slope, aspect and
#' canopy shading (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})
#' as returned by [shortwaveveg()].(ii) `skyviewfact`: a topographic sky view correction factor
#' as returned by [skyviewveg()]. (iii) `canopyfact` a canopy shading factor weighted by the ratio of direct to diffuse radiation.
#' (iv) `whselt`: a ground level wind shelter coefficient as returned by
#' [windcoef()]. (v) `windspeed`: wind speed (m /s) at 1 m, corrected for topographic
#' sheltering. (vi) `slope`: the slope (decimal degrees) of the location specified
#' by `lat` and `long`. Either the same as that input, or determined from `dem`.
#' (vii) `aspect`: the aspect (decimal degrees) of the location specified by `lat` and
#' `long`. Either the same as that input, or determined from `dem`.
#'
#' (3) `tref`: a dataframe with five columns (i) `tref` the reference temperature, derived from
#' NCEP data. (ii) `elev` the elevation of the location. (iii) `elevncep`: the mean elevation
#' of the 2.5º grid cell corresponding to that from which NCEP data were derived. (iv)
#' `telev` the expected difference in temperature due to elevation, calculated by
#' applying [lapserate()]. (v) `tcad` the expected difference in temperature due to cold
#' air drainage effects, determined using [pcad()] and [cadconditions()].
#'
#' (4) `dailyprecip` a vector of daily rainfalls derived from NCEP.
#'
#' @details The package `NicheMapR` (Kearney & Porter 2016), includes a suite of
#' programs for mechanistic modelling of heat and mass exchange and provides a
#' time-series of hourly estimates of above- and below‐ground conditions. It has the
#' significant advantage over `microclima` in that it works entirely from first
#' principles does not require measurements of temperature for calibration, but does
#' not explicitly model elevation and cold-air drainage effects. The slope, aspect
#' and canopy effects must be specified by the user, and their effect on radiation
#' are more simplistically modelled than by `microclima`.
#'
#' @references Kearney MR,  Porter WP (2016) NicheMapR – an R package for biophysical
#' modelling: the microclimate model. Ecography 40: 664-674.
#'
#' @seealso Function [get_dem()] for retrieving digital elevation data and
#' [hourlyNCEP()] for retrieving and processing hourly data from NCEP.

#' @export
#' @import raster
#' @import sp
#' @import zoo
#'
#' @examples
#' mnr <- microclimaforNMR(50, -5.2, '15/01/2015', '15/02/2015', 1, 1)
#' head(mnr$hourlydata)
#' head(mnr$hourlyradwind)
#' head(mnr$tref)
#' head(mnr$dailyprecip)
microclimaforNMR <- function(lat, long, dstart, dfinish, l, x, hourlydata = NA,
                                   dailyprecip = NA, dem = NA, albr =0.15,
                                   resolution = 30, zmin = 0, slope = NA, aspect = NA,
                                   windthresh = 4.5, emthresh = 0.78, reanalysis2 = TRUE) {
  tme <- seq(as.POSIXlt(dstart, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'),
             as.POSIXlt(dfinish, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'),
             by = 'days')
  if (class(dem) == "logical") {
    cat("Downloading digital elevation data \n")
    dem <- get_dem(r = NA, lat = lat, long = long, resolution = resolution, zmin = zmin)
  }
  if (class(hourlydata) == "logical") {
    cat("Extracting climate data from NCEP \n")
    hourlydata <- hourlyNCEP(ncepdata = NA, lat, long, tme, reanalysis2)
  }
  if (class(dailyprecip) == "logical") {
    cat("Extracting rainfall data from NCEP \n")
    dailyprecip <- dailyprecipNCEP(lat, long, tme, reanalysis2)
  }
  cat("Downscaling radiation and wind speed \n")
  radwind <- .pointradwind(hourlydata, dem, lat, long, l, x, albr, zmin, slope, aspect)
  cat("Calculating elevation and cold-air drainage effects \n")
  elev <- .eleveffects(hourlydata, dem, lat, long, windthresh, emthresh)
  return(list(hourlydata = hourlydata, hourlyradwind = radwind, tref = elev, dailyprecip = dailyprecip))
}

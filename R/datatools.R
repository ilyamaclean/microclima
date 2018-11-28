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
#' @param xdims optional numeric value specifying the number of longitudinal grid cells
#' @param ydims optional numeric value specifying the number of latitudinal grid cells
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
#' unspecified, then `lat` and `long` are used to retrieve a xdims x ydims cell raster
#' centred on `lat` and `long`. Bathymetry is also returned automatically form AWS,
#' so `zmin` allows the user to specify a sea-level height. Elevations below `zmin`
#' are set at `zmin`. For use with `microclima` the returned raster has a Mercator
#' Coordinate Reference System if the latitude specified by `lat` is between 80
#' degrees S and 84 degrees N. In polar regions WGS 84 / NSIDC Sea Ice Polar
#' Stereographic projections are used. X, y and Z units are all in metres. If `r` is
#' specified, the Coordinate Reference System of `r` is used.
#'
#' @examples
#' library(raster)
#' dem50m <- get_dem(dtm100m, resolution = 50)
#' plot(dem50m) # 50m raster, OSGB projection system
#'
#' dem5m <- get_dem(lat = 49.97, long = -5.22, resolution = 5)
#' plot(dem5m) # 5m raster, Mercator projection system
get_dem <- function(r = NA, lat, long, resolution = 30, zmin = 0, xdims = 200, ydims = 200) {
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
    e <- extent(c(xy$x - floor(xdims / 2) * resolution, xy$x + ceiling(xdims / 2) * resolution,
                  xy$y - floor(ydims / 2) * resolution, xy$y + ceiling(ydims / 2) * resolution))
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
      latdif <- abs(ll$y - as.numeric(rownames(v)))
      londif <- abs(ll$x%%360 - as.numeric(colnames(v)))
      v <- v[which.min(latdif), which.min(londif),]
    }
    v[sel]
  }
  tme <- as.POSIXlt(tme + 0, tz = "UTC")
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
  dfout <- data.frame(obs_time = tme2[sel], timezone = "UTC", Tk, Tkmin, Tkmax, sh, pr, wu, wv, dlw, ulw, dsw, tcdc)
  rownames(dfout) <- NULL
  return(dfout)
}

#' Calculates the solar index for a flat surface
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the solar altitude is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar altitude is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian Julian day expressed as an integer as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#'
#' @return a numeric value representing the solar altitude (º).
#' @export
#'
#' @examples
#' # solar index at noon GMT on 21 June 2010, Porthleven, Cornwall
#' jd <- julday (2010, 6, 21) # Julian day
#' siflat(12, 50.08, -5.31, jd)
#'
siflat <- function(localtime, lat, long, julian, merid = round(long / 15, 0) * 15, dst = 0){
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
# time sort out
#' @export
.tme.sort <- function(tme, UTC = TRUE) {
  zero <- function(x) ifelse(x < 10, paste0("0", x), paste0("", x))
  pasted <- function(tme, i) {
    paste0(tme$year[i] + 1900, "-", zero(tme$mon[i] + 1), "-", zero(tme$mday[i]))
  }
  tmeseq <- function(dstart, dfinish, tz) {
    seq(as.POSIXlt(dstart, format = "%Y-%m-%d", origin = "1900-01-01", tz = tz),
        as.POSIXlt(dfinish, format = "%Y-%m-%d", origin = "1900-01-01", tz = tz),
        by = 'days')
  }
  if (UTC) {
    tme <- as.POSIXlt(tme + 0, tz = "UTC")
    tme <- tmeseq(pasted(tme, 1), pasted(tme, length(tme)), "UTC")
  } else tme <- tmeseq(pasted(tme, 1), pasted(tme, length(tme)), "")
  tme
}
#' Interpolate NCEP data for running microclima to hourly
#'
#' @description `hourlyNCEP` optionally downloads the required NCEP climate and radiation forcing data
#' required for running microclima and interpolates 4x daily data to hourly.
#'
#' @param ncepdata an optional  data frame of climate variables as returned by [get_NCEP()].
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
#' @details If `ncepdata` is not provided, then [get_NCEP()] is called and data are downloaded from NCEP Atmospheric
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
  tmeout <- tme
  tme <- .tme.sort(tme)
  tme <- c(tme, tme[length(tme)] + 24 * 3600)
  if (class(ncepdata) != "logical") {
    tme <- as.POSIXlt(ncepdata$obs_time)
    tme <- tme[5:(length(tme) - 4)]
  }
  long <- ifelse(long > 180, long - 360, long)
  int <- as.numeric(tme[2]) - as.numeric(tme[1])
  lgth <- (length(tme) * int) / (24 * 3600)
  tme2 <- as.POSIXlt(c(0:(lgth - 1)) * 3600 * 24, origin = min(tme), tz = 'UTC')
  if (class(ncepdata) == "logical") {
    ncepdata <- get_NCEP(lat, long, tme2, reanalysis2)
  }
  tme6 <- as.POSIXlt(ncepdata$obs_time)
  tme6 <- as.POSIXlt(tme6 + 0, tz = "UTC")
  n <- (length(tme6) - 1) * 6 + 1
  h_pr <- spline(tme6, ncepdata$pr, n = n)$y
  h_sh <- spline(tme6, ncepdata$sh, n = n)$y
  h_cc <- bound(spline(tme6, ncepdata$tcdc, n = n)$y, mx = 100)
  tmo  <- spline(tme6, ncepdata$sh, n = n)$x
  tmo <- as.POSIXlt(tmo, origin = "1970-01-01 00:00",
                    tz = "UTC")
  tmo2 <- as.POSIXlt(seq(0,(length(tme6) / 4), by = 1/24) * 3600 * 24, origin = min(tme6), tz = 'UTC')
  jd <- julday(tmo2$year + 1900,  tmo2$mon + 1, tmo2$mday)
  si <- siflat(tmo2$hour, lat, long, jd, merid = 0)
  am <- airmasscoef(tmo2$hour, lat, long, jd, merid = 0)
  si_m <- 0
  am_m <- 0
  for (i in 1:(length(tme6))) {
    st <- (i - 1) * 6 + 1
    ed <- st + 5
    si_m[i] <- mean(si[st:ed], na.rm = T)
    am_m[i] <- mean(am[st:ed], na.rm = T)
  }
  jd <- julday(tme6$year + 1900, tme6$mon + 1, tme6$mday)
  dp <- 0
  for (i in 1:length(jd)) {
    dp[i] <- .difprop2(ncepdata$dsw[i], jd[i], tme6$hour[i], lat, long)
  }
  am_m[am_m > 5] <- 5
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
  edni <- dnir / ((4.87 / 0.0036) * (1 - dp))
  edif <- difr / ((4.87 / 0.0036) * dp)
  odni <- bound((log(edni) / -am_m), mn = 0.24, mx = 1.7)
  odif <- bound((log(edif) / -am_m), mn = 0.24, mx = 1.7)
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
  h_oi <- bound(spline(tme6, odni, n = n)$y, mn = 0.24, mx = 1.7)
  h_od <- bound(spline(tme6, odif, n = n)$y, mn = 0.24, mx = 1.7)
  tmorad <- as.POSIXlt(tmo + 3600 * 3)
  jd <- julday(tmorad$year + 1900, tmorad$mon + 1, tmorad$mday)
  szenith <- 90 - solalt(tmorad$hour, lat, long, jd, merid = 0)
  si <- siflat(tmorad$hour, lat, long, jd, merid = 0)
  am <- airmasscoef(tmorad$hour, lat, long, jd, merid = 0)
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
                                    lat = lat, long = long, merid = 0))
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
                          szenith = szenith[thsel], timezone = "UTC")
  tmetest <- as.POSIXlt(hourlyout$obs_time)
  tmeout <- as.POSIXlt(tmeout + 0, tz = "UTC")
  sel <- which(tmetest >= min(tmeout) & tmetest <= max(tmeout))
  return(hourlyout[sel,])
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
  tmeout <- tme
  tme <- .tme.sort(tme)
  long <- ifelse(long > 180, long - 360, long)
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
    latdif <- abs(ll$y - as.numeric(rownames(pre)))
    londif <- abs(ll$x%%360 - as.numeric(colnames(pre)))
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
  tmeout <- as.POSIXlt(tmeout + 0, tz = "UTC")
  sel <- which(tma >= min(tmeout) & tma <= (max(tmeout)))
  pre <- pre[sel]
  dpre <- t(matrix(pre, nrow = 4))
  dpre <- apply(dpre, 1, sum)
  return(dpre)
}
#' get shortwave radiation (time series)
#' @export
.shortwave.ts <- function(dni, dif, jd, localtime, lat, long, slope, aspect,
                          ha = 0, svv = 1, x = 1, l = 0, albr = 0.15,
                          merid = round(long / 15, 0) * 15, dst = 0, difani = TRUE) {
  sa <- solalt(localtime, lat, long, jd, merid)
  alt <- sa * (pi / 180)
  zen <- pi / 2 - alt
  saz <- solazi(localtime, lat, long, jd, merid)
  azi <- saz * (pi / 180)
  sl <- slope * (pi / 180)
  aspe <- aspect * (pi / 180)
  si <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - aspe)
  si[sa < ha] <- 0
  si[si < 0] <- 0
  dirr <- si * dni
  if (difani) {
    k <- dni / 4.87
  } else k <- 0
  isor <- 0.5 * dif * (1 + cos(sl)) * (1 - k)
  cisr <- k * dif * si
  sdi <- (slope + ha) * (pi / 180)
  refr <- 0.5 * albr * (1 - cos(sdi)) * dif
  fd <- dirr + cisr
  fdf <- isor + refr
  kk <- ((x ^ 2 + 1 / (tan(sa * (pi / 180)) ^ 2)) ^ 0.5) /
    (x + 1.774 * (x + 1.182) ^ (-0.733))
  trd <- exp(-kk * l)
  fr <- as.vector(canopy(array(l, dim = c(1,1)), array(x, dim = c(1,1))))
  trf <- (1 - fr)
  fgd <- fd * trd * (1 - albr)
  fged <- fdf * trf * (1 - albr) * svv
  fgc <- fgd + fged
  cfc <- ((1 - trd) * fd + fr * fdf) / (fd + fdf)
  cfc[is.na(cfc)] <- ((1 - trd[is.na(cfc)]) * 0.5 + fr * 0.5) / (0.5 + 0.5)
  return(xxx=list(swrad = fgc, canopyfact = cfc))
}
#' get radiation and wind for use with NicheMapR
#' @export
.pointradwind <- function(hourlydata, dem, lat, long, l, x, albr = 0.15, zmin = 0,
                          slope = NA, aspect = NA, horizon = NA, difani = TRUE) {
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
    wscr <- windcoef(dem, i*10, hgt = 2, res = reso)
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
  if (class(horizon) == "logical") {
    ha36 <- 0
    for (i in 0:35) {
      har <- horizonangle(dem, i*10, reso)
      ha36[i + 1] <- atan(extract(har, xy)) * (180/pi)
    }
  } else ha36 <- rep(horizon, 36)
  tme <- as.POSIXlt(hourlydata$obs_time)
  tme <- as.POSIXlt(tme + 0, tz = 'UTC')
  ha <- 0
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  for (i in 1:length(tme)) {
    saz <- solazi(tme$hour[i], lat, long, jd[i], merid = 0)
    saz <- round(saz / 10, 0) + 1
    saz <- ifelse(saz > 36, 1, saz)
    ha[i] <- ha36[saz]
  }
  sw <-.shortwave.ts(hourlydata$rad_dni, hourlydata$rad_dif, jd, tme$hour,
                        lat, long, slope, aspect, ha, svf, x, l, albr, merid = 0,
                        difani = difani)
  windsp <- windheight(hourlydata$windspeed, 10, 1)
  hourlyrad <- data.frame(swrad = sw$swrad, skyviewfact = svf, canopyfact = sw$canopyfact,
                          whselt = wsheltatground, windspeed = wshelt *  windsp,
                          slope = slope, aspect = aspect)
  hourlyrad
}
#' Cold air drainage direct from emissivity
#' @export
.cadconditions2 <- function (em, wind, startjul, lat, long, starttime = 0, hourint = 1,
                             windthresh = 4.5, emthresh = 0.725, con = TRUE)
{
  jd <- floor(c(1:length(em)) * hourint/24 - hourint/24 + startjul +
                starttime/24)

  st <- suntimes(jd, lat, long, merid = 0)
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
#' @export
.eleveffects <- function(hourlydata, dem, lat, long, windthresh = 4.5,
                          emthresh = 0.78) {
  xy <- data.frame(x = long, y = lat)
  elevncep <- extract(demworld, xy)
  coordinates(xy) = ~x + y
  proj4string(xy) = "+init=epsg:4326"
  xy <- as.data.frame(spTransform(xy, crs(dem)))
  elev <- extract(dem, xy)
  if (is.na(elev)) elev <- 0
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
  if (is.na(cdif)) cdif <- 0
  cadt <- cdif * lr * cad
  tout <- data.frame(tref = hourlydata$temperature,
                     elev = elev, elevncep = elevncep,
                     telev = lr * elev,
                     tcad = cadt)
  return(list(tout = tout, basins = basins, flowacc = fa))
}
#' Calculates coastal exposure automatically
#' @export
.invls.auto <- function(r, steps = 8, use.raster = T, zmin = 0, plot.progress = TRUE, tidyr = FALSE) {
  tidydems <- function(rfine, rc) {
    rfine[is.na(rfine)] <- zmin
    rc <- trim(rc)
    aggf <- floor(mean(res(rc)[1:2]) / mean(res(rfine)))
    if (aggf > 1) rfine2 <- suppressWarnings(aggregate(rfine, aggf, max))
    rfine2 <- suppressWarnings(resample(rfine2, rc, method = 'ngb'))
    rfine2 <- crop(rfine2, extent(rfine))
    if (dim(rfine2)[1] * dim(rfine2)[2] == 1) {
      xx <- mean(is_raster(rfine), na.rm = T)
      rfine2 <- if_raster(as.matrix(xx), rfine2)
    }
    a <- array(-999, dim = dim(rfine2)[1:2])
    rc2 <- raster(a)
    extent(rc2) <- extent(rfine2)
    rc2 <- mosaic(rc, rc2, fun = min)
    rc2 <- mosaic(rc2, rfine2, fun = max)
    rc2[rc2 == zmin] <- NA
    rc2
  }
  adjust.lsr <- function(lsr, rs) {
    m <- is_raster(lsr)
    m[m < 0] <- 0
    s <- c(0, (8:10000) / 8) ^ 2 * 30
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
  ress <- c(30, 90, 500, 1000, 10000)
  ress <- ress[ress > mean(res(r))]
  ress <- c(mean(res(r)), ress)
  ress <- rev(ress)
  # Create a list of dems
  cat("Downloading land sea data \n")
  dem.list <- list()
  rmet <- projectRaster(rll, crs = crs(r))
  dem <- get_dem(rmet, resolution = ress[1], zmin = zmin)
  dem.list[[1]] <- projectRaster(dem, crs = crs(r))
  dc <- ceiling(max(dim(r)) / 2)
  rres <- mean(res(r))
  for (i in 2:(length(ress)-1)) {
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
  if (tidyr & mean(res(r)) <= 30)
    warning("raster tidying ignored as resolution <= 30")
  if (tidyr & mean(res(r)) > 30) {
    rx <- raster(extent(r))
    res(rx) <- 30
    crs(rx) <- crs(r)
    rfine <- get_dem(rx, resolution = 30, zmin = zmin)
    r <- tidydems(rfine, r)
  }
  dem.list[[i + 1]] <- r
  for (j in 1:i) dem.list[[j]] <- tidydems(r, dem.list[[j]])
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
      lsa[is.na(lsa)] <- 0
      lsm <- is_raster(lsa)
      if (max(lsm) > min(lsm)) {
        if (min(dim(lsm)) < 2) {
          r2 <- aggregate(r, 100)
          xx<- resample(lsa, r2, method = "ngb")
          xx <- resample(xx, r)
          xx <- mask(xx, r)
        } else {
          xx <- resample(lsa, r)
          x <- mask(xx, r)
        }
        mx <- mean(is_raster(xx), na.rm = T)
        if  (is.na(mx)) xx <- xx * 0 + mean(lsm, na.rm = T)
        lsa.list[[i]] <- xx
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
  rone[is.na(rone)] <- zmin
  dem <- dem.list[[1]]
  m <- is_raster(dem)
  m[m == zmin] <- NA
  dem <- if_raster(m, dem)
  e <- extent(r)
  xy <- data.frame(x = (e@xmin + e@xmax) / 2,
                   y = (e@ymin + e@ymax) / 2)
  coordinates(xy) = ~x+y
  proj4string(xy) = crs(r)
  for (dct in 0:(steps - 1)) {
    direction <- dct * (360 / steps)
    lsa <- invls(dem, extent(rone), direction)
    f1 <- extract(lsa, xy)
    if (is.na(f1)) f1 <- mean(is_raster(lsa), na.rm = T)
    mncep <- mean(is_raster(lsa), na.rm = T) / f1
    cncep[dct + 1] <- mncep
  }
  cat("Adjusting coastal exposure by ncep means \n")
  for (i in 1:steps) lsa.array[,,i] <- lsa.array[,,i] * cncep[i]
  return(lsa.array)
}
#' Downloads sea-surface temperature data
#' @export
#' @import rnoaa ncdf4
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
  tme <- as.POSIXlt(tme + 0, tz = "GMT")
  sel <- which(tmeh >=min(tme) & tmeh <= max(tme))
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
#' @param tidyr logical value indicating whether to download 30m resolution digital elevation data
#' and use these data to improve the assignment of land and sea pixels in `landsea`.
#' @details coastalNCEP downloads digital elevation and varying resolutions to calculate
#' coastal effects by applying a variant of [invls()] over a greater extent then that specified by landsea, but the resulting outputs are
#' cropped to the same extent as landsea. Land temperatyres as a function of coastal exposure
#' and sea-surface temperature data are then calculated. Sea surface temperature data are
#' automatically downloaded from the NOAA.
#'
#' @return a three-dimension array of hourly temperature values for each pixel of `landsea`
#' @export
#' @import raster sp zoo rnoaa ncdf4
#' @examples
#' library(raster)
#' # Download NCEP data
#' ll <- latlongfromraster(dtm100m)
#' tme <-as.POSIXlt(c(0:31) * 3600 * 24, origin = "2015-03-15", tz = "GMT")
#' ncephourly<-hourlyNCEP(NA, ll$lat, ll$long, tme)
#' aout <- coastalNCEP(dtm100m, ncephourly)
#' # Calculate mean temperature and convert to raster
#' mtemp <- if_raster(apply(aout, c(1, 2), mean), dtm100m)
#' plot(mtemp, main = "Mean temperature")
coastalNCEP <- function(landsea, ncephourly, steps = 8, use.raster = T, zmin = 0, plot.progress = T, tidyr = F) {
  bound <- function(x, mn = 0.1, mx = 1) {
    x[x < mn] <- mn
    x[x > mx] <- mx
    x
  }
  lsr1 <- .invls.auto(landsea, steps, use.raster, zmin, plot.progress, tidyr)
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
#' @param dstart start date as character of time period required in format DD/MM/YYYY (UTC timezones assumed)
#' @param dfinish end date as character of time period required in format DD/MM/YYYY (UTC timezones assumed)
#' @param l  a single numeric value of the leaf area index for location (see details).
#' @param x  a single numeric value representing the ratio of vertical to horizontal
#' projections of leaf foliage for the location (see details).
#' @param coastal an optional logical value indicating whether to calculate coastal effects using [coastalNCEP()].
#' @param hourlydata an optional data frame of ourly weather and radiation data as returned
#' by function [hourlyNCEP()]. Function called if data not provided.
#' @param dailyprecip an optional data frame of ourly weather and radiation data as returned
#' by function [dailyprecipNCEP()]. Function called if data not provided.
#' @param dem an ptional raster object of elevations covering the location for which point
#' data are required. If not provided, one is downloaded from the registry of Open
#' Data on AWS. Used to calculate coastal, wind and radiation effects
#' @param demmeso an optional raster object of elevations covering the location for which point
#' data are required. Used to calculate elevation and cold air drainage effects.
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
#' @param steps  an optional integer. Coastal effects are calculated in specified directions upwind. Steps defines the total number of directions used. If the default 8 is specified, coastal effects are calculated at 45º intervals.
#' @param use.raster an optional logical value indicating whether to mask the output values by `landsea`.
#' @param plot.progress logical value indicating whether to produce plots to track progress.
#' @param tidyr logical value indicating whether to download 30m resolution digital elevation data.
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
#' (5) `acoast` an array of temperatures that account for coastal effects for each
#' pixel of `dem` and each hour of hourlydata.
#'
#' (6) `basins` a raster object of basins used for calculating cold air drainage.
#'
#' (7) `flowacc` a raster object of accumulated flow.
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
#' [hourlyNCEP()] for retrieving and processing hourly data from NCEP and
#' [coastalNCEP()] for retrieving and processing coastal effects.
#'
#' @export
#' @import raster
#' @import sp
#' @import zoo
#'
#' @examples
#' mnr <- microclimaforNMR(50, -5.2, '15/01/2015', '15/02/2015', 1, 1, coastal = FALSE)
#' head(mnr$hourlydata)
#' head(mnr$hourlyradwind)
#' head(mnr$tref)
#' head(mnr$dailyprecip)
microclimaforNMR <- function(lat, long, dstart, dfinish, l, x, coastal = TRUE, hourlydata = NA,
                             dailyprecip = NA, dem = NA, demmeso = dem, albr =0.15,
                             resolution = 100, zmin = 0, slope = NA, aspect = NA,
                             windthresh = 4.5, emthresh = 0.78, reanalysis2 = TRUE,
                             steps = 8, use.raster = TRUE, plot.progress = TRUE, tidyr = FALSE) {
  tme <- seq(as.POSIXlt(dstart, format = "%d/%m/%Y", origin = "01/01/1900"),
             as.POSIXlt(dfinish, format = "%d/%m/%Y", origin = "01/01/1900")
             + 3600 * 24, by = 'hours')
  tme <- tme[-length(tme)]
  tme <- as.POSIXlt(tme)
  if (class(dem) == "logical") {
    cat("Downloading digital elevation data \n")
    dem <- get_dem(r = NA, lat = lat, long = long, resolution = resolution, zmin = zmin)
  }
  if (class(demmeso) == "logical") demmeso <- dem
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
  info <- .eleveffects(hourlydata, demmeso, lat, long, windthresh, emthresh)
  elev <- info$tout
  if (coastal) {
    m <- is_raster(dem)
    m[m == zmin] <- NA
    dem <- if_raster(m, dem)
    acoast <- coastalNCEP(dem, hourlydata, steps, use.raster, zmin, plot.progress, tidyr = FALSE)
    xy <- SpatialPoints(data.frame(x = long, y = lat))
    proj4string(xy) = crs("+init=epsg:4326")
    pr <- crs(dem, asText = T)
    xy <- spTransform(xy, CRS(pr))
    e <- extent(dem)
    xi <- round((xy$x - e@xmin) / res(dem)[1] + 0.5, 0)
    yi <- round((xy$y - e@ymin) / res(dem)[2] + 0.5, 0)
    ce <- acoast[xi, yi, ]
    if (is.na(ce[1])) ce <- apply(acoast, 3, mean, na.rm = T)
    elev$tref <- ce + elev$telev
  }  else acoast <- NA
  return(list(hourlydata = hourlydata, hourlyradwind = radwind, tref = elev,
              dailyprecip = dailyprecip, acoast = acoast, basins = info$basins,
              flowacc = info$flowacc))
}
#' Function for automatically generating microclimate surfaces for anywhere in the word
#'
#' @description This function generating microclimate temperature surfaces for anywhere
#' in the word. It first downloads coarse-resolution climate and radiationdata from The
#' NCEP-NCAR or NCEP–DOE Atmospheric Model Intercomparison Project (Kanamitso et al
#' 2002) and interpolates these data to hourly. It calculates mesoclimatic effects and derives parameters for fitting the
#' microclimate model using the `NicheMapR` package (Kearny & Porter 2016). Using digital
#' elevation data either downloaded or provided by the user, and canopy
#' characteristics either specified by the user, or derived from habitat,
#' it runs the microclimate model in hourly timesteps to generate an array oftemperatures.
#'
#' @param r a raster object defining the extent and resolution for which microclimate
#' temperature data are required.
#' @param dstart start date as character of time period required in format DD/MM/YYYY
#' @param dfinish end date as character of time period required in format DD/MM/YYYY
#' @param hgt the height (in m) above or below ground for which temperature estimates
#' are required. If negative, below-ground temperatures are derived (range -2 to 2).
#' @param l  an optional single numeric value, vector, raster, matrix or array of leaf area
#' index values. If a single numeric value, the same leaf area is assumed at all locations
#' and times. If a vector, the same leaf area is assumed at all locations, but leaf area is
#' assumed ot vary through time. If a matrix or raster, leaf area is assumed to vary by
#' location, but is assumed static in time. If an array, variable leaf areas in time and
#' space are assumed.
#' @param x  a single numeric value, matrix or raster representing the ratio of vertical to horizontal projections of foliage as,
#' for example, returned by [leaf_geometry()]. if a single numeric value, the same value of `x`
#' is assumed at all locations. Varying `x` through time is not supported.
#' @param habitat a character string, numeric value or matrix or raster of numeric
#' values of habitat type(s). See [habitats()]. Ignored if l or x are provided. If a single
#' numeric value, the same habitat is assumed at all locations. Varying habitats in time are not
#' supported.
#' @param coastal optional logical value indicating whether or not to calculate
#' coastal effects.
#' @param use.raster optional logical value indicating whether to use `r` in the
#' calculation of coastal effects. If used, NAs or values corresponding to `zmin`
#' must represent sea.
#' @param r.is.dem optional logical value indicating whether 'r' is a digital
#' elevation dataset used for calculating microclimatic effects. If FALSE, then a
#' dem is downloaded.
#' @param hourlydata Optional data.frame of hourly climate data as returned by [hourlyNCEP()]
#' @param dailyprecip Optional data.frame of daily precipitation data as returned by [dailyprecipNCEP()]
#' @param albg an optional single value, raster object, two-dimensional array or
#' matrix of values representing the albedo(s) of the ground as returned by [albedo2()].
#' @param albr albr an optional single value, raster object, two-dimensional array
#' or matrix of values representing the albedo(s) of adjacent surfaces as returned
#' by [albedo_reflected()].
#' @param albc an optional single value, raster object, two-dimensional array or
#' matrix of values representing the albedo(s) of the vegetated canopy as returned
#' by [albedo2()].
#' @param mesoresolution optional numeric value indicating resolution of the dem used for
#' modelling mesoclimate.
#' @param zmin assumed sea-level height. Values below this are set to zmin (see details).
#' @param slope an optional slope value for the location in decimal degrees. If not specified, then
#' the slope is calculated from the retrieved dem.
#' @param aspect an optional aspect value for the location in decimal degrees. If not specified, then
#' the slope is calculated from the retrieved dem.
#' @param windthresh an optional threshold wind value (m /s) above which cold air drainage is
#' assumed not to occur.
#' @param emthresh an optional threshold emissivity value above which cold air drainage is
#' assumed not to occur.
#' @param reanalysis2 Optional logical. Should data be obtained from the Reanalysis II dataset (default) or
#' from Reanalysis I (data prior to 1979).
#' @param steps  an optional integer. Coastal effects are calculated in specified
#' directions upwind. Steps defines the total number of directions used. If the
#' default 8 is specified, coastal effects are calculated at 45º intervals.
#' @param plot.progress an optional logical indicating whether to produce plots to track progress.
#' @param continious an optional logical value indicating whether to treat wind speed as a continious variable
#' @param summarydata an optional logical idicating whether to calculate summary data
#' (frost hours and maximum, minimum and mean temperature) for each pixel and return these to the output.
#' @param save.memory An optional logical indicatign whether to save
#' memory by storing temperature x 1000 as an integer values.
#'
#' @references Kearney MR,  Porter WP (2016) NicheMapR – an R package for biophysical
#' modelling: the microclimate model. Ecography 40: 664-674.
#'
#' Kanamitsu M, Ebisuzaki W, Woollen J, Yang SK, Hnilo JJ, Fiorino M, Potter GL (2002)
#'Ncep–doe amip-ii reanalysis (r-2). Bulletin of the American Meteorological Society 83: 1631-1644.
#'
#' @return a list with the following objects:
#' (1) temps: an array of temperatures for each pixel of r and hour of the time sequence.
#' (2) e: an extent object given the extent of `temps`. Generally the same as `raster::extent(r)`
#' though note that edge cells are NA as slopes cannot be calculate for these cells.
#' (3) units: the units of `temps`. Either deg C or dec C * 1000 if `save.memory` is TRUE.
#' (4) tmax: If `summarydata` is TRUE, a matrix of maximum temperatures
#' (5) tmin: If `summarydata` is TRUE, a matrix of minimum temperatures
#' (6) tmean: If `summarydata` is TRUE, a matrix of mean temperatures
#' (7) frosthours: If `summarydata` is TRUE, a matrix of hours below 0 deg C.
#'
#' @import raster NicheMapR
#' @export
#'
#' @examples
#' library(raster)
#' library(NicheMapR)
#' # Get dem for Pico, Azores
#' r <- get_dem(lat = 38.467429, long = -28.398995, resolution = 30)
#' plot(r)
#' # Takes ~ c. 15 minutes to run
#' temps <- runauto.ncep(r, "10/06/2010", "15/06/2010", hgt = 0.1, l = NA, x = NA,
#'                      habitat = "Barren or sparsely vegetated")
#' mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange",
#'                           "red"))(255)
#' meantemp <- temps$tmean
#' extent(meantemp) <- temps$e
#' par(mfrow = c(1, 1))
#' plot(meantemp, main = "Mean temperature", col = mypal)
#' # Interactive 3D plot
#' require(plotly)
#' zrange<-list(range = c(0, 3000))
#' plot_ly(z = ~is_raster(r)) %>%
#' add_surface(surfacecolor = ~is_raster(meantemp)) %>%
#' layout(scene = list(zaxis = zrange))

runauto.ncep <- function(r, dstart, dfinish, hgt = 0.05, l, x, habitat = NA,
                         use.raster = FALSE, coastal = TRUE, r.is.dem = TRUE,
                         hourlydata = NA, dailyprecip = NA, albg = 0.15,
                         albr =0.15, albc = 0.23, mesoresolution = 100,
                         zmin = 0, slope = NA, aspect = NA, windthresh = 4.5,
                         emthresh = 0.78, reanalysis2 = TRUE, steps = 8,
                         plot.progress = TRUE, continious = TRUE,
                         summarydata = TRUE, save.memory = FALSE) {
  longwaveveg2 <- function(le0, lwsky, x, fr, svv, albc) {
    lw1 <- (1 - fr) * lwsky
    lr <- (2 / 3) * log(x + 1)
    r <- 1 / (1 + exp(-1 * lr))
    lw2 <- albc * r * fr * lwsky
    lw3 <- r * (1 - albc) * le0 * fr
    lwd <- lw1 + lw2 + lw3
    lwrad <- (le0 - lwd) * svv
    lwrad
  }
  habfun <- function(habitat, lat, long, sel = NA, hoy, year, la, hgt) {
    if (hgt < 0) hgt <- 0
    xx <- laifromhabitat(habitat[1], lat, long, year = year)
    lv <- xx$lai[hoy]
    hh <- xx$height
    mm <- (1 - (hgt / hh))
    if (hgt > hh) mm <- 0
    if (class(sel) == "logical") {
      la <- round(lv * mm, 3)
    } else la[sel] <- round(lv, 3)
    la
  }
  laitime <- function(l, habitat, lat, long, tme, i, hgt) {
    yr <- tme$year[i] + 1900
    hoy <- (as.numeric(strftime(tme[i], format = "%j")) - 1) * 24 +
      as.numeric(tme$hour[i]) + 1
    if (class(habitat) == "matrix") {
      u <- unique(as.vector(habitat))
      u <- u[is.na(u) == F]
      la <- array(NA, dim = c(dim(r)[1:2]))
      for (ii in 1:length(u)) {
        sel <- which(habitat == u[ii], arr.ind = T)
        la <-  habfun(habitat[sel], lat, long, sel, hoy, yr, la, hgt)
      }
    }
    if  (class(habitat) == "numeric" | class(habitat) == "character") {
      la <- 0
      la <-  habfun(habitat, lat, long, sel = NA, hoy, yr, la, hgt)
      la <- array(la, dim = dim(r)[1:2])
    }
    if (class(l) == "numeric" | class(l) == "integer") {
      if (length(l) == 1) {
        la <- array(l, dim = dim(r)[1:2])
      } else la <- array(l[i], dim = dim(r)[1:2])
    }
    if (class(l) == "matrix") la <- l
    if (class(l) == "array") la <- l[,,i]
    la
  }
  checkfun <- function(l, x, habitat, r, tme) {
    if (class(l) == "logical" & class(habitat) == "logical") {
      stop ("No habitat or leaf area provided")
    }
    if (class(x) == "logical" & class(habitat) == "logical") {
      stop ("No habitat or leaf angle distribution provided")
    }
    if (class(l) == "matrix" | class(l) == "array") {
      if (dim(l)[1] != dim(r)[1] & dim(l)[2] != dim(r)[2]) {
        stop ("l must have same xy dimensions as r")
      }
    }
    if (class(l)== "array") {
      if (dim(l)[3] != length(tme)) {
        stop ("Length of l must equal 1 or number of hours in time sequence")
      }
    }
    if (class(l)== "numeric" | class(l) == "integer") {
      if (length(l) != 1 & length(l)  != length(tme)) {
        stop ("Length of l must equal 1 or number of hours in time sequence")
      }
    }
    if (class(habitat) == "matrix") {
      if (dim(habitat)[1] != dim(r)[1] & dim(habitat)[2] != dim(r)[2]) {
        stop("habitat must have same dimensions as r")
      }
    }
    if (class(x) == "matrix") {
      if (dim(x)[1] != dim(r)[1] & dim(x)[2] != dim(r)[2]) {
        stop("x must have same dimensions as r")
      }
    }
    if (class(x) == "numeric" | class(x) == "integer") {
      if (length(x) > 1) stop("varying x through time not supported")
    }
    if (class(x) == "array") stop("varying x through time not supported")
  }
  frost <- function(x) {
    sel <- which(x <= 0)
    y <- length(sel)
    if (is.na(mean(x, na.rm = TRUE))) y <- NA
    y
  }
  #if (!requireNamespace("NicheMapR", quietly = TRUE)) {
  #  stop("package 'NicheMapR' is needed. Please install it from Github: 'mrke/NicheMapR'",
  #       call. = FALSE)
  #}
  # Lat long and time
  lat <- latlongfromraster(r)$lat
  long <- latlongfromraster(r)$long
  tme <- seq(as.POSIXlt(dstart, format = "%d/%m/%Y", origin = "01/01/1900"),
             as.POSIXlt(dfinish, format = "%d/%m/%Y", origin = "01/01/1900")
             + 3600 * 24, by = 'hours')
  tme <- tme[-length(tme)]
  tme <- as.POSIXlt(tme)
  # l and x
  l <- is_raster(l)
  x <- is_raster(x)
  habitat <- is_raster(habitat)
  checkfun(l, x, habitat, r, tme)
  if (class(x) == "numeric" | class(x) == "integer") x <- array(x, dim = dim(r)[1:2])
  if (class(x) == "logical") {
    if (class(habitat)  == "numeric" | class(habitat) == "integer" |
        class(habitat) == "character") {
      x <- laifromhabitat(habitat[1], lat, long, year = 2015)$x
      x <- array(x, dim = dim(r)[1:2])
    } else {
      u <- unique(as.vector(habitat))
      u <- u[is.na(u) == F]
      x <- array(NA, dim = c(dim(r)[1:2]))
      for (ii in 1:length(u)) {
        sel <- which(habitat == u[ii])
        x[sel] <-  laifromhabitat(u[ii], lat, long, year = 2015)$x
      }
    }
  }
  # Run NicheMapR
  loc <- as.numeric(c(long, lat))
  reso <- mean(res(r))
  if (reso < 100) {
    cat("Downloading DEM for mesoscale modelling \n")
    dem <- get_dem(r = NA, lat = lat, long = long, resolution = mesoresolution, zmin = zmin)
    if (mesoresolution > 30) {
      dm30 <- get_dem(r = dem, lat = lat, long = long, resolution = 30, zmin = zmin)
      dem <- resample(dm30, dem)
    }
  } else {
    if (r.is.dem) {
      dem <- r
    } else dem <- suppressWarnings(get_dem(r = r, lat = lat, long = long, resolution = reso, zmin = zmin))
  }
  if (r.is.dem == F) {
    r <- suppressWarnings(get_dem(r = r, lat = lat, long = long, resolution = reso, zmin = zmin))
  }
  if (hgt < -2) {
    warning("At depths > 2 m temperature assumed to to be stable. Depth set to 2 m")
    hgt <- -2
  }
  if (hgt > 2) {
    warning("User height greater than reerence height. Height set to 2 m")
    hgt <- 2
  }
  dep <- c(2.5, 5, 10, 15, 20, 30, 50, 100, 200)
  if (hgt < 0) {
    shgt <- hgt * -100
    dif <- abs(dep - shgt)
    selsoil <- which(dif == min(dif))
    dep[selsoil] <- shgt
    hgt2 <- 0.05
  } else if (hgt < 0.005) {
    hgt2 <- 0.005
  } else hgt2 <- hgt
  dep <- c(0, dep)
  loc <- c(long, lat)
  dem[dem == zmin] <- NA
  r[r == zmin] <- NA
  micronmr <- .micro_ncep2(dstart = dstart, dfinish = dfinish, dem = r, dem2 = dem, LAI = 0,
                           loc = loc, Usrhyt = hgt2, Refhyt = 2, coastal = coastal,
                           DEP = dep)
  ma <- micronmr$microclima.out
  hourlydata <- ma$hourlydata
  #####################
  tme <- as.POSIXlt(tme + 0, tz = "UTC")
  nmrout <- as.data.frame(micronmr$metout)
  rwind <- ma$hourlyradwind
  netrad = (1 - 0.15) * rwind$swrad - 0.98 * hourlydata$netlong * rwind$skyviewfact * (1 - rwind$canopyfact)
  cat("Paramaterising model using NicheMapR \n")
  if (hgt > 0) {
    mft <- data.frame(temperature = nmrout$TALOC,
                      reftemp = nmrout$TAREF,
                      wind = rwind$windspeed,
                      netrad = netrad)
    params <- fitmicro(mft, alldata = TRUE, continious = T)
  }
  if (hgt <= 0) {
    soiltemps <- as.data.frame(micronmr$soil)
    mft <- data.frame(temperature = soiltemps$D0cm,
                      reftemp = nmrout$TAREF,
                      wind = rwind$windspeed,
                      netrad = netrad)
    params <- fitmicro(mft, alldata = TRUE, continious = T)
  }
  if (hgt < 0) {
    pred0 <- mft$reftemp + params$Estimate[1] + params$Estimate[3] * netrad +
      params$Estimate[2] * log(rwind$windspeed + 1) +
      params$Estimate[4] * netrad * log(rwind$windspeed + 1)
    soiltemp <- soiltemps[,selsoil+3]
    dfsoil <- .fitsoil(pred0, mean(soiltemps$D200cm), soiltemp)
  }
  if (coastal) {
    tref <- ma$acoast
  } else {
    tref <- rep(hourlydata$temperature, each = dim(r)[1] * dim(r)[2])
    tref <- array(tref, dim = c(dim(r)[1:2], length(tme)))
  }
  tout <- array(NA, dim = dim(tref))
  if (hgt < 0) {
    tout <- abind(tout, tout[,,1])
    tout[,,1] <- pred0[1]
  }
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  m <- is_raster(dem)
  m[is.na(m)] <- zmin
  dem2 <- if_raster(m, dem)
  m <- is_raster(r)
  m[is.na(m)] <- zmin
  r2 <- if_raster(m, r)
  cat("Calculating wind shelter coefficients \n")
  wca <- array(NA, dim = c(dim(r)[1:2], 16))
  rff <- raster(r)
  res(rff) <- res(r)[1] / 4
  r3 <- resample(r2, rff)
  hgtw <- ifelse(hgt < 0, 0, hgt)
  for (i in 0:15) {
    wcr <- windcoef(r3, i * 22.5, hgt = hgtw, res = res(r3))
    wcr <- aggregate(wcr, 8)
    wca[,,i + 1] <- is_raster(resample(wcr, r))
  }
  wca2 <- wca
  for (i in 1:16) {
    mn <- i - 1
    mn <- ifelse(mn == 0, 8, mn)
    mx <- i + 1
    mx <- ifelse(mx >8, mx - 8, mx)
    wca2[,,i] <- 0.25 * wca[,,mn]  + 0.5 * wca[,,i] + 0.25 * wca[,,mx]
  }
  wca <- wca2
  wca2 <- NULL
  cad <- ma$tref
  cad <- cad$tcad
  basin <- mask(ma$basins, dem)
  fa <- mask(ma$flowacc, dem)
  la <-  laitime(l, habitat, lat, long, tme, 1, hgt)
  svv <- skyviewveg(is_raster(r2), la, x, steps = 36, res = reso)
  ha <- mean_slope(is_raster(r2), res= reso)
  cat("Running model in hourly time-steps \n")
  if (hgt < 0) tout[,,1] <- soiltemp[1]
  if  (plot.progress) {
    mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange", "red"))(255)
    mypal2 <- rev(colorRampPalette(c("white", "blue", "green"))(255))
  }
  mnwind <- as.numeric(quantile(rwind$windspeed, 0.05))
  mxwind <- as.numeric(quantile(rwind$windspeed, 0.95))
  for (i in 1:length(jd)) {
    # elevation effects
    tr <- mean(tref[,,i], na.rm = T)
    lr <- lapserate(tr, hourlydata$humidity[i], hourlydata$pressure[i])
    tr <- tr + is_raster(r2) * lr
    # cold air drainage effects
    if (cad[i] > 0) {
      cadp <- pcad(dem, basin, fa, hourlydata$temperature[i],
                   hourlydata$humidity[i], hourlydata$pressure[i])
      cadp <- projectRaster(cadp, crs = crs(r))
      cadp <- resample(cadp, r2)
      cadp <- mask(cadp, r)
      tr <- tr + cadp
    }
    if (i%%100 == 0) {
      la <-  laitime(l, habitat, lat, long, tme, i, hgt)
      svv <- skyviewveg(is_raster(r2), la, x, steps = 36, res = reso)
      ha <- mean_slope(is_raster(r2), res = reso)
    }
    fr <- canopy(la, x)
    radsw <- shortwaveveg(hourlydata$rad_dni[i], hourlydata$rad_dif[i], jd[i],
                          tme$hour[i], dtm = r2, svv = svv, albg = albg,
                          fr = fr, albr = albr, ha = ha, res = reso, merid = 0,
                          x = x, l = la)
    radlw <- longwaveveg2(hourlydata$uplong[i], hourlydata$downlong[i], x, fr,
                          svv, is_raster(albc))
    nr <- is_raster(radsw) - 0.98 * radlw
    nr <- if_raster(nr, r2)
    wi <- round(hourlydata$winddir[i] / 22.5, 0) + 1
    wi <- ifelse(wi > 16, wi - 16, wi)
    wi <-ifelse(wi < 1, wi + 16, wi)
    wspd <- windheight(hourlydata$windspeed[i], 2, 1)
    ws <- wca[,,wi] * wspd
    ws <- ifelse(ws < mnwind, mnwind, ws)
    ws <- ifelse(ws > mxwind, mxwind, ws)
    ws <- if_raster(ws, r2)
    tea <- runmicro(params, nr, ws, continious = TRUE)
    tea <- mask(tea, r)
    if (hgt < 0) {
      soil0 <- is_raster(tea) + tr
      heatdown <-  soil0 - tout[,,i]
      heatup <-  mean(soiltemps$D200cm) - tout[,,i]
      xx <- tout[,,i] + dfsoil$x1 * heatdown + dfsoil$x2 * heatup
      if (save.memory) xx <- round(xx * 1000, 0)
      tout[,,i + 1] <- xx
    } else  {
      xx <- is_raster(tr) + is_raster(tea)
      if (save.memory) xx <- round(xx * 1000, 0)
      tout[,,i] <- xx
    }
    if (plot.progress & i%%100 == 0) {
      tempr <- if_raster(tout[,,i], r)
      if (save.memory) tempr <- tempr / 1000
      plot(tempr, main = tme[i], col = mypal)
    }
  }
  if (hgt < 0) {
    dtr <- dim(tout)[3]
    tout <- tout[,,-dtr]
  }
  tmax <- NA
  tmin <- NA
  tmean <- NA
  hfrost <- NA
  if (summarydata) {
    cat("Calculating summary data \n")
    tmax <- if_raster(apply(tout, c(1, 2), max), r)
    tmin <- if_raster(apply(tout, c(1, 2), min), r)
    tmean <- if_raster(apply(tout, c(1, 2), mean), r)
    tfrost <- if_raster(apply(tout, c(1, 2), frost), r)
  }
  if (plot.progress) {
    if (save.memory) {
      tmax <- tmax / 1000
      tmin <- tmin / 1000
      tmean <- tmean / 1000
    }
    par(mfrow=c(2,2))
    plot(tmax, main = "Maximum temperature", col = mypal)
    plot(tmin, main = "Minimum temperature", col = mypal)
    plot(tmean, main = "Mean temperature", col = mypal)
    plot(tfrost, main = "Frost hours", col = mypal2)
  }
  if (save.memory) {
    return(list(temps = tout, e = extent(r), tme = tme, units = "deg C * 1000",
                tmax = tmax, tmin = tmin, tmean = tmean, frosthours = hfrost))
  } else {
    return(list(temps = tout, e = extent(r), tme = tme, units = "deg C",
                tmax = tmax, tmin = tmin, tmean = tmean, frosthours = hfrost))
  }
}
#' Temporary NMR Function
#' @export
.micro_ncep2 <- function (loc = c(-5.3, 50.13), dstart = "01/01/2017", dfinish = "31/12/2017",
                          dem = NA, dem2 = dem, nyears = as.numeric(substr(dfinish, 7, 10)) - as.numeric(substr(dstart, 7, 10)) + 1,
                          REFL = 0.15, slope = NA, aspect = NA, DEP = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
                          Refhyt = 2, Usrhyt = 0.01, Z01 = 0, Z02 = 0, ZH1 = 0, ZH2 = 0,
                          run.gads = 1, write_input = 0, writecsv = 0, reanalysis = TRUE,
                          windfac = 1, warm = 0, ERR = 1.5, RUF = 0.004, EC = 0.0167238,
                          SLE = 0.95, Thcond = 2.5, Density = 2.56, SpecHeat = 870,
                          BulkDensity = 1.3, PCTWET = 0, rainwet = 1.5, cap = 1, CMH2O = 1,
                          hori = rep(0, 24), runmoist = 1, PE = rep(1.1, 19),
                          KS = rep(0.0037, 19), BB = rep(4.5, 19), BD = rep(1.3, 19),
                          DD = rep(2.56, 19), maxpool = 10000, rainmult = 1, evenrain = 0,
                          SoilMoist_Init = c(0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.3, 0.3, 0.3, 0.3),
                          L = c(0, 0, 8.2, 8, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4, 1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0) * 10000,
                          R1 = 0.001, RW = 2.5e+10, RL = 2e+06, PC = -1500, SP = 10,
                          IM = 1e-06, MAXCOUNT = 500, LAI = 0.1, LOR = 1, snowmodel = 1,
                          snowtemp = 1.5, snowdens = 0.375, densfun = c(0.5979, 0.2178, 0.001, 0.0038),
                          snowmelt = 1, undercatch = 1, rainmelt = 0.0125, shore = 0,
                          tides = 0, hourly = 1, rainhourly = 0, rainhour = 0, rainoff = 0,
                          lamb = 0, IUV = 0, soilgrids = 0, message = 0, fail = nyears * 24 * 365,
                          spatial = NA, save = 0, snowcond = 0, intercept = 0/100 * 0.3,
                          grasshade = 0, coastal = T) {
  ystart <- as.numeric(substr(dstart, 7, 10))
  yfinish <- as.numeric(substr(dfinish, 7, 10))
  yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
  if (is.numeric(loc) == FALSE) {
    if (!requireNamespace("dismo", quietly = TRUE)) {
      stop("dismo needed for the place name geocode function to work. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("XML", quietly = TRUE)) {
      stop("XML needed for the place name geocode function to work. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("proj4", quietly = TRUE)) {
      stop("package 'proj4' is needed. Please install it.",
           call. = FALSE)
    }
    longlat <- dismo::geocode(loc)[3:4]
    if (nrow(longlat > 1)) {
      longlat <- longlat[1, ]
    }
    x <- t(as.matrix(as.numeric(c(longlat[1, 1], longlat[1,
                                                         2]))))
  } else {
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1], loc[2]))))
  }
  errors <- 0
  if (as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900") >
      Sys.time() | as.POSIXct(dstart, format = "%d/%m/%Y",
                              origin = "01/01/1900") > Sys.time()) {
    cat("sorry, no NCEP data for these times - please choose a different date range \n")
  }
  if (DEP[2] - DEP[1] > 3 | DEP[3] - DEP[2] > 3) {
    cat("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if (DEP[2] - DEP[1] < 2) {
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if (is.numeric(loc[1])) {
    if (loc[1] > 180 | loc[2] > 90) {
      cat("ERROR: Latitude or longitude (longlat) is out of bounds.\n          Please enter a correct value.",
          "\n")
      errors <- 1
    }
  }
  if (run.gads %in% c(0, 1) == FALSE) {
    cat("ERROR: the variable 'run.gads' be either 0 or 1.\n        Please correct.",
        "\n")
    errors <- 1
  }
  if (write_input %in% c(0, 1) == FALSE) {
    cat("ERROR: the variable 'write_input' be either 0 or 1.\n        Please correct.",
        "\n")
    errors <- 1
  }
  if (EC < 0.0034 | EC > 0.058) {
    cat("ERROR: the eccentricity variable (EC) is out of bounds.\n        Please enter a correct value (0.0034 - 0.058).",
        "\n")
    errors <- 1
  }
  if (RUF < 1e-04) {
    cat("ERROR: The roughness height (RUF) is too small ( < 0.0001).\n        Please enter a larger value.",
        "\n")
    errors <- 1
  }
  if (RUF > 2) {
    cat("ERROR: The roughness height (RUF) is too large ( > 2).\n        Please enter a smaller value.",
        "\n")
    errors <- 1
  }
  if (DEP[1] != 0) {
    cat("ERROR: First soil node (DEP[1]) must = 0 cm.\n        Please correct",
        "\n")
    errors <- 1
  }
  if (length(DEP) != 10) {
    cat("ERROR: You must enter 10 different soil depths.",
        "\n")
    errors <- 1
  }
  for (i in 1:9) {
    if (DEP[i + 1] <= DEP[i]) {
      cat("ERROR: Soil depth (DEP array) is not in ascending size",
          "\n")
      errors <- 1
    }
  }
  if (DEP[10] > 500) {
    cat("ERROR: Deepest soil depth (DEP array) is too large (<=500 cm)",
        "\n")
    errors <- 1
  }
  if (Thcond < 0) {
    cat("ERROR: Thermal variable conductivity (THCOND) is negative.\n        Please input a positive value.",
        "\n")
    errors <- 1
  }
  if (Density < 0) {
    cat("ERROR: Density variable (Density) is negative.\n        Please input a positive value.",
        "\n")
    errors <- 1
  }
  if (SpecHeat < 0) {
    cat("ERROR: Specific heat variable (SpecHeat) is negative.\n        Please input a positive value.",
        "\n")
    errors <- 1
  }
  if (min(BulkDensity) < 0) {
    cat("ERROR: Bulk density value (BulkDensity) is negative.\n        Please input a positive value.",
        "\n")
    errors <- 1
  }
  if (REFL < 0 | REFL > 1) {
    cat("ERROR: Soil reflectivity value (REFL) is out of bounds.\n        Please input a value between 0 and 1.",
        "\n")
    errors <- 1
  }
  if (is.na(slope) == FALSE & slope > 90) {
    cat("ERROR: Slope value (slope) is out of bounds.\n        Please input a value between 0 and 90.",
        "\n")
    errors <- 1
  }
  if (is.na(aspect) == FALSE & aspect > 365) {
    cat("ERROR: Aspect value (aspect) is out of bounds.\n        Please input a value between 0 and 365.",
        "\n")
    errors <- 1
  }
  if (max(hori) > 90 | min(hori) < 0) {
    cat("ERROR: At least one of your horizon angles (hori) is out of bounds.\n        Please input a value between 0 and 90",
        "\n")
    errors <- 1
  }
  if (length(hori) != 24) {
    cat("ERROR: You must enter 24 horizon angle values.",
        "\n")
    errors <- 1
  }
  if (SLE < 0.5 | SLE > 1) {
    cat("ERROR: Emissivity (SLE) is out of bounds.\n        Please enter a correct value (0.05 - 1.00).",
        "\n")
    errors <- 1
  }
  if (ERR < 0) {
    cat("ERROR: Error bound (ERR) is too small.\n        Please enter a correct value (> 0.00).",
        "\n")
    errors <- 1
  }
  if (Usrhyt < RUF) {
    cat("ERROR: Reference height (Usrhyt) smaller than roughness height (RUF).\n        Please use a larger height above the surface.",
        "\n")
    errors <- 1
  }
  if (Usrhyt < 0.005 | Usrhyt > Refhyt) {
    cat("ERROR: Local height (Usrhyt) is out of bounds.\n        Please enter a correct value (0.005 - Refhyt).",
        "\n")
    errors <- 1
  }
  if (CMH2O < 0.5 | CMH2O > 2) {
    cat("ERROR: Preciptable water in air column (CMH2O) is out of bounds.\n        Please enter a correct value (0.1 - 2cm).",
        "\n")
    errors <- 1
  }
  if (errors == 0) {
    moist.lai <- 0.1
    tme <- as.POSIXlt(seq(ISOdate(ystart, 1, 1, tz = "UTC") -
                            3600 * 12, ISOdate((ystart + nyears), 1, 1, tz = "UTC") -
                            3600 * 13, by = "days"))
    doys12 <- c(15, 46, 74, 105, 135, 166, 196, 227, 258,
                288, 319, 349)
    microdaily <- 1
    daystart <- 1
    idayst <- 1
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("package 'raster' is needed. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("RNCEP", quietly = TRUE)) {
      stop("package 'RNCEP' is needed. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.",
           call. = FALSE)
    }
    if (is.numeric(loc) == FALSE) {
      if (!requireNamespace("dismo", quietly = TRUE)) {
        stop("dismo needed for the place name geocode function to work. Please install it.",
             call. = FALSE)
      }
      if (!requireNamespace("XML", quietly = TRUE)) {
        stop("XML needed for the place name geocode function to work. Please install it.",
             call. = FALSE)
      }
      if (!requireNamespace("proj4", quietly = TRUE)) {
        stop("package 'proj4' is needed. Please install it.",
             call. = FALSE)
      }
      longlat <- dismo::geocode(loc)[3:4]
      if (nrow(longlat > 1)) {
        longlat <- longlat[1, ]
      }
      x <- t(as.matrix(as.numeric(c(longlat[1, 1], longlat[1,
                                                           2]))))
    } else {
      longlat <- loc
      x <- t(as.matrix(as.numeric(c(loc[1], loc[2]))))
    }
    ALREF <- abs(trunc(x[1]))
    HEMIS <- ifelse(x[2] < 0, 2, 1)
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2]) - ALAT) * 60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1]) - ALONG) * 60
    azmuth <- aspect
    lat <- as.numeric(longlat[2])
    long <- as.numeric(longlat[1])
    loc <- c(long, lat)
    if (class(dem)[1] == "RasterLayer") {
      cat("using DEM provided to function call \n")
    }
    if (save != 2 & class(dem)[1] != "RasterLayer") {
      cat("downloading DEM via package elevatr \n")
      dem <- get_dem(lat = lat, long = long)
    }
    if (save == 1) {
      save(dem, file = "dem.Rda")
    }
    if (save == 2) {
      load("dem.Rda")
    }
    if (save != 2) {
      if (soilgrids == 1) {
        cat("extracting data from SoilGrids \n")
        require(rjson)
        require(sp)
        require(GSIF)
        pnts <- data.frame(lon = x[1], lat = x[2], id = c("p1"))
        coordinates(pnts) <- ~lon + lat
        proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")
        soilgrids.r <- REST.SoilGrids(c("BLDFIE", "SLTPPT",
                                        "SNDPPT", "CLYPPT"))
        ov <- over(soilgrids.r, pnts)
        if (length(ov) > 3) {
          soilpro <- cbind(c(0, 5, 15, 30, 60, 100, 200),
                           t(ov[3:9])/1000, t(ov[11:17]), t(ov[19:25]),
                           t(ov[27:33]))
          colnames(soilpro) <- c("depth", "blkdens",
                                 "clay", "silt", "sand")
          soil.hydro <- pedotransfer(soilpro = as.data.frame(soilpro),
                                     DEP = DEP)
          PE <- soil.hydro$PE
          BB <- soil.hydro$BB
          BD <- soil.hydro$BD
          KS <- soil.hydro$KS
          BulkDensity <- BD[seq(1, 19, 2)]
        } else {
          cat("no SoilGrids data for this site, using user-input soil properties \n")
        }
      }
    } else {
      if (soilgrids == 1) {
        cat("loading SoilGrids data from previous run \n")
        load("PE.Rda")
        load("BB.Rda")
        load("BD.Rda")
        load("KS.Rda")
        load("BulkDensity.Rda")
      }
    }
    if (save == 1 & soilgrids == 1) {
      cat("saving SoilGrids data for later \n")
      save(PE, file = "PE.Rda")
      save(BB, file = "BB.Rda")
      save(BD, file = "BD.Rda")
      save(KS, file = "KS.Rda")
      save(BulkDensity, file = "BulkDensity.Rda")
    }
    hori <- rep(0, 24)
    VIEWF <- 1
    days <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900"),
                as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900"),
                by = "days")
    alldays <- seq(as.POSIXct("01/01/1900", format = "%d/%m/%Y",
                              origin = "01/01/1900"), Sys.time() - 60 * 60 * 24,
                   by = "days")
    startday <- which(as.character(format(alldays, "%d/%m/%Y")) ==
                        format(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900"),
                               "%d/%m/%Y"))
    endday <- which(as.character(format(alldays, "%d/%m/%Y")) ==
                      format(as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900"),
                             "%d/%m/%Y"))
    countday <- endday - startday + 1
    tt <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UTC"),
              as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UTC") +
                23 * 3600, by = "hours")
    dates2 <- seq(as.POSIXct(dstart, format = "%d/%m/%Y",
                             tz = "UTC"), as.POSIXct(dfinish, format = "%d/%m/%Y",
                                                     tz = "UTC") + 23 * 3600, by = "days")
    if (save == 2) {
      load("tref.Rda")
      load("SLOPE.Rda")
      load("ASPECT.Rda")
      load("HORIZON.Rda")
      elev <- tref$elev[1]
      ALTT <- elev
    }
    if (save != 2) {
      cat("extracting weather data with RNCEP \n")
      ncquery <- function(filename, var, start, count,
                          year) {
        nc <- nc_open(paste(spatial, "/", filename, year,
                            ".nc", sep = ""))
        out <- as.numeric(ncvar_get(nc, varid = var,
                                    start = start, count = count))
        nc_close(nc)
        out
      }
      if (is.na(spatial) == FALSE) {
        sorttimes <- function(tme) {
          tme2 <- as.POSIXct(tme)
          tme2 <- c((tme2 - 24 * 3600), tme2, (tme2 +
                                                 24 * 3600), (tme2 + 48 * 3600))
          tme2 <- as.POSIXlt(tme2)
          yrs <- unique(tme2$year + 1900)
          mths <- unique(tme2$mon + 1)
          yrs <- yrs[order(yrs)]
          mths <- mths[order(mths)]
          tma <- 0
          dm <- c(31, 28, 31, 30, 31, 30, 31, 31, 30,
                  31, 30, 31) * 4 - 1
          for (yr in 1:length(yrs)) {
            dm[2] <- ifelse(yrs[yr]%%4 == 0, 115, 111)
            for (mth in 1:length(mths)) {
              tmym <- as.POSIXct(c(0:dm[mth]) * 3600 *
                                   6, origin = paste0(yrs[yr], "-", mths[mth],
                                                      "-01 00:00"), tz = "GMT")
              tma <- c(tma, tmym)
            }
          }
          tma <- as.POSIXlt(tma[-1], origin = "1970-01-01",
                            tz = "GMT")
          sel <- which(tma >= min(tme2) & tma <= max(tme2))
          sel <- sel[-length(sel)]
          return(list(tme = tma, sel = sel))
        }
        tme2 <- sorttimes(tme)$tme
        sel <- sorttimes(tme)$sel
        cat(paste0("extracting weather data locally from ",
                   spatial, " \n"))
        years <- as.numeric(unique(format(tme, "%Y")))
        nyears <- length(years)
        nc <- nc_open(paste(spatial, "/air.2m.gauss.",
                            years[1], ".nc", sep = ""))
        lon2 <- matrix(ncvar_get(nc, "lon"))
        lat2 <- matrix(ncvar_get(nc, "lat"))
        lon_1 <- long
        if (lon_1 < 0) {
          lon_1 <- 180 - (long * -1) + 180
        }
        lat_1 <- lat
        dist1 <- abs(lon2 - lon_1)
        index1 <- which.min(dist1)
        dist2 <- abs(lat2 - lat_1)
        index2 <- which.min(dist2)
        start <- c(index1, index2, 1, 1)
        count <- c(1, 1, 1, -1)
        start2 <- c(index1, index2, 1, 1460 - 3)
        count2 <- c(1, 1, 1, 4)
        nc_close(nc)
        for (j in 1:(nyears + 2)) {
          if (j == 1) {
            Tkmin <- ncquery("tmin.2m.gauss.", "tmin",
                             start2, count2, years[j] - 1)
            Tkmax <- ncquery("tmax.2m.gauss.", "tmax",
                             start2, count2, years[j] - 1)
            Tk <- ncquery("air.2m.gauss.", "air", start2,
                          count2, years[j] - 1)
            sh <- ncquery("shum.2m.gauss.", "shum", start2,
                          count2, years[j] - 1)
            pr <- ncquery("pres.sfc.gauss.", "pres",
                          start2[c(1, 2, 4)], count2[c(1, 2, 4)],
                          years[j] - 1)
            tcdc <- ncquery("tcdc.eatm.gauss.", "tcdc",
                            start2[c(1, 2, 4)], count2[c(1, 2, 4)],
                            years[j] - 1)
            dsw <- ncquery("dswrf.sfc.gauss.", "dswrf",
                           start2[c(1, 2, 4)], count2[c(1, 2, 4)],
                           years[j] - 1)
            dlw <- ncquery("dlwrf.sfc.gauss.", "dlwrf",
                           start2[c(1, 2, 4)], count2[c(1, 2, 4)],
                           years[j] - 1)
            ulw <- ncquery("ulwrf.sfc.gauss.", "ulwrf",
                           start2[c(1, 2, 4)], count2[c(1, 2, 4)],
                           years[j] - 1)
            wu <- ncquery("uwnd.10m.gauss.", "uwnd",
                          start2, count2, years[j] - 1)
            wv <- ncquery("vwnd.10m.gauss.", "vwnd",
                          start2, count2, years[j] - 1)
            prate <- ncquery("prate.sfc.gauss.", "prate",
                             start2[c(1, 2, 4)], count2[c(1, 2, 4)],
                             years[j] - 1)
          } else {
            if (j <= nyears + 1) {
              cat(paste("reading weather input for ",
                        years[j - 1], " \n", sep = ""))
              Tkmin <- c(Tkmin, ncquery("tmin.2m.gauss.",
                                        "tmin", start, count, years[j - 1]))
              Tkmax <- c(Tkmax, ncquery("tmax.2m.gauss.",
                                        "tmax", start, count, years[j - 1]))
              Tk <- c(Tk, ncquery("air.2m.gauss.", "air",
                                  start, count, years[j - 1]))
              sh <- c(sh, ncquery("shum.2m.gauss.", "shum",
                                  start, count, years[j - 1]))
              pr <- c(pr, ncquery("pres.sfc.gauss.",
                                  "pres", start[c(1, 2, 4)], count[c(1,
                                                                     2, 4)], years[j - 1]))
              tcdc <- c(tcdc, ncquery("tcdc.eatm.gauss.",
                                      "tcdc", start[c(1, 2, 4)], count[c(1,
                                                                         2, 4)], years[j - 1]))
              dsw <- c(dsw, ncquery("dswrf.sfc.gauss.",
                                    "dswrf", start[c(1, 2, 4)], count[c(1,
                                                                        2, 4)], years[j - 1]))
              dlw <- c(dlw, ncquery("dlwrf.sfc.gauss.",
                                    "dlwrf", start[c(1, 2, 4)], count[c(1,
                                                                        2, 4)], years[j - 1]))
              ulw <- c(ulw, ncquery("ulwrf.sfc.gauss.",
                                    "ulwrf", start[c(1, 2, 4)], count[c(1,
                                                                        2, 4)], years[j - 1]))
              wu <- c(wu, ncquery("uwnd.10m.gauss.",
                                  "uwnd", start, count, years[j - 1]))
              wv <- c(wv, ncquery("vwnd.10m.gauss.",
                                  "vwnd", start, count, years[j - 1]))
              prate <- c(prate, ncquery("prate.sfc.gauss.",
                                        "prate", start[c(1, 2, 4)], count[c(1,
                                                                            2, 4)], years[j - 1]))
            } else {
              Tkmin <- c(Tkmin, ncquery("tmin.2m.gauss.",
                                        "tmin", start, count2, years[j - 2] +
                                          1))
              Tkmax <- c(Tkmax, ncquery("tmax.2m.gauss.",
                                        "tmax", start, count2, years[j - 2] +
                                          1))
              Tk <- c(Tk, ncquery("air.2m.gauss.", "air",
                                  start, count2, years[j - 2] + 1))
              sh <- c(sh, ncquery("shum.2m.gauss.", "shum",
                                  start, count2, years[j - 2] + 1))
              pr <- c(pr, ncquery("pres.sfc.gauss.",
                                  "pres", start[c(1, 2, 4)], count2[c(1,
                                                                      2, 4)], years[j - 2] + 1))
              tcdc <- c(tcdc, ncquery("tcdc.eatm.gauss.",
                                      "tcdc", start[c(1, 2, 4)], count2[c(1,
                                                                          2, 4)], years[j - 2] + 1))
              dsw <- c(dsw, ncquery("dswrf.sfc.gauss.",
                                    "dswrf", start[c(1, 2, 4)], count2[c(1,
                                                                         2, 4)], years[j - 2] + 1))
              dlw <- c(dlw, ncquery("dlwrf.sfc.gauss.",
                                    "dlwrf", start[c(1, 2, 4)], count2[c(1,
                                                                         2, 4)], years[j - 2] + 1))
              ulw <- c(ulw, ncquery("ulwrf.sfc.gauss.",
                                    "ulwrf", start[c(1, 2, 4)], count2[c(1,
                                                                         2, 4)], years[j - 2] + 1))
              wu <- c(wu, ncquery("uwnd.10m.gauss.",
                                  "uwnd", start, count2, years[j - 2] +
                                    1))
              wv <- c(wv, ncquery("vwnd.10m.gauss.",
                                  "vwnd", start, count2, years[j - 2] +
                                    1))
              prate <- c(prate, ncquery("prate.sfc.gauss.",
                                        "prate", start[c(1, 2, 4)], count2[c(1,
                                                                             2, 4)], years[j - 2] + 1))
            }
          }
        }
        dsw[dsw < 0] <- 0
        prate[prate < 0] <- 0
        prate <- prate * 3600 * 6
        ncepdata <- data.frame(obs_time = tme2[sel],
                               Tk, Tkmin, Tkmax, sh, pr, wu, wv, dlw, ulw,
                               dsw, tcdc)
        hourlydata <- hourlyNCEP(ncepdata = ncepdata,
                                 lat, long, tme, TRUE)
        microclima.out <- microclimaforNMR(lat = longlat[2],
                                           long = longlat[1], dstart = dstart, dfinish = dfinish,
                                           l = mean(LAI), x = LOR, coastal = coastal, hourlydata = hourlydata,
                                           dailyprecip = prate, dem = dem, demmeso = dem2, albr = REFL,
                                           resolution = 30, zmin = 0, slope = slope, aspect = aspect,
                                           windthresh = 4.5, emthresh = 0.78)
        dailyrain <- microclima.out$dailyprecip[-c(1:4)]
        dailyrain <- dailyrain[1:(length(dailyrain) -
                                    4)]
        dailyrain <- aggregate(dailyrain, by = list(format(hourlydata$obs_time[seq(1,
                                                                                   nrow(hourlydata), 6)], "%Y-%m-%d")), sum)$x
      } else {
        microclima.out <- microclimaforNMR(lat = longlat[2],
                                           long = longlat[1], dstart = dstart, dfinish = dfinish,
                                           l = mean(LAI), x = LOR, coastal = coastal, hourlydata = NA, dailyprecip = NA,
                                           dem = dem, demmeso = dem2, albr = REFL, resolution = 30, zmin = 0,
                                           slope = slope, aspect = aspect, windthresh = 4.5,
                                           emthresh = 0.78)
        hourlydata <- microclima.out$hourlydata
        dailyrain <- microclima.out$dailyprecip
      }
      hourlyradwind <- microclima.out$hourlyradwind
      tref <- microclima.out$tref
      ZENhr <- hourlydata$szenith
      ZENhr[ZENhr > 90] <- 90
      cat("computing radiation and elevation effects with package microclima \n")
      SLOPE <- hourlyradwind$slope[1]
      ASPECT <- hourlyradwind$aspect[1]
      HORIZON <- hori
      if (save == 1) {
        save(SLOPE, file = "SLOPE.Rda")
        save(ASPECT, file = "ASPECT.Rda")
        save(HORIZON, file = "HORIZON.Rda")
      }
      if (save == 1) {
        save(tref, file = "tref.Rda")
      }
      elev <- tref$elev[1]
      ALTT <- elev
      TAIRhr <- tref$tref + tref$telev + tref$tcad
      SOLRhr <- hourlyradwind$swrad/0.0036
      SOLRhr[SOLRhr < 0] <- 0
      CLDhr <- hourlydata$cloudcover
      CLDhr[CLDhr < 0] <- 0
      CLDhr[CLDhr > 100] <- 100
      IRDhr <- hourlydata$downlong/0.0036
      RHhr <- suppressWarnings(humidityconvert(h = hourlydata$humidity,
                                               intype = "specific", p = hourlydata$pressure,
                                               tc = TAIRhr)$relative)
      RHhr[RHhr > 100] <- 100
      RHhr[RHhr < 0] <- 0
      WNhr <- hourlyradwind$windspeed
      WNhr[is.na(WNhr)] <- 0.1
      RAINhr <- WNhr * 0
      PRESShr <- hourlydata$pressure
      RAINFALL <- dailyrain
      RAINFALL[RAINFALL < 0.1] <- 0
      ZENhr2 <- ZENhr
      ZENhr2[ZENhr2 != 90] <- 0
      dmaxmin <- function(x, fun) {
        dx <- t(matrix(x, nrow = 24))
        apply(dx, 1, fun)
      }
      TMAXX <- dmaxmin(TAIRhr, max)
      TMINN <- dmaxmin(TAIRhr, min)
      CCMAXX <- dmaxmin(CLDhr, max)
      CCMINN <- dmaxmin(CLDhr, min)
      RHMAXX <- dmaxmin(RHhr, max)
      RHMINN <- dmaxmin(RHhr, min)
      WNMAXX <- dmaxmin(WNhr, max)
      WNMINN <- dmaxmin(WNhr, min)
      PRESS <- dmaxmin(PRESShr, min)
      if (save == 1) {
        cat("saving met data for later \n")
        save(CCMAXX, file = "CCMAXX.Rda")
        save(CCMINN, file = "CCMINN.Rda")
        save(WNMAXX, file = "WNMAXX.Rda")
        save(WNMINN, file = "WNMINN.Rda")
        save(TMAXX, file = "TMAXX.Rda")
        save(TMINN, file = "TMINN.Rda")
        save(RHMAXX, file = "RHMAXX.Rda")
        save(RHMINN, file = "RHMINN.Rda")
        save(RAINFALL, file = "RAINFALL.Rda")
        save(PRESS, file = "PRESS.Rda")
        save(CLDhr, file = "CLDhr.Rda")
        save(WNhr, file = "WNhr.Rda")
        save(TAIRhr, file = "TAIRhr.Rda")
        save(RHhr, file = "RHhr.Rda")
        save(RAINhr, file = "RAINhr.Rda")
        save(SOLRhr, file = "SOLRhr.Rda")
        save(ZENhr, file = "ZENhr.Rda")
        save(IRDhr, file = "IRDhr.Rda")
        save(microclima.out, file = "microclima.out.Rda")
      }
    } else {
      cat("loading met data from previous run \n")
      load("CCMAXX.Rda")
      load("CCMINN.Rda")
      load("WNMAXX.Rda")
      load("WNMINN.Rda")
      load("TMAXX.Rda")
      load("TMINN.Rda")
      load("RHMAXX.Rda")
      load("RHMINN.Rda")
      load("RAINFALL.Rda")
      load("PRESS.Rda")
      load("CLDhr.Rda")
      load("WNhr.Rda")
      load("TAIRhr.Rda")
      load("RHhr.Rda")
      load("RAINhr.Rda")
      load("SOLRhr.Rda")
      load("ZENhr.Rda")
      load("IRDhr.Rda")
      load("microclima.out.Rda")
    }
    slope <- 0
    azmuth <- 0
    ndays <- length(TMAXX)
    doynum <- ndays
    leapyears <- seq(1900, 2100, 4)
    for (k in 1:nyears) {
      if (k == 1) {
        cyear <- ystart
      } else {
        cyear <- cyear + 1
      }
      if (cyear %in% leapyears) {
        dinyear <- 366
      } else {
        dinyear <- 365
      }
      if (k == 1) {
        doy <- seq(1, dinyear)
      } else {
        doy <- c(doy, seq(1, dinyear))
      }
    }
    ida <- ndays
    idayst <- 1
    maxshades <- rep(0.1, ndays)
    minshades <- rep(0, ndays)
    shademax <- maxshades
    maxshade <- 0.1
    minshade <- 0
    if (run.gads == 1) {
      relhum <- 1
      optdep.summer <- as.data.frame(rungads(longlat[2],
                                             longlat[1], relhum, 0))
      optdep.winter <- as.data.frame(rungads(longlat[2],
                                             longlat[1], relhum, 1))
      optdep <- cbind(optdep.winter[, 1], rowMeans(cbind(optdep.summer[,
                                                                       2], optdep.winter[, 2])))
      optdep <- as.data.frame(optdep)
      colnames(optdep) <- c("LAMBDA", "OPTDEPTH")
      a <- lm(OPTDEPTH ~ poly(LAMBDA, 6, raw = TRUE), data = optdep)
      LAMBDA <- c(290, 295, 300, 305, 310, 315, 320, 330,
                  340, 350, 360, 370, 380, 390, 400, 420, 440,
                  460, 480, 500, 520, 540, 560, 580, 600, 620,
                  640, 660, 680, 700, 720, 740, 760, 780, 800,
                  820, 840, 860, 880, 900, 920, 940, 960, 980,
                  1000, 1020, 1080, 1100, 1120, 1140, 1160, 1180,
                  1200, 1220, 1240, 1260, 1280, 1300, 1320, 1380,
                  1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580,
                  1600, 1620, 1640, 1660, 1700, 1720, 1780, 1800,
                  1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120,
                  2150, 2200, 2260, 2300, 2320, 2350, 2380, 2400,
                  2420, 2450, 2490, 2500, 2600, 2700, 2800, 2900,
                  3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700,
                  3800, 3900, 4000)
      TAI <- predict(a, data.frame(LAMBDA))
    } else {
      TAI <- c(0.0670358341290886, 0.0662612704779235,
               0.065497075238002, 0.0647431301168489, 0.0639993178022531,
               0.0632655219571553, 0.0625416272145492, 0.0611230843885423,
               0.0597427855962549, 0.0583998423063099, 0.0570933810229656,
               0.0558225431259535, 0.0545864847111214, 0.0533843764318805,
               0.0522154033414562, 0.0499736739981675, 0.047855059159556,
               0.0458535417401334, 0.0439633201842001, 0.0421788036108921,
               0.0404946070106968, 0.0389055464934382, 0.0374066345877315,
               0.0359930755919066, 0.0346602609764008, 0.0334037648376212,
               0.0322193394032758, 0.0311029105891739, 0.0300505736074963,
               0.0290585886265337, 0.0281233764818952, 0.0272415144391857,
               0.0264097320081524, 0.0256249068083005, 0.0248840604859789,
               0.0241843546829336, 0.0235230870563317, 0.0228976873502544,
               0.0223057135186581, 0.0217448478998064, 0.0212128934421699,
               0.0207077699817964, 0.0202275105711489, 0.0197702578594144,
               0.0193342605242809, 0.0189178697551836, 0.0177713140039894,
               0.0174187914242432, 0.0170790495503944, 0.0167509836728154,
               0.0164335684174899, 0.0161258546410128, 0.0158269663770596,
               0.0155360978343254, 0.0152525104459325, 0.0149755299703076,
               0.0147045436435285, 0.0144389973831391, 0.0141783930434343,
               0.0134220329447663, 0.0131772403830191, 0.0129356456025128,
               0.0126970313213065, 0.0124612184223418, 0.0122280636204822,
               0.01199745718102, 0.0115436048739351, 0.0110993711778668,
               0.0108808815754663, 0.0106648652077878, 0.0104513876347606,
               0.0102405315676965, 0.00982708969547694, 0.00962473896278535,
               0.00903679230300494, 0.00884767454432418, 0.0083031278398166,
               0.00796072474935954, 0.00755817587626185, 0.00718610751850881,
               0.00704629977586921, 0.00684663903049612, 0.00654155580333479,
               0.00642947339729728, 0.00627223096874308, 0.00603955966866779,
               0.00580920937536261, 0.00568506186880564, 0.00563167068287251,
               0.00556222005081865, 0.00550522989971023, 0.00547395763028062,
               0.0054478983436216, 0.00541823364504573, 0.00539532163908382,
               0.00539239864119488, 0.00541690124712384, 0.00551525885358836,
               0.00564825853509463, 0.00577220185074264, 0.00584222986640171,
               0.00581645238345584, 0.00566088137411449, 0.00535516862329704,
               0.00489914757707667, 0.00432017939770409, 0.0036813032251836,
               0.00309019064543606, 0.00270890436501562, 0.00276446109239711,
               0.00356019862584603)
    }
    if (warm != 0) {
      TMAXX <- TMAXX + seq(0, ndays - 1)/(ndays - 1) *
        warm
      TMINN <- TMINN + seq(0, ndays - 1)/(ndays - 1) *
        warm
      TAIRhr <- TAIRhr + seq(0, ndays - 1)/(ndays - 1) *
        warm
    }
    RAINFALL <- RAINFALL + rainoff
    RAINhr <- RAINhr + rainoff
    ALLMINTEMPS <- TMINN
    ALLMAXTEMPS <- TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)
    WNMAXX <- WNMAXX * windfac
    WNMINN <- WNMINN * windfac
    WNhr <- WNhr * windfac
    MAXSHADES <- maxshades
    MINSHADES <- minshades
    REFLS <- rep(REFL, ndays)
    PCTWET <- rep(PCTWET, ndays)
    soilwet <- RAINFALL
    soilwet[soilwet <= rainwet] = 0
    soilwet[soilwet > 0] = 90
    PCTWET <- pmax(soilwet, PCTWET)
    Intrvls <- rep(0, ndays)
    Intrvls[1] <- 1
    Numtyps <- 1
    Numint <- 1
    Nodes <- matrix(data = 0, nrow = 10, ncol = ndays)
    Nodes[1, 1] <- 10
    ALREF <- abs(trunc(x[1]))
    HEMIS <- ifelse(x[2] < 0, 2, 1)
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2]) - ALAT) * 60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1]) - ALONG) * 60
    avetemp <- (sum(TMAXX) + sum(TMINN))/(length(TMAXX) *
                                            2)
    soilinit <- rep(avetemp, 20)
    tannul <- mean(unlist(ALLTEMPS))
    if (nyears == 1) {
      avetemp <- (sum(TMAXX) + sum(TMINN))/(length(TMAXX) *
                                              2)
      tannulrun <- rep(avetemp, ndays)
    } else {
      avetemp <- rowMeans(cbind(TMAXX, TMINN), na.rm = TRUE)
      if (length(TMAXX) < 365) {
        tannulrun <- rep((sum(TMAXX) + sum(TMINN))/(length(TMAXX) *
                                                      2), length(TMAXX))
      } else {
        tannulrun <- movingFun(avetemp, n = 365, fun = mean,
                               type = "to")
        yearone <- rep((sum(TMAXX[1:365]) + sum(TMINN[1:365]))/(365 *
                                                                  2), 365)
        tannulrun[1:365] <- yearone
      }
    }
    SLES <- matrix(nrow = ndays, data = 0)
    SLES <- SLES + SLE
    moists2 <- matrix(nrow = 10, ncol = ndays, data = 0)
    moists2[1, ndays] <- 0.2
    moists <- moists2
    if (runmoist == 1) {
      moists2 <- matrix(nrow = 10, ncol = ndays, data = 0)
      moists2[1:10, ] <- SoilMoist_Init
      moists <- moists2
    }
    soilprops <- matrix(data = 0, nrow = 10, ncol = 5)
    soilprops[, 1] <- BulkDensity
    soilprops[, 2] <- min(0.26, 1 - BulkDensity/Density)
    soilprops[, 3] <- Thcond
    soilprops[, 4] <- SpecHeat
    soilprops[, 5] <- Density
    if (cap == 1) {
      soilprops[1:2, 3] <- 0.2
      soilprops[1:2, 4] <- 1920
    }
    if (cap == 2) {
      soilprops[1:2, 3] <- 0.1
      soilprops[3:4, 3] <- 0.25
      soilprops[1:4, 4] <- 1920
      soilprops[1:4, 5] <- 1.3
      soilprops[1:4, 1] <- 0.7
    }
    ALTT <- as.numeric(ALTT)
    ALREF <- as.numeric(ALREF)
    ALMINT <- as.numeric(ALMINT)
    ALONG <- as.numeric(ALONG)
    AMINUT <- as.numeric(AMINUT)
    ALAT <- as.numeric(ALAT)
    runshade <- 0
    IR <- 0
    microinput <- c(ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps,
                    Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT,
                    ALONG, ALMINT, ALREF, slope, azmuth, ALTT, CMH2O,
                    microdaily, tannul, EC, VIEWF, snowtemp, snowdens,
                    snowmelt, undercatch, rainmult, runshade, runmoist,
                    maxpool, evenrain, snowmodel, rainmelt, writecsv,
                    densfun, hourly, rainhourly, lamb, IUV, RW, PC, RL,
                    SP, R1, IM, MAXCOUNT, IR, message, fail, snowcond,
                    intercept, grasshade)
    if (hourly == 0) {
      TAIRhr = rep(0, 24 * ndays)
      RHhr = rep(0, 24 * ndays)
      WNhr = rep(0, 24 * ndays)
      CLDhr = rep(0, 24 * ndays)
      SOLRhr = rep(0, 24 * ndays)
      ZENhr = rep(-1, 24 * ndays)
      IRDhr = rep(-1, 24 * ndays)
    } else {
      CLDhr = rep(0, 24 * ndays)
    }
    if (rainhourly == 0) {
      RAINhr = rep(0, 24 * ndays)
    } else {
      RAINhr = RAINhr
    }
    doy1 = matrix(data = 0, nrow = ndays, ncol = 1)
    SLES1 = matrix(data = 0, nrow = ndays, ncol = 1)
    MAXSHADES1 = matrix(data = 0, nrow = ndays, ncol = 1)
    MINSHADES1 = matrix(data = 0, nrow = ndays, ncol = 1)
    TMAXX1 = matrix(data = 0, nrow = ndays, ncol = 1)
    TMINN1 = matrix(data = 0, nrow = ndays, ncol = 1)
    CCMAXX1 = matrix(data = 0, nrow = ndays, ncol = 1)
    CCMINN1 = matrix(data = 0, nrow = ndays, ncol = 1)
    RHMAXX1 = matrix(data = 0, nrow = ndays, ncol = 1)
    RHMINN1 = matrix(data = 0, nrow = ndays, ncol = 1)
    WNMAXX1 = matrix(data = 0, nrow = ndays, ncol = 1)
    WNMINN1 = matrix(data = 0, nrow = ndays, ncol = 1)
    REFLS1 = matrix(data = 0, nrow = ndays, ncol = 1)
    PCTWET1 = matrix(data = 0, nrow = ndays, ncol = 1)
    RAINFALL1 = matrix(data = 0, nrow = ndays, ncol = 1)
    tannul1 = matrix(data = 0, nrow = ndays, ncol = 1)
    moists1 = matrix(data = 0, nrow = 10, ncol = ndays)
    doy1[1:ndays] <- doy
    SLES1[1:ndays] <- SLES
    MAXSHADES1[1:ndays] <- MAXSHADES
    MINSHADES1[1:ndays] <- MINSHADES
    TMAXX1[1:ndays] <- TMAXX
    TMINN1[1:ndays] <- TMINN
    CCMAXX1[1:ndays] <- CCMAXX
    CCMINN1[1:ndays] <- CCMINN
    RHMAXX1[1:ndays] <- RHMAXX
    RHMINN1[1:ndays] <- RHMINN
    WNMAXX1[1:ndays] <- WNMAXX
    WNMINN1[1:ndays] <- WNMINN
    REFLS1[1:ndays] <- REFLS
    PCTWET1[1:ndays] <- PCTWET
    RAINFALL1[1:ndays] <- RAINFALL
    tannul1[1:ndays] <- tannul
    moists1[1:10, 1:ndays] <- moists
    LAI <- moist.lai
    if (length(LAI) < ndays) {
      LAI <- rep(LAI[1], ndays)
      LAI1 <- LAI
    }
    if (shore == 0) {
      tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3)
    }
    TIMAXS <- c(1, 1, 0, 0)
    TIMINS <- c(0, 0, 1, 1)
    micro <- list(tides = tides, microinput = microinput,
                  doy = doy, SLES = SLES1, DEP = DEP, Nodes = Nodes,
                  MAXSHADES = MAXSHADES, MINSHADES = MINSHADES, TIMAXS = TIMAXS,
                  TIMINS = TIMINS, TMAXX = TMAXX1, TMINN = TMINN1,
                  RHMAXX = RHMAXX1, RHMINN = RHMINN1, CCMAXX = CCMAXX1,
                  CCMINN = CCMINN1, WNMAXX = WNMAXX1, WNMINN = WNMINN1,
                  TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr,
                  SOLRhr = SOLRhr, RAINhr = RAINhr, ZENhr = ZENhr,
                  IRDhr = IRDhr, REFLS = REFLS1, PCTWET = PCTWET1,
                  soilinit = soilinit, hori = hori, TAI = TAI, soilprops = soilprops,
                  moists = moists1, RAINFALL = RAINFALL1, tannulrun = tannulrun,
                  PE = PE, KS = KS, BB = BB, BD = BD, DD = DD, L = L,
                  LAI = LAI1)
    if (write_input == 1) {
      if (dir.exists("micro csv input") == FALSE) {
        dir.create("micro csv input")
      }
      write.table(as.matrix(microinput), file = "micro csv input/microinput.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(doy, file = "micro csv input/doy.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(SLES, file = "micro csv input/SLES.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(DEP, file = "micro csv input/DEP.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(Nodes, file = "micro csv input/Nodes.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(MAXSHADES, file = "micro csv input/Maxshades.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(MINSHADES, file = "micro csv input/Minshades.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMAXS, file = "micro csv input/TIMAXS.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMINS, file = "micro csv input/TIMINS.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TMAXX, file = "micro csv input/TMAXX.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TMINN, file = "micro csv input/TMINN.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMAXX, file = "micro csv input/RHMAXX.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMINN, file = "micro csv input/RHMINN.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMAXX, file = "micro csv input/CCMAXX.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMINN, file = "micro csv input/CCMINN.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMAXX, file = "micro csv input/WNMAXX.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMINN, file = "micro csv input/WNMINN.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(REFLS, file = "micro csv input/REFLS.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(PCTWET, file = "micro csv input/PCTWET.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(soilinit, file = "micro csv input/soilinit.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(hori, file = "micro csv input/hori.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TAI, file = "micro csv input/TAI.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(soilprops, file = "micro csv input/soilprop.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(moists, file = "micro csv input/moists.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RAINFALL, file = "micro csv input/rain.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(tannulrun, file = "micro csv input/tannulrun.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(PE, file = "micro csv input/PE.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(BD, file = "micro csv input/BD.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(DD, file = "micro csv input/DD.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(BB, file = "micro csv input/BB.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(KS, file = "micro csv input/KS.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(L, file = "micro csv input/L.csv", sep = ",",
                  col.names = NA, qmethod = "double")
      write.table(LAI, file = "micro csv input/LAI.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(tides, file = "micro csv input/tides.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(TAIRhr, file = "micro csv input/TAIRhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RHhr, file = "micro csv input/RHhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(WNhr, file = "micro csv input/WNhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(CLDhr, file = "micro csv input/CLDhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(SOLRhr, file = "micro csv input/SOLRhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(RAINhr, file = "micro csv input/RAINhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(ZENhr, file = "micro csv input/ZENhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
      write.table(IRDhr, file = "micro csv input/IRDhr.csv",
                  sep = ",", col.names = NA, qmethod = "double")
    }
    if (is.numeric(loc[1])) {
      location <- paste("long", loc[1], "lat", loc[2])
    } else {
      location <- loc
    }
    cat(paste("running microclimate model for ", ndays, " days from ",
              tt[1], " to ", tt[length(tt)], " at site ", location,
              "\n"))
    ptm <- proc.time()
    microut <- microclimate(micro)
    print(proc.time() - ptm)
    metout <- microut$metout
    shadmet <- microut$shadmet
    soil <- microut$soil
    shadsoil <- microut$shadsoil
    if (runmoist == 1) {
      soilmoist <- microut$soilmoist
      shadmoist <- microut$shadmoist
      humid <- microut$humid
      shadhumid <- microut$shadhumid
      soilpot <- microut$soilpot
      shadpot <- microut$shadpot
      plant <- microut$plant
      shadplant <- microut$shadplant
    } else {
      soilpot <- soil
      soilmoist <- soil
      shadpot <- soil
      shadmoist <- soil
      humid <- soil
      shadhumid <- soil
      plant <- cbind(soil, soil[, 3:4])
      shadplant <- cbind(soil, soil[, 3:4])
      soilpot[, 3:12] <- 0
      soilmoist[, 3:12] <- 0.5
      shadpot[, 3:12] <- 0
      shadmoist[, 3:12] <- 0.5
      humid[, 3:12] <- 0.99
      shadhumid[, 3:12] <- 0.99
      plant[, 3:14] <- 0
      shadplant[, 3:14] <- 0
    }
    if (snowmodel == 1) {
      sunsnow <- microut$sunsnow
      shdsnow <- microut$shdsnow
    }
    if (max(metout[, 1] == 0)) {
      cat("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP",
          "\n")
    }
    if (lamb == 1) {
      drlam <- as.data.frame(microut$drlam)
      drrlam <- as.data.frame(microut$drrlam)
      srlam <- as.data.frame(microut$srlam)
      if (snowmodel == 1) {
        return(list(soil = soil, shadsoil = shadsoil,
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist,
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid,
                    soilpot = soilpot, shadpot = shadpot, sunsnow = sunsnow,
                    shdsnow = shdsnow, plant = plant, shadplant = shadplant,
                    RAINFALL = RAINFALL, ndays = ndays, elev = ALTT,
                    REFL = REFL[1], MAXSHADES = MAXSHADES, longlat = longlat,
                    nyears = nyears, minshade = minshade, maxshade = maxshade,
                    DEP = DEP, drlam = drlam, drrlam = drrlam,
                    srlam = srlam, SLOPE = SLOPE, ASPECT = ASPECT,
                    HORIZON = HORIZON, dates = tt, dem = dem, dates2 = dates2,
                    microclima.out = microclima.out))
      } else {
        return(list(soil = soil, shadsoil = shadsoil,
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist,
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid,
                    soilpot = soilpot, shadpot = shadpot, plant = plant,
                    shadplant = shadplant, RAINFALL = RAINFALL,
                    ndays = ndays, elev = ALTT, REFL = REFL[1],
                    MAXSHADES = MAXSHADES, longlat = longlat, nyears = nyears,
                    minshade = minshade, maxshade = maxshade, DEP = DEP,
                    drlam = drlam, drrlam = drrlam, srlam = srlam,
                    SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON,
                    dates = tt, dem = dem, dates2 = dates2, microclima.out = microclima.out))
      }
    } else {
      if (snowmodel == 1) {
        return(list(soil = soil, shadsoil = shadsoil,
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist,
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid,
                    soilpot = soilpot, shadpot = shadpot, sunsnow = sunsnow,
                    shdsnow = shdsnow, plant = plant, shadplant = shadplant,
                    RAINFALL = RAINFALL, ndays = ndays, elev = ALTT,
                    REFL = REFL[1], MAXSHADES = MAXSHADES, longlat = longlat,
                    nyears = nyears, minshade = minshade, maxshade = maxshade,
                    DEP = DEP, SLOPE = SLOPE, ASPECT = ASPECT,
                    HORIZON = HORIZON, dates = tt, dem = dem, dates2 = dates2,
                    microclima.out = microclima.out))
      } else {
        return(list(soil = soil, shadsoil = shadsoil,
                    metout = metout, shadmet = shadmet, soilmoist = soilmoist,
                    shadmoist = shadmoist, humid = humid, shadhumid = shadhumid,
                    soilpot = soilpot, shadpot = shadpot, plant = plant,
                    shadplant = shadplant, RAINFALL = RAINFALL,
                    ndays = ndays, elev = ALTT, REFL = REFL[1],
                    MAXSHADES = MAXSHADES, longlat = longlat, nyears = nyears,
                    minshade = minshade, maxshade = maxshade, DEP = DEP,
                    SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON,
                    dates = tt, dem = dem, dates2 = dates2, microclima.out = microclima.out))
      }
    }
  }
}

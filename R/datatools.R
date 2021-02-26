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
#' (see https://mapzen.com/documentation/terrain-tiles/data-sources/ for details
#' and https://github.com/tilezen/joerd/blob/master/docs/images/footprints-preview.png
#' for a map showing the resolution of data available for different parts of the globe).
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
  if (class(r)[1] != "RasterLayer") {
    xy <- data.frame(x = long, y = lat)
    xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = 4326)
	
    if (lat >= -80 & lat <= 84)
      xy <- sf::st_transform(xy, 3395)
    if (lat > 84)
      xy <- sf::st_transform(xy, 3413)
    if (lat < -80)
      xy <- sf::st_transform(xy, 3976)
	  
      e <- extent(c(sf::st_coordinates(xy)[1] - floor(xdims/2) * resolution,
                    sf::st_coordinates(xy)[1] + ceiling(xdims/2) * resolution,
                    sf::st_coordinates(xy)[2] - floor(ydims/2) * resolution,
                    sf::st_coordinates(xy)[2] + ceiling(ydims/2) * resolution))
					
    r <- raster(e)
    res(r) <- resolution
	crs(r) <- CRS(SRS_string = st_crs(xy)$wkt)
    
  } else {
    lat <- latlongfromraster(r)$lat
    long <- latlongfromraster(r)$long
    res(r) <- resolution
  }
  z <- ceiling(log((cos(lat * pi/180) * 2 * pi * 6378137) / (256 * resolution), 2))
  z <- ifelse(z > 14, 14, z)
  p <- as(extent(r), 'SpatialPoints')
  p <- as.data.frame(p)
  xx <- sp::proj4string(r)
  r2 <- elevatr::get_elev_raster(p, z = z, src = "aws", prj = xx)
  r2 <- resample(r2, r)
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
get_NCEP <- function(lat, long, tme, reanalysis2 = FALSE) {
  ncepget1 <- function(climvar, tme2, ll) {
    yrs <- unique(tme2$year + 1900)
    vv <- list()
    for (yr in 1:length(yrs)) {
      sel <- which(tme2$year + 1900 == yrs[yr])
      tme3 <- tme2[sel]
      mths <- unique(tme3$mon + 1)
      v <- NCEP.gather(climvar, level = 'gaussian',
                       years.minmax = c(yrs[yr],yrs[yr]),
                       months.minmax = c(min(mths):max(mths)),
                       lat.southnorth = c(ll$y,ll$y), lon.westeast = c(ll$x,ll$x),
                       reanalysis2 = reanalysis2, return.units = FALSE, status.bar = FALSE)
      if (is.null(dim(v)) == F) {
        latdif <- abs(ll$y - as.numeric(rownames(v)))
        londif <- abs(ll$x%%360 - as.numeric(colnames(v)))
        vv[[yr]] <- as.numeric(v[which.min(latdif), which.min(londif),])
      } else vv[[yr]] <- as.numeric(v)
    }
    vv <- as.vector(unlist(vv))
  }
  tme <- as.POSIXlt(tme + 0, tz = "UTC")
  ll <- data.frame(x = long, y = lat)
  tme2 <- as.POSIXct(tme)
  tme2 <- c((tme2 - 24 * 3600), tme2, (tme2 + 24 * 3600), (tme2 + 42 * 3600))
  tme2 <- as.POSIXlt(tme2)
  # These variables are forecasts valid 6 hours after the reference time.
  Tk <- ncepget1('air.2m', tme2, ll)
  Tkmin <- ncepget1('tmin.2m', tme2, ll)
  Tkmax <- ncepget1('tmax.2m', tme2, ll)
  sh <- ncepget1('shum.2m', tme2, ll)
  pr <- ncepget1('pres.sfc', tme2, ll)
  wu <- ncepget1('uwnd.10m', tme2, ll)
  wv <- ncepget1('vwnd.10m', tme2, ll)
  # These variables are 6 hour averages starting at the reference time.
  dlw <- ncepget1('dlwrf.sfc', tme2, ll)
  ulw <- ncepget1('ulwrf.sfc', tme2, ll)
  dsw <- ncepget1('dswrf.sfc', tme2, ll)
  tcdc <- ncepget1('tcdc.eatm', tme2, ll)
  # get correct data
  ogn <- paste0(tme2$year[1] + 1900, "-", tme2$mon[1] + 1, "-01 00:00")
  xx <- (c(1:length(Tk)) - 1) * 3600 * 6
  tma <- as.POSIXlt(xx, origin = ogn, tz = "UTC")
  sel <- which(tma >= tme2[1] & tma <= tme2[length(tme2)])
  dfout <- data.frame(obs_time = tma[sel], timezone = "UTC", Tk = Tk[sel],
                      Tkmin = Tkmin[sel], Tkmax = Tkmax[sel], sh = sh[sel],
                      pr = pr[sel], wu = wu[sel], wv = wv[sel], dlw = dlw[sel],
                      ulw = ulw[sel], dsw = dsw[sel], tcdc = tcdc[sel])
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
#' @return a dataframe with the following variables:
#' \describe{
#'   \item{obs_time}{POSIXlt object of times in UTC}
#'   \item{temperature}{emperatures at 2m (ºC)}
#'   \item{humidity}{specific humidity at 2m (Kg / Kg)}
#'   \item{pressure}{surface pressure (Pa)}
#'   \item{windspeed}{wind speed at 2m (metres per second}
#'   \item{winddir}{wind direction (degrees from N)}
#'   \item{emissivity}{emissivity of the atmosphere (0 - 1, downlong / uplong)}
#'   \item{netlong}{Net longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{uplong}{Upward longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{downlong}{Downward longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{rad_dni}{Direct radiation normal to the solar beam (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{rad_dif}{Diffuse radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{szenith}{the zenith angle (degrees)}
#'   \item{cloudcover}{cloud cover (Percentage)}
#' }
#' @export
#' @import zoo
#'
#' @seealso [get_NCEP()]
#' @details If `ncepdata` is not provided, then [get_NCEP()] is called and data are downloaded from NCEP Atmospheric
#' Model Intercomparison Project (Kanamitso et al 2002). Six-hourly data are interpolated as follows.
#' Pressure, humidity and the u and v wind vectors are converted to hourly using spline interpolation. Wind speeed and direction and then
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
hourlyNCEP <- function(ncepdata = NA, lat, long, tme, reanalysis2 = FALSE) {
  bound <- function(x, mn = 0, mx = 1) {
    x[x > mx] <- mx
    x[x < mn] <- mn
    x
  }
  tz <- format(tme, format="%Z")
  if (tz[1] != "UTC" & tz[1]!= "GMT") {
    warning(paste("NCEP data use UTC/GMT. Timezone converted from", tz[1], "to UTC/GMT"))
  }
  tme <- .tme.sort(tme)
  tmeout <- tme
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
  csr <- clearskyrad(tmo2, lat, long, merid = 0, dst = 0)
  csr[is.na(csr)] <- 0
  csr_m <- 0
  for (i in 1:(length(tme6))) {
    st <- (i - 1) * 6 + 1
    ed <- st + 6
    csr_m[i] <- mean(csr[st:ed], na.rm = T)
  }
  od <- bound(ncepdata$dsw / csr_m)
  if (is.na(od)[1]) od[1] <- mean(od, na.rm = T)
  if (is.na(od)[length(od)]) od[length(od)] <- mean(od, na.rm = T)
  h_od <- bound(spline(tme6, od, n = n)$y)
  tmorad <- as.POSIXlt(tmo + 3600 * 3)
  csrh <- clearskyrad(tmorad, lat, long, merid = 0, dst = 0)
  rad_h <- h_od * csrh
  rad_h[is.na(rad_h)] <- 0
  jd <- julday(tmorad$year + 1900, tmorad$mon + 1, tmorad$mday)
  dp <- 0
  for (i in 1:length(jd)) {
    dp[i] <- .difprop2(rad_h[i], jd[i], tmorad$hour[i], lat, long)
  }
  dp[is.na(dp)] <- 1
  si <- siflat(tmorad$hour, lat, long, jd, merid = 0)
  h_dif <- dp *  rad_h
  h_dni <- ((1 - dp) * rad_h) / si
  h_dni[is.na(h_dni)] <- 0
  h_dni[h_dni > 1353] <- 1353
  ###
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
  h_dni <- h_dni * 0.0036
  h_dif <- h_dif * 0.0036
  h_tc<-suppressWarnings(hourlytemp(julian = jd, em = em_h[thsel], dni = h_dni[thsel],
                                    dif = h_dif[thsel], mintemp = tmin, maxtemp = tmax,
                                    lat = lat, long = long, merid = 0))
  hlwu <- 2.043e-10 * (h_tc + 273.15)^4
  hlwd <- em_h[thsel] *  hlwu
  h_nlw <- hlwu - hlwd
  h_uw <- spline(tme6, ncepdata$wu, n = n)$y
  h_vw <- spline(tme6, ncepdata$wv, n = n)$y
  h_ws <- sqrt(h_uw^2 + h_vw^2)
  h_ws <- windheight(h_ws, 10, 2)
  h_wd <- atan2(h_uw, h_vw) * 180/pi + 180
  h_wd <- h_wd%%360
  jd <- julday(tmorad$year + 1900, tmorad$mon + 1, tmorad$mday)
  szenith <- 90 - solalt(tmorad$hour, lat, long, jd, merid = 0)
  hourlyout <- data.frame(obs_time = tmorad[thsel], temperature = h_tc,
                          humidity = h_sh[thsel2], pressure = h_pr[thsel2],
                          windspeed = h_ws[thsel2], winddir = h_wd[thsel2],
                          emissivity = em_h[thsel2], cloudcover = h_cc[thsel],
                          netlong = h_nlw, uplong = hlwu, downlong = hlwd,
                          rad_dni = h_dni[thsel],
                          rad_dif = h_dif[thsel],
                          szenith = szenith[thsel], timezone = "UTC")
  tmetest <- as.POSIXlt(hourlyout$obs_time)
  tmeout <- as.POSIXlt(tmeout + 0, tz = "UTC")
  sel <- which(tmetest >= min(tmeout) & tmetest <= max(tmeout + 23 * 3600))
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
#' @import RNCEP
#'
#' @examples
#' tme <- as.POSIXlt(c(1:15) * 24 * 3600, origin = "2015-01-15", tz = 'UTC')
#' dailyprecipNCEP(50, -5, tme)
dailyprecipNCEP <- function(lat, long, tme, reanalysis2 = FALSE) {
  tmeout <- tme
  tme <- .tme.sort(tme)
  long <- ifelse(long > 180, long - 360, long)
  if (length(tme) > 1) {
    int <- as.numeric(tme[2]) - as.numeric(tme[1])
    lgth <- (length(tme) * int) / (24 * 3600)
    tme <- as.POSIXlt(c(0:(lgth - 1)) * 3600 * 24, origin = min(tme), tz = 'UTC')
  } else tme <- as.POSIXlt(0, origin = min(tme), tz = 'UTC')
  yrs <- unique(tme$year + 1900)
  apre <- list()
  for (y in 1:length(yrs)) {
    sely <- which(tme$year + 1900 == yrs[y])
    mths <- unique(tme$mon[sely] + 1)
    ll <- data.frame(x = long, y = lat)
    pre <- NCEP.gather('prate.sfc', level = 'gaussian',
                       years.minmax = c(min(yrs[y]),max(yrs[y])),
                       months.minmax = c(min(mths):max(mths)),
                       lat.southnorth = c(ll$y,ll$y), lon.westeast = c(ll$x,ll$x),
                       return.units = FALSE, status.bar = FALSE, reanalysis2 = reanalysis2)
    if (is.null(dim(pre)) == F) {
      latdif <- abs(ll$y - as.numeric(rownames(pre)))
      londif <- abs(ll$x%%360 - as.numeric(colnames(pre)))
      pre <- pre[which.min(latdif), which.min(londif),]
    }
    apre[[y]] <-as.numeric(pre)
  }
  pre <- as.vector(unlist(apre))
  pre <- pre * 6 * 3600
  tma <- 0
  dm <- c(31,28,31,30,31,30,31,31,30,31,30,31) * 4 - 1
  for (yr in min(yrs):max(yrs)) {
    dm[2] <- ifelse(yr%%4 == 0, 115, 111)
    sely <- which(tme$year + 1900 == yr)
    mths <- unique(tme$mon[sely] + 1)
    for (mth in min(mths):max(mths)) {
      tmym <- as.POSIXct(c(0:dm[mth]) * 3600 * 6,
                         origin = paste0(yr,"-",mth,"-01 00:00"),
                         tz = 'UTC')
      tma <- c(tma, tmym)
    }
  }
  tma <- as.POSIXlt(tma[-1], origin = '1970-01-01', tz = 'UTC')
  tmeout <- as.POSIXlt(tmeout + 0, tz = "UTC")
  sel <- which(tma >= min(tmeout) & tma < max(tmeout))
  pre <- pre[sel]
  dpre <- t(matrix(pre, nrow = 4))
  dpre <- apply(dpre, 1, sum)
  return(dpre)
}
#' get shortwave radiation (time series)
#' @export
.shortwave.ts <- function(dni, dif, jd, localtime, lat, long, slope, aspect,
                          ha = 0, svv = 1, x = 1, l = 0, albr = 0,
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
  fd <- dirr + cisr  # direct
  fdf <- isor + refr  # diffuse
  kk <- ((x ^ 2 + 1 / (tan(sa * (pi / 180)) ^ 2)) ^ 0.5) /
    (x + 1.774 * (x + 1.182) ^ (-0.733))
  trd <- exp(-kk * l)
  fr <- as.vector(canopy(array(l, dim = c(1,1))))
  trf <- (1 - fr)
  fgd <- fd * trd
  fged <- fdf * trf * svv
  fgc <- fgd + fged
  cfc <- ((1 - trd) * fd + fr * fdf) / (fd + fdf)
  cfc[is.na(cfc)] <- ((1 - trd[is.na(cfc)]) * 0.5 + fr * 0.5) / (0.5 + 0.5)
  return(xxx=list(swrad = fgc, canopyfact = cfc))
}
#' get radiation and wind for use with NicheMapR
#' @export
.pointradwind <- function(hourlydata, dem, lat, long, l, x, albr = 0.15, zmin = 0,
                          slope = NA, aspect = NA, horizon = NA, svf = NA, difani = TRUE) {
  m <- is_raster(dem)
  m[is.na(m)] <- zmin
  m[m < zmin] <- zmin
  dem <- if_raster(m, dem)
  xy <- data.frame(x = long, y = lat)
  xy <- sf::st_as_sf(xy, coords = c("x","y"), crs = 4326)
  xy <- st_transform(xy, st_crs(dem)$wkt)
  reso <- res(dem)[1]
  wsc36 <- 0
  wsc36atground <- 0
  for (i in 0:35) {
    wscr <- windcoef(dem, i*10, hgt = 2, res = reso)
    wscr2 <- windcoef(dem, i*10, hgt = 0, res = reso)
    wsc36[i + 1] <- raster::extract(wscr, xy)
    wsc36atground[i + 1] <- raster::extract(wscr2, xy)
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
  if (class(slope)[1] == "logical") {
    slope <- terrain(dem, unit = 'degrees')
    slope <- raster::extract(slope, xy)
  }
  if (class(aspect)[1] == "logical") {
    aspect <- terrain(dem, opt = 'aspect', unit = 'degrees')
    aspect <- raster::extract(aspect, xy)
  }
  fr <- canopy(array(l, dim = dim(dem)[1:2]))
  if (class(svf)[1] == "logical") {
    svf <- skyviewveg(dem, array(l, dim = dim(dem)[1:2]),
                      array(x, dim = dim(dem)[1:2]), res = reso)
    svf <- raster::extract(svf, xy)
  }
  if (class(horizon)[1] == "logical") {
    ha36 <- 0
    for (i in 0:35) {
      har <- horizonangle(dem, i*10, reso)
      ha36[i + 1] <- atan(raster::extract(har, xy)) * (180/pi)
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
  windsp <- hourlydata$windspeed
  hourlyrad <- data.frame(swrad = sw$swrad, skyviewfact = svf, canopyfact = sw$canopyfact,
                          whselt = wsheltatground, windspeed = wshelt *  windsp,
                          slope = slope, aspect = aspect)
  return(hourlyrad)
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
                         emthresh = 0.78,  weather.elev, cad.effects) {
  xy <- data.frame(x = long, y = lat)
  if (weather.elev == 'ncep') {
    elevncep <- extract(demworld, xy)
  } else if (weather.elev == 'era5') {
    lar<-round(lat*4,0)/4
    lor<-round(long*4,0)/4
    e<-extent(lor-0.875,lor+0.875,lar-0.875,lar+0.875)
    e5d<-get_dem(r=NA,lat,long,resolution = 10000, xdims = 10, ydims = 10)
    e5d<-projectRaster(e5d,crs="+init=epsg:4326")
    rte<-raster(e)
    res(rte)<-0.25
    e5d<-resample(e5d,rte)
    xy <- data.frame(x = long, y = lat)
    elevncep <- extract(e5d, xy)
  } else elevncep<-as.numeric(weather.elev)
  if (is.na(elevncep)) {
    warnings("elevation of input weather data NA. Setting to zero")
    elevncep<-0
  }
  coordinates(xy) = ~x + y
  proj4string(xy) = "+init=epsg:4326"
  xy <- as.data.frame(spTransform(xy, crs(dem)))
  elev <- extract(dem, xy)
  if (is.na(elev)) elev <- 0
  lr <- lapserate(hourlydata$temperature, hourlydata$humidity, hourlydata$pressure)
  elevt <- lr * (elev - elevncep) + hourlydata$temperature
  tme <- as.POSIXlt(hourlydata$obs_time)
  ## cad effects
  if (cad.effects) {
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
  } else {
    basins<-dem*0+1
    fa<-basins
    cadt<-rep(0,length(tme))
  }
  tout <- data.frame(tref = hourlydata$temperature,
                     elev = elev, elevncep = elevncep,
                     telev = elevt, tcad = cadt, lapserate = lr)
  return(list(tout = tout, basins = basins, flowacc = fa))
}
#' Calculates coastal exposure automatically
#' @export
.invls.auto <- function(r, steps = 8, use.raster = T, zmin = 0, plot.progress = TRUE, tidyr = FALSE) {
  tidydems <- function(rfine, rc) {
    rfine[is.na(rfine)] <- zmin
    rc <- trim(rc)
    aggf <- floor(mean(res(rc)[1:2]) / mean(res(rfine)))
    if (aggf > 1) {
      rfine2 <- suppressWarnings(aggregate(rfine, aggf, max))
    } else rfine2 <- rfine
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
  yr <-  unique(tme$year) + 1900
  if (length(yr) > 1) stop("SST intepolation only works on data from one year. Select dates in yearly blocks")
  if(yr%%4 == 0) dms[2] <- 29
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
#' @param landsea a raster object with NAs (representing sea) or any non-NA value and a projection system defined.
#' @param ncephourly a dataframe of hourly climate data as returned by [hourlyNCEP()].
#' @param steps an optional integer. Coastal effects are calculated in specified directions upwind. Steps defines the total number of directions used. If the default 8 is specified, coastal effects are calculated at 45º intervals.
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
    xx <- p1[i] * d1^b1[i] + p2[i] * d2^b2[i]
    xx[xx > 6] <- 6
    xx[xx < -6] <- -6
    pdT <- dT[i] + xx
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
#' @param weather.elev optional value indicating the elevation of values in `hourlydata`. Either a numeric value, corresponding to
#' the elevation in (m) of the location from which `hourlydata` were obtained, or one of `ncep` (default, data derive
#' from NOAA-NCEP reanalysis) project or `era5` (derived from Copernicus ERA5 climate reanalysis).
#' @param cad.effects optional logical indicating whether to calaculate cold air drainage effects
#' (TRUE = Yes, slower. FALSE =  No, quicker)
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
#' [windcoef()]. (v) `windspeed`: wind speed (m /s) at 2 m, corrected for topographic
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
                             resolution = 100, zmin = 0, slope = NA, aspect = NA, horizon = NA,
                             svf = NA, difani = TRUE, windthresh = 4.5, emthresh = 0.78, reanalysis2 = FALSE,
                             steps = 8, use.raster = TRUE, plot.progress = TRUE, tidyr = FALSE,
                             weather.elev = 'ncep', cad.effects = TRUE) {
  tme <- seq(as.POSIXlt(dstart, format = "%d/%m/%Y", origin = "01/01/1900", tz = "UTC"),
             as.POSIXlt(dfinish, format = "%d/%m/%Y", origin = "01/01/1900", tz = "UTC")
             + 3600 * 24, by = 'hours')
  tme2 <- seq(as.POSIXlt(dstart, format = "%d/%m/%Y", origin = "01/01/1900"),
              as.POSIXlt(dfinish, format = "%d/%m/%Y", origin = "01/01/1900")
              + 3600 * 24, by = 'hours')
  tz1 <- format(as.POSIXlt(tme), format="%Z")
  tz2 <- format(as.POSIXlt(tme2), format="%Z")
  xx <- as.numeric(tme) == as.numeric(tme2)
  sel <- which(xx == FALSE)
  if (length(sel) > 1) {
    warning(paste("Data sequence in UTC/GMT. Some or all dates using system timezone are in", tz2[sel[1]]))
  }
  tme <- tme[-length(tme)]
  tme <- as.POSIXlt(tme)
  if (class(dem)[1] == "logical") {
    cat("Downloading digital elevation data \n")
    dem <- get_dem(r = NA, lat = lat, long = long, resolution = resolution, zmin = zmin)
  }
  if (class(demmeso)[1] == "logical") demmeso <- dem
  if (class(hourlydata) == "logical") {
    cat("Extracting climate data from NCEP \n")
    hourlydata <- hourlyNCEP(ncepdata = NA, lat, long, tme, reanalysis2)
  }
  if (class(dailyprecip) == "logical") {
    cat("Extracting rainfall data from NCEP \n")
    dailyprecip <- dailyprecipNCEP(lat, long, tme, reanalysis2)
  }
  cat("Downscaling radiation and wind speed \n")
  radwind <- .pointradwind(hourlydata, dem, lat, long, l, x, albr, zmin, slope, aspect,
                           horizon, svf, difani)
  cat("Calculating meso-scale terrain effects \n")
  info <- .eleveffects(hourlydata, demmeso, lat, long, windthresh, emthresh, weather.elev, cad.effects)
  elev <- info$tout
  if (coastal) {
    m <- is_raster(dem)
    m[m == zmin] <- NA
    dem <- if_raster(m, dem)
    acoast <- coastalNCEP(dem, hourlydata, steps, use.raster, zmin, plot.progress, tidyr = FALSE)
    xi <- floor(dim(dem)[1] / 2)
    yi <- floor(dim(dem)[2] / 2)
    ce <- acoast[xi, yi, ]
    if (is.na(ce[1])) ce <- apply(acoast, 3, mean, na.rm = T)
    elev$tref <- ce - elev$tref + elev$telev
  }  else acoast <- NA
  return(list(hourlydata = hourlydata, hourlyradwind = radwind, tref = elev,
              dailyprecip = dailyprecip, acoast = acoast, basins = info$basins,
              flowacc = info$flowacc))
}
#' Function for automatically generating microclimate surfaces for anywhere in the word
#'
#' @description This function generating microclimate temperature surfaces for anywhere
#' in the word. If hourly weather data are not provided, it first downloads coarse-resolution climate and radiation data from the
#' NCEP-NCAR or NCEP–DOE Atmospheric Model Intercomparison Project (Kanamitso et al
#' 2002) and interpolates these data to hourly. It calculates mesoclimatic effects and derives parameters for fitting the
#' microclimate model using the `NicheMapR` package (Kearny & Porter 2016). Using digital
#' elevation data either downloaded or provided by the user, and canopy
#' characteristics either specified by the user, or derived from habitat,
#' it runs the microclimate model in hourly timesteps to generate an array of temperatures.
#'
#' @param r a raster object defining the extent and resolution for which microclimate
#' temperature data are required. Supplied raster must have a projection such that the units of
#' x, y and z are identical, and grid cells are square. NA values are assumed to be sea, which
#' is important in the calculation of coastal and cold air drainage effects.
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
#' @param hourlydata an optional dataframe of hourly climate forcing variables. If not supplied
#' downloaded using [hourlyNCEP()]
#' \describe{
#'   \item{obs_time}{POSIXlt object of times in UTC}
#'   \item{temperature}{Temperatures at 2 m (ºC)}
#'   \item{humidity}{Specific humidity at 2m (Kg / Kg)}
#'   \item{pressure}{Surface pressure (Pa)}
#'   \item{windspeed}{Wind speed at 2 m (metres per second)}
#'   \item{winddir}{Wind direction (degrees from N)}
#'   \item{emissivity}{Emissivity of the atmosphere (0 - 1, downlong / uplong)}
#'   \item{netlong}{Net longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{uplong}{Upward longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{downlong}{Downward longwave radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{rad_dni}{Direct radiation normal to the solar beam (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{rad_dif}{Diffuse radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})}
#'   \item{szenith}{The zenith angle (degrees)}
#'   \item{cloudcover}{Cloud cover (Percentage)}
#' }
#' @param dailyprecip a vector of daily rainfall (mm / day) for the period specified.
#' If not supplied downloaded using [dailyprecipNCEP()].
#' @param coastal optional logical value indicating whether or not to calculate
#' coastal effects.
#' @param use.raster optional logical value indicating whether to use `r` in the
#' calculation of coastal effects. If used, NAs or values corresponding to `zmin`
#' must represent sea.
#' @param r.is.dem optional logical value indicating whether 'r' is a digital
#' elevation dataset used for calculating microclimatic effects. If FALSE, then a
#' dem is downloaded.
#' @param save optional integer: don't save forcing data (0), save the forcing data (1)
#' or read previously saved data (2).
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
#' @param continuous an optional logical value indicating whether to treat wind speed as a continuous variable
#' @param summarydata an optional logical indicating whether to calculate summary data
#' (frost hours and maximum, minimum and mean temperature) for each pixel and return these to the output.
#' @param save.memory An optional logical indicatign whether to save
#' @param weather.elev optional value indicating the elevation of values in `hourlydata`. Either a numeric value, corresponding to
#' the elevation in (m) of the location from which `hourlydata` were obtained, or one of `ncep` (default, data derive
#' from NOAA-NCEP reanalysis) project or `era5` (derived from Copernicus ERA5 climate reanalysis).
#' @param cad.effects optional logical indicating whether to calaculate cold air drainage effects
#' (TRUE = Yes, slower. FALSE =  No, quicker)
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
#' though note that edge cells are NA as slopes cannot be calculated for these cells.
#' (3) units: the units of `temps`. Either deg C or dec C * 1000 if `save.memory` is TRUE.
#' (4) tmax: If `summarydata` is TRUE, a matrix of maximum temperatures
#' (5) tmin: If `summarydata` is TRUE, a matrix of minimum temperatures
#' (6) tmean: If `summarydata` is TRUE, a matrix of mean temperatures
#' (7) frosthours: If `summarydata` is TRUE, a matrix of hours below 0 deg C.
#'
#' @import raster
#' @export
#'
#' @examples
#' library(raster)
#' require(NicheMapR)
#' # Get DEM for Pico, Azores
#' r <- get_dem(lat = 38.467429, long = -28.398995, resolution = 30)
#' plot(r)
#' # Takes ~ c. 5 minutes to run
#' temps <- runauto(r, "10/06/2010", "15/06/2010", hgt = 0.1, l = NA, x = NA,
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

runauto <- function(r, dstart, dfinish, hgt = 0.05, l, x, habitat = NA,
                    hourlydata = NA, dailyprecip = NA, use.raster = FALSE,
                    coastal = TRUE, r.is.dem = TRUE, save = 0, albg = 0.15, albr =0.15,
                    albc = 0.23, mesoresolution = 100, zmin = 0, slope = NA,
                    aspect = NA, windthresh = 4.5, emthresh = 0.78, reanalysis2 = FALSE,
                    steps = 8, plot.progress = TRUE, continuous = TRUE,
                    summarydata = TRUE, save.memory = FALSE, weather.elev = 'ncep',
                    cad.effects = TRUE) {
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
    if (class(habitat)[1] == "matrix") {
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
    if (class(l)[1] == "matrix") la <- l
    if (class(l)[1] == "array") la <- l[,,i]
    la
  }
  checkfun <- function(l, x, habitat, r, tme) {
    if (class(l)[1] == "logical" & class(habitat)[1] == "logical") {
      stop ("No habitat or leaf area provided")
    }
    if (class(x)[1] == "logical" & class(habitat)[1] == "logical") {
      stop ("No habitat or leaf angle distribution provided")
    }
    if (class(l)[1] == "matrix" | class(l)[1] == "array") {
      if (dim(l)[1] != dim(r)[1] & dim(l)[2] != dim(r)[2]) {
        stop ("l must have same xy dimensions as r")
      }
    }
    if (class(l)[1] == "array") {
      if (dim(l)[3] != length(tme)) {
        stop ("Length of l must equal 1 or number of hours in time sequence")
      }
    }
    if (class(l)[1] == "numeric" | class(l)[1] == "integer") {
      if (length(l) != 1 & length(l)  != length(tme)) {
        stop ("Length of l must equal 1 or number of hours in time sequence")
      }
    }
    if (class(habitat)[1] == "matrix") {
      if (dim(habitat)[1] != dim(r)[1] & dim(habitat)[2] != dim(r)[2]) {
        stop("habitat must have same dimensions as r")
      }
    }
    if (class(x)[1] == "matrix") {
      if (dim(x)[1] != dim(r)[1] & dim(x)[2] != dim(r)[2]) {
        stop("x must have same dimensions as r")
      }
    }
    if (class(x)[1] == "numeric" | class(x) == "integer") {
      if (length(x) > 1) stop("varying x through time not supported")
    }
    if (class(x)[1] == "array") stop("varying x through time not supported")
  }
  frost <- function(x) {
    sel <- which(x <= 0)
    y <- length(sel)
    if (is.na(mean(x, na.rm = TRUE))) y <- NA
    y
  }
  if (!requireNamespace("NicheMapR", quietly = TRUE)) {
    #  stop("package 'NicheMapR' is needed. Please install it from Github: 'mrke/NicheMapR'",
    #       call. = FALSE)
  }
  # Lat long and time
  lat <- latlongfromraster(r)$lat
  long <- latlongfromraster(r)$long
  tme <- seq(as.POSIXlt(dstart, format = "%d/%m/%Y", origin = "01/01/1900", tz = "UTC"),
             as.POSIXlt(dfinish, format = "%d/%m/%Y", origin = "01/01/1900", tz = "UTC")
             + 3600 * 24, by = 'hours')
  tme2 <- seq(as.POSIXlt(dstart, format = "%d/%m/%Y", origin = "01/01/1900"),
              as.POSIXlt(dfinish, format = "%d/%m/%Y", origin = "01/01/1900")
              + 3600 * 24, by = 'hours')
  tz1 <- format(as.POSIXlt(tme), format="%Z")
  tz2 <- format(as.POSIXlt(tme2), format="%Z")
  xx <- as.numeric(tme) == as.numeric(tme2)
  sel <- which(xx == FALSE)
  if (length(sel) > 1) {
    cat(paste("Note: Model run in UTC/GMT. Some or all dates using system timezone are in", tz2[sel[1]]), "\n")
  }
  tme <- tme[-length(tme)]
  tme <- as.POSIXlt(tme)
  # l and x
  l <- is_raster(l)
  x <- is_raster(x)
  habitat <- is_raster(habitat)
  checkfun(l, x, habitat, r, tme)
  if (class(x)[1] == "numeric" | class(x)[1] == "integer") x <- array(x, dim = dim(r)[1:2])
  if (class(x)[1] == "logical") {
    if (class(habitat)[1]  == "numeric" | class(habitat)[1] == "integer" |
        class(habitat)[1] == "character") {
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
  micronmr <- micro_ncep(dstart = dstart, dfinish = dfinish, dem = r, dem2 = dem, LAI = 0,
                         loc = loc, Usrhyt = hgt2, Refhyt = 2, coastal = coastal, reanalysis = reanalysis2,
                         DEP = dep, save = save, hourlydata = hourlydata, dailyprecip = dailyprecip,
                         weather.elev = weather.elev, cad.effects = cad.effects)
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
    params <- fitmicro(mft, alldata = TRUE, continuous = continuous)
  }
  if (hgt <= 0) {
    soiltemps <- as.data.frame(micronmr$soil)
    mft <- data.frame(temperature = soiltemps$D0cm,
                      reftemp = nmrout$TAREF,
                      wind = rwind$windspeed,
                      netrad = netrad)
    params <- fitmicro(mft, alldata = TRUE, continuous = continuous)
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
    fr <- canopy(la, 0.23)
    alb <- fr * albc + (1 - fr) * albg
    radsw <- shortwaveveg(hourlydata$rad_dni[i], hourlydata$rad_dif[i], jd[i],
                          tme$hour[i], dtm = r2, svv = svv, alb = alb,
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
    tea <- runmicro(params, nr, ws, continuous = continuous)
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

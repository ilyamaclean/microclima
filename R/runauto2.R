#' Zero-plane displacement height
.zeroplanedis <- function(hgt, l) {
  r <- hgt
  hgt <- is_raster(hgt)
  l <- is_raster(l)
  l[l < 0.1] <- 0.1
  l[l > 20] <- 20
  m <- 0.0609 * log(l) + 0.5894
  d <- m * hgt
  d <- if_raster(d, r)
  d
}
#' Roughness length
.roughlength <- function(hgt, l) {
  r <- hgt
  hgt <- is_raster(hgt)
  l <- is_raster(l)
  m <- 0.1 * l + 0.08
  sel <- which(l > 0.6)
  m[sel] <- -0.0239 * log(l[sel]) + 0.1275
  zm <- m * hgt
  zm <- ifelse(zm < 0.004, 0.004, zm)
  zm <- if_raster(zm, r)
  zm
}
#' Crop reference evaporanspiration
.evapref <- function(Rnet, tc, u2, relhum, pk) {
  Rnet <- is_raster(Rnet)
  Delta<-4098*(0.6108*exp(17.27*tc/(tc+237.3)))/(tc+237.3)^2
  es <- 0.6108*exp(17.27*tc/(tc+237.3))
  ea <- (relhum/100) * es
  G <- ifelse(Rnet > 0, 0.1 * Rnet, 0.5 * Rnet)
  Gamma <-  (0.001013 * pk) / (0.622 * 2.45)
  lhs <- 0.408 * Delta * (Rnet - G)
  rhs <- Gamma * (37 / (tc + 273)) * u2 * (es - ea)
  btm <- Delta + Gamma * (1 + 0.34 * u2)
  cre <- (lhs + rhs) / btm
  cre[cre < 0] <- 0
  cre
}
#' Diabatic correction factor
.diabatic_cor <- function(tc, pk = 101.3, H = 0, uf, zi = 2, d) {
  Tk <- tc + 273.15
  ph <- 42.35
  cp <-  29.3
  st <- -(0.4 * 9.81 * (zi - d) * H) / (ph * cp * Tk * uf^3)
  # Stable flow
  sel <- which(st < 0) # unstable
  # Stable
  psi_h <- suppressWarnings(6 * log(1 + st))
  psi_m <- psi_h
  # Unstable
  psi_h[sel] <-   -2 * log((1 + (1 - 16 * st[sel])^0.5) / 2)
  psi_m[sel] <- 0.6 * psi_h[sel]
  psi_m <-ifelse(psi_m > 4, 4, psi_m)
  psi_h <-ifelse(psi_h > 4, 4, psi_h)
  return(list(psi_m = psi_m, psi_h = psi_h))
}
#' Mixing length
.mixinglength <- function(veghgt, l, x) {
  Ld <- l / veghgt
  lmg <- (4 * 0.07 / (pi * Ld))^0.5
  lms <- (6 * 0.07^2 * veghgt / (pi * l))^(1 / 3)
  wgt<-x/2
  sel <- which(x > 1)
  wgt[sel] <- 1-(0.5/x[sel])
  l_m <- wgt * lms + (1 - wgt) * lmg
  l_m
}
#' Wind attenuation coefficient
.attencoef <- function(veghgt, l, x) {
  r<-veghgt
  veghgt<-is_raster(veghgt)
  l<-is_raster(l)
  x<-is_raster(x)
  l_m <- .mixinglength(veghgt, l, x)
  a <- ((0.2 * l * veghgt) / (l_m))^0.5
  a <- if_raster(a, r)
  a
}
#' Downscaled shortwave radiation
.shortwaveveg2 <- function(dni, dif, julian, localtime, lat = NA, long = NA,
                           dtm = array(0, dim = c(1, 1)), slope = NA, aspect = NA,
                           svv = 1, albg = 0.15, albc = 0.23, fr, albr = 0.23, ha = 0,
                           res = 1, merid = NA, dst = 0, shadow = TRUE,
                           x, l, difani = TRUE) {
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
  if (class(merid) == "logical") merid <- round(long / 15, 0) * 15
  slope <- is_raster(slope)
  aspect <- is_raster(aspect)
  dtm <- is_raster(dtm)
  dtm[is.na(dtm)] <- 0
  si <- solarindex(slope, aspect, localtime, lat, long, julian, dtm, res,
                   merid, dst, shadow)
  dirr <- si * dni
  a <- slope * (pi / 180)
  if (difani) {
    k <- dni / 4.87
  } else k <- 0
  k <- ifelse(k > 1, 1, k)
  isor <- 0.5 * dif * (1 + cos(a)) * (1 - k)
  cisr <- k * dif * si
  sdi <- (slope + ha) * (pi / 180)
  refr <- 0.5 * albr * (1 - cos(sdi)) * dif
  fd <- dirr + cisr
  fdf <- isor + refr
  saltitude <- solalt(localtime, lat, long, julian, merid, dst)
  s <- 1 - albc
  zen <- 90 - saltitude
  kk <- sqrt((x^2 + (tan(zen * (pi/180))^2)))/(x + 1.774 * (x + 1.182)^(-0.733))
  trd <- exp(-kk * s * l)
  trf <- exp(-s * l)
  fgd <- fd * trd * (1 - albg)
  fged <- fdf * trf * (1 - albg) * svv
  fgc <- fgd + fged
  fgd2 <- fd * trd
  fged2 <- fdf * trf * svv
  fgc2 <- fgd2 + fged2
  fgc <- if_raster(fgc, r)
  fgc2 <- if_raster(fgc2, r)
  return(list(r1 = fgc, r2 = fgc2))
}
# Converts rasters or matrices to arrays
.rastertoarray <- function(r, n) {
  r <-is_raster(r)
  if (length(dim(r)) == 2) {
    a <- array(r, dim = c(dim(r), n))
  }
  if (length(dim(r)) == 3 & dim(r)[3] != n) {
    a <- array(r, dim = c(dim(r)[1:2], n))
    for (i in 1:dim(r)[1]) {
      for (j in 1:dim(r)[2]) {
        tst <- mean(r[i,j,], na.rm = T)
        if (is.na(tst) == F) {
          y <- r[i,j,]
          x <- c(1:length(y))
          a[i,j,] <- spline(y~x,n=n)$y
          a[i,j,] <- ifelse(a[i,j,] < 0, 0, a[i,j,])
        }
      }
    }
  }
  a
}
#' Ground and canopy albedo etc.
.albedosort <- function(alb, l) {
  alb <- is_raster(alb)
  l <- is_raster(l)
  if (length(dim(alb)) == 2) {
    dm <-1
  } else dm <- dim(alb)[3]
  if (length(dim(l)) == 2) {
    dm[2] <-1
  } else dm[2] <- dim(l)[3]
  if (dm[1] > dm[2]) l <- .rastertoarray(l, dim(alb)[3])
  if (dm[2] > dm[1]) alb <- .rastertoarray(alb, dim(l)[3])
  if (length(dim(l)) == 2) {
    fr <- canopy(l)
    albg <- albedo2(alb, fr)
    albc <- albedo2(alb, fr, ground = FALSE)
  } else {
    fr <- alb; albg <- alb; albc <- alb
    for (i in 1:dim(l)[3]) {
      fr[,,i] <- canopy(l[,,i])
      albg[,,i] <- albedo2(alb[,,i], fr[,,i])
      albc[,,i] <- albedo2(alb[,,i], fr[,,i], ground = FALSE)
    }
  }
  return(list(fr = fr, albg = albg, albc = albc, l = l, alb = alb))
}
# Calculates roughness lengths etc as arrays
.roughnessort <- function(veghgt, l) {
  veghgt <- is_raster(veghgt)
  l <- is_raster(l)
  if (length(dim(veghgt)) == 2) {
    dm <-1
  } else dm <- dim(veghgt)[3]
  if (length(dim(l)) == 2) {
    dm[2] <-1
  } else dm[2] <- dim(l)[3]
  if (dm[1] > dm[2]) stop("dims of l must be greater than or equal to veghgt")
  if (dm[2] > dm[1]) veghgt <- .rastertoarray(veghgt, dim(l)[3])
  if (length(dim(l)) == 2) {
    zm <- .roughlength(veghgt, l)
    d <- .zeroplanedis(veghgt, l)
  } else {
    zm <- alb; d <- d
    for (i in 1:dim(l)[3]) {
      zm[,,i] <- .roughlength(veghgt[,,i], l[,,i])
      d <- .zeroplanedis(veghgt[,,i], l[,,i])
    }
  }
  return(list(zm = zm, d = d, veghgt = veghgt))
}
#' Function for automatically generating microclimate surfaces for areas with variable habitat
#'
#' @description This function generating microclimate temperature surfaces from
#' first principles. It can be used in place of [runauto()] when the region for which
#' microclimate surfaces are wanted contains variable habitat. No mesoclimate adjustments
#' are performed using this function. In contrast to [runauto()] bespoke parameterisations are performed for each pixel though the methods for computing soil and latent
#' heat fluxes are more simplistic. If hourly weather data are not provided, it downloads
#' coarse-resolution climate and radiation data from the NCEP-NCAR or NCEP–DOE Atmospheric Model Intercomparison Project (Kanamitso et al
#' 2002) and interpolates these data to hourly. It runs the microclimate model in hourly
#' time increments to generate an array of temperatures, relative humidities and wind speeds.
#' Incoming solar and downward longwave radiation are also returned.
#'
#' @param dem a raster object of elevations. Also defines extent and resolution of returned
#' microclimate data. Supplied raster must have a projection such that the units of
#' x, y and z are identical, and grid cells are square. NA values are assumed to be sea,
#' with elevation zero
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
#' @param dstart start date as character of time period required in format DD/MM/YYYY. Ignored
#' if hourly weather data supplied
#' @param dfinish end date as character of time period required in format DD/MM/YYYY. Ignored
#' if hourly weather data supplied
#' @param hgt the height (in m) above or below ground for which temperature estimates
#' are required. Below-ground not supported. Must be less than two metres.
#' @param l a raster, matrix or array of leaf area index values, as for example returned by [lai()] (see details)
#' @param x a raster, matrix or array of representing the ratio of vertical to horizontal projections of foliage as,
#' for example, returned by [leaf_geometry()]. (see details)
#' @param veghgt a raster, matrix or array of vegetation heights in metres (see details)
#' @param alb a raster, matrix or array of surface albedo in range 0 to 1 (see details)
#' @param wind.agg number of cells over which to average wind in calculation of vertical wind
#' profiles must be lower than dimensions of `dem`
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param plot.progress optional logial indicating whether to plot temperatures every 100 hours.
#'
#' @return a list with the following objects:
#' (1) `temperature`: an array of temperatures for each pixel of `dem` and hour of the time sequence (degrees C.
#' (2) `relhum`: an array of temperatures for each pixel of `dem` and hour of the time sequence (Percentage).
#' (3) `windspeed` an array of temperatures for each pixel of `dem` and hour of the time sequence (m / s).
#' (4) `swdown` an array of incoming solar radiation for each pixel of `dem` and hour of the time sequence (MJ / m^2 / hr).
#' (5) `lwdown` an array of downward longwave radiation for each pixel of `dem` and hour of the time sequence (MJ / m^2 / hr).
#' @import raster
#' @export
#'
#' @details This function uses K-theory to determine wind and temperature profiles above ground
#' as a function of sensible heat flux, a roughness length and zero-plane displacement height. Sensible
#' heat is determined from net radiation and simple approximations of latent heat and ground heat
#' fluxes. Roughness lengths and zero-plane displacement heights are determined from vegetation height and
#' leaf area. If `hgt` is below `veghgt`, the above vegetation temperature profile is extrapolated to a
#' point 0.8 x the height of the vegetation as here extrapolated wind speeds, which ultimately govern
#' the temperature profile, are approximately equivelent to the average wind speed within the canopy.
#' Incoming radiation and outgoing radiation are  adjusted by the proportion of `l` above `hgt` assuming
#' uniform vertical density in vegtation. The returned wind profile values for below canopy are calculated
#' using the method describedin Campbell & Norman (2012) Environmental Biophysics, Springer. Since
#' wind profiles are determined by fetch as well as surface roughness, the vertical wind profiles are spatially
#' averaged by an amount specified by `wind.agg`. if `l`, `x`, `veghgt` or `alb` are arrays and the third dimension of the array does not correspond
#' to the number of hours in the time sequence, hourly values are derived by interpolation.
#'
#' @examples
#' library(raster)
#' # =========================================================================== #
#' # ~~~~~~~~~~~~~~~~~ Create spatial datasets needed for input ~~~~~~~~~~~~~~~~ #
#' # =========================================================================== #
#' e<-extent(c(169450,169550,12450,12550))
#' dem <- crop(dtm1m,e) # elevation
#' alb <- albedo(aerial_image[,,1], aerial_image[,,2], aerial_image[,,3],
#'               aerial_image[,,4]) # albedo
#' plot(if_raster(alb, dtm1m), col = gray(0:255/255))
#' # Adjust albedo using modis input
#' img <- if_raster(alb, dtm1m)
#' mds <- raster(modis, xmn = 169000, xmx = 170000, ymn = 12000, ymx = 13000)
#' alb <- albedo_adjust(img, mds)
#' alb <- crop(alb, e)
#' plot(alb, col = gray(0:255/255))
#' # Leaf area index
#' leaf <- lai(aerial_image[,,3], aerial_image[,,4])
#' leaf<- raster(leaf,template=dtm1m)
#' l<-crop(leaf,dem)
#' # Leaf angle coefficient
#' x <- leaf_geometry(veg_hgt)
#' x <-crop(x,e)
#' veghgt <- crop(veg_hgt, e) # vegetation height
#' # =========================================================================== #
#' # ~~~~~~~~~~~ Run model in hourly time-steps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#' # ~~~~~~~~~~~ Weather data automatically downloaded as hourlydata = NA ~~~~~~ #
#' # =========================================================================== #
#' microout <- runauto2(dem, hourlydata = NA, dstart = "01/05/2019",
#'                      dfinish = "31/05/2019", hgt = 0.05, l = l, x = x,
#'                      veghgt = veghgt, alb = alb)
#' # =========================================================================== #
#' # ~~~~~~~~~~~ Extract temperature data and plot results ~~~~~~~~~~~~~~~~~~~~~ #
#' # =========================================================================== #
#' temps <- microout$temperature
#' tmin <- apply(temps,c(1,2),min) # spatial min
#' tmax <- apply(temps,c(1,2),max) # spatial max
#' tmin2 <- apply(temps,3,min, na.rm = T) # time series min
#' tmax2 <- apply(temps,3,max, na.rm = T) # time series max
#' # Create palette:
#' mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange", "red"))(255)
#' # Raster plots
#' par(mfrow=(c(1,2)))
#' plot(if_raster(tmin,dem), col = mypal)
#' plot(if_raster(tmax,dem), col = mypal)
#' # Time trend plots
#' ymn<-min(tmin2); ymx <- max(tmax2) # limits for plot
#' tme <- as.POSIXct(c(0:743)*3600,origin="2019-05-01 00:00", tz = "GMT")
#' par(mfrow = c(1,1))
#' plot(tmin2 ~ tme, type = "l", col = "blue", ylim = c(ymn, ymx),
#'      xlab = "", ylab = "")
#' par(new=T)
#' plot(tmax2 ~ tme, type = "l", col = "red", ylim = c(ymn, ymx),
#'      xlab = "Day", ylab = "Temperature")
runauto2 <- function(dem, hourlydata = NA, dstart, dfinish, lat = NA, long = NA, hgt = 0.05,
                     l, x, veghgt, alb, wind.agg = 10, merid = 0, dst = 0, plot.progress = TRUE) {
  if (class(hourlydata) == "logical") {
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
      warnings(paste("Note: Model run in UTC/GMT. Some or all dates using system timezone are in", tz2[sel[1]]), "\n")
    }
    tme <- tme[-length(tme)]
    tme <- as.POSIXlt(tme)
    hourlydata <- hourlyNCEP(ncepdata = NA, 50, -5, tme)
  }
  mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange", "red"))(255)
  if (plot.progress) cat("Processing input data\n")
  # get values
  r <- dem; r2 <- r; r2[is.na(r2)] <- 0
  dem<-is_raster(dem)
  # sort l and xfor sky view
  l <- is_raster(l)
  x <- is_raster(x)
  if (length(dim(l)) == 2) {
    l2 <- l
  } else l2 <- apply(l,c(1,2), mean, na.rm = T)
  if (length(dim(x)) == 2) {
    x2 <- x
  } else x2 <- apply(x,c(1,2), mean, na.rm = T)
  # sort albedos
  albx <- .albedosort(alb, l)
  fr <- albx$fr; albc <- albx$albc; albg <- albx$albg
  l <- albx$l; alb <- albx$alb
  # sort xzero plane displacement and roughness length
  zmd <- .roughnessort(veghgt, l)
  zm <- zmd$zm; d <- zmd$d; veghgt <- zmd$veghgt
  # Calculate smoothed wind
  dm <- min(dim(dem))
  if  (wind.agg > dm) {
    warning("Wind aggregation factor is larger than extent\n")
    wind.agg <- dm
  }
  if (length(dim(l)) == 2) {
    a <- .attencoef(veghgt, l, x)
    hgt3<-ifelse(hgt>veghgt,veghgt,hgt)
    mm<-(hgt3/veghgt)-1
    mm[is.na(mm)]<-1
    wc2<-exp(a*mm)
    ln1s <- suppressWarnings(log((2-d)/zm))
    ln1s[is.na(ln1s)] <- 1
    ln1s <- if_raster(ln1s, r)
    ln1s <- aggregate(ln1s, wind.agg)
    ln1s <- resample(ln1s, r)
    ln1s[is.na(r)] <- NA
    ln1s <- is_raster(ln1s)
  } else {
    ln1s <- l; wc2 <- l; a <- l
    for (i in 1:dim(l)[3]) {
      a <- .attencoef(veghgt[,,i], l[,,i], x[,,i])
      hgt3<-ifelse(hgt>veghgt[,,i],veghgt[,,i],hgt)
      mm<-(hgt3/veghgt[,,i])-1
      mm[is.na(mm)]<-1
      wc2[,,i]<-exp(a[,,i]*mm)
      ln1 <- suppressWarnings(log((2-d[,,i])/zm[,,i]))
      ln1[is.na(ln1s)] <- 1
      ln1 <- if_raster(ln1, r)
      ln1 <- aggregate(ln1, wind.agg)
      ln1 <- resample(ln1, r)
      ln1[is.na(r)] <- NA
      ln1s[,,i] <- is_raster(ln1)
    }
  }
  # Convert to arrays
  l <- .rastertoarray(l, length(tme))
  x <- .rastertoarray(x, length(tme))
  veghgt <- .rastertoarray(veghgt, length(tme))
  alb <- .rastertoarray(alb, length(tme))
  fr <- .rastertoarray(fr, length(tme))
  albc <- .rastertoarray(albc, length(tme))
  albg <- .rastertoarray(albg, length(tme))
  zm <- .rastertoarray(zm, length(tme))
  d <-  .rastertoarray(d, length(tme))
  ln1s <- .rastertoarray(ln1s, length(tme))
  wc2 <- .rastertoarray(wc2, length(tme))
  a <- .rastertoarray(a, length(tme))
  # Skyviews
  if (class(lat) == "logical") lat <- latlongfromraster(r)$lat
  if (class(long) == "logical") long <- latlongfromraster(r)$long
  tme<-as.POSIXlt(hourlydata$obs_time)
  reso<- res(r)[1]
  svv <- skyviewveg(dem, l2, x2, steps = 36, res = reso)
  sv <- skyviewtopo(dem, res = reso)
  ha <- mean_slope(dem)
  # Do height adjustments
  hgtm <- (veghgt - hgt) / veghgt
  hgtm[hgtm < 0] <- 0
  l2 <- hgtm * l
  l3 <- l; l3[l3>8] <- 8; l3[l3<0.1] <- 0.1
  # Sort out where height is needed for
  hgt2 <- ifelse(hgt < veghgt, 0.8 * veghgt, hgt)
  # Create wind shelter coefficient array
  wc <- array(NA, dim = c(dim(r)[1:2],8))
  for (i in 0:7) {
    dct <- i * 45
    w <- windcoef(r2, dct, hgt = 2, res = reso)
    w2 <- aggregate(w, wind.agg)
    w <- resample(w2, w)
    wc[,,i+1] <- is_raster(w)
  }
  wc[wc>1] <- 1; wc[wc<0] <- 0
  # juldian day and local time
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  lt <- tme$hour
  # Create arrays
  swin <- array(NA,dim = c(dim(dem)[1:2],length(tme)))
  lwin <- swin
  tempout <- swin
  rhout <- swin
  wsout <- swin
  tr <- exp(-0.9848858 * l2)
  albr <- apply(alb,3,mean)
  # Run model in hourly timesteps
  if (plot.progress) cat("Running model in hourly time steps \n")
  for (hr in 1:length(tme)) {
    # Calculate radiation
    rad <- hourlydata$rad_dni + hourlydata$rad_dif
    if (rad[hr] > 0) {
      swout <- .shortwaveveg2(hourlydata$rad_dni[hr], hourlydata$rad_dif[hr], jd[hr],
                              lt[hr], dtm = r2, svv = svv, albg = albg[,,hr], albc = albc[,,hr], fr = fr[,,hr],
                              albr = albr[hr], ha = ha, res = reso, merid = merid,
                              dst = dst, x = x[,,hr], l = l2[,,hr])
      swout1 <- is_raster(swout$r1)
      swout2 <- is_raster(swout$r2)
    } else {
      swout1 <- 0
      swout2 <- 0
    }
    up <- hourlydata$uplong[hr]*tr[,,hr]*svv
    dw <- 0.97*hourlydata$downlong[hr]*tr[,,hr]*svv
    lwout <- is_raster(up-dw)
    swin[,,hr] <- swout2
    lwin[,,hr] <-dw/0.97
    # Net radiation (topo) for latent heat
    sw <- shortwavetopo(hourlydata$rad_dni[hr], hourlydata$rad_dif[hr], jd[hr],
                        lt[hr], dtm = r2, svf = sv, alb = alb[,,hr], albr = albr[hr], ha = ha,
                        res = reso, merid = merid, dst = dst, shadow = TRUE)
    lw <- sv*(hourlydata$uplong[hr]-hourlydata$downlong[hr])
    netrad <- is_raster(sw) - lw
    # Latent heat flux
    rh <- suppressWarnings(humidityconvert(hourlydata$humidity[hr], intype = "specific",
                                           hourlydata$temperature[hr], hourlydata$pressure[hr])$relative)
    rh[rh>100]<-100
    cre <- .evapref(netrad, hourlydata$temperature[hr], hourlydata$windspeed[hr], rh,
                    hourlydata$pressure[hr]/1000)
    kc <- 0.2693*log(l3[,,hr]) + 0.9543
    L <- cre * kc * 2.45
    # Ground heat flux
    G <- ifelse(netrad>0,0.0383*exp(0.9839*netrad),0.5*netrad)
    H <- swout1 - lwout - L - G
    H <- H / 0.0036
    H2 <- netrad - L - G
    H2[is.na(H2)] <- mean(H2,na.rm=T)
    H2 <- if_raster(H2, r)
    H2 <- aggregate(H2, wind.agg)
    H2 <- resample(H2, r)
    H2 <- is_raster(H2 / 0.0036)
    # Calculate friction velocity
    dct<-round(hourlydata$winddir[hr]/45,0) + 1
    dct <- ifelse(dct > 8, dct -8, dct)
    windspeed <- hourlydata$windspeed[hr] * wc[,,dct]
    uf <- (0.4 * windspeed) / ln1s[,,hr]
    uf[uf<0.1]<-0.1
    db <- .diabatic_cor(hourlydata$temperature[hr],hourlydata$pressure[hr]/1000,
                        H2,uf,zi=2,d[,,hr])
    ln1 <- ln1s[,,hr] + db$psi_m
    ln1[ln1 < 0] <- 0
    uf <- (0.4 * windspeed) / ln1
    uf[uf<0.1]<-0.1
    # Temperature
    ln2 <- log((2 - d[,,hr])/(zm[,,hr] * 0.2)) + db$psi_h
    ln2[ln2<0] <- 0
    mf <- ln2/(496.4*uf)
    mf[mf < 0] <- 0
    mf[mf > 0.09] <- 0.09
    dT0 <- H*mf
    dT0[dT0 < 4.5] <- -4.5
    dT0[dT0 > 20] <- 20
    gr <- -dT0/ln2
    dcm <- log((hgt2[,,hr]-d[,,hr])/(0.2*zm[,,hr])) / log((2 - d[,,hr])/(0.2*zm[,,hr]))
    ln2b <- log((hgt2[,,hr]-d[,,hr])/(0.2*zm[,,hr])) + dcm * db$psi_m
    ln2b[ln2b<0]<-0
    mr <- gr*ln2b
    mr<-ifelse(mr > abs(dT0), abs(dT0), mr)
    mr<-ifelse(mr < -abs(dT0), -abs(dT0), mr)
    tempout[,,hr] <- (mr+dT0)+hourlydata$temperature[hr]
    # Relative humidity
    ea <- suppressWarnings(humidityconvert(hourlydata$humidity[hr],intype="specific",
                                           hourlydata$temperature[hr],hourlydata$pressure[hr])$vapour_pressure)
    es<-0.6108*exp(17.27*tempout[,,hr]/(tempout[,,hr]+237.3))
    rh<-(ea/es)*100
    rh[rh>100]<-100
    # Wind speed
    uf <- (0.4 * windspeed) / ln1s[,,hr]
    xx<-(veghgt[,,hr]-d[,,hr])/zm[,,hr]
    xx[xx<1]<-1; xx[is.na(xx)]<-1
    uh <- (uf/0.4)*log(xx)
    wsout[,,hr]<-uh*wc2[,,hr]
    # Save arrays
    rhout[,,hr] <- rh
    if (plot.progress & hr%%100 == 0) {
      ro <- if_raster(tempout[,,hr],r)
      plot(ro,main=tme[hr], col = mypal)
    }
  }
  return(list(temperature = tempout, relhum = rhout, windspeed = wsout,
              swdown = swin, lwdown = lwin))

}

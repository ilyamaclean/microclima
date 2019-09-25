#' Calculates surface albedo
#'
#' @description `albedo` is used to calculate surface albedo.
#'
#' @param blue a raster object, two-dimensional array or matrix of reflectance values in the blue spectral band (0 to `max.val`).
#' @param green a raster object, two-dimensional array or matrix of reflectance values in the green spectral band (0 to `max.val`).
#' @param red a raster object, two-dimensional array or matrix of reflectance values in the red spectral band (0 to `max.val`).
#' @param nir a raster object, two-dimensional array or matrix of reflectance values in the near-infrared spectral band (0 to `max.val`).
#' @param maxval a single numeric value representing the maximum reflectance in any of the spectral bands.
#' @param bluerange an optional numeric vector of length 2 giving the range (minimum, maximum) of wavelength values captured by the blue spectral band sensor (nm).
#' @param greenrange an optional numeric vector of length 2 giving the range (minimum, maximum) of wavelength values captured by the green spectral band sensor (nm).
#' @param redrange an optional numeric vector of length 2 giving the range (minimum, maximum) of wavelength values captured by the red spectral band sensor (nm).
#' @param nirrange an optional numeric vector of length 2 giving the range (minimum, maximum) of wavelength values captured by the near-infrared spectral band sensor (nm).
#'
#' @details
#' The function assumes that image reflectance has been captured using four spectral bands.
#' If `blue` is a raster object, then a raster object is returned.
#'
#' @seealso Function [albedo_adjust()] for adjusted albedo for image brightness and contrast.
#'
#' @return a raster object or two-dimensional array of numeric values representing surface albedo (range 0 to 1).
#' @import raster
#' @export
#'
#' @examples
#' library(raster)
#' alb <- albedo(aerial_image[,,1], aerial_image[,,2], aerial_image[,,3],
#'               aerial_image[,,4])
#' plot(if_raster(alb, dtm1m), main = "Albedo", col = gray(0:255/255))
albedo <- function(blue, green, red, nir, maxval = 255,
                   bluerange = c(430, 490), greenrange = c(535, 585),
                   redrange = c(610, 660), nirrange = c(835, 885)) {
  planck <- function(d) {
    d <- d * 1e-09
    h <- 6.62606957e-34
    cc <- 299792458
    tt <- 5250 + 273.16
    k <- 1.3806488e-23
    b <- 2 * h * cc ^ 2 / d ^ 5 * (1 / (exp((h * cc) / (k * tt * d)) - 1))
    b
  }
  r <- blue
  blue <- is_raster(blue)
  green <- is_raster(green)
  red <- is_raster(red)
  nir <- is_raster(nir)
  d <- c(1:3000)
  p <- planck(d)
  weight1 <- mean(p[bluerange[1]:bluerange[2]])
  weight2 <- mean(p[greenrange[1]:greenrange[2]])
  weight3 <- mean(p[redrange[1]:redrange[2]])
  weight4 <- mean(p[nirrange[1]:nirrange[2]])
  weightall <- weight1 + weight2 + weight3 + weight4
  rat1 <- weight1 / weightall
  rat2 <- weight2 / weightall
  rat3 <- weight3 / weightall
  rat4 <- weight4 / weightall
  albedo <- blue * rat1 + green * rat2 + red * rat3 + nir * rat4
  albedo <- albedo / maxval
  albedo <- if_raster(albedo, r)
  albedo
}
#' Adjusts albedo to correct for image brightness and contrast
#'
#' @description `albedo_adjust` is used to correct albedo derived from aerial imagery for image brightness and contrast using MODIS data.
#'
#' @param alb_image a raster object with numeric albedo values derived from high-resolution aerial imagery, as derived by [albedo()] and converted to a raster object.
#' @param alb_modis a raster object with numeric albedo values derived from MODIS imagery covering the same extent as `alb.image` (range 0 to 1). If the extent of the `alb.modis` is greater than that of `alb.image`, `alb.modis` is cropped.
#'
#' @return a raster object with numeric values representing the adjusted surface albedo values (range 0 to 1).
#' @import raster
#' @export
#'
#' @examples
#' library(raster)
#' alb <- albedo(aerial_image[,,1], aerial_image[,,2], aerial_image[,,3],
#'               aerial_image[,,4])
#' img <- if_raster(alb, dtm1m)
#' mds <- raster(modis, xmn = 169000, xmx = 170000, ymn = 12000, ymx = 13000)
#' alb2 <- albedo_adjust(img, mds)
#' par(mfrow=c(2, 1))
#' plot(img, main = "Raw albedo", col = gray(0:255/255))
#' plot(alb2, main = "Adjusted albedo", col = gray(0:255/255))
albedo_adjust <- function(alb_image, alb_modis) {
  alb_modis <- crop(alb_modis, extent(alb_image))
  alb_imagec <- resample(alb_image, alb_modis)
  logit_imagec <- log(alb_imagec / (1 - alb_imagec))
  logit_image <- log(alb_image / (1 - alb_image))
  logit_modis <- log(alb_modis / (1 - alb_modis))
  sd_adj <- sd(getValues(logit_modis)) / sd(getValues(logit_imagec))
  mn_adj <- mean(getValues(logit_modis)) - mean(getValues(logit_imagec))
  logit_image <- sd_adj * (logit_image + mn_adj)
  alb_adj <- 1 / (1 + exp(-1 * logit_image))
  alb_adj
}
#' Calculates Leaf Area Index
#'
#' @description `lai` is used to calculate the total one-sided area of leaf tissue per unit ground surface area from the Normalized Difference Vegetation Index.
#'
#' @param red a raster object, two-dimensional array or matrix of reflectance values in the red spectral band.
#' @param nir a raster object, two-dimensional array or matrix of reflectance values in the near-infrared spectral band.
#' @param maxlai an optional single numeric value representing the likely upper limit of Leaf Area Index to which all values are capped.
#'
#' @seealso function [lai_adjust()] for calculated leaf area index values at specified heights above the ground.
#'
#' @return a raster object or two-dimensional array of numeric values representing the total one-sided area of leaf tissue per unit ground surface area.
#' @import raster
#' @export
#'
#' @details
#' If `red` is a raster object, a raster object is returned. This function has been calibrated
#' using data derived from a small area of Cornwall only. It is strongly recommended that
#' locally calibrated values are obtained.

#' @examples
#' library(raster)
#' leaf <- lai(aerial_image[,,3], aerial_image[,,4])
#' plot(if_raster(leaf, dtm1m), main = "Leaf area index")
lai <- function(red, nir, maxlai = 20) {
  r <- red
  red <- is_raster(red)
  nir <- is_raster(nir)
  ndvi <- (nir - red) / (nir + red)
  logl <- 3.4838 * ndvi - 0.148
  l <- 10 ^ logl
  l[l > maxlai] <- maxlai
  if_raster(l, r)
}
#' Calculates Leaf Area Index for specified height above ground
#'
#' @description `lai_adjust` is used to adjust the total one-sided area of leaf tissue per unit ground surface area to derive values at at a specified height above the ground.

#' @param l raster object, two-dimensional array or matrix of leaf area index values as returned by [lai()].
#' @param veghgt a raster object, two-dimensional array or matrix of vegetation heights (m).
#' @param hgt a numeric value representing the height above the ground for which Leaf Area Index is required (m).
#' @import raster
#' @export
#'
#' @details
#' If `l` is a raster object, a raster object is returned. Temperatures are often required for a specified height above the ground,
#' and in short vegetation, the leaf area can be substantially less at this
#' height than at gorund level. This function enables the user to estimate
#' leaf area for a specified height about the ground.
#'
#' @return a raster object or a two-dimensional area of numeric values representing the Leaf Area Index values for a specified height above the ground
#'
#' @examples
#' library(raster)
#' l <- lai(aerial_image[,,3], aerial_image[,,4])
#' la<-lai_adjust(l, veg_hgt)
#' par(mfrow=c(2, 1))
#' plot(if_raster(l, dtm1m), main = "Leaf area index")
#' plot(if_raster(la, dtm1m), main = "Adjusted leaf area index")
lai_adjust <- function(l, veghgt, hgt = 0.05) {
  r <- l
  l <- is_raster(l)
  veghgt <- is_raster(veghgt)
  cf <- (veghgt - hgt) / veghgt
  cf[cf < 0] <- 0
  l <- l * cf
  if_raster(l, r)
}
#' Calculates leaf orientation
#'
#' @description `leaf_geometry` is used to calculates the ratio of vertical to horizontal projections of leaf foliage.
#'
#' @param veghgt a raster object, two-dimensional array or matrix of vegetation heights (m).
#' @param maxx a theoretical upper limit for the ratio of vertical to horizontal projections of leaf foliage, to which all values are capped.
#' @import raster
#' @export
#'
#' @details
#' Under vegetated canopies, canopy transmission not only decreases with canopy cover but is
#' also affected by leaf structure. At low solar angles, radiation is lower when leaves are
#' more vertically oriented and `leaf_geometry` is hence used to calculate an approximate  factor
#' indicating the degree to which vegetation is vertically orientated based on the premise that
#' shorter vegetation is more likely to have more vertically orientated leaves.
#' If `veghgt` is a raster object, a raster object is returned.
#' This function has been calibrated using data derived from a small area of
#' Cornwall only. It is strongly recommended that locally calibrated values
#  are obtained. The projection system associated with `veghgt` must be such that
#' units of x, y and z are identical. Use [projectRaster()] to convert the
#' projection to a Universal Transverse Mercator type projection system.

#' @return a raster object or a two-dimensional array of numeric values representing the
#' ratio of vertical to horizontal projections of leaf foliage. The output tends towards
#' zero as vegetation is more vertically orientated and maxx as it is more horizontally
#' orientated.
#'
#' @examples
#' library(raster)
#' x <- leaf_geometry(veg_hgt)
#' plot(x, main = "Leaf geometry")
leaf_geometry <- function(veghgt, maxx = 20) {
  r <- veghgt
  veghgt <- is_raster(veghgt)
  logx <- 5.5246 * log(veghgt + 1) - 1.4384
  x <- 10 ^ logx
  x[x > maxx] <- maxx
  if_raster(x, r)
}
#' Calculates canopy cover
#'
#' @description `canopy` is used to calculate fractional canopy cover.
#'
#' @param l a raster object, two-dimensional array or matrix of leaf area index values as returned by [lai()].
#' @param x a raster object, two-dimensional array of numeric values representing the ratio of vertical to horizontal projections of leaf foliage as returned by [leaf_geometry()].
#' @import raster
#' @export
#'
#' @return a raster object or a two-dimensional array of numeric values representing fractional canopy cover estimated as the proportion of isotropic radiation transmitted through the canopy.
#' @details
#' Canopy cover calculated by this function is defined as 1 - the proportion of isotropic
#' radiation transmitted through the canopy. If `l` is a raster object, a raster object is
#' returned.
#'
#' @examples
#' library(raster)
#' l <- lai(aerial_image[,,3], aerial_image[,,4])
#' l <- if_raster(l, dtm1m) # convert to raster
#' x <- leaf_geometry(veg_hgt)
#' fr <- canopy(l, x)
#' plot(fr, main = "Fractional canopy cover")
canopy <- function(l, x) {
  r <- l
  l <- is_raster(l)
  x <- is_raster(x)
  xs <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
  p1s <- c(-25.287, -25.254, -25.196, -24.932, -24.278, -22.088, -19.097,
           -15.255, -10.159, -7.105)
  p2s <- c(86.439, 86.42, 86.399, 86.306, 86.078, 85.517, 85.168, 85.228,
           85.963, 86.708)
  xp1 <- spline(xs, p1s, n = 1000)
  xp2 <- spline(xs, p2s, n = 1000)
  xl <- round(x * 100, 0)
  xl <- ifelse(xl < 1, 1, xl)
  xl <- ifelse(xl > 1000, 1000, xl)
  p1 <- array(xp1$y[xl], dim = dim(x))
  p2 <- array(xp2$y[xl], dim = dim(x))
  a <- p1 * l ^ (1 / 3) + p2
  arad <- a * (pi / 180)
  k <- sqrt(x ^ 2 + tan(arad) ^ 2) / (x + 1.774 * (x + 1.182) ^ -0.733)
  fr <- 1 - exp(-k * l)
  if_raster(fr, r)
}
#' Partitions surface albedo between ground and canopy albedo
#'
#' @description `albedo2` is used to calculate either bare ground or canopy albedo from surface albedo.
#'
#' @param alb a raster object, two-dimensional array or matrix of surface albedo values (range 0 - 1) derived using [albedo()] or [albedo_adjust()].
#' @param fr a raster object, two-dimensional array or matrix of fractional canopy cover as returned by [canopy()].
#' @param ground a logical value indicating whether to return ground albedo (TRUE) or canopy albedo (FALSE).
#' @import raster
#' @export
#'
#' @return If ground is `TRUE`, a raster object or a two-dimensional array of numeric values representing ground albedo (range 0 to 1).
#' @return If ground is `FALSE`, a raster object or two-dimensional array of numeric values representing canopy albedo (range 0 to 1).
#'
#' @details
#' If alb is a raster object, a raster object is returned. For calculation of net radiation, both
#' ground and canopy albedo may be are needed. Areas with high canopy cover typically have lower
#' albedo values than areas with low canopy cover and mean values for these can be derived from
#' parts of the image with very high or very low canopy cover. It is assumed that albedo of the
#' image as returned by [albedo()] is a function of both, weighted by canopy cover, such that
#' canopy albedos are closer to the image-derived albedos in areas of high canopy cover and ground
#' albedos closer to the image-derived albedo in areas with low canopy cover.
#'
#' @examples
#' library(raster)
#' # ======================
#' # Calculate image albedo
#' # ======================
#' alb <- albedo(aerial_image[,,1], aerial_image[,,2], aerial_image[,,3],
#'             aerial_image[,,4])
#' # ======================
#' # Calculate canopy cover
#' # ======================
#' l <- lai(aerial_image[,,3], aerial_image[,,4])
#' x <- leaf_geometry(veg_hgt)
#' fr <- canopy(l, x)
#' # ===========================================
#' # Calculate and plot ground and canopy albedo
#' # ===========================================
#' ag <- albedo2(alb, fr)
#' ac <- albedo2(alb, fr, ground = FALSE)
#' par(mfrow=c(2, 1))
#' plot(if_raster(ag, dtm1m), main = "Ground albedo", col = gray(0:255/255))
#' plot(if_raster(ac, dtm1m), main = "Canopy albedo", col = gray(0:255/255))
albedo2 <- function(alb, fr, ground = TRUE) {
  r <- alb
  alb <- is_raster(alb)
  fr <- is_raster(fr)
  mg <- mean(alb[fr < 0.03], na.rm = TRUE)
  mc <- mean(alb[fr > 0.97], na.rm = TRUE)
  alb2 <- alb * fr + mc * (1 - fr)
  if (ground) alb2 <- alb * (1 - fr) + mg * fr
  if_raster(alb2, r)
}
#' Calculates the mean albedo of surfaces surrounding each location
#'
#' @description `albedo_reflected` is used to calculate mean albedo of surfaces surrounding a location from which radiation is reflected.
#'
#' @param alb a raster object, two-dimensional array or matrix of surface albedo values (range 0 - 1) derived using [albedo()] or [albedo_adjust()].
#' @param e an optional extent object indicating the geographic extent of `alb`.
#' @import raster rgdal
#' @importFrom sp coordinates
#' @export
#'
#' @details
#' A small proportion of radiation received at any given location is in the
#' form of reflected radiation, and this function permits the albedo of
#' surrounding surfaces to be calculated. An inverse distance-weighting is
#' applied (range 0 to 1). If `alb` is a raster object, then a raster object is returned
#' and `e` can be determined from `alb`.
#'
#'
#' @return a raster object, two-dimensional array or matrix of values representing the mean albedo of surfaces surrounding each pixel of a two-dimension albedo array (range 0 - 1).
#'
#' @examples
#' library(raster)
#' alb <- albedo(aerial_image[,,1], aerial_image[,,2], aerial_image[,,3],
#'               aerial_image[,,4])
#' alb <- alb[901:1000, 901:1000]
#' e <- extent(c(169900, 170000, 12900, 13000))
#' rt <-raster(e, res = 1)
#' r_alb <- albedo_reflected(alb, e)
#' par(mfrow = c(2, 1))
#' plot(if_raster(alb, rt), main = "Surface albedo", col= gray(0:255/255))
#' plot(if_raster(r_alb, rt), main = "Albedo of surrounding surfaces",
#'      col= gray(0:255/255))
albedo_reflected <- function(alb, e = extent(alb)) {
  rr <- alb
  alb <- is_raster(rr)
  albr <- array(NA, dim = dim(alb))
  r <- raster(alb)
  extent(r) <- e
  res <- xres(r)
  maxdist <- max(e@xmax - e@xmin, e@ymax - e@ymin)
  s <- c((4:2000) / 4) ^ 2 * res
  s <- c(0, s)
  s <- s[s <= maxdist]
  dct <- runif(length(s), 0, 360)
  for (yy in 1:dim(alb)[1]) {
    for (xx in 1:dim(alb)[2]) {
      if (is.na(alb[yy, xx]) == FALSE) {
        x <- xx * res + e@xmin - res / 2
        y <- e@ymax + res / 2 - yy * res
        xdist <- round(s * sin(dct * pi / 180), 0)
        ydist <- round(s * cos(dct * pi / 180), 0)
        xy <- data.frame(x = x + xdist, y = y + ydist)
        coordinates(xy) <- ~x + y
        ar <- extract(r, xy)
        albr[yy, xx] <- mean(ar, na.rm = T)
      }
    }
  }
  if_raster(albr, rr)
}
#' Calculates the mean slope to the horizon
#'
#' @description `mean_slope` is used to calculates the mean slope to the horizon in all directions.
#'
#' @param dtm a raster object, two-dimensional array or matrix of elevations (m).
#' @param steps an optional integer. The mean slope is calculated from the horizon angle in specified directions. Steps defines the total number of directions used. If the default 36 is specified, then the horizon angle is calculated at 10º intervals.
#' @param res a single numeric value representing the spatial resolution of `dtm` (m).
#'
#' @return a raster object or a two-dimensional array of the mean slope angle to the horizon in all directions (º).
#' @import raster
#' @export
#'
#' @details
#' If `dtm` is a raster object, a raster object is returned.
#' The projection system associated with `dtm` must be such that
#' units of x, y and z are identical. Use [projectRaster()] to convert the
#' projection to a Universal Transverse Mercator type projection system.
#'
#' @examples
#' library(raster)
#' ms <- mean_slope(dtm100m, res = 100)
#' plot(ms, main = "Mean slope to horizon")
mean_slope <- function(dtm, steps = 36, res = 1) {
  r <- dtm
  dtm <- is_raster(dtm)
  dtm[is.na(dtm)] <- 0
  ha <- array(0, dim(dtm))
  for (s in 1:steps) {
    ha <- ha + atan(horizonangle(dtm, s * 360 / steps, res))
  }
  ha <- ha / steps
  ha <- tan(ha) * (180 / pi)
  if_raster(ha, r)
}
#' calculates the proportion of sky in view
#'
#' @description `skyviewtopo` is used to calculate a coefficient to correct for the proportion of sky in view when calculating net shortwave or longwave radiation above the canopy.
#'
#' @param dtm dtm a raster object, two-dimensional array or matrix of elevations (m).
#' @param steps an optional integer. The sky view is calculated from the horizon angle in specified directions. Steps defines the total number of directions used. If the default 36 is specified, then the horizon angle is calculated at 10º intervals.
#' @param res a single numeric value representing the spatial resolution of `dtm` (m).
#' @import raster
#' @export
#'
#' @details If a proportion of the sky of partially obscured, then the isotropic radiation flux received by a surface can be determined by integrating the single direction radiation flux over the proportion of sky in view. This function returns the integrated flux over the proportion of sky in view expressed as a proportion of the integrated flux over the entire hemisphere.
#'
#' @seealso The function [skyviewveg()] calculates a sky view correction factor underneath vegetation.
#'
#' @return a raster object or two-dimension array of values representing the proportion of isotropic radiation received by a partially obscured surface relative to the full hemisphere.
#'
#' @details
#' If `dtm` is a raster object, a raster object is returned.
#' The projection system associated with `dtm` must be such that
#' units of x, y and z are identical. Use [projectRaster()] to convert the
#' projection to a Universal Transverse Mercator type projection system.
#
#' @examples
#' library(raster)
#' sv <- skyviewtopo(dtm100m)
#' plot(sv, main = "Sky view factor")
skyviewtopo <- function(dtm, steps = 36, res = 100) {
  r <- dtm
  dtm <- is_raster(dtm)
  ha <- mean_slope(dtm, steps, res)
  ha <- ha * (pi / 180)
  svf <- 0.5 * cos(2 * ha) + 0.5
  if_raster(svf, r)
}
#' Calculates a sky view correction factor underneath vegetation
#'
#' @description
#' `skyviewveg` is used to calculate a coefficient to correct for the
#' proportion of sky obscured by topography when calculating net shortwave or
#' longwave radiation above the canopy.
#'
#' @param dtm a raster object, two-dimensional array or matrix of elevations (m).
#' @param l a raster object, two-dimensional array or matrix of leaf area index values as returned by [lai()].
#' @param x a raster object, two-dimensional array of numeric values representing the ratio of vertical to horizontal projections of leaf foliage as returned by [leaf_geometry()].
#' @param steps an optional integer. The sky view is calculated from the horizon angle in specified directions. Steps defines the total number of directions used. If the default 36 is specified, then the horizon angle is calculated at 10º intervals.
#' @param res a single numeric value representing the spatial resolution of `dtm` (m).
#' @import raster
#' @importFrom  stats runif sd spline
#' @export
#'
#' @details
#' If `dtm` is a raster object, a raster object is returned.
#' The projection system associated with `dtm` must be such that
#' units of x, y and z are identical. Use [projectRaster()] to convert the
#' projection to a Universal Transverse Mercator type projection system.
#' If a proportion of the sky of partially obscured, then the isotropic
#' radiation flux received by a surface underneath canopy can be determined by
#' integrating the single direction radiation transmission over the proportion
#' of sky in view. This function returns a computationally efficient
#' approximation of the integrated transmission over the proportion of sky in
#' view expressed as a proportion of the integrated transmission over the
#' entire hemisphere.
#'
#' @seealso The function [skyviewtopo()] calculates a sky view correction factor above vegetation.
#'
#' @return a raster object or a two-dimensional array of numeric values representing the proportion of isotropic radiation received by a surface partially obscured by topography relative to the full hemisphere underneath vegetation.
#'
#' @examples
#' library(raster)
#' l <- lai(aerial_image[,, 3], aerial_image[,, 4])
#' x <- leaf_geometry(veg_hgt)
#' sv <- skyviewveg(dtm1m, l, x)
#' plot(sv, main = "Sky view factor")
skyviewveg <- function(dtm, l, x, steps = 36, res = 1) {
  r <- dtm
  dtm <- is_raster(dtm)
  l <- is_raster(l)
  x <- is_raster(x)
  xs <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
  p1s <- c(0.5421, 0.542, 0.5409, 0.5375, 0.5265, 0.4775, 0.3941, 0.2805,
           0.1462, 0.0811)
  p2s <- c(0.6959, 0.6962, 0.6963, 0.6945, 0.689, 0.6692, 0.6422, 0.6169,
           0.6029, 0.6075)
  xp1 <- spline(xs, p1s, n = 1000)
  xp2 <- spline(xs, p2s, n = 1000)
  xl <- round(x * 100, 0)
  xl <- ifelse(xl < 1, 1, xl)
  xl <- ifelse(xl > 1000, 1000, xl)
  p1 <- array(xp1$y[xl], dim = dim(x))
  p2 <- array(xp2$y[xl], dim = dim(x))
  ha <- mean_slope(dtm, steps, res)
  ha <- ha * (pi / 180)
  a <- p1 * l ^ p2 + 0.564
  H2 <- (pi / 2) * (ha ^ a) / ((pi / 2) ^ a)
  svv <- 0.5 * cos(2 * H2) + 0.5
  if_raster(svv, r)
}
#' calculates net longwave radiation above canopy
#'
#' @description
#' `longwavetopo` is used to calculate a high-resolution dataset of the net
#' longwave radiation flux density emmited from the Earth, ignoring canopy
#' effects.
#'
#' @param h a single numeric value, raster object, two-dimensional array or matrix of specific humidities (\ifelse{html}{\out{kg kg<sup>{-1}</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc a single numeric value, raster object, two-dimensional array or matrix of temperatures (ºC).
#' @param p an optional single numeric value, raster object, two-dimensional array or matrix of sea-level pressures (Pa).
#' @param n a single numeric value, raster object, two-dimensional array or matrix of fractional cloud cover values (range 0 - 1).
#' Technically assumed to be 1 - ratio of measured to clear-sky radiation.
#' @param svf an optional single value, raster object, two-dimensional array or matrix of values representing the proportion of isotropic radiation received by a partially obscured surface relative to the full hemisphere, as returned by [skyviewtopo()].
#' @param co parameter relationship between vapor pressure and temperature near the ground Brutsaert (1975).
#' @import raster
#' @export
#'
#' @seealso The function [longwaveveg()] returns the net longwave radiation under vegetation.
#' The function [humidityconvert()] can be used to derive specific humidy from other meaures of humidity.
#' [cloudfromrad()] can be used to derive derive cloud cover from radiation.
#'
#' @details
#' if `svf` is a raster object, a raster object is returned.
#' If no values for `p` are provided, a default value of 101300 Pa, typical of
#' sea-level pressure, is assumed. If single values of `h`, `tc`, `p` and `n`
#' are given, and `svf` is an array or matrix, then the entire area is assumed to
#' have the same values of `h`, `tc`, `p` and `n`. If no value for `svf` is
#' provided then the entire hemisphere is assumed to be in view. If single
#' values of `h`, `tc`, `p` and `n` are given, and no value of `svf' is
#' provided, a single value is returned, and it is assumed that the entire
#' hemisphere is in view.
#'
#' @return a single numeric value, raster object, two-dimensional array pr matrix of values representing net longwave radiation (MJ per metre squared per hour).
#'
#' @examples
#' library(raster)
#' # =================================
#' # Extract data for 2010-05-24 11:00
#' # =================================
#' h <- huss[,,144]
#' p <- pres[,,144]
#' tc <- tas[,,144] + dtr[,,144]
#' n <-cfc[,,3444]
#' sv <- skyviewtopo(dtm100m)
#' # ===========================
#' # Resample to 100m resolution
#' # ===========================
#' hr <- if_raster(h, dtm1km)
#' tr <- if_raster(tc, dtm1km)
#' pr <- if_raster(p, dtm1km)
#' nr <- raster(n, xmn = -5.40, xmx = -5.00, ymn = 49.90, ymx = 50.15)
#' crs(nr) <- '+init=epsg:4326'
#' nr <- projectRaster(nr, crs = '+init=epsg:27700')
#' hr <- resample(hr, dtm100m)
#' tr <- resample(tr, dtm100m)
#' pr <- resample(pr, dtm100m)
#' nr <- resample(nr, dtm100m)
#' # =========================================
#' # Calculate and plot net longwave radiation
#' # =========================================
#' netlong100m <- longwavetopo(hr, tr, pr, nr, sv)
#' netlong100m <- mask(netlong100m, dtm100m)
#' plot(netlong100m, main = "Net longwave radiation")
longwavetopo <- function(h, tc, p = 101300, n, svf = 1, co = 1.24) {
  r <- svf
  h <- is_raster(h)
  tc <- is_raster(tc)
  p <- is_raster(p)
  n <- is_raster(n)
  svf <- is_raster(svf)
  pk <- p/1000
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  ws <- 0.622 * e0/pk
  rh <- (h/ws) * 100
  rh <- ifelse(rh > 100, 100, rh)
  ea <- e0 * (rh/100)
  eo <- co * (10 * ea /(tc + 273.15))^(1/7)
  em <- n + (1 - n) * eo
  Ln <- (1 - em) * 2.043e-10 * (tc + 273.15)^4
  lwr <- Ln * svf
  if_raster(lwr, r)
}
#' Calculates net longwave radiation below canopy
#'
#' @description `longwaveveg` is used to calculate a high-resolution dataset of the net longwave radiation flux density emmited from the Earth, accounting for canopy effects.
#'
#' @param h a single numeric value, raster object, two-dimensional array or matrix of specific humidities (\ifelse{html}{\out{kg kg<sup>{-1}</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc a single numeric value, raster object, two-dimensional array or matrix of temperatures (ºC).
#' @param p an optional single numeric value, raster object, two-dimensional array or matrix of sea-level pressures (Pa).
#' @param n a single numeric value, raster object, two-dimensional array or matrix of fractional cloud cover (range 0 - 1).
#' Technically assumed to be 1 - ratio of measured to clear-sky radiation.
#' @param x a raster object, two-dimensional array or matrix of numeric values representing the ratio of vertical to horizontal projections of leaf foliage as returned by [leaf_geometry()].
#' @param fr a raster object, two-dimensional array or matrix of fractional canopy cover as returned by [canopy()].
#' @param svv an optional raster object, two-dimensional array or matrix of values representing the proportion of isotropic radiation received by a surface partially obscured by topography relative to the full hemisphere underneath vegetation as returned by [skyviewveg()].
#' @param albc an optional single value, raster object, two-dimensional array or matrix of values representing the albedo(s) of the vegetated canopy as returned by [albedo2()].
#' @param co parameter relationship between vapor pressure and temperature near the ground Brutsaert (1975).
#' @import raster
#' @export
#'
#' @seealso The function [longwavetopo()] returns the net longwave radiation above vegetation.
#' The function [humidityconvert()] can be used to derive specific humidy from other meaures of humidity.
#' [cloudfromrad()] can be used to derive derive cloud cover from radiation.
#'
#' @return a single numeric value, raster object or two-dimensional array of values representing net longwave radiation (MJ per metre squared per hour).
#'
#' @details
#' If `svv` is a raster object, a raster object is returned.
#' If no values for `p` are provided, a default value of 101300 Pa, typical of
#' sea-level pressure, is assumed. If no value for `albc` is provided, then
#' the entire area is assumed to have a default value of 0.23, typical of
#' well-watered grass. If single values of `h`, `tc`, `p`, `n` or
#' `albc` are given, then the entire area is assumed to have the same values.
#' If no value for `svv` is provided then the entire hemisphere is assumed to
#' be in view.
#'
#' @examples
#' library(raster)
#' # =================================
#' # Extract data for 2010-05-24 11:00
#' # =================================
#' h <- microvars$humidity[564]
#' p <- microvars$pressure[564]
#' n <- microvars$cloudcover[264]
#' tcr <- raster(temp100[,,564], xmn = 169000, xmx = 170000, ymn = 12000, ymx = 13000)
#' tc <- resample(tcr, dtm1m) # Resample temperature raster to 1m
#' # =========================
#' # calculate input variables
#' # =========================
#' x <- leaf_geometry(veg_hgt)
#' l <- lai(aerial_image[,,3], aerial_image[,,4])
#' l <- lai_adjust(l, veg_hgt)
#' svv <- skyviewveg(dtm1m, l, x)
#' fr <- canopy(l, x)
#' alb <- albedo(aerial_image[,,1], aerial_image[,,2], aerial_image[,,3],
#'              aerial_image[,,4])
#' albc <- albedo2(alb, fr, ground = FALSE)
#' # =====================================
#' # calculate and plot longwave radiation
#' # =====================================
#' netlong1m <-longwaveveg(h, tc, p, n, x, fr, svv, albc)
#' nlr <- mask(netlong1m, dtm1m)
#' plot(nlr, main = "Net longwave radiation")
longwaveveg <- function(h, tc, p = 101300, n, x, fr, svv = 1, albc = 0.23, co = 1.24) {
  rr <- svv
  h <- is_raster(h)
  tc <- is_raster(tc)
  p <- is_raster(p)
  n <- is_raster(n)
  x <- is_raster(x)
  fr <- is_raster(fr)
  svv <- is_raster(svv)
  albc <- is_raster(albc)
  pk <- p / 1000
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ws <- 0.622 * e0 / pk
  rh <- (h / ws) * 100
  rh <- ifelse(rh > 100, 100, rh)
  ea <- e0 * rh / 100
  eo <- co * (10 * ea / (tc + 273.15)) ^ (1/7)
  em <- n + (1 - n) * eo
  le0 <- 2.043e-10 * (tc + 273.15) ^ 4
  lwsky <- em * le0
  lw1 <- (1 - fr) * lwsky
  lr <- (2 / 3) * log(x + 1)
  r <- 1 / (1 + exp(-1 * lr))
  lw2 <- albc * r * fr * lwsky
  lw3 <- r * (1 - albc) * fr * 2.043e-10 * (tc + 273.15) ^ 4
  lwd <- lw1 + lw2 + lw3
  lwrad <- (le0 - lwd) * svv
  if_raster(lwrad, rr)
}
#' Downscales shortwave radiation accounting for topographic effects
#'
#' @description
#' `shortwavetopo` is used to downscale components of the flux density of
#' shortwave radiation received at the surface of the Earth using a
#' high-resolution digital elevation dataset, ignoring canopy effects.
#'
#' @param dni a single numeric value, raster object, two-dimensional array or matrix of coarse-resolution direct radiation perpendicular to the solar beam (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}}).
#' @param dif a single numeric value, raster object, two-dimensional array or matrix of diffuse radiation horizontal ot the surface (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}}).
#' @param julian a single integer representing the Julian as returned by [julday()].
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock).
#' @param lat a single numeric value representing the mean latitude of the location for which downscaled radiation is required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the mean longitude of the location for which downscaled radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param dtm an optional raster object, two-dimensional array or matrix of elevations (m), orientated as if derived using [is_raster()]. I.e. `[1, 1]` is the NW corner.
#' @param slope a single value, raster object, two-dimensional array or matrix of slopes (º). If an array or matrix, then orientated as if derived using [is_raster()]. I.e. `[1, 1]` is the NW corner.
#' @param aspect a single value, raster object, two-dimensional array or matrix of aspects (º). If an array or matrix, then orientated as if derived using [is_raster()]. I.e. `[1, 1]` is the NW corner.
#' @param svf an optional single value, raster object, two-dimensional array or matrix of values representing the proportion of isotropic radiation received by a partially obscured surface relative to the full hemisphere as returned by [skyviewtopo()].
#' @param alb an optional single value, raster object, two-dimensional array or matrix of surface albedo(s) (range 0 - 1) derived using [albedo()] or [albedo_adjust()].
#' @param albr an optional single value, raster object, two-dimensional array or matrix of values of albedo(s) of adjacent surfaces (range 0 - 1) as returned by [albedo_reflected()].
#' @param ha an optional raster object, two-dimensional array or matrix of values representing the mean slope to the horizon (decimal degrees) of surrounding surfaces from which radiation is reflected for each cell of `dtm` as returned by [mean_slope()].
#' @param res a single numeric value representing the spatial resolution of `dtm` (m).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param shadow an optional logical value indicating whether topographic shading should be considered (False = No, True = Yes).
#' @param component an optional character string of the component of radiation to be returned. One of "sw" (net shortwave radiation, i.e. accounting for albedo), "sw2" (total incoming shortwave radiation), "dir" (direct), "dif" (diffuse), "iso" (isotropic diffuse), "ani" (anistopic diffuse), "ref" (reflected).
#' @param difani an optinional logical indicating whether to treat a proportion of the diffuse radiation as anistropic (see details).
#' @import raster
#' @export
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
#' object is returned. If `dtm` is a raster object, a raster object is returned. If `dni` or `dif` are raster
#' objects, two-dimensional arrays or matrices, then it is assumed that they have been
#' derived from coarse-resolution data by interpolation, and have the same extent as `dtm`.
#' If no value for `ha` is provided, the mean slope to the horizon is assumed
#' to be 0. If no value for `svv` is provided, then the entire hemisphere is
#' assumed to be in view. If no value for `svf` is provided, then the entire
#' hemisphere is assumed to be in view. If values of `alb` and `albr` are not
#' specified, then a default value of 0.23, typical of well-watered grass is
#' assumed. If single values of `alb` and `albr` are given, then the entire
#' area is assumed to have the same albedo. If `dtm` is specified, then the
#' projection system used must be such that the units of x, y and z are
#' identical. Use [projectRaster()] to convert the projection to a Universal
#' Transverse Mercator type projection system. If no value for `dtm` is
#' provided, radiation is downscaled by deriving values on the inclined
#' surfaces specified in `slope` and `aspect` and topographic shadowing is
#' ignored. If single values are provided for `slope` and `aspect` single
#' values of components of shortwave radiation for an inclined surface are
#' returned. Only single values of `lat` and `long` are taken as inputs. Under partially
#' cloudy conditions, a proportion of diffuse radiation is typically anistropic
#' (directional). If `difani` is TRUE (the default), then the assumption is made that
#' hourly direct radiation transmission can define the portions of the diffuse
#' radiation to be treated as anistropic and isotropic.
#' If `difani` is FALSE, all diffuse radiation is treated as isotropic.
#' If `dtm` covers a large extent, the `dtm` is best divided into blocks and
#' seperate calculations performed on each block. Since horizon angles,
#' topographic shading and sky view correction factors may be influenced by
#' locations beyond the extent of `dtm`, it is best to ensure `dtm` covers a
#' larger extent than that for which radiation values are needed, and to
#' ensure sub-divided blocks overlap in extent. Calculations are faster if values
#' for all inputs are provided.
#'
#' @seealso Function [shortwaveveg()] returns net shortwave radiation below a canopy.
#'
#' @return If component is "sw", a raster object or two-dimensional array of numeric values representing net shortwave radiation (MJ m^-2 hr^-1).
#' @return If component is "sw2", a raster object or two-dimensional array of numeric values representing total incoming shortwave radiation (MJ m^-2 hr^-1).
#' @return If component is "dir", a raster object or two-dimensional array of numeric values representing direct shortwave radiation (MJ m^-2 hr^-1).
#' @return If component is "dif", a raster object or two-dimensional array of numeric values representing diffuse shortwave radiation (MJ m^-2 hr^-1).
#' @return If component is "iso", a raster object or two-dimensional array of numeric values representing isotropic diffuse shortwave radiation (MJ m^-2 hr^-1).
#' @return If component is unspecified, then the default "sw" is returned.
#'
#' @examples
#' library(raster)
#' # =================================
#' # Extract data for 2010-05-24 11:00
#' # =================================
#' dni <-dnirad[,,3444]
#' dif <-difrad[,,3444]
#' # ===========================
#' # Resample to 100m resolution
#' # ===========================
#' dnir <- raster(dni, xmn = -5.40, xmx = -5.00, ymn = 49.90, ymx = 50.15)
#' difr <- raster(dif, xmn = -5.40, xmx = -5.00, ymn = 49.90, ymx = 50.15)
#' crs(dnir) <- '+init=epsg:4326'
#' crs(difr) <- '+init=epsg:4326'
#' dnir <- projectRaster(dnir, crs = "+init=epsg:27700")
#' difr <- projectRaster(difr, crs = "+init=epsg:27700")
#' dni <- resample(dnir, dtm100m)
#' dif <- resample(difr, dtm100m)
#' sv <- skyviewtopo(dtm100m)
#' jd <- julday(2010, 5, 24)
#' ha <- mean_slope(dtm100m)
#' # ================================================================
#' # Calculate and plot net shortwave radiation for 2010-05-24 11:00
#' # ================================================================
#' netshort100m <- shortwavetopo(dni, dif, jd, 11, dtm = dtm100m,
#'                               svf = sv, ha = ha)
#' plot(mask(netshort100m, dtm100m),
#'      main = "Net shortwave radiation")
shortwavetopo <- function(dni, dif, julian, localtime, lat = NA, long = NA,
                          dtm = array(0, dim = c(1, 1)), slope = NA,
                          aspect = NA, svf = 1, alb = 0.23, albr = 0.23,
                          ha = 0, res = 100, merid = round(long / 15, 0) * 15, dst = 0,
                          shadow = TRUE, component = "sw", difani = TRUE) {
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
  dni <- is_raster(dni)
  dif <- is_raster(dif)
  slope <- is_raster(slope)
  aspect <- is_raster(aspect)
  svf <- is_raster(svf)
  alb <- is_raster(alb)
  albr <- is_raster(albr)
  ha <- is_raster(ha)
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
  isor <- 0.5 * dif * (1 + cos(a)) * (1 - k) * svf
  cisr <- k * dif * si
  sdi <- (slope + ha) * (pi / 180)
  refr <- 0.5 * albr * (1 - cos(sdi)) * dif
  difr <- isor + cisr + refr
  sw2r <- difr + dirr
  swr <- (1 - alb) * sw2r
  rad <- NA
  if (component == "sw") rad <- swr
  if (component == "sw2") rad <- sw2r
  if (component == "dir") rad <- dirr
  if (component == "dif") rad <- difr
  if (component == "iso") rad <- isor
  if (component == "ani") rad <- cisr
  if (component == "ref") rad <- refr
  if_raster(rad, r)
}
#' Downscales net shortwave radiation accounting for topography and vegetation
#'
#' @description
#' `shortwaveveg` is used to downscale the flux density of shortwave radiation
#' received at the surface of the Earth, accounting for both topographic and
#' canopy effects.
#'
#' @param dni a single numeric value, raster object, two-dimensional array or matrix of coarse-resolution direct radiation perpendicular to the solar beam (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}}).
#' @param dif a single numeric value, raster object, two-dimensional array or matrix of coarse-resolution diffuse radiation horizontal ot the surface (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}}).
#' @param julian a single integer representing the Julian as returned by [julday()].
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock).
#' @param lat an optional single numeric value representing the mean latitude of the location for which downscaled radiation is required (decimal degrees, -ve south of equator).
#' @param long an optional single numeric value representing the mean longitude of the location for which downscaled radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param dtm an optional raster object, two-dimensional array or matrix of elevations (m), orientated as if derived using [is_raster()]. I.e. `[1, 1]` is the NW corner.
#' @param slope an optional single value, raster object, two-dimensional array or matrix of slopes (º). If an array or matrix, then orientated as if derived using [is_raster()]. I.e. `[1, 1]` is the NW corner.
#' @param aspect an optional single value, raster object, two-dimensional array or matrix of aspects (º). If an array or matrix, then orientated as if derived using [is_raster()]. I.e. `[1, 1]` is the NW corner.
#' @param svv an optional raster object, two-dimensional array or matrix of values representing the proportion of isotropic radiation received by a surface partially obscured by topography relative to the full hemisphere underneath vegetation as returned by [skyviewveg()].
#' @param albg an optional single value, raster object, two-dimensional array or matrix of values representing the albedo(s) of the ground as returned by [albedo2()].
#' @param fr a raster object, two-dimensional array or matrix of fractional canopy cover as returned by [canopy()].
#' @param albr an optional single value, raster object, two-dimensional array or matrix of values representing the albedo(s) of adjacent surfaces as returned by [albedo_reflected()].
#' @param ha an optional raster object, two-dimensional array or matrix of values representing the mean slope to the horizon (decimal degrees) of surrounding surfaces from which radiation is reflected for each cell of `dtm` as returned by [mean_slope()].
#' @param res a single numeric value representing the spatial resolution of `dtm` (m).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param shadow an optional logical value indicating whether topographic shading should be considered (False = No, True = Yes).
#' @param x a raster object, two-dimensional array or matrix of numeric values representing the ratio of vertical to horizontal projections of leaf foliage as returned by [leaf_geometry()].
#' @param l a raster object, two-dimensional array or matrix of leaf area index values as returned by [lai()].
#' @param difani an optinional logical indicating whether to treat a proportion of the diffuse radiation as anistropic (see details).
#' @import raster
#' @export
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
#' object is returned. If `dtm` is a raster object, a raster object is returned. If `dni` or `dif` are raster
#' objects, two-dimensional arrays or matrices, then it is assumed that they have been
#' derived from coarse-resolution data by interpolation, and have the same extent as `dtm`.
#' If no value for `ha` is provided, the mean slope to the horizon is assumed
#' to be 0. If no value for `svv` is provided, then the entire hemisphere is
#' assumed to be in view. If values of `albg` and `albr` are not specified,
#' then a default value of 0.23, typical of well-watered grass is assumed. If
#' single values of `albg` and `albr` are given, then the entire area is
#' assumed to have the same albedo. If `dtm` is specified, then the projection
#' system used must be such that the units of x, y and z are identical. Use
#' [projectRaster()] to convert the projection to a Universal Transverse
#' Mercator type projection system. If no value for `dtm` is provided,
#' radiation is downscaled by deriving values on the inclined surfaces
#' specified in `slope` and `aspect` and topographic shadowing is ignored. If
#' single values are provided for `slope` and `aspect`, the entire extent
#' covered by `fr` is assumed to have the same slope and aspect. Only single
#' values of `lat` and `long` are taken as inputs. Under partially
#' cloudy conditions, a proportion of diffuse radiation is typically anistropic
#' (directional). If `difani` is TRUE (the default), then the assumption is made that
#' hourly direct radiation transmission can define the portions of the diffuse
#' radiation to be treated as anistropic and isotropic. If `dtm` covers a large
#' extent, the `dtm` is best divided into blocks and seperate calculations
#' performed on each block. Since horizon angles, topographic shading and
#' sky view correction factors may be influenced by locations beyond the extent of `dtm`, it is best to ensure
#' `dtm` covers a larger extent than that for which radiation values are
#' needed, and to ensure sub-divided blocks overlap in extent. Calculations are faster
#' if values for all inputs are provided.
#'
#' @seealso Function [shortwavetopo()] returns net shortwave radiation, or components thereof, above the canopy.
#'
#' @return a raster object, two-dimensional array of numeric values representing net shortwave radiation (MJ per metre squared per hour).
#' The raster package function [terrain()] can be used to derive slopes and aspects from `dtm` (see example).
#'
#' @examples
#' library(raster)
#' # =================================
#' # Extract data for 2010-05-24 11:00
#' # =================================
#' dni <- microvars$dni[564]
#' dif <- microvars$dif[564]
#' # ==========================
#' # Calculate input paramaters
#' # ==========================
#' x <- leaf_geometry(veg_hgt)
#' l <- lai(aerial_image[,,3], aerial_image[,,4])
#' l <- lai_adjust(l, veg_hgt)
#' fr <- canopy(l, x)
#' alb <- albedo(aerial_image[,,1], aerial_image[,,2], aerial_image[,,3],
#'              aerial_image[,,4])
#' albg <- albedo2(alb, fr)
#' sv <- skyviewveg(dtm1m, l, x)
#' jd <- julday(2010, 5, 24)
#' ha <- mean_slope(dtm1m)
#' # ===============================================================
#' # Calculate and plot net shortwave radiation for 2010-05-24 11:00
#' # ===============================================================
#' netshort1m <- shortwaveveg(dni, dif, jd, 11, dtm = dtm1m, svv = sv, albg = albg,
#'                            fr = fr, ha = ha, x = x, l = l)
#' plot(mask(netshort1m, dtm1m), main = "Net shortwave radiation")
shortwaveveg <- function(dni, dif, julian, localtime, lat = NA, long = NA,
                         dtm = array(0, dim = c(1, 1)), slope = NA, aspect = NA,
                         svv = 1, albg = 0.23, fr, albr = 0.23, ha = 0,
                         res = 1, merid = round(long / 15, 0) * 15, dst = 0, shadow = TRUE,
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
  dni <- is_raster(dni)
  dif <- is_raster(dif)
  slope <- is_raster(slope)
  aspect <- is_raster(aspect)
  svv <- is_raster(svv)
  albg <- is_raster(albg)
  fr <- is_raster(fr)
  albr <- is_raster(albr)
  ha <- is_raster(ha)
  x <- is_raster(x)
  l <- is_raster(l)
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
  kk <- ((x ^ 2 + 1 / (tan(saltitude * (pi / 180)) ^ 2)) ^ 0.5) /
        (x + 1.774 * (x + 1.182) ^ (-0.733))
  trd <- exp(-kk * l)
  trf <- (1 - fr)
  fgd <- fd * trd * (1 - albg)
  fged <- fdf * trf * (1 - albg) * svv
  fgc <- fgd + fged
  if_raster(fgc, r)
}
#' Calculates the diffuse fraction from incoming shortwave radiation
#'
#' @description `difprop` calculates proportion of incoming shortwave radiation that is diffuse radiation using the method of Skartveit et al. (1998) Solar Energy, 63: 173-183.
#'
#' @param rad a vector of incoming shortwave radiation values (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}})
#' @param julian the Julian day as returned by [julday()]
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock)
#' @param lat a single numeric value representing the latitude of the location for which partitioned radiation is required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the longitude of the location for which partitioned radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param hourly specifies whether values of `rad` are hourly (see details).
#' @param watts a logical value indicating  whether the units of `rad` are \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}} (TRUE) or \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} (FALSE).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param corr an optional numeric value representing a correction to account for over- or under-estimated diffuse proportions. Values > 1 will apportion a greater ammount of total radiation as diffuse than originally calculated by the formula.
#'
#' @return a vector of diffuse fractions (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}}).
#' @export
#'
#' @details
#' The method assumes the environment is snow free. Both overall cloud cover and heterogeneity in
#' cloud cover affect the diffuse fraction. Breaks in an extensive cloud deck may primarily
#' enhance the beam irradiance, whereas scattered clouds may enhance the diffuse irradiance and
#' leave the beam irradiance unaffected.  In consequence, if hourly data are available, an index
#' is applied to detect the presence of such variable/inhomogeneous clouds, based on variability
#' in radiation for each hour in question and values in the preceding and deciding hour.  If
#' hourly data are unavailable, an average variability is determined from radiation intensity.
#'
#' @examples
#' rad <- c(5:42) / 0.036 # typical values of radiation in W/m^2
#' jd <- julday(2017, 6, 21) # julian day
#' dfr <- difprop(rad, jd, 12, 50, -5)
#' plot(dfr ~ rad, type = "l", lwd = 2, xlab = "Incoming shortwave radiation",
#'      ylab = "Diffuse fraction")
difprop <- function(rad, julian, localtime, lat, long, hourly = FALSE,
                  watts = TRUE, merid = round(long / 15, 0) * 15, dst = 0,
                  corr = 1) {
  if (watts) rad <- rad * 0.0036
  sa <- solalt(localtime, lat, long, julian, merid, dst)
  alt <- sa * (pi / 180)
  k1 <- 0.83 - 0.56 * exp(- 0.06 * sa)
  si <- cos(pi / 2 - alt)
  si[si < 0] <- 0
  k <- rad / (4.87 * si)
  k[is.na(k)] <- 0
  k <- ifelse(k > k1, k1, k)
  k[k < 0] <- 0
  rho <- k / k1
  if (hourly) {
    rho <- c(rho[1], rho, rho[length(rho)])
    sigma3  <- 0
    for (i in 1:length(rad)) {
      sigma3[i] <- (((rho[i + 1] - rho[i]) ^ 2 + (rho[i + 1] - rho[i + 2]) ^ 2)
                    / 2) ^ 0.5
    }
  } else {
    sigma3a <- 0.021 + 0.397 * rho - 0.231 * rho ^ 2 - 0.13 *
               exp(-1 * (((rho - 0.931) / 0.134) ^ 2) ^ 0.834)
    sigma3b <- 0.12 + 0.65 * (rho - 1.04)
    sigma3 <- ifelse(rho <= 1.04, sigma3a, sigma3b)
  }
  k2 <- 0.95 * k1
  d1 <- ifelse(sa > 1.4, 0.07 + 0.046 * (90 - sa) / (sa + 3), 1)
  K <- 0.5 * (1 + sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2))
  d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K ^ 2))
  d3 <- (d2 * k2) * (1 - k) / (k * (1 - k2))
  alpha <- (1 / sin(alt)) ^ 0.6
  kbmax <- 0.81 ^ alpha
  kmax <- (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2))
  dmax <- (d2 * k2) * (1 - kmax) / (kmax * (1 - k2))
  d4 <- 1 - kmax * (1 - dmax) / k
  d <- ifelse(k <= kmax, d3, d4)
  d <- ifelse(k <= k2, d2, d)
  d <- ifelse(k <= 0.22, 1, d)
  kX <- 0.56 - 0.32 * exp(-0.06 * sa)
  kL <- (k - 0.14) / (kX - 0.14)
  kR <- (k - kX) / 0.71
  delta <- ifelse(k >= 0.14 & k < kX, -3 * kL ^ 2 *(1 - kL) * sigma3 ^ 1.3, 0)
  delta <- ifelse(k >= kX & k < (kX + 0.71), 3 * kR * (1 - kR) ^ 2 * sigma3 ^
                  0.6, delta)
  d[sigma3 > 0.01] <- d[sigma3 > 0.01] + delta[sigma3 > 0.01]
  d[rad == 0] <- 0.5
  d[sa < 0] <- 1
  # apply correction
  dif_val <- rad * d
  dif_val_adj <- dif_val * corr
  d <- dif_val_adj /rad
  d[d > 1] <- 1
  d[d < 0] <- 1
  d[is.na(d)] <- 0.5
  d
}
#' calculates clearsky radiation
#'
#' @description
#' `clearskyrad` is used to calculate clear-sky shortwave irradiance using the
#'  Crawford & Duchon (1999) method
#' @param tme a single value or vector of POSIXlt objects indicating the time(s)
#' for which clearksy radiation is required.
#' @param lat latitude in decimal degrees
#' @param long longitude in decimal degrees
#' @param h an optional single value or vector of specific humidities (\ifelse{html}{\out{kg kg<sup>{-1}</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc an optional single value or vector of temperatures (ºC).
#' @param p an optional single value or vector of pressures (Pa).
#' @param G an optional single value or vector describing he moisture profile
#' in the atmosphere (per Smith 1966).
#' @param Ie an optional single value for extra-terrestrail radiation to permit adjustment for
#' sun-earth distances (see details).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @seealso The function [cloudfromrad()] uses this function to return 1 - ratio of measured to clearsky
#' radiation for input when computing longwave radiation when using [longwavetopo()] or [longwaveveg()]
#' @export
#'
#' @details
#' The units returned are the same as for `Ie`, with the default option in W / M^2.
#' If no values for `p` are provided, a default value of 101300 Pa, typical of
#' sea-level pressure, is assumed. The method used is that detailed in
#' Crawford & Duchon (1999) Quarterly Journal of the Royal Meteorological
#' Society 122: 1127-1151. The method is not greatly sensitive to humidity,
#' temperature and pressure so approximate values can be provided, or the defaults
#' chosen, if these data are unavailable.
#'
#' @return a single value or vector of clearksy radiation.
#'
#' @examples
#' tme <- as.POSIXlt(c(0:23) * 3600, origin = "2010-05-23 00:00", tz = "GMT")
#' Io <- clearskyrad(tme, 50, -5, 0.007953766 , 11)
#' plot(Io ~ as.POSIXct(tme), type = "l")
clearskyrad <- function(tme, lat, long, h = 0.00697, tc = 15, p = 101300, G = 2.78, Ie = 1352.778,
                        merid = round(long/15, 0) * 15, dst = 0) {
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  lt <- tme$hour + tme$min / 60 + tme$sec / 3600
  sa <- solalt(lt, lat, long, jd, merid, dst)
  sa[sa < 0] <- NA
  z <- (90 - sa) * (pi / 180)
  m <- 35 * cos(z) * ((1224 * cos(z)^2 + 1)^(-0.5))
  TrTpg <- 1.021 - 0.084 * (m * 0.000949 * 0.01 * p + 0.051)^0.5
  pk <- p / 1000
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  ws <- 0.622 * e0/pk
  rh <- (h/ws) * 100
  rh <- ifelse(rh > 100, 100, rh)
  xx <- log(rh / 100) + ((17.27 * tc) / (237.3 + tc))
  Td <- (237.3 * xx) / (17.27 - xx)
  u <- exp(0.1133 - log(G + 1) + 0.0393 * Td)
  Tw <- 1 - 0.077 * (u * m) ^ 0.3
  Ta <- 0.935 * m
  od <- TrTpg * Tw * Ta
  Ic <- Ie * (cos(z)) * TrTpg * Tw * Ta
  Ic[Ic > Ie] <- NA
  Ic[Ic < 0] <- NA
  Ic
}
#' calculates cloud cover form shortwave radiation
#'
#' @description
#' `cloudfromrad` is used to derive a cloud cover index as 1 - the ratio of measured to clearksy radiation
#' missing values (e.g. at night) are interpolated
#' @param rad a single numeric value or vector of solar irradiance(s). Units must be the same
#' as those for `Ie`. Default is Watts /m^2.
#' @param tme a single value or vector of POSIXlt objects indicating the time(s)
#' for which clearksy radiation is required.
#' @param lat latitude in decimal degrees
#' @param long longitude in decimal degrees
#' @param h a single value or vector of specific humidities (\ifelse{html}{\out{kg kg<sup>{-1}</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc a single value or vector of temperatures (ºC).
#' @param p an optional single value or vector of pressures (Pa).
#' @param G an optional single value or vector describing he moisture profile
#' in the atmosphere (per Smith 1966).
#' @param Ie an optional single value for extra-terrestrail radiation to permit adjustment for
#' sun-earth distances (see details).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @import zoo
#' @export
#'
#' @seealso The function [clearskyrad()] is uses to derive clear-sky irradiance. Can be used
#' to derive cloud covers for computing longwave radiation when using [longwavetopo()] or [longwaveveg()]
#'
#' @details
#' If no values for `p` are provided, a default value of 101300 Pa, typical of
#' sea-level pressure, is assumed. The method used is that detailed in
#' Crawford & Duchon (1999) Quarterly Journal of the Royal Meteorological
#' Society 122: 1127-1151. The method is not greatly sensitive to humidity,
#' temperature and pressure so approximate values can be provided, or the defaults
#' chosen, if these data are unavailable.
#'
#' @return a single value or vector of clearksy radiation.
#'
#' @examples
#' tme <- as.POSIXlt(c(0:23) * 3600, origin = "2010-05-23 00:00", tz = "GMT")
#' rad <- clearskyrad(tme, 50, -5, 0.007953766 , 11) * 0.75
#' cfc <- cloudfromrad(rad, tme, 50, -5, 0.007953766 , 11)
#' plot(cfc ~ as.POSIXct(tme), type = "l") # should be 0.25
cloudfromrad <- function(rad, tme, lat, long, h = 0.00697, tc = 15, p = 101300, G = 2.78,
                         Ie = 1352.778, merid = round(long/15, 0) * 15, dst = 0) {
  Ic <- clearskyrad(tme, lat, long, h, tc, p, G, Ie, merid, dst)
  s <- rad / Ic
  s[s > 1] <- 1
  s[s < 0] <- 0
  if (is.na(s[1])) s[1] <- mean(s, na.rm = T)
  if (is.na(s[length(s)])) s[length(s)] <- mean(s, na.rm = T)
  s <- na.approx(s)
  cfc <- 1 - s
  cfc
}

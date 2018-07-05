#' A 1 m resolution aerial image.
#'
#' A dataset containing image reflectance values for the area bounded by 169000, 170000, 12000,
#' 13000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: +init=epsg:27700).
#'
#' @format An array with 1000 rows, 1000 columns and 4 layers:
#' \describe{
#'   \item{1}{blue band reflectance values (430 to 490 nm), in the range 0 to 255}
#'   \item{2}{green band reflectance values (535 to 585 nm), in the range 0 to 255}
#'   \item{3}{red band reflectance values (610 to 660 nm), in the range 0 to 255}
#'   \item{4}{near-infrared band reflectance values (835 to 885 nm), in the range 0 to 255}
#'}
#' @source \url{https://www.bluesky-world.com/}
"aerial_image"
#' A raster of basins numbered as integers
#'
#' A raster object of basins numbered as integers for the area bounded by 160000, 181400,
#' 11300, 30000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system
#' (CRS: +init=epsg:27700) as returned by [basindelin()].
#'
#' @format A raster object with 187 rows and 214 columns.
"basins100m"
#' A 0.05º resolution dataset of hourly fractional cloud cover
#'
#' A dataset containing hourly fractional cloud cover values in 2010 for the area bounded by
#' -5.40, -5.00, 49.90, 50.15 (xmin, xmax, ymin, ymax) (CRS: +init=epsg:4326).
#'
#' @format An array with 5 rows, 8 columns and 8670 hourly values
#'
#' @source \url{http://www.cmsaf.eu/}
"cfc"
#' A 0.05º resolution dataset of hourly diffuse radiation
#'
#' A dataset containing hourly diffuse radiation values in 2010 (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})
#' for the area bounded by -5.40, -5.00, 49.90, 50.15 (xmin, xmax, ymin, ymax) (CRS: +init=epsg:4326).
#'
#' @format An array with 5 rows, 8 columns and 8670 hourly values
#'
#' @source \url{http://www.cmsaf.eu/}
"difrad"
#' A 0.05º resolution dataset of hourly direct radiation normal to the direction of the solar beam
#'
#' A dataset containing hourly direct radiation values in 2010 (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})
#' for the area bounded by -5.40, -5.00, 49.90, 50.15 (xmin, xmax, ymin, ymax) (CRS: +init=epsg:4326).
#'
#' @format An array with 5 rows, 8 columns and 8670 hourly values
#'
#' @source \url{http://www.cmsaf.eu/}
"dnirad"
#' A 1 m resolution raster object of elevation for part of the Lizard Peninsula, Cornwall, UK.
#'
#' A raster object containing elevation in metres with sea coded as NA for the area bounded by
#' 169000, 170000, 12000, 13000  (xmin, xmax, ymin, ymax) using the Ordance
#' Survey GB Grid Reference system (CRS: +init=epsg:27700).
#'
#' @format A raster object with 1000 rows and 1000 columns.

#' @source \url{http://www.tellusgb.ac.uk/}
"dtm1m"
#' A 100 m resolution raster object of elevation for the Lizard Peninsula, Cornwall, UK.
#'
#' A raster object containing elevation in metres with sea coded as NA for the area bounded by
#' 160000, 181400, 11300, 30000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid
#' Reference system (CRS: +init=epsg:27700).
#'
#' @format A raster object with 187 rows and 214 columns.

#' @source \url{http://www.tellusgb.ac.uk/}
"dtm100m"
#' A 1 km resolution raster object of elevation for the Lizard Peninsula, Cornwall, UK.
#'
#' A raster object containing elevation in metres with sea coded as NA for the area bounded by
#' 160000, 182000, 11000, 30000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid
#' Reference system (CRS: +init=epsg:27700).
#'
#' @format A raster object with 19 rows and 22 columns.

#' @source \url{http://www.tellusgb.ac.uk/}
"dtm1km"
#' A 1 km resolution dataset of diurnal temperature ranges
#'
#' A spatially interpolated dataset containing diurnal temperature ranges (ºC) in 2010 for the
#' area bounded by 160000, 182000, 11000, 30000  (xmin, xmax, ymin, ymax) using the Ordance
#' Survey GB Grid Reference system (CRS: +init=epsg:27700).
#'
#' @format An array with 19 rows, 22 columns and 365 daily values
#'
#' @source \url{https://www.metoffice.gov.uk/}
"dtr"
#' A 1 km resolution dataset of interpolated daily specific humidity values
#'
#' A dataset containing daily specific humidity values (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}})
#' in 2010 for the area bounded by 160000, 182000, 11000, 30000  (xmin, xmax, ymin, ymax)
#' using the Ordance Survey GB Grid Reference system (CRS: +init=epsg:27700).
#'
#' @format An array with 19 rows, 22 columns and 365 daily values.
#'
#' @source \url{https://eip.ceh.ac.uk/chess/}
"huss"
#' Inverse-distance weighted land-sea ratios in 36 directions
#'
#' A 100m resolution dataset of nverse-distance weighted land-sea ratios, as
#' produced by function [invls()] in each of 36 directions (0º, 10º..350º) for
#' the area bounded by 160000, 181400, 11300, 30000  (xmin, xmax, ymin, ymax)
#' using the Ordance Survey GB Grid Reference system (CRS: +init=epsg:27700).
#'
#' @format An array with 187 rows, 214 columns and 36 directional values.
#'
"landsearatios"
#' 2010 Data for fitting mesoclimate model.
#'
#' A dataset containing the hourly temperature logger data and net
#' radiation and wind data at each logger location in May 2010.
#'
#' @format A data frame with 65772 rows and 7 variables:
#' \describe{
#'   \item{obs_time}{The time of the logger recording, in DD/MM/YY HH:MM:SS format}
#'   \item{temperature}{The mean hourly logger temperature reading (ºC)}
#'   \item{reftemp}{The predicted reference temperature at the logger location as derived from coarse-scale data after adjusted for elevation, cold air drainage and coastal effects (ºC)}
#'   \item{wind}{The predicted wind speed at the logger locations (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{m s^{-1}}})}
#'   \item{netrad}{The predicted net radiation at the logger locations (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})}
#' }
"mesofitdata"
#' May 2010 Data for fitting microclimate model.
#'
#' A dataset containing the hourly temperature logger data and net
#' radiation and wind data at each logger location  in May 2010.
#'
#' @format A data frame with 11761 rows and 5 variables:
#' \describe{
#'   \item{obs_time}{The time of the logger recording, in DD/MM/YY HH:MM:SS format}
#'   \item{temperature}{The mean hourly logger temperature reading (ºC)}
#'   \item{reftemp}{The predicted reference temperature at the logger location as output by the mesoclimate model (ºC)}
#'   \item{wind}{The predicted wind speed at the logger locations (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{m s^{-1}}})}
#'   \item{netrad}{The predicted net radiation at the logger locations (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})}
#' }
"microfitdata"
#' Climate variables for  May 2010.
#'
#' A dataset containing hourly coarse-resolution climate variables in May 2010.
#'
#' @format A data frame with 744 rows and 9 variables:
#' \describe{
#'   \item{dni}{Direct Normal Radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})}
#'   \item{dif}{Diffuse Radiation (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})}
#'   \item{humidity}{Specific humidity (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}})}
#'   \item{pressure}{Sea-level pressure (Pa)}
#'   \item{cloudcover}{Fractional cloud cover}
#'   \item{year}{Year (AD)}
#'   \item{month}{Numeric month}
#'   \item{day}{Day of month}
#'   \item{hour}{Hour of day}
#' }
"microvars"
#' A ~500 m resolution dataset of MODIS-derived albedo values
#'
#' A dataset ground surface albedo for the area bounded by
#' 169000, 170000, 12000, 13000  (xmin, xmax, ymin, ymax) using the Ordance
#' Survey GB Grid Reference system (CRS: +init=epsg:27700) derived from
#' Moderate Resolution Imaging Spectroradiometer (MODIS) imagery
#'
#' @format A matrix with 18 rows and 36 columns.

#' @source \url{https://lpdaac.usgs.gov/}
"modis"
#' A one metre resolution matrix of net longwave radiation
#'
#' A dataset containing net longwave radiation values (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^-2 hr^{-1}}})
#' for 2010-05-24 11:00 across the area bounded by 169000, 170000, 12000, 13000
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: +init=epsg:27700), as produced by [longwaveveg()]
#'
#' @format A matrix with 1000 rows and 1000 columns.
"netlong1m"
#' A 100 m resolution matrix of net longwave radiation
#'
#' A dataset containing net longwave radiation values (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})
#' for 2010-05-24 11:00 GMT across the area bounded by 160000, 181400, 11300, 30000
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: +init=epsg:27700), , as produced by [longwavetopo()]
#'
#' @format A matrix with 187 rows and 214 columns.
"netlong100m"
#' A one metre resolution matrix of net shortwave radiation
#'
#' A dataset containing net shortwave radiation values (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})
#' for 2010-05-24 11:00 GMT across the area bounded by 169000, 170000, 12000, 13000
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: +init=epsg:27700), as produced by [shortwaveveg()]
#'
#' @format A matrix with 1000 rows and 1000 columns.
"netshort1m"
#' A 100 metre resolution matrix of net shortwave radiation
#'
#' A dataset containing net shortwave radiation values (\ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}})
#' for 2010-05-24 11:00 GMT across the area bounded by 160000, 181400, 11300, 30000
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: +init=epsg:27700), as produced by [shortwavetopo()]
#'
#' @format A matrix with 187 rows and 214 columns.
"netshort100m"
#' A 1 km resolution dataset of interpolated daily sea-level pressure values
#'
#' A dataset containing daily sea-level pressure values (Pa) in 2010 for the area bounded by
#' 160000, 182000, 11000, 30000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid
#' Reference system (CRS: +init=epsg:27700).
#'
#' @format An array with 19 rows, 22 columns and 365 daily values
#'
#' @source \url{https://eip.ceh.ac.uk/chess/}
"pres"
#' A 1km resolution dataset of daily sea-level temperature
#'
#' A dataset containing daily sea-level temperature (ºC) in 2010 for the area bounded by
#' 160000, 182000, 11000, 30000  (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid
#' Reference system (CRS: +init=epsg:27700).
#'
#' @format An array with 19 rows, 22 columns and 365 daily values
#'
#' @source \url{https://eip.ceh.ac.uk/chess/}
"tas"
#' A 100 m resolution array of hourly reference temperatures.
#'
#' A 100 m resolution three-dimensional array of hourly reference temperatures
#' (ºC) for the area bounded by 169000, 170000, 12000, 13000  (xmin, xmax, ymin, ymax)
#' using the Ordance Survey GB Grid Reference system (CRS: +init=epsg:27700),
#' for May 2010 as output by [runmicro()].
#' @format An array with 10 rows, 10 columns and 744 hourly values.
"temp100"
#' A 1 m resolution raster object of vegetation height.
#'
#' A dataset containing the vegetation height (m) for the area bounded by
#' 169000, 170000, 12000, 13000  (xmin, xmax, ymin, ymax) using the Ordance
#' Survey GB Grid Reference system (CRS: +init=epsg:27700). Data were derived
#' from a digital terrain and digital surface model.
#'
#' @format A matrix with 1000 rows and 1000 columns.

#' @source \url{http://www.tellusgb.ac.uk/}
"veg_hgt"
#' A one metre resolution matrix of windspeed
#'
#' A dataset containing wind speed values (\ifelse{html}{\out{m s<sup>-1</sup>}}{\eqn{m s^{-1}}})
#' for 2010-05-24 11:00 GMT across the area bounded by 169000, 170000, 12000, 13000
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: +init=epsg:27700), as produced with [windcoef()]
#'
#' @format A matrix with 1000 rows and 1000 columns.
"wind1m"
#' Six-hourly wind speed at 10 m for 2010.
#'
#' A dataset containing six-hourly wind speed in 2010 estimated 10 m above the ground.
#'
#' @format A data frame with 1460 rows and 2 variables:
#' \describe{
#'   \item{obs_time}{Time at date in yyyy-mm-dd HH:MM format}
#'   \item{wind10m}{Wind speed 10m above the ground (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{m s^{-1}}})}
#' }
#' @source \url{http://www.ncep.noaa.gov/}
"wind2010"

#' Orders drainage basins by elevation
#'
#' @description
#' `basinsort` is an internal function used by [basindelin()] and [basinmerge()] to number
#' drainage basins sequentially from lowest elevation (of lowest point) to highest.
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#' @param basins a raster object, two-dimensional array or matrix of basins numbered as integers.
#'
#' @return a raster object, two-dimensional array or matrix of sequentially numbered basins
#' @import raster
#' @importFrom dplyr left_join
#' @export
#' @keywords internal
#' @examples
#' basins <- matrix(c(1:4), nrow = 2, ncol = 2)
#' dem <- matrix(c(4:1), nrow = 2, ncol = 2)
#' basinsort(dem, basins)
basinsort <- function(dem, basins) {
  mdf <- function(u) {
    sel <- which(basins == u)
    min(dem[sel], na.rm = TRUE)
  }
  r <- dem
  dem <- is_raster(dem)
  basins <- is_raster(basins)
  u <- unique(as.vector(basins))
  u <- u[is.na(u) == F]
  u2 <- c(1:length(u))
  df1 <- data.frame(old = as.vector(basins))
  df2 <- data.frame(old = c(NA, u), new = c(NA, u2))
  df3 <- left_join(df1, df2, by = "old")
  basins <- array(df3$new, dim = dim(basins))
  u <- unique(as.vector(basins))
  u <- u[is.na(u) == F]
  mnd <- sapply(u, mdf)
  o <- order(mnd)
  u2 <- u[o]
  df1 <- data.frame(old = as.vector(basins))
  df2 <- data.frame(old = c(NA, u), new = c(NA, u2))
  df3 <- left_join(df1, df2, by = "old")
  bm2 <- array(df3$new, dim = dim(basins))
  if_raster(bm2, r)
}
#' Internal funcion to calculate whether neighbouring cell is higher or lower
#' @export
.hgttongbr <- function(m) {
  hgt <- rep(m[2, 2], 8)
  neighbours <- c(m[1, 2], m[1, 3], m[2, 3], m[3, 3], m[3, 2], m[3, 1],
                  m[2, 1], m[1, 1])
  htn <- ifelse(hgt > neighbours, 1, 0)
  intcode <- sum(2 ^ (which(rev(unlist(strsplit(as.character(htn), "")) ==
                                  1)) - 1))
  intcode
}
#' Internal funcion to convert integer to binary
#' @import stringr
#' @export
.integertobinary8 <- function(i) {
  a <- 2 ^ (0:9)
  b <- 2 * a
  binc <- format(sapply(i, function(x) sum(10 ^ (0:9)[(x %% b) >= a])),
                 scientific = FALSE)
  if (nchar(binc) > 8) warning("Integer > 8 bit binary")
  binc <- str_pad(binc, 8, pad = "0")
  binc
}
#' Internal funcion to check whether neighbouring cell is higher or lower
#' @export
.updown <- function(dem) {
  m <- dem
  m2 <- array(9999, dim = c(dim(dem)[1] + 2, dim(dem)[2] + 2))
  m2[2:(dim(dem)[1] + 1), 2:(dim(dem)[2] + 1)] <- m
  updownm <- matrix(rep(NA, length(dem)), nrow = dim(dem)[1])
  for (y in 2:(dim(m2)[2] - 1)) {
    for (x in 2:(dim(m2)[1] - 1)) {
      focalcell <- m2[x, y]
      if (is.na(focalcell) == FALSE) {
        m9 <- matrix(c(m2[x - 1, y - 1], m2[x, y - 1], m2[x + 1, y - 1],
                       m2[x - 1, y], focalcell, m2[x + 1, y],
                       m2[x - 1, y + 1], m2[x, y + 1], m2[x + 1, y + 1]),
                     nrow = 3)
        updownm[(x - 1), (y - 1)] <- .hgttongbr(m9)
      }
    }
  }
  updownm
}
#' Internal function to merge single basin
#' @export
.onebasin_merge <- function(md, mb, b) {
  mb2 <- array(NA, dim = c(dim(mb)[1] + 2, dim(mb)[2] + 2))
  mb2[2:(dim(mb)[1] + 1), 2:(dim(mb)[2] + 1)] <- mb
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  bm <-  b
  sel <- which(mb == b)
  x <- arrayInd(sel, dim(md))[, 1]
  y <- arrayInd(sel, dim(md))[, 2]
  for (i in 1:length(x)) {
    mb9 <- mb2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
    md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
    sel2 <- which(mb9 > b)
    if (length(sel2) > 0) {
      for (j in 1:length(sel2)) {
        xi <- arrayInd(sel2[j], dim(mb9))[, 1]
        yi <- arrayInd(sel2[j], dim(mb9))[, 2]
        x2 <- c(xi - 1, xi, xi + 1)
        y2 <- c(yi - 1, yi, yi + 1)
        x2 <- x2[x2 > 0 & x2 < 4]
        y2 <- y2[y2 > 0 & y2 < 4]
        mb3 <- mb9[x2, y2]
        md3 <- md9[x2, y2]
        md3[mb3 == b] <- NA
        u <- unique(mb3[mb3 > b])
        u <- u[is.na(u) == F]
        for (k in 1:length(u)) {
          vrs <- md3[mb3  == u[k]]
          if (max(vrs, na.rm = T) > md9[2, 2]) bm <- c(bm, u[k])
        }
      }
    }
  }
  unique(bm)
}
#' Internal function used by basinmerge
#' @export
.basinchars <- function(md, mb) {
  mb2 <- array(NA, dim = c(dim(mb)[1] + 2, dim(mb)[2] + 2))
  mb2[2:(dim(mb)[1] + 1), 2:(dim(mb)[2] + 1)] <- mb
  sel <- which(is.na(mb2) == T)
  mb2[sel] <- (-999)
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  md2[sel] <- 9999
  dfm <- data.frame(basin = NA, basinmindepth = NA,
                    pourpointbasin = NA, pourpointhgt = NA)
  for (b in 1:max(mb, na.rm = T)) {
    sel <- which(mb == b)
    x <- arrayInd(sel, dim(md))[, 1]
    y <- arrayInd(sel, dim(md))[, 2]
    df1 <- data.frame(basin = b,
                      basinmindepth = min(md[sel], na.rm = T),
                      pourpointbasin = NA, pourpointhgt = 0)
    b2 <- 0
    pph <- 0
    for (i in 1:length(x)) {
      mb9 <- mb2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
      md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
      sel <- which(as.vector(mb9) != b)
      b2[i] <- NA
      pph[i] <- NA
      if (length(sel) > 0) {
        mb9 <- mb9[sel]
        md9 <- md9[sel]
        sel2 <- which(as.vector(md9) == min(as.vector(md9)))
        b2[i] <- mb9[sel2][1]
        pph[i] <- md9[sel2][1]
      }
    }
    sel <- which(pph == min(pph, na.rm = T))
    df1$pourpointbasin <- b2[sel[1]]
    df1$pourpointhgt <- min(pph, na.rm = T)
    dfm <- rbind(dfm, df1)
  }
  dfm <- dfm[which(is.na(dfm$basin) == F), ]
  dfm
}
#' Internal function used by basinmerge
#' @export
.bvarsrem <- function(bvars, boundary) {
  sel <- which(bvars$pourpointbasin != -999)
  if (length(sel) == 0) warning("all basins flow into sea")
  bvars <- bvars[sel, ]
  bvars$basindepth <- bvars$pourpointhgt - bvars$basinmindepth
  bvars <- bvars[bvars$basindepth < boundary, ]
  o <- order(bvars$pourpointhgt, decreasing = TRUE)
  bvars <- bvars[o, ]
  bvars
}
#' Internal function used to calculate flow direction
#' @export
.flowdir <- function(md) {
  fd <- md * 0
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  v <- c(1:length(md))
  v <- v[is.na(md) == F]
  x <- arrayInd(v, dim(md))[, 1]
  y <- arrayInd(v, dim(md))[, 2]
  for (i in 1:length(x)) {
    md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
    fd[x[i], y[i]] <- round(mean(which(md9 == min(md9, na.rm = T))), 0)
  }
  fd
}
#' Delineates hydrological basins
#'
#' @description `basindelin` uses digital elevation data to delineate hydrological or cold-air drainage basins.
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#'
#' @return a raster object, two-dimensional array or matrix  with individual basins numbered sequentially as integers.
#' @import raster
#' @export
#'
#' @seealso [basindelin_big()] for working with large datasets.
#'
#' @details
#' If `dem` is a raster object, a raster onbject is returned.
#' This function is used to delineate cold-air drainage basins.
#' It iteratively identifies the lowest elevation pixel of `dem`, and
#' assigns any of the eight adjoining pixels to the same basin if higher.
#' The process is repeated on all assigned adjoining pixels until no further
#' higher pixels are found. The next lowest unassigned pixel is then
#' identified, a new basin identity assigned and the processes repeated until
#' all pixels are assigned to a basin. Relative to heuristic algorithms, it is
#' slow and run time increases exponentially with size of `dtm`. However, in
#' contrast to many such algorithms, all basins are correctly seperated by
#' boundaries >0. With large datasets, with > 160,000 pixels, calculations will
#' be slow and [basindelin_big()] should be used instead
#'
#' @examples
#' library(raster)
#' dem <- aggregate(dtm1m, 20)
#' basins <- basindelin (dem)
#' plot(basins, main = "Basins")
basindelin <- function(dem) {
  r <- dem
  dem <- is_raster(dem)
  ngbrow <- c(-1, -1, 0, 1, 1, 1, 0, -1)
  ngbcol <- c(0, 1, 1, 1, 0, -1, -1, -1)
  updownm <- .updown(dem)
  basinsm <- ifelse(is.na(dem), NA, 0)
  donem <- dem
  basincell <- order(dem)[1]
  basin <- 1
  while (basincell != -999) {
    while (basincell != -999) {
      i <- basincell
      irow <- arrayInd(i, dim(dem))[1]
      icol <- arrayInd(i, dim(dem))[2]
      basinsm[irow, icol] <- basin
      neighbours <- .integertobinary8(updownm[i])
      for (n in 1:8) {
        if ((irow + ngbrow[n]) > 0 & (irow + ngbrow[n]) <= dim(basinsm)[1] &
            (icol + ngbcol[n]) > 0 & (icol + ngbcol[n]) <= dim(basinsm)[2]) {
          if (!is.na(basinsm[(irow + ngbrow[n]), (icol + ngbcol[n])])) {
            if (substr(neighbours, n, n) == "0" &
                basinsm[(irow + ngbrow[n]), (icol + ngbcol[n])] == 0) {
              basinsm[(irow + ngbrow[n]), (icol + ngbcol[n])] <- basin
            }
          }
        }
      }
      donem[i] <- NA
      if (length(which(basinsm == basin & !is.na(donem))) > 0) {
        basincell <- min(which(basinsm == basin & !is.na(donem) ))
      }
      else basincell <- -999
    }
    if (length(which(!is.na(donem)) > 0)) {
      basincell <- which(donem == min(donem, na.rm = TRUE))[1]
    }
    else basincell <- -999
    basin <- basin + 1
  }
  basinsm <- if_raster(basinsm, r)
  basinsort(r, basinsm)
}
#' Delineates hydrological basins for large datasets
#'
#' @description
#' `basindelin_big` is for use with large digital elevation datasets, to
#' delineate hydrological or cold-air drainage basins.
#'
#' @param dem a raster object of elevations.
#' @param dirout an optional character vector containing a single path directory for temporarily storing tiles. Deleted after use. Tilde expansion (see [path.expand()]) is done.
#' @param trace a logical value indicating whether to plot and report on progress.
#'
#' @return a raster object with individual basins numbered sequentially as integers.
#' @import raster
#' @export
#' @seealso [basindelin()] for working with smaller datasets.
#'
#' @details
#' The function `basindelin_big` divides the large dataset into tiles and then
#' uses [basindelin()] to delineate basins for each tile before mosaicing back
#' together and merging basins along tile edges if not seperated by a boundary
#' > 0. If `dirout` is unspecified, then a directory `basinsout` is
#' temporarily created within the working directory. If `trace` is TRUE (the
#' default) then progress is tracked during three stages: (1) the basins
#' of each tile are plotted, (2) basins after mosaicing, but prior
#' to merging are plotted and (3) on each merge iteration, the number of basins
#' to merge is printed and processed basin is plotted.
#'
#' @examples
#' library(raster)
#' basins <- basindelin_big(dtm1m)
#' plot(basins, main = "Basins")
basindelin_big <- function(dem, dirout = NA, trace = TRUE) {
  dem <- trim(dem)
  dmsx <- ceiling(dim(dem)[2] / 200) - 1
  dmsy <- ceiling(dim(dem)[1] / 200) - 1
  if (dmsx < 1 & dmsy < 1) {
    basins <- basindelin(dem)
  }
  else {
    xres <- xres(dem)
    yres <- yres(dem)
    ed <- extent(dem)
    if (is.na(dirout)) dirout <- "basinsout/"
    dir.create(dirout)
    fol <- ""
    ii <- 1
    for (i in 0:dmsx) {
      for (j in 0:dmsy) {
        xmn <- ed@xmin + i * 200 * xres
        xmx <- min((xmn + 200 * xres), ed@xmax)
        ymn <- ed@ymin + j * 200 * yres
        ymx <- min((ymn + 200 * yres), ed@ymax)
        e <- extent(c(xmn, xmx, ymn, ymx))
        r <- crop(dem, e)
        v <- getValues(r)
        if (is.na(mean(v, na.rm = TRUE)) == FALSE) {
          b <- basindelin(r)
          if (trace) {
            progress <- paste0(round(ii / ((dmsx + 1) * (dmsy + 1)) * 100, 1),
                               "%")
            plot(b, main = paste0("basin slice progress: ", progress))
          }
          fo <- paste0(dirout, "b", i, "_", j, ".tif")
          fol <- c(fol, fo)
          writeRaster(b, filename = fo, overwrite = TRUE)
        }
        ii <- ii + 1
      }
    }
    fol <- fol[fol != ""]
    basins <- raster(fol[1])
    ta <- max(getValues(basins), na.rm = TRUE)
    if (length(fol) > 1) {
      for (i in 2:length(fol)) {
        r <- raster(fol[i]) + ta
        basins <- mosaic(basins, r, fun = mean)
        ta <- max(getValues(basins), na.rm = TRUE)
      }
    }
    if (trace) {
      plot(basins, main = "Basin mosaic complete")
    }
    unlink(dirout)
    iter <- 1
    test <- 0
    basins[is.na(basins)] <- 9999
    dem[is.na(dem)] <- 9999
    while (test != 1) {
      basins <- basinsort(dem, basins)
      bmm <- getValues(basins, format = "matrix")
      lst <- as.list(c(1:max(getValues(basins), na.rm = TRUE)))
      ii <- 1
      if (dmsx > 0) {
        for (i in 1:dmsx) {
          xmn <- ed@xmin + i * 200 * xres - 1
          xmx <- min((xmn + 2 * xres), ed@xmax)
          e <- extent(c(xmn, xmx, ed@ymin, ed@ymax))
          r <- crop(basins, e)
          ds <- crop(dem, e)
          md <- getValues(ds, format = "matrix")
          mb <- getValues(r, format = "matrix")
          u <- unique(getValues(r))
          u <- u[is.na(u) == F]
          u <- u[order(u)]
          if (length(u) > 0) {
            for (k in 1:length(u)) {
              lst[[ii]] <- .onebasin_merge(md, mb, u[k])
              ii <- ii + 1
            }
          }
        }
      }
      if (dmsy > 0) {
        for (j in 1:dmsy) {
          ymn <- ed@ymin + j * 200 * yres - 1
          ymx <- ymn + 2 * yres
          e <- extent(c(ed@xmin, ed@xmax, ymn, ymx))
          r <- crop(basins, e)
          ds <- crop(dem, e)
          md <- getValues(ds, format = "matrix")
          mb <- getValues(r, format = "matrix")
          u <- unique(getValues(r))
          u <- u[is.na(u) == F]
          u <- u[order(u)]
          if (length(u) > 0) {
            for (k in 1:length(u)) {
              lst[[ii]] <- .onebasin_merge(md, mb, u[k])
              ii <- ii + 1
            }
          }
        }
      }
      lst <- lst[1:ii]
      ub <- unique(unlist(lst))
      ub <- ub[order(ub)]
      ub2 <- ub
      for (i in 1:length(ub)) {
        for (j in 1:ii) {
          v <- lst[[j]]
          tst <- which(v == ub[i])
          if (length(tst) > 0) ub2[i] <- min(ub2[i], min(v))
        }
      }
      sel <- which(ub != ub2)
      if (trace) {
        print(paste0("Basins to merge: ", length(sel)))
      }
      if (length(sel) == 0) {
        test <- 1
        if (trace) plot(basins, main = "Merge complete")
      }
      else {
        mbout <- bmm
        for (i in 1:length(sel)) {
          mbout[bmm == ub[sel[i]]] <- ub2[sel[i]]
        }
        basins <- raster(mbout, template = basins)
        if (trace) {
          plot(basins, main = paste0("Basin merge iteration: ", iter))
        }
        iter <- iter + 1
      }
    }
  }
  basins[dem == 9999] <- NA
  if (trace) plot (basins, main = "Merge complete")
  basins <- basinsort(dem, basins)
  basins
}
#' Merges adjoining basins
#'
#' @description
#' `basinmerge` merges adjoining basins if the height differences between the
#' bottom of the basin and the pour point is less than than that
#' specified by `boundary`.
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#' @param basins a raster object, two-dimensional array or matrix with basins numbered as integers as returned by [basindelin()].
#' @param boundary a single numeric value. Basins seperated by boundaries below this height are merged (should have same units as `dtm`).
#'
#' @return a raster object, two-dimensional array or matrix with basins numbered as integers.
#' @import raster
#' @export
#'
#' @details
#' If `dem` is a raster object, then a raster object is returned.
#' If the differences in height between the pour-point and bottom of the basin is
#' less than that specified by `boundary` the basin is merged with basin to which
#' water or air would pour.
#'
#' @examples
#' library(raster)
#' basins2 <- basinmerge(dtm100m, basins100m, 1)
#' par(mfrow=c(1, 2))
#' plot(basins100m, main = "Basins")
#' plot(basins2, main = "Merged basins")
basinmerge <- function(dem, basins, boundary) {
  if (all.equal(dim(basins)[1:2], dim(dem)[1:2]) == FALSE)
    stop ("basins and dem have different dimensions")
  r <- basins
  dem <- is_raster(dem)
  basins <- is_raster(basins)
  test <- F
  while (test == F) {
    mb2 <- basins
    bvars <- .basinchars(dem, basins)
    bkeep <- .bvarsrem(bvars, boundary)
    if (dim(bkeep)[1] == 0) test <- T
    for (b in 1:(dim(bkeep)[1])) {
      sel <- which(bkeep$pourpointbasin == bkeep$basin[b])
      if (length(sel) > 0) {
        bkeep$pourpointbasin[sel] <- bkeep$pourpointbasin[b]
      }
      sel <- which(basins == bkeep$basin[b])
      mb2[sel] <- bkeep$pourpointbasin[b]
    }
    u <- unique(as.vector(mb2))
    sel <- which(is.na(u) == F)
    u <- u[sel]
    for (i in 1:length(u)) {
      sel <- which(mb2 == u[i])
      basins[sel] <- i
    }
  }
  basins <- basinsort(dem, basins)
  if_raster(basins, r)
}
#' Calculates accumulated flow
#'
#' @description
#' `flowacc` is used by [pcad()] to calculate accumulated flow to each cold air drainage
#' basin
#'
#' @param dem a raster object, two-dimensional array or matrix of elevations.
#'
#' @return a raster object, two-dimensional array or matrix of accumulated flow.
#' @details Accumulated flow is expressed in terms of number of cells.
#' @export
#'
#' @examples
#' library(raster)
#' fa <- flowacc(dtm100m)
#' plot(fa, main = 'Accumulated flow')
flowacc <- function (dem)
{
  dm <- is_raster(dem)
  fd <- .flowdir(dm)
  fa <- fd * 0 + 1
  o <- order(dm, decreasing = T, na.last = NA)
  for (i in 1:length(o)) {
    x <- arrayInd(o[i], dim(dm))[1]
    y <- arrayInd(o[i], dim(dm))[2]
    f <- fd[x, y]
    x2 <- x + (f - 1)%%3 - 1
    y2 <- y + (f - 1)%/%3 - 1
    if (x2 > 0 & x2 < dim(dm)[1] & y2 > 0 & y2 < dim(dm)[2])
      fa[x2, y2] <- fa[x, y] + 1
  }
  if_raster(fa, dem)
}
#' Calculates whether conditions are right for cold air drainage
#'
#' @description
#' `cadconditions` determines whether wind speed, humidity and cloud cover  are such that
#' cold air drainage is likely to occur
#'
#' @param h a vector of hourly specific humidities (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc a single numeric value, raster object, two-dimensional array or matrix of temperatures (ºC).
#' @param n a vector of hourly fractional cloud cover values (range 0 - 1).
#' @param p an optional vector of hourly atmospheric pressure values (Pa).
#' @param wind a vector of wind speed values at one metre above the ground (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{m s^{-1}}})
#' @param startjul julian day of first observation as returned by [julday()]
#' @param lat latitude of the location for which cold air drainage conditions are required (decimal degrees, -ve south of equator).
#' @param long longitude of the location for which cold air drainage conditions are required (decimal degrees, -ve west of Greenwich meridian).
#' @param starttime the hour of the first observation (decimal, 0-23).
#' @param hourint the interval (in hours) between successive observations
#' @param windthresh an optional threshold value of wind speed below which cold air conditions can occur (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{m s^{-1}}})
#' @param emthresh an optional threshold value of emissivity below which cold air conditions can occur (range 0 - 1)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param con an optional logical value indicating whether or not to allow cold air drainage conditions to occur only if conditions are right for three or more consecutive hours. Ignored if hourint != 1.
#'
#' @return a vector of binary values indicating whether cold air drainage conditions occur (1) or not (0)
#' @export
#'
#' @details
#' `cadconditions` uses a time series of wind and emissivity data to determine
#' whether cold air drainage conditions are likely to occur and returns a binary
#' vector of the same length as `em` and `wind` indicating whether conditions
#' occur (1) or not (0). They are assumed to occur at night or within three
#' hours of dawn only and when both `em` and `wind` are below the values
#' specified by `windthresh` and `emthresh`. If no start time is specified it
#' is assumed that first index of `em` and `wind` occurs at midnight on the
#' date specified by `startjul`. If no `hourint` is provided, the time
#' interval between indices of `em` and `wind` are assumed to be hourly.
#' If `con` is TRUE  and `hourint` is 1 (the default), cold air drainage conditions are
#' assumed to persist only when conditions are right for three consecutive hours or more.
#' The first index of the output, for which prior conditions cannot be assessed
#' is set to 1 if conditions are right, irrespective of prior conditions. The
#' second index is set to one only if conditions are right in both that hour
#' and the preceding hour. If `con` is FALSE or `hourint != 1` prior conditions
#' are ignored.
#'
#' @examples
#' # ===============================================
#' # Mean daily climate for Lizard, Cornwall in 2010
#' # ===============================================
#' h <- apply(huss, 3, mean, na.rm = TRUE)
#' p <- apply(pres, 3, mean, na.rm = TRUE)
#' tmin <- apply((tas - dtr), 3, mean, na.rm = TRUE)[2:364]
#' tmax <- apply((tas + dtr), 3, mean, na.rm = TRUE)[2:364]
#' # =====================================
#' # hourly climate 2nd Jan to 30 Dec 2010
#' # =====================================
#' h <- spline(h, n = 8737)$y[13:8724]
#' p <- spline(p, n = 8737)$y[13:8724]
#' n <- apply(cfc[,,13:8724], 3, mean)
#' rdni <- apply(dnirad[,,13:8724], 3, mean)
#' rdif <- apply(difrad[,,13:8724], 3, mean)
#' jd <- julday(2010, 2, 1)
#' jd <- c(jd:(jd+362))
#' tc <- hourlytemp(jd, h, n, p, rdni, rdif, tmin, tmax, 50.05, -5.19)
#' ws10m <- spline(wind2010$wind10m, n = 8755)$y[25:8736]
#' ws1m <- windheight(ws10m, 10, 1)
#' # =================================================================
#' # Calculate whether cold air drainage, persists and plot proportion
#' # =================================================================
#' startjul <- julday(2010,1,2)
#' hist(cadconditions(h, tc, n, p, ws1m, startjul, 50.05, -5.19),
#'      main = "", xlab = "Cold air drainage conditions (1 = Y)")
cadconditions <- function(h, tc, n, p = 100346.13, wind, startjul, lat, long,
                          starttime = 0, hourint = 1, windthresh = 4.5, emthresh = 0.5,
                          merid = round(long / 15, 0) * 15, dst = 0, con = TRUE) {
  pk <- p / 1000
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ws <- 0.622 * e0 / pk
  rh <- (h / ws) * 100
  rh[rh > 100] <- 100
  ea <- e0 * rh / 100
  emcs <- 0.23 + 0.433 * (ea / (tc + 273.15)) ^ (1 / 8)
  em <- emcs * (1 - n ^ 2) + 0.976 * n ^ 2
  jd <- floor(c(1:length(em)) * hourint / 24 - hourint / 24 + startjul +
                starttime / 24)
  st <- suntimes(jd, lat, long, merid, dst)
  hrs <- (c(1:length(em)) * hourint - hourint + starttime)%%24
  dn <- ifelse(hrs > (st$sunrise + 3) & hrs < st$sunset, 1, 0)
  cad <- ifelse(dn < 1 & wind < windthresh & em < emthresh, 1, 0)
  if (hourint == 1 & con) {
    cad[2] <- ifelse(cad[1] & cad[2] == 1, 1, 0)
    cad <- c(cad[1:2], cad[3:length(cad)] * cad[2:(length(cad) - 1)] *
             cad[1:(length(cad) - 2)])
  }
  cad
}
#' Calculates cold air drainage potential
#'
#' @description
#' ` pcad` calculates the expected temperature differences resulting from cold air
#' drainage.
#'
#' @param dem a raster object of elevation (m).
#' @param basins a raster object, two-dimensional array or matrix with basins numbered as integers as returned by [basindelin()].
#' @param fa a raster object of accumulated flow, as returned by [flowacc()]
#' @param tc a single numeric value, raster object, two-dimensional array or matrix of values with the dimensions as `dem` of temperature (ºC).
#' @param h a single numeric value, raster object, two-dimensional array or matrix of values with the dimensions as `dem` of specific humidity (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}}).
#' @param p an optional single numeric value, raster object, two-dimensional array or matrix of values with the dimensions as `dem` of atmospheric pressure (Pa).
#' @param out specifies the type of output to provide (see details). Possible values are "cadp", "tempdif" and "pflow".
#' @return a raster object:
#'  \describe{
#'   \item{if `out = "tempdif"`}{the expected temperature difference (ºC).}
#'   \item{if `out = "pflow"`}{the accumulated flow logairthmically transformed expressed as a proportion of the maximum in each basin.}
#'   \item{if `out = "cadp"` (the default)}{the expected temperature difference x the proportion of accumulated flow.}
#' }
#' @import raster
#' @export
#'
#' @details To derive expected temperature differences, `pcad` calculates the
#' difference difference in elevation between each point and the highest in the
#' basin and multiplies this by the lapse rate. Accumulated flow is calculated
#' using [flowacc()], logarithmically transformed and expressed as proportion of
#' the maximum in each basin as an indication of each locations potential to
#' experience cold air drainage. Cold air flow is possible over shallow
#' boundaries, `basins` is best derived using [basinmerge()].
#' Warning: function is quite slow on large datasets.
#'
#' @examples
#' library(raster)
#' basins <- basinmerge(dtm100m, basins100m, 2)
#' h <- humidityconvert(50, intype = "relative", 20)$specific
#' fa <- flowacc(dtm100m, basins)
#' cp1 <- pcad(dtm100m, basins, fa, 20, h)
#' cp2 <- pcad(dtm100m, basins, fa, 20, h, out = "tempdif")
#' cp3 <- pcad(dtm100m, basins, fa, 20, h, out = "pflow")
#' par(mfrow=c(1, 3))
#' plot(cp3, main = "Accumulated flow proportion")
#' plot(cp2, main = "Expected temperature difference")
#' plot(cp1, main = "Cold air drainage potential")
pcad <- function(dem, basins, fa, tc, h, p = 101300, out = "cadp") {
  h <- is_raster(h)
  tc <- is_raster(tc)
  p <- is_raster(p)
  dm <- is_raster(dem)
  bm <- is_raster(basins)
  fm <- is_raster(fa) * xres(dem) * yres(dem)
  pfa <- fm * 0
  td <- fm * 0
  lr <- lapserate(tc, h, p)
  if (is.array(lr) == FALSE) {
    lr <- array(lr, dim = dim(dm))
  }
  for (b in 1:max(bm, na.rm = TRUE)) {
    sel <- which(bm == b)
    fao <- log(fm[sel])
    pfa[sel] <- fao / max(fao, na.rm = TRUE)
    ed <- max(dm[sel], na.rm = TRUE) - dm[sel]
    td[sel] <- ed * lr[sel]
  }
  cdp <- NA
  if (out == "cadp") cdp <- pfa * td
  if (out == "tempdif") cdp <- td
  if (out == "pflow") cdp <- pfa
  if (is.na(max(cdp, na.rm = T))) {
    stop(paste0(out, " not recognized"))
  }
  if_raster(cdp, dem)
}

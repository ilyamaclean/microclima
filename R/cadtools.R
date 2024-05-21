# ============================================================================ #
# ~~~~~~~~~ Basin delineation worker functions here  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' Check if input is a SpatRaster or PackedSpatRaster and convert to matrix or array
#' if it is
#'@noRd
.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-terra::rast(r)
  if (class(r)[1] == "SpatRaster") {
    if (dim(r)[3] == 1) {
      m<-as.matrix(r,wide=TRUE)
    } else m<-as.array(r)
  } else {
    m<-r
  }
  return(m)
}
#' Create SpatRaster object using a template
#' @import terra
#' @export
#' @keywords internal
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' Wrapper for C++ function
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclima, .registration = TRUE
#' @export
.basindelinCpp<-function(dtm) {
  dm<-.is(dtm)
  dm[is.na(dm)]<-9999
  dm2<-array(9999,dim=c(dim(dm)[1]+2,dim(dm)[2]+2))
  dm2[2:(dim(dm)[1]+1),2:(dim(dm)[2]+1)]<-dm
  # (2) create blank basin file
  bsn<-dm2*NA
  dun<-array(0,dim=dim(bsn))
  bsn<-basinCpp(dm2, bsn, dun)
  dd<-dim(bsn)
  bsn<-bsn[2:(dd[1]-1),2:(dd[2]-1)]
  r<-.rast(bsn,dtm)
  return(r)
}
#' Internal function for delineating basins with option for boundary > 0
#'@noRd
.basindelin<-function(dtm, boundary = 0) {
  # Delineate basins
  dm<-dim(dtm)
  me<-mean(as.vector(dtm),na.rm=TRUE)
  if (is.na(me) == FALSE) {
    bsn<-.basindelinCpp(dtm)
    # Merge basins if boundary > 0
    if (boundary > 0) {
      mx<-max(as.vector(bsn),na.rm=T)
      tst<-1
      while (tst == 1) {
        u<-unique(as.vector(bsn))
        u<-u[is.na(u) == F]
        if (length(u) > 1) {
          bsn<-.basinmerge(dtm,bsn,boundary)
          u<-unique(as.vector(bsn))
          u<-u[is.na(u) == F]
          if (length(u) > 1) {
            bsn<-.basinmerge(dtm,bsn,boundary)
            u<-unique(as.vector(bsn))
            u<-u[is.na(u) == F]
          }
          u<-unique(as.vector(bsn))
          u<-u[is.na(u) == F]
        }
        if (length(u) == 1) tst<-0
        mx2<-max(as.vector(bsn),na.rm=T)
        if (mx2 ==  mx) {
          tst<-0
        } else mx<-mx2
      } # end while
    } # end if boundary
  } else bsn<-dtm # end if boundary
  return(bsn)
}
#' function to identify which basin edge cells are less or equal to boundary
#'@noRd
.edge<-function(v) {
  o<-0
  if (is.na(v[1]) == FALSE) {
    if (max(v,na.rm=TRUE) > v[1]) o<-1
  }
  o
}
#' function to assign which surrounding cells should be merged
#'@noRd
.edgec<-function(v) {
  o<-v*0
  if (is.na(v[1]) == FALSE) {
    s<-which(v>v[1])
    o[s]<-1
  }
  o
}
#' function to grab neighbouring cells and reassign basin number
#'@noRd
.asign3<-function(bm2,bea,rw,cl) {
  b3<-bm2[rw:(rw+2),cl:(cl+2)]
  v<-bea[rw,cl,]
  if (is.na(v[2])==FALSE & v[2] > 0) b3[2,1]<-b3[2,2]
  if (is.na(v[3])==FALSE & v[3] > 0)  b3[2,3]<-b3[2,2]
  if (is.na(v[4])==FALSE & v[4] > 0)  b3[1,2]<-b3[2,2]
  if (is.na(v[5])==FALSE & v[5] > 0)  b3[1,1]<-b3[2,2]
  if (is.na(v[6])==FALSE & v[6] > 0)  b3[1,3]<-b3[2,2]
  if (is.na(v[7])==FALSE & v[7] > 0) b3[3,2]<-b3[2,2]
  if (is.na(v[8])==FALSE & v[8] > 0)  b3[3,1]<-b3[2,2]
  if (is.na(v[9])==FALSE & v[9] > 0)  b3[3,3]<-b3[2,2]
  b3
}
#' Merge basins based on specified boundary
#'@noRd
.basinmerge<-function(dtm,bsn,boundary=0.25) {
  # Put buffer around basin and dtn
  bm<-.is(bsn)
  bm2<-array(NA,dim=c(dim(bm)[1]+2,dim(bm)[2]+2))
  bm2[2:(dim(bm)[1]+1),2:(dim(bm)[2]+1)]<-bm
  dm<-.is(dtm)
  dm2<-array(NA,dim=c(dim(dm)[1]+2,dim(dm)[2]+2))
  dm2[2:(dim(dm)[1]+1),2:(dim(dm)[2]+1)]<-dm
  # Create 3D array of  basin numbers  with adjoining cells
  bma<-array(NA,dim=c(dim(bm),9))
  bma[,,1]<-bm # rw, cl
  bma[,,2]<-bm2[2:(dim(bm)[1]+1),1:dim(bm)[2]] # rw, cl-1
  bma[,,3]<-bm2[2:(dim(bm)[1]+1),3:(dim(bm)[2]+2)] # rw, cl+1
  bma[,,4]<-bm2[1:dim(bm)[1],2:(dim(bm)[2]+1)] # rw-1, cl
  bma[,,5]<-bm2[1:dim(bm)[1],1:dim(bm)[2]] # rw-1, cl-1
  bma[,,6]<-bm2[1:dim(bm)[1],3:(dim(bm)[2]+2)] # rw-1, cl+1
  bma[,,7]<-bm2[3:(dim(bm)[1]+2),2:(dim(bm)[2]+1)] # rw+1, cl
  bma[,,8]<-bm2[3:(dim(bm)[1]+2),1:dim(bm)[2]] # rw+1, cl-1
  bma[,,9]<-bm2[3:(dim(bm)[1]+2),3:(dim(bm)[2]+2)] # rw+1, cl+1
  # Create 3D array of elevation differences with adjoining cells
  dma<-array(NA,dim=c(dim(dm),9))
  dma[,,1]<-dm # rw, cl
  dma[,,2]<-dm2[2:(dim(dm)[1]+1),1:dim(dm)[2]]-dm # rw, cl-1
  dma[,,3]<-dm2[2:(dim(dm)[1]+1),3:(dim(dm)[2]+2)]-dm  # rw, cl+1
  dma[,,4]<-dm2[1:dim(dm)[1],2:(dim(dm)[2]+1)]-dm  # rw-1, cl
  dma[,,5]<-dm2[1:dim(dm)[1],1:dim(dm)[2]]-dm  # rw-1, cl-1
  dma[,,6]<-dm2[1:dim(dm)[1],3:(dim(dm)[2]+2)]-dm  # rw-1, cl+1
  dma[,,7]<-dm2[3:(dim(dm)[1]+2),2:(dim(dm)[2]+1)]-dm  # rw+1, cl
  dma[,,8]<-dm2[3:(dim(dm)[1]+2),1:dim(dm)[2]]-dm  # rw+1, cl-1
  dma[,,9]<-dm2[3:(dim(dm)[1]+2),3:(dim(dm)[2]+2)]-dm  # rw+1, cl+1
  dma2<-dma*0
  dma2[abs(dma)<boundary]<-1
  bma<-bma*dma2
  bma[,,1]<-bm
  # identify edge and basin merge cells
  be<-apply(bma,c(1,2),.edge)
  bea<-aperm(apply(bma,c(1,2),.edgec),c(2,3,1))
  s<-which(be>0,arr.ind=TRUE)
  for (i in 1:dim(s)[1]) {
    rw<-as.numeric(s[i,1])
    cl<-as.numeric(s[i,2])
    b3<-.asign3(bm2,bea,rw,cl)
    bm2[rw:(rw+2),cl:(cl+2)]<-b3
  }
  # reassign basin number
  u<-unique(as.vector(bm2))
  u<-u[is.na(u)==FALSE]
  u<-u[order(u)]
  bm3<-bm2
  for (i in 1:length(u)) {
    s<-which(bm2==u[i])
    bm3[s]<-i
  }
  dd<-dim(bm3)
  bsn<-bm3[2:(dd[1]-1),2:(dd[2]-1)]
  r<-.rast(bsn,dtm)
}

#' Mosaic tiled basins merging common joins
#'@noRd
.basinmosaic<-function(b1,b2) {
  e1<-ext(b1)
  e2<-ext(b2)
  reso<-res(b1)
  # *********** Do this if the tiles are vertically adjoined  *************** #
  if (abs(e1$ymax-e2$ymax) > reso[1]) {
    if (e2$ymax > e1$ymax) {  # b2 above b1
      m1<-.is(b1)
      m2<-.is(b2)
    } else {  # b1 above b2
      m1<-.is(b2)
      m2<-.is(b1)
    }
    for (itr in 1:3) {
      # merge based on top row of b1
      v1<-m1[1,] # top row of b1
      n<-dim(m2)[1] # mumber of rows
      v2<-m2[n,] # bottom row of b2
      # Create unique pairs matrix
      mup<-as.matrix(cbind(v1,v2))
      mup<-unique(mup)
      s<-which(is.na(mup[,1])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      s<-which(is.na(mup[,2])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      # Create vector of unique v1s
      u1<-unique(v1)
      u1<-u1[is.na(u1)==FALSE]
      u1<-u1[order(u1)]
      ras2<-list() # list of basins in m2 that should be re-asigned for each basin in u1
      ras1<-list() # list of basins in m1 that should be re-asigned for each basin in u1
      if (length(u1) > 0) {
        for (i in 1:length(u1)) {
          s<-which(mup[,1]==u1[i])
          u2<-mup[s,2] # list of basins in m2 that need reassigned
          u2<-u2[order(u2)]
          ras2[[i]]<-u2
          # list of basins in m1 that need reassinged
          s<-which(mup[,2]==u2[1])
          u1n<-mup[s,1] # list of basins in m2 that need reassigned
          if (length(u2) > 1) {
            for (j in 2:length(u2)) {
              s<-which(mup[,2]==u2[j])
              u1n<-c(u1n,mup[s,1]) # list of basins in m2 that need reassigned
            }
          }
          u1n<-unique(u1n)
          u1n<-u1n[u1n>u1[i]]
          ras1[[i]]<-u1n[order(u1n)]
          u2<-ras2[[i]]
          # Reassign basins in m2
          if (length(u2) > 0) for (j in 1:length(u2)) m2[m2==u2[j]]<-u1[i]
          u1n<-ras1[[i]]
          # Reassign basins in m1
          if (length(u1n) > 0) for (j in 1:length(u1n)) m1[m1==u1n[j]]<-u1[i]
        } # end for u1
      } # end if u1
    } # end iter
    # Convert back to SpatRasts
    if (e2$ymax > e1$ymax) {  # b2 above b1
      b1n<-.rast(m1,b1)
      b2n<-.rast(m2,b2)
    } else {  # b1 above b2
      b1n<-.rast(m2,b1)
      b2n<-.rast(m1,b2)
    }
  } else {# end do this if the tiles are vertically adjoined
    # *********** Do this if the tiles are horizontally adjoined  ************** #
    if (e2$xmax > e1$xmax) {  # b2 right of b1
      m1<-.is(b1)
      m2<-.is(b2)
    } else {  # b2 left of b1
      m1<-.is(b2)
      m2<-.is(b1)
    }
    for (itr in 1:3) {
      # merge based on right hand column of b1
      n<-dim(m1)[2]
      v1<-m1[,n] # right-hand column of b1
      v2<-m2[,1] # left-hand column of b2
      # Create unique pairs matrix
      mup<-as.matrix(cbind(v1,v2))
      mup<-unique(mup)
      s<-which(is.na(mup[,1])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      s<-which(is.na(mup[,2])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      # Create vector of unique v1s
      u1<-unique(v1)
      u1<-u1[is.na(u1)==FALSE]
      u1<-u1[order(u1)]
      ras2<-list() # list of basins in m2 that should be re-asigned for each basin in u1
      ras1<-list() # list of basins in m1 that should be re-asigned for each basin in u1
      if (length(u1) > 0) {
        for (i in 1:length(u1)) {
          s<-which(mup[,1]==u1[i])
          u2<-mup[s,2] # list of basins in m2 that need reassigned
          u2<-u2[order(u2)]
          ras2[[i]]<-u2
          # list of basins in m1 that need reassinged
          s<-which(mup[,2]==u2[1])
          u1n<-mup[s,1] # list of basins in m2 that need reassigned
          if (length(u2) > 1) {
            for (j in 2:length(u2)) {
              s<-which(mup[,2]==u2[j])
              u1n<-c(u1n,mup[s,1]) # list of basins in m2 that need reassigned
            }
          }
          u1n<-unique(u1n)
          u1n<-u1n[u1n>u1[i]]
          ras1[[i]]<-u1n[order(u1n)]
          u2<-ras2[[i]]
          # Reassign basins in m2
          if (length(u2) > 0) for (j in 1:length(u2)) m2[m2==u2[j]]<-u1[i]
          u1n<-ras1[[i]]
          # Reassign basins in m1
          if (length(u1n) > 0) for (j in 1:length(u1n)) m1[m1==u1n[j]]<-u1[i]
        } # end for u1
      } # end if u1
    } # end iter
    # Convert back to SpatRasts
    if (e2$xmax > e1$xmax) {  # b2 above b1
      b1n<-.rast(m1,b1)
      b2n<-.rast(m2,b2)
    } else {  # b1 above b2
      b1n<-.rast(m2,b1)
      b2n<-.rast(m1,b2)
    }
  }
  # ********************************** Mosaic ******************************* #
  bout<-mosaic(b1n,b2n)
  return(bout)
}
#' Do an entire column of tiled basins
#'@noRd
.docolumn<-function(dtm,tilesize,boundary,x) {
  e<-ext(dtm)
  reso<-res(dtm)
  ymxs<-as.numeric(ceiling((e$ymax-e$ymin)/reso[2]/tilesize))-1
  xmn<-as.numeric(e$xmin)+reso[1]*tilesize*x
  xmx<-xmn+reso[1]*tilesize
  ymn<-as.numeric(e$ymin)+reso[2]*tilesize*0
  ymx<-ymn+reso[2]*tilesize
  if (xmx > e$xmax) xmx<-e$xmax
  if (ymx > e$ymax) ymx<-e$ymax
  ec<-ext(xmn,xmx,ymn,ymx)
  dc<-crop(dtm,ec)
  bma<-basindelin(dc,boundary)
  # delineate basins for columns
  for (y in 1:ymxs) {
    xmn<-as.numeric(e$xmin)+reso[1]*tilesize*x
    xmx<-xmn+reso[1]*tilesize
    ymn<-as.numeric(e$ymin)+reso[2]*tilesize*y
    ymx<-ymn+reso[2]*tilesize
    if (xmx > e$xmax) xmx<-e$xmax
    if (ymx > e$ymax) ymx<-e$ymax
    ec<-ext(xmn,xmx,ymn,ymx)
    dc<-crop(dtm,ec)
    ta<-suppressWarnings(max(as.vector(bma),na.rm=T))
    if (is.infinite(ta)) ta<-0
    bo<-basindelin(dc,boundary)+ta
    bma<-.basinmosaic(bma,bo)
  } # end y
  return(bma)
}
#' Function used for delineating basins with big dtms
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclima, .registration = TRUE
#'@noRd
.basindelin_big<-function(dtm, boundary = 0, tilesize = 100, plotprogress = FALSE) {
  # chop into tiles
  e<-ext(dtm)
  reso<-res(dtm)
  xmxs<-as.numeric(ceiling((e$xmax-e$xmin)/reso[1]/tilesize))-1
  bma<-.docolumn(dtm,tilesize,boundary,0)
  for (x in 1:xmxs) {
    ta<-suppressWarnings(max(as.vector(bma),na.rm=T))
    if (is.infinite(ta)) ta<-0
    bo<-.docolumn(dtm,tilesize,boundary,x)+ta
    ed<-Sys.time()
    bma<-.basinmosaic(bma,bo)
    if (plotprogress) plot(bma,main=x)
  }
  m<-.is(bma)
  # renumber basins
  u<-unique(as.vector(m))
  u<-u[is.na(u) == FALSE]
  u<-u[order(u)]
  m<-array(renumberbasin(m,u),dim=dim(m))
  bout<-.rast(m,bout)
  return(bma)
}
#' Calculate flow direction
#'@noRd
.flowdir <- function(md) {
  fd<-md*0
  md2<-array(NA,dim=c(dim(md)[1]+2,dim(md)[2]+2))
  md2[2:(dim(md)[1]+1),2:(dim(md)[2]+1)]<-md
  v<-c(1:length(md))
  v<-v[is.na(md) == F]
  x<-arrayInd(v,dim(md))[,1]
  y<-arrayInd(v,dim(md))[,2]
  for (i in 1:length(x)) {
    md9<-md2[x[i]:(x[i]+2),y[i]:(y[i]+2)]
    fd[x[i],y[i]]<-round(mean(which(md9==min(md9,na.rm=TRUE))),0)
  }
  fd
}
#' @title Delineates hydrological or cold-air drainage basins
#' @description The function `basindelin` uses a digital elevation dataset to delineate
#' hydrological basins, merging adjoining basis separated by a low boundary if specified.
#' @param dtm a SpatRast object of elevations
#' @param boundary optional numeric value. If greater than 0, adjoining basins
#' separated by elevation differences < boundary are merged (see details.
#' @return a SpatRast of basins sequentially numbered as integers.
#' @details This function searches for the lowest grid cell in `dtm` and assigns it
#' as basin 1. All immediately adjacent pixels (in 8 directions) not previously assigned
#' are then assigned as being part of this basin if higher than the focal cell. The process is repeated
#' until no higher further cells are found. The next lowest unassigned grid cell is identified
#' and assigned as basin 2 and the process repeated until all grid cells are assigned a basin number.
#' If `boundary > 0`, edge grid cells are identified and the height difference from all
#' surrounding cells calculated. If the height difference is less than `boundary`, basins
#' are merged and the basins renumbered sequentially.
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclima, .registration = TRUE
#' @export
#' @rdname basindelin
#' @examples
#' library(terra)
#' basins <- basindelin(rast(dtm100m))
#' plot(basins, main = "Basins")
basindelin<-function(dtm, boundary = 0) {
  dm<-dim(dtm)
  if (sqrt(dm[1]*dm[2]) > 250) {
    bsn<-.basindelin_big(dtm, boundary)
  } else bsn<-.basindelin(dtm, boundary)
  return(bsn)
}
#' Calculates accumulated flow
#' @description
#' `flowacc` calculates accumulated flow(used to model cold air drainage)
#' @param dtm a SpatRast elevations (m).
#' @param basins optionally a SpatRast of basins numbered as integers (see details).
#' @return a SpatRast of accumulated flow (number of cells)
#' @details Accumulated flow is expressed in terms of number of cells. If `basins`
#' is provided, accumulated flow to any cell within a basin can only occur from
#' other cells within that basin.
#' @import terra
#' @export
#' @rdname flowacc
#' @examples
#' library(terra)
#' fa <- flowacc(rast(dtm100m))
#' plot(fa, main = 'Accumulated flow')
flowacc <- function (dtm, basins = NA) {
  dm<-.is(dtm)
  fd<-.flowdir(dm)
  fa<-fd*0+1
  if (class(basins) != "logical") ba<-.is(basins)
  o<-order(dm,decreasing=T,na.last=NA)
  for (i in 1:length(o)) {
    x<-arrayInd(o[i],dim(dm))[1]
    y<-arrayInd(o[i],dim(dm))[2]
    f<-fd[x,y]
    x2<-x+(f-1)%%3-1
    y2<-y+(f-1)%/%3-1
    # If basin file provided only add flow accumulation of from same basin
    if (class(basins) != "logical" & x2>0 & y2>0) {
      b1<-ba[x,y]
      b2<-ba[x2,y2]
      if(!is.na(b1) && !is.na(b2)){ # Check if both b1 and b2 valid basins
        if(b1==b2 & x2>0 & x2<dim(dm)[1] & y2>0 & y2<dim(dm)[2]) fa[x2,y2]<-fa[x,y]+1 }
    } else if (x2>0 & x2<dim(dm)[1] & y2>0 & y2<dim(dm)[2]) fa[x2,y2]<-fa[x,y]+1
  }
  fa<-.rast(fa,dtm)
  return(fa)
}
#' Calculates whether conditions are right for cold air drainage
#'
#' @description
#' `cadconditions` determines whether wind speed, humidity and cloud cover  are such that
#' cold air drainage is likely to occur
#'
#' @param h a vector of hourly specific humidities (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc a single numeric value, SpatRaster object, two-dimensional array or matrix of temperatures (ºC).
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
#' @param dem a SpatRaster object of elevation (m).
#' @param basins a SpatRaster object, two-dimensional array or matrix with basins numbered as integers as returned by [basindelin()].
#' @param fa a SpatRaster object of accumulated flow, as returned by [flowacc()]
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
#' @import terra
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
#' library(terra)
#' basins <- basindelin(rast(dtm100m), boundary = 2)
#' h <- humidityconvert(50, intype = "relative", 20)$specific
#' fa <- flowacc(rast(dtm100m), basins)
#' cp1 <- pcad(rast(dtm100m), basins, fa, 20, h)
#' cp2 <- pcad(rast(dtm100m), basins, fa, 20, h, out = "tempdif")
#' cp3 <- pcad(rast(dtm100m), basins, fa, 20, h, out = "pflow")
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

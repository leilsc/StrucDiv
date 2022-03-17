############################################################################################################################
#
# getValuesWindow
#
############################################################################################################################
#
# # Dependent Function from Raster package - modified
#
############################################################################################################################

# GET NEIGHBOURHOOD VALUES FOR EACH PIXEL

# Modified R Code from raster package - getValuesFocal
# to retrieve values of defined neighbourhood
# returns one row per pixel (row containing pixel neighbourhood (wsl)
# currently programmed with NA padding

############################################################################################################################

.getValuesWindow <- function(x, wsl, padValue, aroundTheGlobe, ...) {
  
  message("Extracting values from raster file")

  wsl <- c(wsl, wsl)
  nl <- raster::nlayers(x)
  if (nl == 0) {
    stop("x has no values")
  }

  xx <- raster::raster(x)
  nc <- ncol(xx)
  nr <- nrow(xx)

  stopifnot(is.atomic(padValue))
  # geo <- raster::couldBeLonLat(xx) ## returns FALSE when projected

  ngbr <- floor(wsl[1] / 2)
  ngbc <- floor(wsl[2] / 2)

  startrow <- 1 - ngbr
  endrow <- nr + ngbr

  vv <- matrix(raster::getValues(x, 1, nr), ncol = 1) # vv: one column contains the values returned for one row
  v <- matrix(vv, ncol = nc, byrow = TRUE)

  if (ngbr > 0) {
    v <- rbind(matrix(padValue, nrow = ngbr, ncol = ncol(v)), v, matrix(padValue, nrow = ngbr, ncol = ncol(v)))
  }

  ## If the image goes around the globe, add stripes of the image itself on the left and the right
  res <- FALSE
  tolerance <- 0.1
  scale <- raster::xres(x)
  if (isTRUE(all.equal(raster::xmin(xx), -180, tolerance=tolerance, scale=scale)) &
      isTRUE(all.equal(raster::xmax(xx),  180, tolerance=tolerance, scale=scale))) {
    if (raster::couldBeLonLat(xx, warnings=FALSE)) {
      res <- TRUE
    }
  }
  gll <- as.integer(res)

  if (gll & aroundTheGlobe) {
    nv <- ncol(v)
    if (ngbc < nv) {
      v <- cbind(v[, (nv - ngbc + 1):nv], v, v[, 1:ngbc])
    } else {
      stop("neighbourhood is too big")
    }
  } else {
    add <- matrix(padValue, ncol = ngbc, nrow = nrow(v))
    v <- cbind(add, v, add)
  }

  v <- .focal_get(v, as.integer(dim(v)), as.integer(wsl))

  m <- matrix(v, nrow = nrow(x) * nc, byrow = TRUE)

  return(m)
}



.isGlobalLonLat <- function(x) {
  res <- FALSE
  tolerance <- 0.1
  scale <- raster::xres(x)
  if (isTRUE(all.equal(raster::xmin(x), -180, tolerance=tolerance, scale=scale)) & 
      isTRUE(all.equal(raster::xmax(x),  180, tolerance=tolerance, scale=scale))) {
    if (raster::couldBeLonLat(x, warnings=FALSE)) {
      res <- TRUE
    }
  }
  res
}


.G <- function(Mat, overlap, edge){
  for (i in overlap:(1 + edge)){
    Mat[ ,i] <- (i - edge) / (overlap - edge)
    Mat[ ,(ncol(Mat) - i + 1)] <- (i - edge) / (overlap - edge)
    Mat[i, ] <- (i - edge) / (overlap - edge)
    Mat[(nrow(Mat) - i + 1), ] <- (i - edge) / (overlap - edge)
  }
  Mat[ ,0:edge] <- 0 
  Mat[ ,(ncol(Mat) - edge + 1):ncol(Mat)] <- 0 
  Mat[0:edge, ] <- 0 
  Mat[(nrow(Mat) - edge + 1):nrow(Mat), ] <- 0 
  return(Mat)
}


############################################################################################################################

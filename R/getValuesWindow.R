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

#' @name getValuesWindow
#' @rdname getValuesWindow
#' @title Get the values from a raster object centered at each pixel
#' @description Modified R Code from raster package - getValuesFocal,
#' to retrieve values of defined neighbourhood. Returns one row per pixel 
#' (row containing pixel neighbourhood (wsl)
#' @param x raster layer. Input raster layer for which 
#' @param wsl integer. size of neighbourhood
#' @param padValue atomic. If a pixel is on the edges of an image, what value should it be padded with.
#' @param aroundTheGlobe logical. Does the image go around the globe.
#' @return A matrix
#' @export

getValuesWindow <- function(x, wsl, padValue, aroundTheGlobe, ...) {
  
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


.WM <- function(nrow, ncol, ul, nNA){
  ul <- ul - 2*nNA
  ul1 <- ul + 1
  w <- matrix(1, nrow, ncol)
  w[nNA+(1:ul),] <- 1:ul/ul1
  w[nrow-nNA-(0:(ul-1)),] <- 1:ul/ul1
  w[, nNA+1:ul] <- w[,nNA+1:ul] * (rep(1:ul, each=ncol))/ul1
  w[, ncol-nNA-(0:(ul-1))] <- w[,ncol-nNA-(0:(ul-1))] * (rep(1:ul, each=ncol))/ul1
  w[1:nNA,] <- 0
  w[,1:nNA] <- 0
  w[nrow-(0:(nNA-1)),] <- 0
  w[,ncol-(0:(nNA-1))] <- 0
  w
}


############################################################################################################################

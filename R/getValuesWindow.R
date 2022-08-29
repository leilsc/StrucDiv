#============================================================================================
#
# StrucDiv::getValuesWindow
# Modified function from raster package, raster:: getValuesFocal()
#
#============================================================================================

#' @name getValuesWindow
#' @rdname getValuesWindow
#' @title Retrieve pixel values of a defined area. 
#' The area is defined by the size of a window, which is centered on one pixel.
#' @description Modified R Code from raster package raster::getValuesFocal.
#' Returns one row per pixel, which contains the values
#' of the pixel neighborhood that is defined by the size of the window.
#' The size of the window is defined by the window side length (wsl).
#' The window is centered on one specific pixel.
#' @param x raster layer. The input raster layer.
#' @param wsl integer. The window side length. The window is defined by \code{wsl x wsl}.
#' @param padValue atomic. If a pixel is on the edge of an image, padding should be used?
#' Can be NA or a value.
#' @param aroundTheGlobe logical. Does the image go around the globe?
#' @return Returns a matrix.
#' The matrix contains the values of the defined window centered on the respective pixel.
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

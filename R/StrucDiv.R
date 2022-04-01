#' @name StrucDiv2
#' @rdname StrucDiv2
#' @title Calculate spatial (horizontal) structural diversity for an arbitrary raster layer. 
#' @description
#' This is a wrapper function that returns a spatial structural diversity map 
#' as a raster layer. 
#' Pixel pairs are considered at user defined distances and angles within a moving window, 
#' the resulting frequency table is turned into a gray level co-occurrence matrix (GLCM).
#' Based on the GLCM, spatial structural diversity metrics are calculated, 
#' the metric value is assigned to the center pixel of the moving window. 
#' The output is a diversity map, which is returned as a raster layer.
#' @param x raster layer. Input raster layer for which 
#' horizontal structural diversity should be calculated.
#' @param wsl uneven integer. The window side length, 
#' wsl x wsl defines the size of the moving window.
#' The window must be smaller than the dimensions of the input raster.
#' @param dist integer. The distance between two pixels that should be considered as a pair, 
#' defaults to dist = 1 (direct neighbors).
#' @param angle string. The angle on which pixels should be considered as pairs. 
#' Takes 5 options: "horizontal", "vertical", "diagonal45", "diagonal135", "all". 
#' "all" considers all of the 4 angles. Defaults to "all".
#' @param rank logical. Should pixel values be replaced with ranks in each GLCM? Defaults to FALSE.
#' @param fun function, the diversity metric. Takes one of the following: entropy,
#' entropyNorm, contrast, dissimilarity, or homogeneity.
#' @param delta integer, takes 3 options: 0, 1, or 2. 
#' Delta is the difference weight parameter, 
#' it defines how the differences between pixel values within a pixel pair should be weighted.  
#' If rank = TRUE, delta defines how the differences between ranks should be weighted.  
#' The default value is 0 (no weight). Set delta = 1 for absolute difference weight, 
#' or delta = 2 for squared difference weight. 
#' The delta parameter can only be set when the metric entropy is used. 
#' Dissimilarity automatically employs delta = 1, and contrast employs delta = 2.
#' @param na.handling na.omit or na.pass. 
#' If na.handling = na.omit, NAs are ignored, diversity metrics are calculated with less values. 
#' In this case the GLCM does not sum to 1.
#' If na.handling = na.pass and if there is at least one missing value inside the moving window,
#' an NA is assigned to the center pixel. Therefore, the diversity map will contain more 
#' NAs than the input raster layer.
#' Defaults to na.pass.
#' @param padValue numeric. The value of the padded cells at the edges of the input raster. 
#' Defaults to NA.
#' @param aroundTheGlobe logical. If the input raster goes around the whole globe, 
#' set aroundTheGlobe = TRUE, and the input raster will be "glued together" from both sides
#' to calculate diversity without edge effects on the sides.
#' Defaults to FALSE.
#' @param filename character. If the output raster should be written to a file, define file name (optional).
#' @param ... possible further arguments.
#' @return A raster layer
#' @examples 
#' 
#' a <- raster::raster(matrix(rnorm(648), 18, 36))
#' raster::plot(a)
#' b <- StrucDiv(a, 3, fun = contrast, na.handling = na.omit, rank = FALSE)
#' raster::plot(b)
#' 
#' @export


StrucDiv2 <- function(x, wsl, dist = 1, angle = "all",
                     rank = FALSE, fun, delta = 0, 
                     na.handling = na.pass, padValue = NA, 
                     aroundTheGlobe = FALSE, filename = "", display_progress = TRUE, ...) {
  
  dotArgs <- list(...)
  
  if(isTRUE(aroundTheGlobe) & !isTRUE(.isGlobalLonLat(x))){
    warning("The raster image does not go around the globe.")
  }
  
  suppressWarnings(
    if ( identical(na.handling, na.pass) && anyNA(raster::values(x)) ) {
      warning("Raster layer contains missing values. Wherever there are missing values,
             an NA will be returned. if you want to proceed without NAs, 
             set na.handling = na.omit.")
    }
  )
  
  out <- raster::raster(x)
  
  stopifnot(raster::hasValues(x))
  
  if (wsl == 0) {
    stop("Window side length must be > 0.")
  }
  if (wsl %% 2 == 0) {
    stop("Window side length must be an odd number.")
  }
  if (wsl > min(dim(out)[1:2])) {
    stop("Window must be smaller than the raster dimensions.")
  }
  
  if (dist > min(dim(out)[1:2]) / 2 - 1) {
    stop("Distance value is too big.")
  }
  
  if (!(angle %in% c("horizontal", "vertical", "diagonal45", "diagonal135", "all"))) {
    stop('Angle must be one of "horizontal", "vertical", "diagonal45", "diagonal135", or "all".')
  }
  
  if(!is.function(fun)){
    stop("This diversity metric is not available.")
  }
  
  if ( !(delta %in% c(0,1,2)) ) {
    stop("Delta must be 0, 1, or 2.")
  }
  
  filename <- glue::trim(filename)
  
  if (raster::canProcessInMemory(out)) {
    
    vMat <- .getValuesWindow(x, wsl = wsl, padValue = padValue, 
                                  aroundTheGlobe = aroundTheGlobe)
    Hetx <- vMat
    
    suppressWarnings(
      if( identical(na.handling, na.omit) && anyNA(raster::values(x)) ) {
        
        narm <- 1
        
      } 
      else{
        
        narm <- 0
      
      }
    )
    
    
    switch_angle <- function(angle) {
      
      switch(angle,
             "horizontal" = .ProbabilityMatrixHorizontalDynamic(vMat = vMat, d = dist, narm = narm),
             "vertical" = .ProbabilityMatrixVerticalDynamic(vMat = vMat, d = dist, narm = narm),
             "diagonal45" = .ProbabilityMatrixDiagonal45Dynamic(vMat = vMat, d = dist, narm = narm),
             "diagonal135" = .ProbabilityMatrixDiagonal135Dynamic(vMat = vMat, d = dist, narm = narm),
             "all" = .ProbabilityMatrixAllDynamic(vMat = vMat, d = dist, narm = narm),
             .ProbabilityMatrixHorizontalDynamic(vMat = vMat, d = dist, narm = narm)
      )
      
    }
    
    SpatMat <- switch_angle(angle)
    
    if (angle %in% c("horizontal", "vertical")) {
      nrp <- 2*wsl*(wsl - dist)
    }
    
    if (angle %in% c("diagonal45", "diagonal135")) {
      nrp <- (wsl - dist) * 2 *(wsl - dist)
    }
    
    if (angle == "all") {
      nrp <- 4 * (wsl - dist) * (2 * wsl - dist)
    }
    
    v <- do.call(fun, list(rank = rank, Hetx = Hetx, SpatMat = SpatMat, delta = delta,
                           nrp = nrp, narm = narm, display_progress = display_progress))
    
    out <- raster::setValues(out, v)
    return(out)
    
    if (filename != "") {
      out <- raster::writeRaster(out, filename)
    }
    
  } else {
    stop("Cannot proccess in memory.")
  }
  
}

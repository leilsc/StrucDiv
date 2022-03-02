#' @name Block_StrucDiv
#' @rdname Block_StrucDiv
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
#' @param window is a vector of length 2. This is the block size.
#' @param wsl uneven integer. The window side length, 
#' wsl x wsl defines the size of the moving window.
#' The window must be smaller than the dimensions of the input raster.
#' @param WSLw an uneven integer smaller than the dimensions of the input raster 
#' and larger than wsl.
#' @param overlap an integer that defines the overlap between the blocks.
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


Block_StrucDiv <- function(x, window, wsl, WSLw, 
                           overlap,
                           dist = 1, angle = "all",
                           rank = FALSE, fun, delta = 0, 
                           na.handling = na.pass, padValue = NA, 
                           aroundTheGlobe = FALSE, filename = "", ...) {
  
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
  
  x <- raster::crop(x, raster::extent(x, 1, (nrow(x) - nrow(x) %% window[1]), 1, 
                                      (ncol(x) - ncol(x) %% window[2])))
  
  out <- raster::crop(x, raster::extent(x, 1, 1, 1, 1))

 
  # if (window[1] %% 2 != 0 | window[2] %% 2 != 0) {
  #   stop("window must have even sides")
  # }  ## Is is required to be an even number?


  filename <- raster::trim(filename)

  if (raster::canProcessInMemory(x)) {
    
    
    switch_angle <- function(angle) {
      
      switch(angle,
             "horizontal" = .ProbabilityMatrixHorizontalNested(vMat = vMat, 
                                                               vMat_big = vMat_big,
                                                               d = dist),
             "vertical" = .ProbabilityMatrixVerticalNested(vMat = vMat, 
                                                           vMat_big = vMat_big, 
                                                           d = dist),
             "diagonal45" = .ProbabilityMatrixDiagonal45Nested(vMat = vMat, 
                                                               vMat_big = vMat_big,
                                                               d = dist),
             "diagonal135" = .ProbabilityMatrixDiagonal135Nested(vMat = vMat, 
                                                                 vMat_big = vMat_big,
                                                                 d = dist),
             "all" = .ProbabilityMatrixAllNested(vMat = vMat, 
                                                 vMat_big = vMat_big,
                                                 d = dist),
             .ProbabilityMatrixHorizontalNested(vMat = vMat,
                                                vMat_big = vMat_big,
                                                d = dist)
      )
      
    }
    
    
    if (angle %in% c("horizontal", "vertical")) {
      nrp <- 2*wsl*(wsl - dist)
      # nrp_big <- 2*WSLw*(WSLw - dist)
      
    }
    
    if (angle %in% c("diagonal45", "diagonal135")) {
      nrp <- (wsl - dist) * 2 *(wsl - dist)
      # nrp_big <- (WSLw - dist) * 2 *(WSLw - dist)
    }
    
    if (angle == "all") {
      nrp <- 4 * (wsl - dist) * (2 * wsl - dist)
      # nrp_big <- 4 * (WSLw - dist) * (2 * WSLw - dist)
    }
    
    wmx <- matrix(1, ncol = window[2], nrow = window[1])
    
    G <- function(Mat, overlap, edge){
      #Browser()
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
    
    wmx <- G(Mat = wmx, overlap = overlap, edge = WSLw)
    
    
      #Browser()
      # num <- raster(nrow=nrow(x), ncol = ncol(x)) # numerator file
      # denom <- raster(nrow=nrow(x), ncol = ncol(x)) # denominator file
      # 
      # raster::crs(num) <- raster::crs(x)
      # raster::res(num) <- raster::res(x)
      # 
      # raster::crs(denom) <- raster::crs(x)
      # raster::res(denom) <- raster::res(x)
    
    num <- out
    denom <- out
      
      #create row and col indexes - modification from FocalRcppHybrid
      stepsize <- window[1]-overlap
      rowtimes <- (nrow(x) - window[1]) / (window[1] - overlap) + 1
      times <- rep(stepsize, rowtimes - 1)
      RowIndex <- c(1, cumsum(times) + 1)
      
      stepsize <- window[2]-overlap
      coltimes <- (ncol(x) - window[2]) / (window[2] - overlap) + 1
      times <- rep(stepsize, coltimes - 1)
      ColIndex <- c(1, cumsum(times) + 1)
      
      for (i in RowIndex) {
        for (j in ColIndex) {
          block <- raster::getValuesBlock(x = x, row = i, nrows = window[1], 
                                          col = j, ncols = window[2], format = "matrix")
          blockra <- raster::raster(block)
          raster::extent(blockra) <- raster::extent(x, i, i + window[1] - 1, j, j + window[2] - 1)
          raster::crs(blockra) <- raster::crs(x)
          raster::res(blockra) <- raster::res(x)
          
          vMat <- .getValuesWindow(blockra, wsl = wsl, padValue = padValue, 
                                   aroundTheGlobe = aroundTheGlobe)
          
          suppressMessages( vMat_big <- .getValuesWindow(blockra, wsl = WSLw, padValue = padValue, 
                                                         aroundTheGlobe = aroundTheGlobe)  )
          
          Hetx <- vMat
          
          suppressWarnings(
            if( identical(na.handling, na.omit) && anyNA(raster::values(x)) ) {
              
              narm <- 1
              
            } 
            else{
              
              narm <- 0
              
            }
          )
          
          # switch_angle <- function(angle) {
          #   
          #   switch(angle,
          #          "horizontal" = .ProbabilityMatrixHorizontalNested(vMat = vMat, 
          #                                                            vMat_big = vMat_big,
          #                                                            d = dist),
          #          "vertical" = .ProbabilityMatrixVerticalNested(vMat = vMat, 
          #                                                        vMat_big = vMat_big, 
          #                                                        d = dist),
          #          "diagonal45" = .ProbabilityMatrixDiagonal45Nested(vMat = vMat, 
          #                                                            vMat_big = vMat_big,
          #                                                            d = dist),
          #          "diagonal135" = .ProbabilityMatrixDiagonal135Nested(vMat = vMat, 
          #                                                              vMat_big = vMat_big,
          #                                                              d = dist),
          #          "all" = .ProbabilityMatrixAllNested(vMat = vMat, 
          #                                              vMat_big = vMat_big,
          #                                              d = dist),
          #          .ProbabilityMatrixHorizontalNested(vMat = vMat,
          #                                             vMat_big = vMat_big,
          #                                             d = dist)
          #   )
          #   
          # }
          
          SpatMat <- switch_angle(angle)
          
          if (angle %in% c("horizontal", "vertical")) {
            nrp <- 2*wsl*(wsl - dist)
            # nrp_big <- 2*WSLw*(WSLw - dist)
            
          }
          
          if (angle %in% c("diagonal45", "diagonal135")) {
            nrp <- (wsl - dist) * 2 *(wsl - dist)
            # nrp_big <- (WSLw - dist) * 2 *(WSLw - dist)
          }
          
          if (angle == "all") {
            nrp <- 4 * (wsl - dist) * (2 * wsl - dist)
            # nrp_big <- 4 * (WSLw - dist) * (2 * WSLw - dist)
          }
          
          sdiv <- do.call(fun, list(rank = rank, Hetx = Hetx, SpatMat = SpatMat, delta = delta,
                                 nrp = nrp, narm = narm))
          
          # multiply each block with the spatial weights matrix 
          wdiv <- sdiv * wmx
          # rasterize weighted blocks
          wdiv <- raster::raster(wdiv)
          raster::extent(wdiv) <- raster::extent(x, i, i + window[1] - 1, j, j + window[2] - 1)
          raster::crs(wdiv) <- raster::crs(x)
          raster::res(wdiv) <- raster::res(x)
          
          wmr <- raster(wmx) # rasterize weights matrix so we can create denominator layer
          WMRext <- setValues(wdiv, values = values(wmr)) # each block gets a weights matrix layer 
          # that has the same extent, i.e. is in the same place as the block
          
          # create numerator - take the sum of overlapping weighted blocks (i.e. take the sum where they overlap)
          num <- raster::mosaic(num, wdiv, fun = sum) 
          # create denominator - take the sum of overlapping weights (i.e. take the sum where they overlap)
          denom <- raster::mosaic(denom, WMRext, fun = sum)
          # create final raster layer
          #out <- num/denom
        }
      }
      out <- num/denom
      
      # out <- raster::setValues(out, v)

      return(out)
    
    if (filename != "") {
      out <- raster::writeRaster(out, filename)
    }
    
  } else {
    stop("Cannot proccess in memory.")
  }
  }
  

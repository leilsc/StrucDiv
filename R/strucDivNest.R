#' @name strucDivNest
#' @rdname strucDivNest
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
#' @param wslI uneven integer. The window side length, 
#' wslI x wslI defines the size of the moving window.
#' The window must be smaller than the dimensions of the input raster. Default is NULL
#' in which case no nesting is used.
#' @param wslO an uneven integer smaller than the dimensions of the input raster 
#' and larger than wslI.
#' @param dimB is a vector of length 2 or a logical. This is the block size.
#' Default is FALSE in which case the image is not cut into blocks.
#' @param oLap an integer that defines the overlap between the blocks. It must
#' be at least 2*(wslO-1) or bigger. Blocks can overlap by a maximum of half the 
#' rows of blocks in row-direction, and by half the columns of blocks in column-direction.
#' @param priorB logical Should the blocks be used for prior estimates? 
#' If priorB = TRUE, then the spatial structure in a block serves as prior for 
#' spatial structural diversity estimates of the moving windows inside the block. 
#' If priorB = FALSE, then the blocks are only used to increase speed through 
#' parallelization, not as prior estimates.
#' @param domain logical. Should the spatial structure in the domain 
#' (i.e. the input raster) be used as prior estimate for the moving windows? 
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
#' @param parallelize logical. Set to TRUE to parallelize the computation on multiple cores.
#' Default is FALSE. Please set to TRUE only for bigger problems.
#' @param display_progress logical
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


strucDivNest <- function(x, wslI = NULL, wslO = NULL, dimB = FALSE, oLap = NULL, 
                         priorB = FALSE, domain = FALSE, dist = 1, angle = "all",
                         rank = FALSE, fun, delta = 0, 
                         na.handling = na.pass, padValue = NA, 
                         aroundTheGlobe = FALSE, 
                         parallelize = FALSE, 
                         display_progress = FALSE, 
                         filename = "", ...) {
  
  dotArgs <- list(...)
  
  ## General warnings and errors!
  
  if(!is.logical(aroundTheGlobe)){
    stop("aroundTheGlobe must be either TRUE, or FALSE.")
  }
  
  if(!is.logical(parallelize)){
    stop("parallelize must be either TRUE, or FALSE.")
  }
  
  if(isTRUE(aroundTheGlobe) & !isTRUE(.isGlobalLonLat(x))){
    warning("The raster image does not go around the globe.")
  }
  
  if ( identical(na.handling, na.pass) && anyNA(raster::values(x)) ) {
      warning("Raster layer contains missing values. Wherever there are missing values,
              an NA will be returned. if you want to proceed without NAs, 
              set na.handling = na.omit.")
    }
  
  out <- raster::raster(x)
  
  ## Check if the raster has values inside.
  
  stopifnot(raster::hasValues(x))
  
  ## Walk through all the argument checks
  
  if(is.null(wslI)) {
    stop("wslI must be specified")
  }
  if (wslI %% 2 == 0) {
    stop("Window side length must be an odd number.")
  }
  if (wslI > min(dim(out)[1:2])) {
    stop("Window must be smaller than the raster dimensions.")
  }
  if (wslI == 0) {
    stop("Window side length must be > 0.")
  }
  
  if(!is.null(wslO)){
    if (wslO == 0) {
      stop("Window side length must be > 0.")
    }
    if (wslO %% 2 == 0) {
      stop("Window side length must be an odd number.")
    }
    if (wslO > min(dim(out)[1:2])) {
      stop("Window must be smaller than the raster dimensions.")
    }
    if ( wslO < wslI) {
      stop("Large window must be larger than the small window.")
    }
  }
  
  if((dimB == FALSE)[1] & !is.null(oLap)){
    oLap <- NULL
    warning("oLap is ignored because dimB = FALSE")
  }
  
  if((dimB != FALSE)[1] & is.null(oLap)){
    
    minOL <- 2*(wslI-1)
    maxOL <- 0.5* min(dimB)
    oLap <- floor((maxOL-minOL)/2 + minOL)
    warning("oLap is calculated to lie in the middle of its min and max values 
            because it was not specified.")
    
  }
  
  if(!is.null(oLap) & is.numeric(dimB)){
    
    if (oLap < 2*(wslI-1)) {
      stop("Overlap is too small. Please input an overlap of at least 2*(wslI-1).")
    }
    
    if ( oLap > floor(dimB[1]/2) ) {
      stop("Overlap is too big. Please input an overlap of at most half the window size.")
    }
    
  }
  
  if(is.numeric(dimB) & is.null(oLap)){
    stop("You must specify overlap.")
  }
  
  if(!is.logical(priorB)){
    stop("priorB must be logical.")
  }
  
  if(!is.logical(domain)){
    stop("domain must be logical.")
  }
  
  if(priorB == FALSE & (dimB == TRUE)[1]) {
  
  warning("Blocks are only used for parallelization. 
          If you want to use blocks for prior information, set priorB = TRUE.")

  }
  
  if(priorB == TRUE & (dimB == FALSE)[1]) {
    
    stop("Specify dimB and oLap to use blocks for prior information.")
    
  }
  
  if( (!is.null(wslO) & priorB == TRUE) | (!is.null(wslO) & domain == TRUE) | (priorB == TRUE & domain == TRUE) ) {
    stop("No nested scales approach is implemented because more than one scale is 
         specified for prior information. If you want to use a nested scales approach, 
         please specify wslO, or specify dimB and oLap and set priorB = TRUE, or set domain = TRUE."
    )
  }
   
  if( is.null(wslO) & (dimB == FALSE)[1] & domain == FALSE ){
    
    warning("No prior information is provided, hence you are not using a nested-scales approach. 
            The output reflects spatial structural diversity approximation on a scale of wslI $\times$ wslI. 
            If you want to use a nested scales approach, please specify wslO, 
            or specify dimB and oLap and set priorB = TRUE, or set domain = TRUE.")
    
    out <- strucDiv(x = x, wsl = wslI, dist = dist, angle = angle,
            rank = rank, fun = fun, delta = delta, 
            na.handling = na.handling, padValue = padValue, 
            aroundTheGlobe = aroundTheGlobe, filename = filename, 
            display_progress = display_progress)
    
  }
  
  else{ ## START checks of other arguments and do nesting
    
  if (dist > min(dim(out)[1:2]) / 2 - 1) {
    stop("Distance value is too big.")
  }
  
  if (!(angle %in% c("horizontal", "vertical", "diagonal45", "diagonal135", "all"))) {
    stop('Angle must be one of "horizontal", "vertical", "diagonal45", "diagonal135", or "all".')
  }
    
  if(!is.logical(rank)){
    stop("rank mjust be either TRUE or FALSE")
  }
  
  if(!is.function(fun)){
    stop("This diversity metric is not available.")
  }
  
  if ( !(delta %in% c(0,1,2)) ) {
    stop("Delta must be 0, 1, or 2.")
  }
    
  if(!(is.numeric(padValue) | is.na(padValue))){
    stop("padValue must be a number or NA.")
  }
    
  filename <- glue::trim(filename)

  if (raster::canProcessInMemory(x)) {
    
    if(!is.null(wslO)) {  ## SMALLER WINDOW NESTED IN A LARGER WINDOW
    switch_angle <- function(angle) {
      
      switch(angle,
             "horizontal" = .ProbabilityMatrixHorizontalNested(vMat = vMat, 
                                                               vMat_big = vMat_big,
                                                               d = dist,
                                                               display_progress=FALSE),
             "vertical" = .ProbabilityMatrixVerticalNested(vMat = vMat, 
                                                           vMat_big = vMat_big, 
                                                           d = dist,
                                                           display_progress=FALSE),
             "diagonal45" = .ProbabilityMatrixDiagonal45Nested(vMat = vMat, 
                                                               vMat_big = vMat_big,
                                                               d = dist,
                                                               display_progress=FALSE),
             "diagonal135" = .ProbabilityMatrixDiagonal135Nested(vMat = vMat, 
                                                                 vMat_big = vMat_big,
                                                                 d = dist,
                                                                 display_progress=FALSE),
             "all" = .ProbabilityMatrixAllNested(vMat = vMat, 
                                                 vMat_big = vMat_big,
                                                 d = dist,
                                                 display_progress=FALSE),
             .ProbabilityMatrixAllNested(vMat = vMat,
                                                vMat_big = vMat_big,
                                                d = dist,
                                                display_progress=FALSE)
      )
      
    }
    
    vMat <- .getValuesWindow(x, wsl = wslI, padValue = padValue, 
                             aroundTheGlobe = aroundTheGlobe)
    
    suppressMessages( vMat_big <- .getValuesWindow(x, wsl = wslO, padValue = padValue, 
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
    
    if (angle %in% c("horizontal", "vertical")) {
      nrp <- 2*wslI*(wslI - dist)
    }
    if (angle %in% c("diagonal45", "diagonal135")) {
      nrp <- (wslI - dist) * 2 *(wslI - dist)
    }
    if (angle == "all") {
      nrp <- 4 * (wslI - dist) * (2 * wslI - dist)
    }
    
    SpatMat <- switch_angle(angle)
    
    v <- do.call(fun, list(rank = rank, Hetx = Hetx, SpatMat = SpatMat, delta = delta,
                           nrp = nrp, narm = narm, display_progress = display_progress, 
                           parallelize = parallelize))
    
    out <- raster::setValues(out, v)
    
    }
    
    if(domain == TRUE){  ## USE THE WHOLE IMAGE AS A PRIOR
      
      vMat <- .getValuesWindow(x, wsl = wslI, padValue = padValue, 
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
               "horizontal" = .ProbabilityMatrixHorizontalPost(vMat = vMat, 
                                                               x = xMat, d = dist, 
                                                               nrp = nrp, nrp_big = nrp_big,
                                                               display_progress=FALSE),
               "vertical" = .ProbabilityMatrixVerticalPost(vMat = vMat, 
                                                           x = xMat, d = dist, 
                                                           nrp = nrp, nrp_big = nrp_big,
                                                           display_progress=FALSE),
               "diagonal45" = .ProbabilityMatrixDiagonal45Post(vMat = vMat, 
                                                               x = xMat, d = dist, 
                                                               nrp = nrp, nrp_big = nrp_big,
                                                               display_progress=FALSE),
               "diagonal135" = .ProbabilityMatrixDiagonal135Post(vMat = vMat, 
                                                                 x = xMat, d = dist, 
                                                                 nrp = nrp, nrp_big = nrp_big,
                                                                 display_progress=FALSE),
               "all" = .ProbabilityMatrixAllPost(vMat = vMat, 
                                                 x = xMat, d = dist, 
                                                 nrp = nrp, nrp_big = nrp_big,
                                                 display_progress=FALSE),
               .ProbabilityMatrixAllPost(vMat = vMat, 
                                         x = xMat, d = dist, 
                                         nrp = nrp, nrp_big = nrp_big,
                                         display_progress=FALSE)
        )
        
      }
      
      xMat <- matrix(raster::values(x), nrow(x), ncol(x), byrow = TRUE)
      if (angle == "horizontal") {
        nrp <- 2*wslI*(wslI - dist)
        nrp_big <- 2*nrow(xMat)*(ncol(xMat) - dist)
      }
      if (angle == "vertical") {
        nrp <- 2*wslI*(wslI - dist)
        nrp_big <- 2*ncol(xMat)*(nrow(xMat) - dist)
      }
      if (angle %in% c("diagonal45", "diagonal135")) {
        nrp <- (wslI - dist) * 2 *(wslI - dist)
        nrp_big <- (nrow(xMat) - dist) * 2 *(ncol(xMat) - dist)
      }
      if (angle == "all") {
        nrp <- 4 * (wslI - dist) * (2 * wslI - dist)
        nrp_big <- 4 * (nrow(xMat) - dist) * (2 * ncol(xMat) - dist)
      }
      
      SpatMat <- switch_angle(angle)

      v <- do.call(fun, list(rank = rank, Hetx = Hetx, SpatMat = SpatMat, delta = delta,
                             nrp = nrp, narm = narm, display_progress = display_progress, 
                             parallelize = parallelize))
      
      out <- raster::setValues(out, v)
      
    }
    
    if( (dimB != FALSE)[1] ){
      
      if( priorB == FALSE & is.null(wslO)){
        
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
        
      }
      if( priorB == FALSE & !is.null(wslO)) {
        
        switch_angle <- function(angle) {
          
          switch(angle,
                 "horizontal" = .ProbabilityMatrixHorizontalNested(vMat = vMat, 
                                                                   vMat_big = vMat_big,
                                                                   d = dist,
                                                                   display_progress=FALSE),
                 "vertical" = .ProbabilityMatrixVerticalNested(vMat = vMat, 
                                                               vMat_big = vMat_big, 
                                                               d = dist,
                                                               display_progress=FALSE),
                 "diagonal45" = .ProbabilityMatrixDiagonal45Nested(vMat = vMat, 
                                                                   vMat_big = vMat_big,
                                                                   d = dist,
                                                                   display_progress=FALSE),
                 "diagonal135" = .ProbabilityMatrixDiagonal135Nested(vMat = vMat, 
                                                                     vMat_big = vMat_big,
                                                                     d = dist,
                                                                     display_progress=FALSE),
                 "all" = .ProbabilityMatrixAllNested(vMat = vMat, 
                                                     vMat_big = vMat_big,
                                                     d = dist,
                                                     display_progress=FALSE),
                 .ProbabilityMatrixAllNested(vMat = vMat,
                                             vMat_big = vMat_big,
                                             d = dist,
                                             display_progress=FALSE)
          )
          
        }
        
      }
      
      if( priorB == TRUE ) {
        
        switch_angle <- function(angle) {
          
          switch(angle,
                 "horizontal" = .ProbabilityMatrixHorizontalPost(vMat = vMat, 
                                                                 x = blockrama, d = dist, 
                                                                 nrp = nrp, nrp_big = nrp_big,
                                                                 display_progress=FALSE),
                 "vertical" = .ProbabilityMatrixVerticalPost(vMat = vMat, 
                                                             x = blockrama, d = dist, 
                                                             nrp = nrp, nrp_big = nrp_big,
                                                             display_progress=FALSE),
                 "diagonal45" = .ProbabilityMatrixDiagonal45Post(vMat = vMat, 
                                                                 x = blockrama, d = dist, 
                                                                 nrp = nrp, nrp_big = nrp_big,
                                                                 display_progress=FALSE),
                 "diagonal135" = .ProbabilityMatrixDiagonal135Post(vMat = vMat, 
                                                                   x = blockrama, d = dist, 
                                                                   nrp = nrp, nrp_big = nrp_big,
                                                                   display_progress=FALSE),
                 "all" = .ProbabilityMatrixAllPost(vMat = vMat, 
                                                   x = blockrama, d = dist, 
                                                   nrp = nrp, nrp_big = nrp_big,
                                                   display_progress=FALSE),
                 .ProbabilityMatrixAllPost(vMat = vMat, 
                                           x = blockrama, d = dist, 
                                           nrp = nrp, nrp_big = nrp_big,
                                           display_progress=FALSE)
          )
          
        }
        
        
      }
      
      out <- raster::crop(x, raster::extent(x, 1, 1, 1, 1))
      
      if (angle %in% c("horizontal", "vertical")) {
        nrp <- 2*wslI*(wslI - dist)
      }
      
      if (angle %in% c("diagonal45", "diagonal135")) {
        nrp <- (wslI - dist) * 2 *(wslI - dist)
      }
      
      if (angle == "all") {
        nrp <- 4 * (wslI - dist) * (2 * wslI - dist)
      }
      
      wmx <- matrix(1, ncol = dimB[2], nrow = dimB[1])
      
      wsl <- ifelse(is.null(wslO), wslI, wslO)
      
      wmx <- .G(Mat = wmx, overlap = oLap, edge = wsl-1)
      
      num <- out
      denom <- out
      
      stepsize <- dimB[1]-oLap
      rowtimes <- (nrow(x) - dimB[1]) / (dimB[1] - oLap) + 1
      times <- rep(stepsize, rowtimes - 1)
      RowIndex <- c(1, cumsum(times) + 1)
      
      stepsize <- dimB[2]-oLap
      coltimes <- (ncol(x) - dimB[2]) / (dimB[2] - oLap) + 1
      times <- rep(stepsize, coltimes - 1)
      ColIndex <- c(1, cumsum(times) + 1)
      
      # pb = txtProgressBar(min = 0, max = length(RowIndex), initial = 0, style = 3)
      
      for (i in 1:length(RowIndex)) {
        for (j in 1:length(ColIndex)) {
          block <- raster::getValuesBlock(x = x, row = RowIndex[i], nrows = dimB[1],
                                          col = ColIndex[j], ncols = dimB[2], format = "matrix")
          blockra <- raster::raster(block)
          
        if(priorB == TRUE){
          
          if (angle == "horizontal") {
            nrp_big <- 2*nrow(blockra)*(ncol(blockra) - dist)
          }
          if (angle == "vertical") {
            nrp_big <- 2*ncol(blockra)*(nrow(blockra) - dist)
          }
          if (angle %in% c("diagonal45", "diagonal135")) {
            nrp_big <- (nrow(blockra) - dist) * 2 *(ncol(blockra) - dist)
          }
          if (angle == "all") {
            nrp_big <- 4 * (nrow(blockra) - dist) * (2 * ncol(blockra) - dist)
          }
          
        }
          raster::extent(blockra) <- raster::extent(x, RowIndex[i], RowIndex[i] + dimB[1] - 1,
                                                    ColIndex[j], ColIndex[j] + dimB[2] - 1)
          raster::crs(blockra) <- raster::crs(x)
          raster::res(blockra) <- raster::res(x)
          
          if(!is.null(wslO)) {
            suppressMessages( vMat_big <- .getValuesWindow(blockra, wsl = wslO, padValue = padValue,
                                                           aroundTheGlobe = aroundTheGlobe) )
          }
          
          suppressMessages( vMat <- .getValuesWindow(blockra, wsl = wslI, padValue = padValue,
                                                     aroundTheGlobe = aroundTheGlobe) )
          Hetx <- vMat
          
          suppressWarnings(
            if( identical(na.handling, na.omit) && anyNA(raster::values(x)) ) {
              narm <- 1
            }
            else{
              narm <- 0
            }
          )
          
          blockrama <- matrix(raster::values(blockra), nrow(blockra), ncol(blockra), byrow = TRUE)
          
          SpatMat <- switch_angle(angle)
          
          sdiv <- do.call(fun, list(rank = rank, Hetx = Hetx, SpatMat = SpatMat, delta = delta,
                                    nrp = nrp, narm = narm, display_progress = FALSE))
          
          # multiply each block with the spatial weights matrix
          wdiv <- sdiv * wmx
          # rasterize weighted blocks
          wdiv <- raster::raster(wdiv)
          raster::extent(wdiv) <- raster::extent(x, RowIndex[i], RowIndex[i] + dimB[1] - 1,
                                                 ColIndex[j], ColIndex[j] + dimB[2] - 1)
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
        }
        # setTxtProgressBar(pb,i)
      }
      
      # close(pb)
      
      out <- num/denom
      
      return(out)
      
    } ## END of if statement for BLOCK
    
    
  } else {
    stop("Cannot proccess in memory.")
  }
  
  }
  
  if (filename != "") {
    out <- raster::writeRaster(out, filename)
  }
  
  return(out)
  
}
  

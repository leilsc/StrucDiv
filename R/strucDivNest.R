#' @name strucDivNest
#' @rdname strucDivNest
#' @title Quantify Spatial Structural Diversity Across Scales in an Arbitrary Raster Layer 
#' @description
#' This is a wrapper function that returns a 'spatial structural diversity map' 
#' as a raster layer. 
#' 'Spatial structural diversity' will hereafter be used synonymous to 'structural diversity'.
#' Pixels are considered as pairs in user-specified distances and angles. 
#' Angles include horizontal and vertical direction, and the diagonals at 45° and 135°. 
#' The direction-invariant version considers all angles. 
#' Spatial structural diversity is quantified based on the probabilities that pixel values 
#' are arranged in the specified way (distance and angle). The \code{\link{strucDiv}} function employs empirical probabilities of pixel value co-occurrence. 
#' The \code{\link{strucDivNest}} function combines information from two different scales with an empirical Bayesian approach and a Beta-Binomial model.
#' Two scales are nested inside each other - a larger, outer scale and a smaller, inner scale. 
#' Three different nesting schemes are available, whereby the inner scale is always a moving window.
#' The outer scale can either be another mowing window, a block, or the domain (i.e. the input raster).
#' The outer scale is used as prior information, and the inner scale serves as likelihood to estimate posterior probabilities
#' of pixel value co-occurrences. 
#' In the Beta-Binomial model both, the prior and the posterior follow a beta distribution, and the likelihood follows a conditional 
#' binomial distribution.
#' Posterior probabilities are estimated with mean estimates.
#' Structural diversity is quantified based on these posterior probabilities.
#' Structural diversity metrics are calculated on every element of the GLCM, 
#' and their sum is assigned to the center pixel of the moving window. 
#' The final map is called a '(spatial) structural diversity map' and is returned as a raster layer. 
#' The output map represents structural diversity, quantified across different spatial scales, which are defined
#' by the outer scale and the inner scale.
#' @param x raster layer. Input raster layer for which 
#' horizontal structural diversity should be calculated.
#' @param wslI uneven integer. The window side length of the inner scale, 
#' \code{wslI} x \code{wslI} defines the size of the inner moving window.
#' The window must be smaller than the dimensions of the input raster and smaller than the outer scale. 
#' Default is NULL, in which case no nesting is used.
#' @param wslO uneven integer.  The window side length of the outer scale, 
#' \code{wslO} x \code{wslO} defines the size of the outer moving window.
#' The window must be smaller than the dimensions of the input raster and larger than the inner scale (i.e. \code{wslI}). 
#' Defaults to \code{NULL}, in which case no nesting is used.
#' @param dimB a vector of length 2 or logical. This defines the block size (number of rows, number of columns). 
#' The domain (i.e. the input raster) is divided into equal size, overlapping blocks. 
#' Each block provides prior information for the inner window, which moves inside each block. 
#' Structural diversity is quantified in each block.
#' Blocks are merged together in a spatially weighted manner, using linear weights.
#' Defaults to \code{FALSE}, in which case no blocks are used.
#' @param oLap integer. This defines the size of overlap between the blocks. The overlap must
#' be at least \code{wslI-1} or bigger. Blocks can overlap by a maximum of half the 
#' rows of blocks in row-direction, and by half the columns of blocks in column-direction.
#' If oLap is not specified, the minimum overlap is used.
#' Defaults to \code{NULL} in which case no blocks are used.
#' @param priorB logical. Should blocks be used for prior information? 
#' If \code{priorB = TRUE}, then the spatial structure in a block serves 
#' as prior information for the inner scale. 
#' If \code{priorB = FALSE}, then the blocks are only used to increase speed through 
#' parallelization, not for prior information.
#' Defaults to \code{FALSE}.
#' @param domain logical. Should the domain (i.e. the input raster) be used for prior information?
#' If  \code{domain = TRUE}, then it is used as prior for all inner moving windows.
#' Defaults to \code{FALSE}.
#' @param dist integer. The distance between two pixels that should be considered as a pair, 
#' defaults to \code{dist = 1} (direct neighbors).
#' @param angle string. The angle on which pixels should be considered as pairs. 
#' Takes 5 options: \code{"horizontal"}, \code{"vertical"}, \code{"diagonal45"}, \code{"diagonal135"}, \code{"all"}. 
#' The direction-invariant version is \code{"all"}, which considers all of the 4 angles. Defaults to \code{"all"}.
#' @param rank logical. Should pixel values be replaced with ranks in each GLCM? Defaults to \code{FALSE}.
#' @param fun function, the structural diversity metric. Takes one of the following: \code{entropy},
#' \code{entropyNorm}, \code{contrast}, \code{dissimilarity}, or \code{homogeneity}. 
#' Structural diversity entropy is \code{entropy} with different \code{delta} parameters. Shannon entropy is employed, when \code{delta = 0}. 
#' Shannon entropy has a scale-dependent maximum when \code{\link{strucDiv}} is used, but this maximum may be violated in \code{\link{strucDivNest}}, 
#' when information from different scales is combined, depending on the posterior probabilities of pixel value co-occurrences.
#' Additionally, the value gradient is considered when \code{delta = 1} or \code{delta = 2}. 
#' The values of structural diversity entropy with \code{delta = 1} or \code{delta = 2} are not restricted and depend on the values of the input raster.
#' the metric \code{entropyNorm} is Shannon entropy normalized over maximum entropy, which depends on the size of the moving window when no scales are nested. 
#' When information from different scales is combined in \code{\link{strucDivNest}}, the metric \code{entropyNorm} may be larger than 1, 
#' depending on the posterior probabilities of pixel value co-occurrences.
#' The metrics \code{contrast} and \code{dissimilarity} consider the value gradient, their values are not restricted and depend on the values of the input raster.
#' The metric \code{homogeneity} quantifies the closeness of empirical probabilities to the diagonal and ranges between 0 and 1 when scales are not nested. 
#' When information from different scales is combined in \code{\link{strucDivNest}}, the metric \code{homogeneity} may be larger than 1, 
#' depending on the posterior probabilities of pixel value co-occurrences.
#' @param delta numeric, takes three options: \code{0}, \code{1}, or \code{2}. 
#' The parameter \code{delta} is the difference weight, 
#' it defines how the differences between pixel values within a pixel pair should be weighted.  
#' If \code{rank = TRUE}, delta defines how the differences between ranks should be weighted.  
#' Defaults to \code{0} (no weight). Set \code{delta = 1} for absolute weights, 
#' or \code{delta = 2} for square weights. 
#' The \code{delta} parameter can only be set when the metric \code{entropy} is used. 
#' the metric \code{dissimilarity} automatically employs \code{delta = 1}, and \code{contrast} employs \code{delta = 2}.
#' @param na.handling \code{na.omit} or \code{na.pass}. 
#' If \code{na.handling = na.omit}, NAs are ignored and structural diversity metrics are calculated with less values. 
#' If \code{na.handling = na.pass} and if there is at least one missing value inside the moving window,
#' an NA is assigned to the center pixel. Therefore, the diversity map will contain more 
#' NAs than the input raster layer.
#' Defaults to \code{na.pass}.
#' @param padValue numeric or \code{NA}. The value of the padded cells at the edges of the input raster. 
#' Defaults to \code{NA}.
#' @param aroundTheGlobe logical. If the input raster goes around the whole globe, 
#' set \code{aroundTheGlobe = TRUE}, and the input raster will be 'glued together' from both sides
#' to calculate structural diversity without edge effects.
#' Defaults to \code{FALSE}.
#' @param ncores integer. The number of cores the computation will be parallelized on.
#' Parallelization is only available with blocks - both when they are used as prior, 
#' and when they are simply used to cut the image.
#' @param display_progress logical. If \code{display_progress = TRUE}, a progress bar will be visible.
#' @param filename character. If the output raster should be written to a file, define file name (optional).
#' @param ... possible further arguments.
#' @importFrom raster raster
#' @importFrom raster setValues
#' @importFrom raster values
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit na.pass
#' @importFrom glue trim
#' @details The memory requirement of the function is determined 
#' by \code{raster::canProcessInMemory()}. 
#' If the raster file cannot be processed in memory, its size needs to be reduced before \code{\link{strucDivNest}} can be used. 
#' @return The output is a (spatial) structural diversity map, returned as a raster layer.
#' If the outer scale is a moving window or the domain, then the output raster has the same dimensions as the input raster.
#' If the outer scale is a block, then the output raster may be smaller than the input raster 
#' because if there are edges that do not fit inside the blocks, they are cut off.
#' When \code{na.handling = na.pass}, then the output map will have an NA-edge of 0.5*(\code{wslO}-1), 
#' and it will contain more missing values than the input raster.
#' The output represents spatial structural diversity quantified across different spatial scale, which are defined by the 
#' size of the inner and the outer scale.
#' @examples
#' # Construct a small raster file containing realizations of normal random variables:
#' a <- raster::raster(matrix(rnorm(648), 18, 36))
#' raster::plot(a)
#' # Calculate structural diversity entropy with delta = 2
#' sde_a <- strucDivNest(a, wslI = 3, wslO = 7, fun = entropy, delta = 2, na.handling = na.omit, 
#'     rank = FALSE)
#' raster::plot(sde_a)
#' 
#' # Calculate structural diversity entropy with delta = 1
#' b <- raster::raster(matrix(rnorm(2500), 50, 50))
#' raster::plot(b)
#' sde_b <- strucDivNest(b, wslI = 3, dimB = c(10, 10), oLap = 4, priorB = TRUE, fun = entropy, 
#'     delta = 1, na.handling = na.pass, rank = FALSE)
#' raster::plot(sde_b)
#' 
#' # Calculate contrast on NDVI data, use block nesting scheme 
#' ndvi <- raster::raster(ndvi)
#' contrastNest_ndvi <- strucDivNest(ndvi, wslI = 9, dimB = c(50, 50), oLap = 20, priorB = TRUE,
#'      fun = contrast, na.handling = na.pass, rank = FALSE)
#' raster::plot(contrastNest_ndvi)
#' 
#' # Calculate entropy on NDVI data binned to 15 gray levels, use domain nesting scheme 
#' ndvi.15gl <- raster::raster(ndvi.15gl)
#' entropyNest_ndvi15 <- strucDivNest(ndvi.15gl, wslI = 5, domain = TRUE, fun = entropy, 
#'     na.handling = na.pass, rank = FALSE)
#' raster::plot(entropyNest_ndvi15)
#' 
#' @export

strucDivNest <- function(x, wslI = NULL, wslO = NULL, dimB = FALSE, oLap = NULL, 
                         priorB = FALSE, domain = FALSE, dist = 1, angle = "all",
                         rank = FALSE, fun, delta = 0, 
                         na.handling = na.pass, padValue = NA, 
                         aroundTheGlobe = FALSE, 
                         ncores = 1,
                         display_progress = TRUE, 
                         filename = "", ...) {
  
  #browser()
  
  dotArgs <- list(...)
  
  ## General warnings and errors!
  
  if(!is.logical(aroundTheGlobe)){
    stop("aroundTheGlobe must be either TRUE, or FALSE.")
  }
  
  if(isTRUE(aroundTheGlobe) & !isTRUE(.isGlobalLonLat(x))){
    warning("The raster image does not go around the globe.")
  }
  
  if ( identical(na.handling, na.pass) && anyNA(raster::values(x)) ) {
    warning("Raster layer contains missing values. Wherever there are missing values,
              an NA will be returned. if you want to proceed without NAs, 
              set na.handling = na.omit.")
  }
  
  if ( identical(na.handling, na.pass) && anyNA(raster::values(x)) && domain == TRUE ) {
    return(NA)
    warning("The domain contains missing values. You may want to consider the argument 'na.handling = na.ignore'.")
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
    
    minOL <- (wslI-1)
    maxOL <- 0.5* min(dimB)
    oLap <- floor((maxOL-minOL)/2 + minOL)
    warning("oLap is calculated to lie in the middle of its min and max values 
            because it was not specified.")
    
  }
  
  if(!is.null(oLap) & is.numeric(dimB)){
    
    if (oLap < (wslI-1)) {
      stop("Overlap is too small. The overlap must be oLap >= wslI-1.")
      
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
  
  if(priorB == FALSE & is.numeric(dimB[1])) {
    
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
                                                                   display_progress=display_progress),
                 "vertical" = .ProbabilityMatrixVerticalNested(vMat = vMat, 
                                                               vMat_big = vMat_big, 
                                                               d = dist,
                                                               display_progress=display_progress),
                 "diagonal45" = .ProbabilityMatrixDiagonal45Nested(vMat = vMat, 
                                                                   vMat_big = vMat_big,
                                                                   d = dist,
                                                                   display_progress=display_progress),
                 "diagonal135" = .ProbabilityMatrixDiagonal135Nested(vMat = vMat, 
                                                                     vMat_big = vMat_big,
                                                                     d = dist,
                                                                     display_progress=display_progress),
                 "all" = .ProbabilityMatrixAllNested(vMat = vMat, 
                                                     vMat_big = vMat_big,
                                                     d = dist,
                                                     display_progress=display_progress),
                 .ProbabilityMatrixAllNested(vMat = vMat,
                                             vMat_big = vMat_big,
                                             d = dist,
                                             display_progress=display_progress)
          )
          
        }
        
        vMat <- getValuesWindow(x, wsl = wslI, padValue = padValue, 
                                aroundTheGlobe = aroundTheGlobe)
        
        suppressMessages( vMat_big <- getValuesWindow(x, wsl = wslO, padValue = padValue, 
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
        
        v <- do.call(fun, list(rank = rank, Hetx = Hetx, vMat_big = vMat_big, SpatMat = SpatMat, delta = delta,
                               nrp = nrp, narm = narm, display_progress = display_progress))
        
        out <- raster::setValues(out, v)
        
      }
      
      if(domain == TRUE){  ## USE THE WHOLE IMAGE AS A PRIOR
        
        vMat <- getValuesWindow(x, wsl = wslI, padValue = padValue, 
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
                                                                 display_progress=display_progress),
                 "vertical" = .ProbabilityMatrixVerticalPost(vMat = vMat, 
                                                             x = xMat, d = dist, 
                                                             nrp = nrp, nrp_big = nrp_big,
                                                             display_progress=display_progress),
                 "diagonal45" = .ProbabilityMatrixDiagonal45Post(vMat = vMat, 
                                                                 x = xMat, d = dist, 
                                                                 nrp = nrp, nrp_big = nrp_big,
                                                                 display_progress=display_progress),
                 "diagonal135" = .ProbabilityMatrixDiagonal135Post(vMat = vMat, 
                                                                   x = xMat, d = dist, 
                                                                   nrp = nrp, nrp_big = nrp_big,
                                                                   display_progress=display_progress),
                 "all" = .ProbabilityMatrixAllPost(vMat = vMat, 
                                                   x = xMat, d = dist, 
                                                   nrp = nrp, nrp_big = nrp_big,
                                                   display_progress=display_progress),
                 .ProbabilityMatrixAllPost(vMat = vMat, 
                                           x = xMat, d = dist, 
                                           nrp = nrp, nrp_big = nrp_big,
                                           display_progress=display_progress)
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
          nrp_big <- 2*((nrow(xMat)-dist)*(2*ncol(xMat)-dist)+(ncol(xMat)-dist)*(2*nrow(xMat)-dist))
        }
        
        SpatMat <- switch_angle(angle)
        
        v <- do.call(fun, list(rank = rank, Hetx = Hetx, SpatMat = SpatMat, delta = delta,
                               nrp = nrp, narm = narm, display_progress = display_progress))
        
        out <- raster::setValues(out, v)
        
      }
      
      if( (dimB != FALSE)[1] ){
        
        if( priorB == FALSE & is.null(wslO)){
          
          switch_angle <- function(angle) {
            
            switch(angle,
                   "horizontal" = .ProbabilityMatrixHorizontalDynamic(vMat = vMat, d = dist, narm = narm,
                                                                      display_progress = display_progress),
                   "vertical" = .ProbabilityMatrixVerticalDynamic(vMat = vMat, d = dist, narm = narm,
                                                                  display_progress = display_progress),
                   "diagonal45" = .ProbabilityMatrixDiagonal45Dynamic(vMat = vMat, d = dist, narm = narm,
                                                                      display_progress = display_progress),
                   "diagonal135" = .ProbabilityMatrixDiagonal135Dynamic(vMat = vMat, d = dist, narm = narm,
                                                                        display_progress = display_progress),
                   "all" = .ProbabilityMatrixAllDynamic(vMat = vMat, d = dist, narm = narm,
                                                        display_progress = display_progress),
                   .ProbabilityMatrixAllDynamic(vMat = vMat, d = dist, narm = narm,
                                                       display_progress = display_progress)
            )
            
          }
          
        }
        if( priorB == FALSE & !is.null(wslO)) {
          
          switch_angle <- function(angle) {
            
            switch(angle,
                   "horizontal" = .ProbabilityMatrixHorizontalNested(vMat = vMat, 
                                                                     vMat_big = vMat_big,
                                                                     d = dist,
                                                                     display_progress=display_progress),
                   "vertical" = .ProbabilityMatrixVerticalNested(vMat = vMat, 
                                                                 vMat_big = vMat_big, 
                                                                 d = dist,
                                                                 display_progress=display_progress),
                   "diagonal45" = .ProbabilityMatrixDiagonal45Nested(vMat = vMat, 
                                                                     vMat_big = vMat_big,
                                                                     d = dist,
                                                                     display_progress=display_progress),
                   "diagonal135" = .ProbabilityMatrixDiagonal135Nested(vMat = vMat, 
                                                                       vMat_big = vMat_big,
                                                                       d = dist,
                                                                       display_progress=display_progress),
                   "all" = .ProbabilityMatrixAllNested(vMat = vMat, 
                                                       vMat_big = vMat_big,
                                                       d = dist,
                                                       display_progress=display_progress),
                   .ProbabilityMatrixAllNested(vMat = vMat,
                                               vMat_big = vMat_big,
                                               d = dist,
                                               display_progress=display_progress)
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
        
        wsl <- ifelse(is.null(wslO), wslI, wslO)
        
        wmx <- .WM(nrow = dimB[1], ncol = dimB[2], ul = oLap, nNA = floor(0.5*wsl))
        
        num <- setValues(out, values = NA)
        denom <- setValues(out, values = NA)
        
        stepsize <- dimB[1]-oLap
        rowtimes <- (nrow(x) - dimB[1]) / (dimB[1] - oLap) + 1
        times <- rep(stepsize, rowtimes - 1)
        RowIndex <- c(1, cumsum(times) + 1)
        
        stepsize <- dimB[2]-oLap
        coltimes <- (ncol(x) - dimB[2]) / (dimB[2] - oLap) + 1
        times <- rep(stepsize, coltimes - 1)
        ColIndex <- c(1, cumsum(times) + 1)
        
        rows <- length(RowIndex)
        
        cl <- parallel::makePSOCKcluster(ncores)
        registerDoSNOW(cl)
        
        if(display_progress){
          pb <- txtProgressBar(max=rows, style=3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress=progress)
          suppressWarnings(
            out <- foreach( row.num = 1:rows, .packages = c("raster", "StrucDiv2"),
                            .options.snow=opts) %dopar% {
                              
                              for (j in 1:length(ColIndex)) {
                                block <- raster::getValuesBlock(x = x, row = RowIndex[row.num], nrows = dimB[1],
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
                                    nrp_big <- 2*((nrow(blockra)-dist)*(2*ncol(blockra)-dist)+(ncol(blockra)-dist)*(2*nrow(blockra)-dist))
                                  }
                                  
                                }
                                raster::extent(blockra) <- raster::extent(x, RowIndex[row.num], RowIndex[row.num] + dimB[1] - 1,
                                                                          ColIndex[j], ColIndex[j] + dimB[2] - 1)
                                raster::crs(blockra) <- raster::crs(x)
                                raster::res(blockra) <- raster::res(x)
                                
                                vMat_big = NULL
                                
                                if(!is.null(wslO)) {
                                  suppressMessages( vMat_big <- StrucDiv2::getValuesWindow(blockra, wsl = wslO, padValue = NA,
                                                                                           aroundTheGlobe = aroundTheGlobe) )
                                }
                                
                                vMat <- StrucDiv2::getValuesWindow(blockra, wsl = wslI, padValue = NA, 
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
                                
                                blockrama <- matrix(raster::values(blockra), nrow(blockra), ncol(blockra), byrow = TRUE)
                                
                                SpatMat <- switch_angle(angle)
                                
                                sdiv <- do.call(fun, list(rank = rank, Hetx = Hetx, vMat_big = vMat_big, SpatMat = SpatMat, delta = delta,
                                                          nrp = nrp, narm = narm, display_progress = FALSE))
                                
                                if(anyNA(raster::values(blockra)) && priorB == TRUE && narm == 0) {
                                  sdiv <- NA
                                  warning("At least one block contains missing values. You may want to consider the argument 'na.handling = na.omit'.")
                                }
                                
                                # multiply each block with the spatial weights matrix
                                sdiv <- matrix(sdiv, dimB[1], dimB[2], byrow = TRUE)
                                wdiv <- sdiv * wmx
                                
                                # rasterize weighted blocks
                                wdiv <- raster::raster(wdiv)
                                wdiv <- setValues(blockra, values = values(wdiv))
                                
                                wmr <- raster(wmx) # rasterize weights matrix so we can create denominator layer
                                WMRext <- setValues(blockra, values = values(wmr)) # each block gets a weights matrix layer
                                # that has the same extent, i.e. is in the same place as the block
                                
                                # create numerator - take the sum of overlapping weighted blocks (i.e. take the sum where they overlap)
                                num <- raster::mosaic(num, wdiv, fun = sum)
                                # create denominator - take the sum of overlapping weights (i.e. take the sum where they overlap)
                                denom <- raster::mosaic(denom, WMRext, fun = sum)
                                
                                # create final raster layer
                              }
                              
                              list(num = num, denom = denom)
                              
                            }
          )
          close(pb)
        }
        
        else{
          suppressWarnings(
            out <- foreach( row.num = 1:rows, .packages = c("raster", "StrucDiv2") ) %dopar% {
                              
                              for (j in 1:length(ColIndex)) {
                                block <- raster::getValuesBlock(x = x, row = RowIndex[row.num], nrows = dimB[1],
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
                                    nrp_big <- 2*((nrow(blockra)-dist)*(2*ncol(blockra)-dist)+(ncol(blockra)-dist)*(2*nrow(blockra)-dist))
                                  }
                                  
                                }
                                raster::extent(blockra) <- raster::extent(x, RowIndex[row.num], RowIndex[row.num] + dimB[1] - 1,
                                                                          ColIndex[j], ColIndex[j] + dimB[2] - 1)
                                raster::crs(blockra) <- raster::crs(x)
                                raster::res(blockra) <- raster::res(x)
                                
                                vMat_big = NULL
                                
                                if(!is.null(wslO)) {
                                  suppressMessages( vMat_big <- StrucDiv2::getValuesWindow(blockra, wsl = wslO, padValue = NA,
                                                                                           aroundTheGlobe = aroundTheGlobe) )
                                }
                                
                                vMat <- StrucDiv2::getValuesWindow(blockra, wsl = wslI, padValue = NA, 
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
                                
                                blockrama <- matrix(raster::values(blockra), nrow(blockra), ncol(blockra), byrow = TRUE)
                                
                                SpatMat <- switch_angle(angle)
                                
                                sdiv <- do.call(fun, list(rank = rank, Hetx = Hetx, vMat_big = vMat_big, SpatMat = SpatMat, delta = delta,
                                                          nrp = nrp, narm = narm, display_progress = FALSE))
                                
                                if(anyNA(raster::values(blockra)) && priorB == TRUE && narm == 0) {
                                  sdiv <- NA
                                  warning("At least one block contains missing values. You may want to consider the argument 'na.handling = na.omit'.")
                                }
                                
                                # multiply each block with the spatial weights matrix
                                sdiv <- matrix(sdiv, dimB[1], dimB[2], byrow = TRUE)
                                wdiv <- sdiv * wmx
                                
                                # rasterize weighted blocks
                                wdiv <- raster::raster(wdiv)
                                wdiv <- setValues(blockra, values = values(wdiv))
                                
                                wmr <- raster(wmx) # rasterize weights matrix so we can create denominator layer
                                WMRext <- setValues(blockra, values = values(wmr)) # each block gets a weights matrix layer
                                # that has the same extent, i.e. is in the same place as the block
                                
                                # create numerator - take the sum of overlapping weighted blocks (i.e. take the sum where they overlap)
                                num <- raster::mosaic(num, wdiv, fun = sum)
                                # create denominator - take the sum of overlapping weights (i.e. take the sum where they overlap)
                                denom <- raster::mosaic(denom, WMRext, fun = sum)
                                
                                # create final raster layer
                              }
                              
                              list(num = num, denom = denom)
                              
                            }
          )
          
        }
        
        parallel::stopCluster(cl)
        
        num <- raster::mosaic(out[[1]]$num, out[[2]]$num, fun = sum)
        denom <- raster::mosaic(out[[1]]$denom, out[[2]]$denom, fun = sum)
        
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


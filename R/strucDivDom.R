#' @name strucDivDom
#' @rdname strucDivDom
#' @title Returns the structural diversity value, the gray level co-occurrence matrix (GLCM) and the structural diversity 
#' matrix of the domain. 
#' @description The function \code{strucDivDom} returns the spatial, i.e. horizontal, structural diversity value for the domain (i.e. the input raster). 
#' 'Spatial structural diversity' will hereafter be used synonymous to 'structural diversity'.
#' The function also returns the gray level co-occurrence matrix (GLCM) and the structural diversity matrix of 
#' the domain. Structural diversity is calculated on every element of the GLCM,
#' which generates the structural diversity matrix.
#' @param x raster layer. Input raster layer for which 
#' structural diversity should be calculated.
#' @param dist integer. The distance between two pixels that should be considered as a pair, 
#' defaults to \code{dist = 1} (direct neighbors).
#' @param angle string. The angle on which pixels should be considered as pairs. 
#' Takes 5 options: \code{"horizontal"}, \code{"vertical"}, \code{"diagonal45"}, \code{"diagonal135"}, \code{"all"}. 
#' The direction-invariant version is \code{"all"}, which considers all of the 4 angles. Defaults to \code{"all"}.
#' @param rank logical. Should pixel values be replaced with ranks in each GLCM? Defaults to \code{FALSE}.
#' @param fun function, the structural diversity metric. Takes one of the following: \code{entropyDom},
#' \code{entropyNormDom}, \code{contrastDom}, \code{dissimilarityDom}, or \code{homogeneityDom}. 
#' Structural diversity entropy is \code{entropyDom} with different \code{delta} parameters. 
#' Shannon entropy is employed when \code{delta = 0}. 
#' Shannon entropy has a scale-dependent maximum. Scale-dependent means dependent on the extent of he area within which
#' structural diversity is quantified, because this area defines he total number of pixel pairs.
#' The metric \code{entropyNormDom} is Shannon entropy normalized over scale-dependent maximum entropy. 
#' Additionally, the value gradient is considered with \code{delta = 1} and \code{delta = 2}. 
#' The values of structural diversity entropy with \code{delta = 1} or \code{delta = 2} are not restricted 
#' and depend on the values of the input raster.
#' The metrics \code{contrastDom} and \code{dissimilarityDom} consider the value gradient,
#'  their values are not restricted and depend on the values of the input raster.
#' The metric \code{homogeneityDom} quantifies the closeness of empirical probabilities to the diagonal and ranges between 0 and 1.
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
#' If \code{na.handling = na.pass} and if there is at least one missing value in the domain, an NA will be returned.
#' Defaults to \code{na.pass}.
#' @importFrom raster raster
#' @importFrom raster setValues
#' @importFrom raster values
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit na.pass
#' @importFrom glue trim
#' @details The memory requirement of the function is determined 
#' by \code{raster::canProcessInMemory()}. 
#' If the raster file cannot be processed in memory, its size needs to be 
#' reduced before \code{\link{strucDivDom}} can be used. 
#' @return The output is a list containing the structural diversity value of the domain, which can be accessed with $div. 
#' the list also contains the gray level co-occurrence matrix ($GLCM) and the structural diversity 
#' matrix ($divMat) of the domain. 
#' @examples
#' # Calculate entropy on NDVI data binned to 15 gray levels
#' ndvi15 <- raster::raster(ndvi.15gl)
#' ndvi15Dom <- strucDivDom(ndvi15, fun = entropyDom)
#' ndvi15GLCM <- ndvi15Dom$GLCM
#' ndvi15Div <- ndvi15Dom$div
#' 
#' @export

strucDivDom <- function(x, dist = 1, angle = "all",
                        rank = FALSE, fun, delta = 0, ...) {
  
  #browser()
  
  dotArgs <- list(...)
  
  ## General warnings and errors!
  
  out <- raster::raster(x)
  
  ## Check if the raster has values inside.
  
  stopifnot(raster::hasValues(x))
  
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
  
  if (raster::canProcessInMemory(x)) {
    
    vMat <- as.matrix(values(x), nrow(x), ncol(x), byrow=TRUE)
    Hetx <- vMat
    values <- sort(unique(vMat))
    
    switch_angle <- function(angle) {
      
      switch(angle,
             "horizontal" = .ProbabilityMatrixHorizontal(xMat = vMat, d = dist, 
                                                         Values = values),
             "vertical" = .ProbabilityMatrixVertical(xMat = vMat, d = dist, 
                                                     Values = values),
             "diagonal45" = .ProbabilityMatrixDiagonal45(xMat = vMat, d = dist, 
                                                        Values = values),
             "diagonal135" = .ProbabilityMatrixDiagonal135(xMat = vMat, d = dist, 
                                                           Values = values),
             "all" = .ProbabilityMatrixAll(xMat = vMat, d = dist, 
                                           Values = values),
             .ProbabilityMatrixAll(xMat = vMat, d = dist, 
                                          Values = values)
      )
      
    }
    
    SpatMat <- switch_angle(angle)
    rownames(SpatMat) = values
    colnames(SpatMat) = values
    
    v <- do.call(fun, list(rank = rank, xVal = values, 
                           PMat = SpatMat, delta = delta))

  }
  
  out <- list(GLCM = SpatMat, metric_matrix = v, div = sum(v))

  return(out)
  
}


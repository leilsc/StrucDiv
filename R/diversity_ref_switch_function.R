#' @name DiversityDom
#' @rdname DiversityDom
#' @title Spatial Structural Diversity Metrics
#' @description
#' The functions \code{entropyDom} , \code{entropyNormDom}, \code{contrastDom}, \code{dissimilarityDom} and \code{homogeneityDom}
#' are the spatial structural diversity metrics used in the default configurations of \code{\link{strucDiv}} and \code{\link{strucDivNest}}. 
#' Structural diversity entropy is \code{entropyDom} with different \code{delta} parameters. Shannon entropy is employed, when \code{delta = 0}. 
#' Shannon entropy has a window-dependent maximum when \code{\link{strucDiv}} is used, which may be violated when \code{\link{strucDivNest}} is used, 
#' depending on the posterior probabilities of pixel value co-occurrences.
#' Additionally, the value gradient is considered when \code{delta = 1} or \code{delta = 2}. 
#' The values of structural diversity entropy with \code{delta = 1} or \code{delta = 2} are not restricted and depend on the values of the input raster.
#' the metric \code{entropyNormDom} is Shannon entropy normalized over maximum entropy, which depends on the size of the moving window when no nesting is used. 
#' The metric \code{entropyNormDom} ranges between 0 and 1, when \code{\link{strucDiv}} is used, but may be larger than 1 when \code{\link{strucDivNest}} is used, 
#' depending on the posterior probabilities of pixel value co-occurrences.
#' The metrics \code{contrastDom} and \code{dissimilarityDom} consider the value gradient, their values are not restricted and depend on the values of the input raster.
#' The metric \code{homogeneityDom} quantifies the closeness of empirical probabilities to the diagonal and ranges between 0 and 1 when \code{\link{strucDiv}} is used, 
#' but may be larger than 1 when \code{\link{strucDivNest}} is used, depending on the posterior probabilities of pixel value co-occurrences.
#' @param rank logical. Should values be replaced with ranks in each co-occurrence 
#' matrix (GLCM)? Defaults to \code{FALSE}.
#' @param delta numeric, takes 3 options: \code{0}, \code{1}, or \code{2}. 
#' The parameter \code{delta} is the difference weight parameter, 
#' it defines how the differences between pixel values within a pixel pair should be weighted.  
#' If \code{rank = TRUE}, delta defines how the differences between ranks should be weighted.  
#' The default value is \code{0} (no weight). Set \code{delta = 1} for absolute weights, 
#' or \code{delta = 2} for squared weights. 
#' The \code{delta} parameter can only be set when the metric \code{entropy} is used. 
#' The metric \code{dissimilarity} automatically employs \code{delta = 1}, and \code{contrast} employs \code{delta = 2}.
#' @param PMat the GLCM that is returned by an internal function
#' to the \code{\link{strucDiv}} and \code{\link{strucDivNest}} functions. 
#' @param ... possible further arguments.
#' @details This function is used internally and is called 
#' as an argument to the \code{\link{strucDiv}} and \code{\link{strucDivNest}} functions.
#' @importFrom raster raster
#' @export

homogeneityDom <- function(rank, delta, PMat, xVal) {
  
  switch_function <- function(rank) {
    
    switch(EXPR = paste(rank),
           "TRUE" = .HomogeneityRankRef(PMat = PMat),
           "FALSE" = .HomogeneityValueRef(PMat = PMat, xVal = xVal)
    )
  }
  
  v <- switch_function(rank)
  return(v)
  
}

#' @rdname DiversityDom
#' @export


dissimilarityDom <- function(rank, delta, PMat, xVal) {
  
  switch_function <- function(rank) {
    
    switch(EXPR = paste(rank),
           "TRUE" = .DissimilarityRankRef(PMat = PMat),
           "FALSE" = .DissimilarityValueRef(PMat = PMat, xVal = xVal)
    )
  }
  
  v <- switch_function(rank)
  return(v)
  
}

#' @rdname DiversityDom
#' @export

contrastDom <- function(rank, delta, PMat, xVal) {
  
  switch_function <- function(rank) {
    
    switch(EXPR = paste(rank),
           "TRUE" = .ContrastRankRef(PMat = PMat),
           "FALSE" = .ContrastValueRef(PMat = PMat, xVal = xVal)    )
  }
  
  v <- switch_function(rank)
  return(v)
  
}

#' @rdname DiversityDom
#' @export

entropyDom <- function(rank, delta, PMat, xVal) {
  
  rank_delta <- paste(rank, delta)
  
  
  switch_function <- function(rank_delta) {
    
    switch(EXPR = rank_delta,
           "FALSE 0" = .EntropyRef(PMat = PMat, xVal = xVal),
           "TRUE 0" = .EntropyRef(PMat = PMat, xVal = xVal),
           "TRUE 1" = .WeightedEntropyAbsRankRef(PMat = PMat),
           "FALSE 1" = .WeightedEntropyAbsValueRef(PMat = PMat, xVal = xVal),
           "TRUE 2" = .WeightedEntropySqrRankRef(PMat = PMat),
           "FALSE 2" = .WeightedEntropySqrValueRef(PMat = PMat, xVal = xVal)
            )
  }
  
  v <- switch_function(rank_delta)
  return(v)
  
}








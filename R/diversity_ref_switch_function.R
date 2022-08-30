#' @name DiversityDom
#' @rdname DiversityDom
#' @title Spatial Structural Diversity Metrics
#' @description
#' The functions \code{entropyDom} , \code{entropyNormDom}, \code{contrastDom}, \code{dissimilarityDom} 
#' and \code{homogeneityDom} are the spatial structural diversity metrics used in the default 
#' configurations of \code{\link{strucDivDom}}. 
#' For programming reasons, these metrics have different name endings than the metrics used in the functions
#' \code{\link{strucDiv}} and \code{\link{strucDivNest}}, but they have the same mathematical formulation. 
#' Hence, \code{entropyDom} is specified by the same equation as \code{entropy}, and so forth.
#' Structural diversity entropy is \code{entropyDom} with different \code{delta} parameters. 
#' Shannon entropy is employed, when \code{delta = 0}. 
#' The metric \code{entropyDom} has a scale-dependent maximum. Scale, here, refers to the extent of the domain.
#' The metric \code{entropyNormDom} is Shannon entropy normalized over maximum entropy. 
#' The metric \code{entropyNormDom} ranges between 0 and 1.
#' Additionally, the value gradient is considered with \code{delta = 1} and \code{delta = 2}. 
#' The values of structural diversity entropy with \code{delta = 1} or \code{delta = 2} are not restricted and depend on the values of the input raster.
#' The metric \code{dissimilarityDom} employs \code{delta = 1}, \code{contrastDom} employs \code{delta = 2}. 
#' The values of \code{dissimilarityDom} and \code{contrastDom} are not restricted and depend on the values of the input raster. 
#' The metric \code{homogeneityDom} quantifies the closeness of empirical probabilities to the diagonal and ranges between 0 and 1.
#' @param rank logical. Should values be replaced with ranks in the co-occurrence 
#' matrix (GLCM)? Defaults to \code{FALSE}.
#' @param delta numeric, takes 3 options: \code{0}, \code{1}, or \code{2}. 
#' The parameter \code{delta} is the difference weight parameter, 
#' it defines how the differences between pixel values within a pixel pair should be weighted.  
#' If \code{rank = TRUE}, delta defines how the differences between ranks should be weighted.  
#' The default value is \code{0} (no weight). Set \code{delta = 1} for absolute weights, 
#' or \code{delta = 2} for square weights. 
#' The \code{delta} parameter can only be set when the metric \code{entropyDom} is used. 
#' The metric \code{dissimilarityDom} automatically employs \code{delta = 1}, and \code{contrastDom} employs \code{delta = 2}.
#' @param PMat the GLCM that is returned by an internal function
#' to the \code{\link{strucDivDom}} function. 
#' @param xVal what is xval?
#' @details These functions are used internally and are called 
#' as an argument to the \code{\link{strucDivDom}}.
#' @importFrom raster raster
#' @export

homogeneityDom <- function(rank, delta, PMat, xVal, nrp) {
  
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


dissimilarityDom <- function(rank, delta, PMat, xVal, nrp) {
  
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

contrastDom <- function(rank, delta, PMat, xVal, nrp) {
  
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

entropyDom <- function(rank, delta, PMat, xVal, nrp) {
  
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

#' @rdname Diversity
#' @export

entropyNormDom <- function(rank, delta, PMat, xVal, nrp) {

  v <- .NormalizedEntropyRef( PMat = PMat, xVal = xVal, nrp = nrp )
  return(v)

}






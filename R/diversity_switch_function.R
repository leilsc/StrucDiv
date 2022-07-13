#' @name Diversity
#' @rdname Diversity
#' @title Spatial Structural Diversity Metrics
#' @description
#' The functions \code{entropy} , \code{entropyNorm}, \code{contrast}, \code{dissimilarity} and \code{homogeneity}
#' are the spatial structural diversity metrics used in the default configurations of \code{\link{strucDiv}} and \code{\link{strucDivNest}}. 
#' Structural diversity entropy is \code{entropy} with different \code{delta} parameters. Shannon entropy is employed, when \code{delta = 0}. 
#' Shannon entropy has a window-dependent maximum when \code{\link{strucDiv}} is used, which may be violated when \code{\link{strucDivNest}} is used, 
#' depending on the posterior probabilities of pixel value co-occurrences.
#' Additionally, the value gradient is considered when \code{delta = 1} or \code{delta = 2}. 
#' The values of structural diversity entropy with \code{delta = 1} or \code{delta = 2} are not restricted and depend on the values of the input raster.
#' the metric \code{entropyNorm} is Shannon entropy normalized over maximum entropy, which depends on the size of the moving window when no nesting is used. 
#' The metric \code{entropyNorm} ranges between 0 and 1, when \code{\link{strucDiv}} is used, but may be larger than 1 when \code{\link{strucDivNest}} is used, 
#' depending on the posterior probabilities of pixel value co-occurrences.
#' The metrics \code{contrast} and \code{dissimilarity} consider the value gradient, their values are not restricted and depend on the values of the input raster.
#' The metric \code{homogeneity} quantifies the closeness of empirical probabilities to the diagonal and ranges between 0 and 1 when \code{\link{strucDiv}} is used, 
#' but may be larger than 1 when \code{\link{strucDivNest}} is used, depending on the posterior probabilities of pixel value co-occurrences.
#' @param rank logical. Should values be replaced with ranks in each co-occurrence 
#' matrix (GLCM)? Defaults to \code{FALSE}.
#' @param vMat_big matrix. The value matrix of the outer scale. Defaults to \code{NULL}, in which case no nesting is used.
#' @param SpatMat is the probability matrix that is returned by an internal function
#' to the \code{\link{strucDiv}} and \code{\link{strucDivNest}} functions. 
#' @param delta numeric, takes 3 options: \code{0}, \code{1}, or \code{2}. 
#' The parameter \code{delta} is the difference weight parameter, 
#' it defines how the differences between pixel values within a pixel pair should be weighted.  
#' If \code{rank = TRUE}, delta defines how the differences between ranks should be weighted.  
#' The default value is \code{0} (no weight). Set \code{delta = 1} for absolute weights, 
#' or \code{delta = 2} for squared weights. 
#' The \code{delta} parameter can only be set when the metric \code{entropy} is used. 
#' The metric \code{dissimilarity} automatically employs \code{delta = 1}, and \code{contrast} employs \code{delta = 2}.
#' @param nrp integer. The total number of pixel pairs. \code{nrp} is calculated internally by the 
#' functions \code{\link{strucDiv}} and \code{\link{strucDivNest}} and passed to the diversity functions.
#' @param SpatMat the GLCM that is returned by an internal function
#' to the \code{\link{strucDiv}} and \code{\link{strucDivNest}} functions. 
#' @param Hetx the structural diversity matrix that is returned by an internal function
#' to the \code{\link{strucDiv}} and \code{\link{strucDivNest}} functions. 
#' The structural diversity metric is calculated on every element
#' of the GLCM, which generates the structural diversity matrix \code{Hetx}. The sum of this
#' matrix is assigned to the center pixel of the moving window. 
#' @param narm logical. Should NAs be removed? It is automatically set to 0 if \code{na.handling = na.pass}, 
#' and to 1 if \code{na.handling = na.omit}.
#' @param display_progress logical. Should a progress bar be displayed?
#' @param parallelize logical. Should the computation be parallelized on multiple cores?
#' @param ... possible further arguments.
#' @details This function is used internally and is called 
#' as an argument to the \code{\link{strucDiv}} and \code{\link{strucDivNest}} functions.
#' @importFrom raster raster
#' @export

homogeneity <- function(rank, delta, Hetx, vMat_big = NULL, SpatMat, nrp, narm, display_progress = TRUE, ...) {
  
  switch_function <- function(rank, vMat_big) {
    
    switch(EXPR = paste(rank, is.null(vMat_big)),
           "TRUE TRUE" = .HomogeneityRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE FALSE" = .HomogeneityRankNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE TRUE" = .HomogeneityValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE FALSE" = .HomogeneityValueNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress)
           )
  }
  
  v <- switch_function(rank, vMat_big)
  return(v)
  
}
  
#' @rdname Diversity
#' @export


dissimilarity <- function(rank, delta, Hetx, vMat_big = NULL, SpatMat, nrp, narm, display_progress, ...) {
  
  switch_function <- function(rank, vMat_big) {
    
    switch(EXPR = paste(rank, is.null(vMat_big)),
           "TRUE TRUE" = .DissimilarityRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE FALSE" = .DissimilarityRankNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE TRUE" = .DissimilarityValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE FALSE" = .DissimilarityValueNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress)
           
    )
  }
  
  v <- switch_function(rank, vMat_big)
  return(v)
  
}

#' @rdname Diversity
#' @export

contrast <- function(rank, delta, Hetx, vMat_big = NULL, SpatMat, nrp, narm, display_progress, ...) {
  
  switch_function <- function(rank, vMat_big) {
    
    switch(EXPR = paste(rank, is.null(vMat_big)),
           "TRUE TRUE" = .ContrastRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE FALSE" = .ContrastRankNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE TRUE" = .ContrastValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE FALSE" = .ContrastValueNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress)
    )
  }
  
  v <- switch_function(rank, vMat_big)
  return(v)
  
}

#' @rdname Diversity
#' @export

entropy <- function(rank, delta, Hetx, vMat_big = NULL, SpatMat, nrp, narm, display_progress, 
                    parallelize = FALSE, ...) {
  
  rank_delta <- paste(rank, delta, is.null(vMat_big))
  
  if(parallelize == FALSE){
  
  switch_function <- function(rank_delta) {
    
    switch(EXPR = rank_delta,
           "FALSE 0 TRUE" = .Entropy(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 0 TRUE" = .Entropy(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 1 TRUE" = .WeightedEntropyAbsRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 1 TRUE" = .WeightedEntropyAbsValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 2 TRUE" = .WeightedEntropySqrRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 2 TRUE" = .WeightedEntropySqrValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 0 FALSE" = .EntropyNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 0 FALSE" = .EntropyNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 1 FALSE" = .WeightedEntropyAbsRankNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 1 FALSE" = .WeightedEntropyAbsValueNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 2 FALSE" = .WeightedEntropySqrRankNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 2 FALSE" = .WeightedEntropySqrValueNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, narm = narm, display_progress = display_progress)
    )
  }
  }
  
  else{
  
  switch_function <- function(rank_delta) {
    
    switch(EXPR = rank_delta,
           "FALSE 0" = .EntropyParallel(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 0" = .Entropy(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 1" = .WeightedEntropyAbsRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 1" = .WeightedEntropyAbsValueParallel(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 2" = .WeightedEntropySqrRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 2" = .WeightedEntropySqrValueParallel(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress)
    )
  }
  }
  
  v <- switch_function(rank_delta)
  return(v)
  
}

#' @rdname Diversity
#' @export

entropyNorm <- function(rank, delta, Hetx, vMat_big = NULL, SpatMat, nrp, narm, display_progress, ...) {
  
  if(is.null(vMat_big)){
  
  v <- .NormalizedEntropy(Hetx = Hetx, PMat = SpatMat, nrp = nrp, narm = narm, 
                          display_progress = display_progress)
  }
  
  else{
    
    v <- .NormalizedEntropyNested(Hetx = Hetx, vMat_big = vMat_big, PMat = SpatMat, 
                                  nrp = nrp, narm = narm, display_progress = display_progress)
    
  }
  
  return(v)
  
}







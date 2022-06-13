#' @name Diversity
#' @rdname Diversity
#' @title Diversity functions
#' @description
#' The \code{homogeneity}, \code{contrast}, \code{dissimilarity} and \code{entropy} 
#' functions are used in the default configuration of \code{\link{StrucDiv}}.
#' @param rank logical. Should values be replaced with ranks in each co-occurrence 
#' matrix. Defaults to FALSE.
#' @param Hetx is the diversity matrix that is returned by an internal function
#' to the \code{StrucDiv} function. 
#' @param vMat_big is a matrix. Default is NULL.
#' @param SpatMat is the probability matrix that is returned by an internal function
#' to the \code{StrucDiv} function. 
#' @param delta difference weight. The delta parameter can only be set when the metric entropy is used. 
#' The default value is 0. It can also take the values 1 for absolute value and
#' 2 for squared. Contrast automatically employs delta = 2 and dissimilarity employs delta = 1.
#' @param nrp integer. The number of possible pixel pairs. StrucDiv calculates it 
#' internally and passes it to the diversity functions.
#' @param narm logical. It is automatically set to 0 if na.handling = na.pass, and to 1 if na.handling = na.omit.
#' @param display_progress logical. Shows if the progress bar should be displayed.
#' @param ... possible further arguments.
#' @details This function is used internally and is called as an argument to the \code{StrucDiv}.
#' @examples 
#' \dontrun{
#' a <- raster::raster(matrix(rnorm(25), 5, 5))
#' hetx <- matrix(rnorm(25*9), 5*5, 3*3)
#' z <- matrix(runif(9, 0,1), 3, 3)
#' spatmat <- list(z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z)
#' }
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







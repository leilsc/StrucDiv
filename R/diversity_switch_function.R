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

homogeneity <- function(rank, delta, Hetx, SpatMat, nrp, narm, display_progress = TRUE, ...) {
  
  switch_function <- function(rank) {
    
    switch(EXPR = as.character(rank),
           "TRUE" = .HomogeneityRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE" = .HomogeneityValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress)
           )
  }
  
  v <- switch_function(rank)
  return(v)
  
}
  
#' @rdname Diversity
#' @export


dissimilarity <- function(rank, delta, Hetx, SpatMat, nrp, narm, display_progress, ...) {
  
  switch_function <- function(rank) {
    
    switch(EXPR = as.character(rank),
           "TRUE" = .DissimilarityRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE" = .DissimilarityValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress)
    )
  }
  
  v <- switch_function(rank)
  return(v)
  
}

#' @rdname Diversity
#' @export

contrast <- function(rank, delta, Hetx, SpatMat, nrp, narm, display_progress, ...) {
  
  switch_function <- function(rank) {
    
    switch(EXPR = as.character(rank),
           "TRUE" = .ContrastRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE" = .ContrastValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress)
    )
  }
  
  v <- switch_function(rank)
  return(v)
  
}

#' @rdname Diversity
#' @export

entropy <- function(rank, delta, Hetx, SpatMat, nrp, narm, display_progress, 
                    parallelize = FALSE, ...) {
  
  rank_delta <- paste(rank, delta)
  
  if(parallelize == FALSE){
  
  switch_function <- function(rank_delta) {
    
    switch(EXPR = rank_delta,
           "FALSE 0" = .Entropy(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 0" = .Entropy(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 1" = .WeightedEntropyAbsRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 1" = .WeightedEntropyAbsValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "TRUE 2" = .WeightedEntropySqrRank(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress),
           "FALSE 2" = .WeightedEntropySqrValue(Hetx = Hetx, PMat = SpatMat, narm = narm, display_progress = display_progress)
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

entropyNorm <- function(rank, delta, Hetx, SpatMat, nrp, narm, display_progress, ...) {
  
  v <- .NormalizedEntropy(Hetx = Hetx, PMat = SpatMat, nrp = nrp, narm = narm, 
                          display_progress = display_progress)
  
  return(v)
  
}







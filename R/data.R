#' NDVI
#'
#' @format A matrix with 221 rows and 1092 columns. 
#' Mean Normalized Difference Vegetation Index (NDVI).
#'
#' \describe{
#' \item{\code{Modified remote sensing product}}{MOD13A1v006}
#' \item{\code{Device}}{MODIS sensor}
#' \item{\code{Year}}{2018}
#' \item{\code{Aggregation}}{Mean aggregation over the growing season 2018}
#' \item{\code{Location}}{Study region in North East Eurasia}
#' \item{\code{Data quality}}{Only pixels with sufficient quality flags were used.}
#' \item{\code{NA handling}}{NA gaps were filled with a local neighborhood average.}
#' \item{\code{Value range}}{NDVI values below zero were excluded. NDVI values range between 0 and 1.}
#' \item{\code{Data retrieval}}{Data was pre-processed and downloaded from Google Earth Engine.}
#' }
#' 
#' 
#' For further details, see \url{https://lpdaac.usgs.gov/products/mod13q1v006/}
#' and \url{https://earthengine.google.com/}
#'
"ndvi"


#' NDVI, 15 gray levels
#'
#' @format A matrix with 221 rows and 1092 columns.
#' Mean Normalized Difference Vegetation Index (NDVI), with reduced
#' number of gray levels (15).
#' 
#' \describe{
#' \item{\code{Modified remote sensing product}}{MOD13A2v006}
#' \item{\code{Device}}{MODIS sensor}
#' \item{\code{Year}}{2018}
#' \item{\code{Aggregation}}{Mean aggregation over the growing season 2018}
#' \item{\code{Gray level reduction}}{Data was binned into 15 bins of equal size.}
#' \item{\code{Location}}{Study region in North East Eurasia}
#' \item{\code{Data quality}}{Only pixels with sufficient quality flags were used.}
#' \item{\code{NA handling}}{NA gaps were filled with a local neighborhood average.}
#' \item{\code{Value range}}{NDVI values below zero were excluded. NDVI values range between 0 and 1.}
#' \item{\code{Data retrieval}}{Data was pre-processed and downloaded from Google Earth Engine.}
#' }
#' 
#' 
#' For further details, see \url{https://lpdaac.usgs.gov/products/mod13q1v006/}
#' and \url{https://earthengine.google.com/}
#'
#' @examples
#' # This dataset is essentially constructed via:
#' nGrayLevels <- 15
#' require(raster)
#' ndvi <- raster(StrucDiv::ndvi)
#' ndvi15 <- cut(ndvi, breaks=seq(minValue(ndvi), maxValue(ndvi), len=nGrayLevels + 1), 
#'               include.lowest=TRUE, right=FALSE)
#' 

"ndvi.15gl"

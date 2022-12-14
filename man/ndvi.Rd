% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ndvi}
\alias{ndvi}
\title{NDVI}
\format{
A matrix with 221 rows and 1092 columns. 
Mean Normalized Difference Vegetation Index (NDVI).

\describe{
\item{\code{Modified remote sensing product}}{MOD13A1v006}
\item{\code{Device}}{MODIS sensor}
\item{\code{Year}}{2018}
\item{\code{Aggregation}}{Mean aggregation over the growing season 2018}
\item{\code{Location}}{Study region in North East Eurasia}
\item{\code{Data quality}}{Only pixels with sufficient quality flags were used.}
\item{\code{NA handling}}{NA gaps were filled with a local neighborhood average.}
\item{\code{Value range}}{NDVI values below zero were excluded. NDVI values range between 0 and 1.}
\item{\code{Data retrieval}}{Data was pre-processed and downloaded from Google Earth Engine.}
}

For further details, see \url{https://lpdaac.usgs.gov/products/mod13q1v006/}
and \url{https://earthengine.google.com/}
}
\usage{
ndvi
}
\description{
NDVI
}
\keyword{datasets}

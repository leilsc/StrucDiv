% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getValuesWindow.R
\name{getValuesWindow}
\alias{getValuesWindow}
\title{Retrieve pixel values of a defined area. 
The area is defined by the size of a window, which is centered on one pixel.}
\usage{
getValuesWindow(x, wsl, padValue, aroundTheGlobe, ...)
}
\arguments{
\item{x}{raster layer. The input raster layer.}

\item{wsl}{integer. The window side length. The window is defined by \code{wsl x wsl}.}

\item{padValue}{atomic. If a pixel is on the edge of an image, padding should be used?
Can be NA or a value.}

\item{aroundTheGlobe}{logical. Does the image go around the globe?}

\item{...}{possible further arguments.}
}
\value{
Returns a matrix.
The matrix contains the values of the defined window centered on the respective pixel.
}
\description{
Modified R Code from raster package raster::getValuesFocal.
Returns one row per pixel, which contains the values
of the pixel neighborhood that is defined by the size of the window.
The size of the window is defined by the window side length (wsl).
The window is centered on one specific pixel.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/envelope_2d.R
\name{envelope_2d}
\alias{envelope_2d}
\title{Derive and visualize a species' bi-dimensional envelope}
\usage{
envelope_2d(
  shape = NULL,
  point1 = NULL,
  point2 = NULL,
  env1 = NULL,
  env2 = NULL,
  samp = NULL,
  env1_lab = NULL,
  env2_lab = NULL,
  figname = "envelope_2d",
  scale_pos = "bottomright",
  lon_min = NULL,
  lon_max = NULL,
  lat_min = NULL,
  lat_max = NULL,
  density = FALSE,
  a_lab = "(a)",
  b_lab = "(b)",
  cell_cex = 0.1,
  col_rev = FALSE,
  env2_lab_line = 2.8
)
}
\arguments{
\item{shape}{Optional path to an ESRI shapefile representing the species distribution (\eqn{SpatialPolygonsDataFrame}).}

\item{point1}{Optional path to an ESRI shapefile representing isolated populations not included in \eqn{shape} (\eqn{SpatialPointsDataFrame}).}

\item{point2}{Optional path to an additional ESRI shapefile representing other populations neither included into \eqn{shape} not \eqn{point1} (**SpatialPointsDataFrame**).}

\item{env1}{Raster object describing the first environmental variable to be used to characterize the bi-dimensional envelope (X-axis).}

\item{env2}{Raster object describing the second environmental variable to be used to characterize the bi-dimensional envelope (Y-axis).}

\item{samp}{Optional path to a \eqn{.txt} file reporting longitude (first column) and latitude (second column) of sampled sites.}

\item{env1_lab}{Optional label for the X-axis (e.g., variable name and unit of measure).}

\item{env2_lab}{Optional label for the Y-axis (e.g., variable name and unit of measure).}

\item{figname}{Optional name for the figure produced by the function.}

\item{scale_pos}{Position of the bar reporting the geographic scale of the figure (default: "bottomright"; can be also "bottomleft", "topleft" or "topright").}

\item{lon_min}{An optional number reporting the left bound of the figure reporting the species distribution.}

\item{lon_max}{An optional number reporting the right bound of the figure reporting the species distribution.}

\item{lat_min}{An optional number reporting the lower bound of the figure reporting the species distribution.}

\item{lat_max}{An optional number reporting the upper bound of the figure reporting the species distribution.}

\item{density}{Logical indicating whether or not density is to be estimated for the bi-dimensional envelope.}

\item{a_lab}{Label for the left figure (default: "(a)").}

\item{b_lab}{Label for the right figure (default: "(b)").}

\item{cell_cex}{The size of the points in figure (a) and (b).}

\item{col_rev}{Logical indicating whether the color gradient is to be reversed or not (default: FALSE).}

\item{env2_lab_line}{A number stating the distance \eqn{env2_lab} should be placed from the Y-axis in figure (b) (default: 2.8; this number should not exceed 3.9).}
}
\value{
A \eqn{data.frame} object reporting the environmental information for the species distribution. Each line reports a georeferenced cell with resolution derived from \eqn{env1} and \eqn{env2}. Furthermore, a figure is printed in the working directory reporting the species distribution on the left and the bi-dimensional envelope on the right.
}
\description{
Derive and visualize a species' bi-dimensional envelope
}
\details{
All geographic objects need to have the same CRS. Moreover, \eqn{env1} and \eqn{env2} need an equal resolution and origin. The geographic distribution is colored following a \eqn{ModifiedSpectralScheme11Steps} scheme from the \eqn{colorBlindness} R package. To obtain the color gradient, a PCA is carried out on the values from \eqn{env1} and \eqn{env2} and colors are calibrated on PC1. The obtained gradient is used to represent habitat conditions in the geographic and environmental space (left and right figure, respectively).
}

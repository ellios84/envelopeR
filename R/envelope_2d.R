#' Derive and visualize a species' bi-dimensional envelope
#'
#' @param shape  Optional path to an ESRI shapefile representing the species distribution (\eqn{SpatialPolygonsDataFrame}).
#' @param point1 Optional path to an ESRI shapefile representing isolated populations not included in \eqn{shape} (\eqn{SpatialPointsDataFrame}).
#' @param point2 Optional path to an additional ESRI shapefile representing other populations neither included into \eqn{shape} not \eqn{point1} (**SpatialPointsDataFrame**).
#' @param env1 Raster object describing the first environmental variable to be used to characterize the bi-dimensional envelope (X-axis).
#' @param env2 Raster object describing the second environmental variable to be used to characterize the bi-dimensional envelope (Y-axis).
#' @param samp Optional path to a \eqn{.txt} file reporting longitude (first column) and latitude (second column) of sampled sites.
#' @param env1_lab Optional label for the X-axis (e.g., variable name and unit of measure).
#' @param env2_lab Optional label for the Y-axis (e.g., variable name and unit of measure).
#' @param figname Optional name for the figure produced by the function.
#' @param scale_pos Position of the bar reporting the geographic scale of the figure (default: "bottomright"; can be also "bottomleft", "topleft" or "topright").
#' @param lon_min An optional number reporting the left bound of the figure reporting the species distribution.
#' @param lon_max An optional number reporting the right bound of the figure reporting the species distribution.
#' @param lat_min An optional number reporting the lower bound of the figure reporting the species distribution.
#' @param lat_max An optional number reporting the upper bound of the figure reporting the species distribution.
#' @param density Logical indicating whether or not density is to be estimated for the bi-dimensional envelope.
#' @param a_lab Label for the left figure (default: "(a)").
#' @param b_lab Label for the right figure (default: "(b)").
#' @param cell_cex The size of the points in figure (a) and (b).
#' @param col_rev Logical indicating whether the color gradient is to be reversed or not (default: FALSE).
#' @param env2_lab_line A number stating the distance \eqn{env2_lab} should be placed from the Y-axis in figure (b) (default: 2.8; this number should not exceed 5).
#'
#' @details All geographic objects need to have the same CRS. Moreover, \eqn{env1} and \eqn{env2} need an equal resolution and origin. The geographic distribution is colored following a \eqn{ModifiedSpectralScheme11Steps} scheme from the \eqn{colorBlindness} R package. To obtain the color gradient, a PCA is carried out on the values from \eqn{env1} and \eqn{env2} and colors are calibrated on PC1. The obtained gradient is used to represent habitat conditions in the geographic and environmental space (left and right figure, respectively).
#'
#' @return A \eqn{data.frame} object reporting the environmental information for the species distribution. Each line reports a georeferenced cell with resolution derived from \eqn{env1} and \eqn{env2}. Furthermore, a figure is printed in the working directory reporting the species distribution on the left and the bi-dimensional envelope on the right.
#' @export
#'
#' @examples
envelope_2d <- function(shape=NULL, point1=NULL, point2=NULL,
                        env1=NULL, env2=NULL, samp=NULL,
                        env1_lab=NULL, env2_lab=NULL,
                        figname="envelope_2d", scale_pos="bottomright",
                        lon_min=NULL, lon_max=NULL,
                        lat_min=NULL, lat_max=NULL, density=FALSE,
                        a_lab="(a)", b_lab="(b)", cell_cex=0.1,
                        col_rev=FALSE, env2_lab_line=2.8) {

  cat("\n"); cat("... loading the required R packages .................. ")
  if(!require(rgdal, quietly = T)) install.packages("rgdal")
  if(!require(raster, quietly = T)) install.packages("raster")
  if(!require(rnaturalearth, quietly = T)) install.packages("rnaturalearth")
  if(!require(prettymapr, quietly = T)) install.packages("prettymapr")
  if(!require(colorBlindness, quietly = T)) install.packages("colorBlindness")
  if(!require(MASS, quietly = T)) install.packages("MASS")
  cat("DONE"); cat("\n")

  # loading shapefiles
  cat("\n"); cat("... loading the species' spatial distribution ........ ")
  if(!is.null(shape)) {
    shape <- rgdal::readOGR(shape, verbose = FALSE)
    extent1 <- raster::extent(shape)
    extent1 <- as.data.frame(as.vector(extent1))
    rownames(extent1) <- c("xmin", "xmax", "ymin", "ymax")
  } else extent1 <- data.frame(rep(NA, 4), row.names = c("xmin", "xmax", "ymin", "ymax"))

  if(!is.null(point1)) {
    point1 <- rgdal::readOGR(point1, verbose = FALSE)
    extent2 <- raster::extent(point1)
    extent2 <- as.data.frame(as.vector(extent2))
    rownames(extent2) <- c("xmin", "xmax", "ymin", "ymax")
  } else extent2 <- data.frame(rep(NA, 4), row.names = c("xmin", "xmax", "ymin", "ymax"))

  if(!is.null(point2)) {
    point2 <- rgdal::readOGR(point2, verbose = FALSE)
    extent3 <- raster::extent(point2)
    extent3 <- as.data.frame(as.vector(extent3))
    rownames(extent3) <- c("xmin", "xmax", "ymin", "ymax")
  } else extent3 <- data.frame(rep(NA, 4), row.names = c("xmin", "xmax", "ymin", "ymax"))

  extent_all <- cbind.data.frame(extent1, extent2, extent3)
  xmin <- apply(extent_all[1, ], 1, min, na.rm=T)
  xmax <- apply(extent_all[2, ], 1, max, na.rm=T)
  ymin <- apply(extent_all[3, ], 1, min, na.rm=T)
  ymax <- apply(extent_all[4, ], 1, max, na.rm=T)

  extent_all <- extent(matrix(c(xmin, xmax, ymin, ymax), nrow = 2, byrow = T))
  if (!is.null(lon_min)) extent_all[1] <- lon_min
  if (!is.null(lon_max)) extent_all[2] <- lon_max
  if (!is.null(lat_min)) extent_all[3] <- lat_min
  if (!is.null(lat_max)) extent_all[4] <- lat_max
  cat("DONE"); cat("\n")

  # loading env. variables
  env1 <- raster::raster(env1); env2 <- raster::raster(env2)

  # deriving the data frame with the species coordinates, env. information and colors
  cat("\n"); cat("... deriving 'env1' .................................. ")

  if(!is.null(shape)) {
    cell_shape <- raster::cellFromPolygon(env1, shape, weights=F)
    cell_shape <- unique(unlist(cell_shape))
    coo_shape <- as.data.frame(raster::xyFromCell(env1, cell_shape))
    coo_shape$Sam <- "N"
  } else coo_shape <- data.frame("x"=NA,"y"=NA,"Sam"=NA)

  if(!is.null(point1)) {
    coo_point1 <- as.data.frame(point1@coords)
    colnames(coo_point1) <- c("x","y")
    coo_point1$Sam<-"N"
  } else coo_point1 <- data.frame("x"=NA,"y"=NA,"Sam"=NA)

  if(!is.null(point2)) {
    coo_point2 <- as.data.frame(point2@coords)
    colnames(coo_point2) <- c("x","y")
    coo_point2$Sam<-"N"
  } else coo_point2 <- data.frame("x"=NA,"y"=NA,"Sam"=NA)

  if(!is.null(samp)) {
    samp <- read.table(samp, h=T)
    colnames(samp) <- c("x","y")
    samp$Sam<-"Y"
  } else samp <- data.frame("x"=NA,"y"=NA,"Sam"=NA)

  coo <- rbind.data.frame(coo_shape, coo_point1, coo_point2, samp)

  coo$env1 <- raster::extract(env1, coo[,c(1,2)])
  cat("DONE"); cat("\n")

  cat("\n"); cat("... deriving 'env2' .................................. ")
  coo$env2 <- raster::extract(env2, coo[,c(1,2)])
  cat("DONE"); cat("\n")
  coo <- na.omit(coo)

  cat("\n"); cat("... principal component anaysis ...................... ")
  pca <- prcomp(coo[,c("env1","env2")], center = T, scale. = T)
  cat("DONE"); cat("\n")

  cat("\n"); print(summary(pca)); cat("\n")
  pca <- as.data.frame(pca$x)
  coo <- cbind.data.frame(coo, pca)

  rbPal <- colorRampPalette(colorBlindness::ModifiedSpectralScheme11Steps)

  if (isTRUE(col_rev)) rbPal <- rev(colorRampPalette(colorBlindness::ModifiedSpectralScheme11Steps))

  coo$Col <- rbPal(1000)[as.numeric(cut(coo$PC1, breaks = 1000))]
  cat("\n"); cat("... characterization of the envelope  ................ DONE ")
  cat("\n")

  if(!is.null(samp)) {
    coo_Y <- subset(coo, coo$Sam == "Y")
  }
  coast<-rnaturalearth::ne_coastline(scale = "medium", returnclass = c("sp", "sf"))

  cat("\n"); cat("... plotting niche_2d figure ......................... ")
  postscript(paste0(figname, ".eps"), width = 12, height = 6)
  par(mar=c(1, 0, 1, 0), mfrow = c(1, 2), oma=c(3, 5, 3, 5), mgp=c(1.9, 0.9, 0))
  plot(coo$x, coo$y, col=coo$Col, xlab="", ylab="", pch=16, cex=cell_cex, las=1,
       cex.axis=1.1, xlim=c(extent_all[1], extent_all[2]),
       ylim=c(extent_all[3], extent_all[4]))
  raster::plot(coast, add=T, lwd=.5, col="gray23")
  if(!is.null(samp)) {
    points(coo_Y$x, coo_Y$y, cex=1)
    }
  prettymapr::addscalebar(widthhint = 0.25, unitcategory = "metric", htin = 0.1, padin = c(0.1, 0.1),
						  style = "bar", bar.cols = c("black", "black"), lwd = 1,
						  linecol = "black", tick.cex = 0.7, labelpadin = 0.08, label.cex = 1,
						  label.col = "black", pos = scale_pos)
  box(lwd=2)

  plot(coo$env1, coo$env2, pch=16, cex=cell_cex, col=coo$Col, yaxt="n",
       xlab="", ylab="")
  if (density) {
    dens <- MASS::kde2d(coo$env1, coo$env2)
    levels <- quantile(as.vector(dens$z))
    contour(dens, levels=levels, labels=NULL, add=T, labcex=1, col=rev(gray.colors(length(levels))))
    }
  abline(v=median(coo$env1), lty=5, col="gray23")
  abline(h=median(coo$env2), lty=5, col="gray23")
  if(!is.null(samp)) {
    points(coo_Y$env1, coo_Y$env2, cex=1)
    }
  axis(side = 4, las=1, cex.axis=1.1)
  box(lwd=2)

  mtext(text = "Longitude", side = 1, line = 1.8, outer = T, at = .25, cex = 1.5)
  mtext(text = "Latitude", side = 2, line = 2.8, outer = T, at = .5, cex = 1.5)
  mtext(text = env1_lab, side = 1, line = 1.8, outer = T, at = .75, cex = 1.5)
  mtext(text = env2_lab, side = 4, line = env2_lab_line, outer = T, at = .5, cex = 1.5)
  mtext(text = a_lab, side = 3, line = -.5, outer = T, at = 0.02, cex = 1.5)
  mtext(text = b_lab, side = 3, line = -.5, outer = T, at = 0.52, cex = 1.5)

  dev.off()
  cat("DONE"); cat("\n")

  return(coo)

}

#' Visualizing the bidimensional envelope of a species
#'
#' @param shape  Optional path to an ESRI shapefile (\eqn{.shp}) representing the species distribution as a polygon (\eqn{SpatialPolygonsDataFrame}).
#' @param env1 Raster object describing the first environmental variable to be used to characterize the bi-dimensional envelope (X-axis).
#' @param env2 Raster object describing the second environmental variable to be used to characterize the bi-dimensional envelope (Y-axis).
#' @param samp Optional path to a \eqn{.csv} file reporting longitude (first column) and latitude (second column) of the sampled sites.
#' @param env1_lab Optional label for the X-axis (e.g., variable name and unit of measure).
#' @param env2_lab Optional label for the Y-axis (e.g., variable name and unit of measure).
#' @param figname Optional name for the figure produced by the function.
#' @param scale_pos Position of the bar reporting the geographic scale of the figure (default: "bottomright"; can be also "bottomleft", "topleft" or "topright").
#' @param extent A numeric vector for the left, right, lower and upper bounds of the map.
#' @param density Logical indicating whether or not density is to be estimated for the bi-dimensional envelope.
#' @param a_lab Label for the left figure (default: "(a)").
#' @param b_lab Label for the right figure (default: "(b)").
#' @param cell_cex The size of the points in figure a and b.
#' @param col_rev Logical indicating whether the color gradient is to be reversed or not (default: FALSE).
#' @param env2_lab_line A number stating the distance \eqn{env2_lab} should be placed from the Y-axis in figure b (default: 2.8; this number should not exceed 3.9).
#'
#' @details All the geographic objects used need to be expressed in the WGS84 system (EPSG:4326). Moreover, \eqn{env1} and \eqn{env2} need an equal resolution and origin. The geographic distribution is colored following a \eqn{ModifiedSpectralScheme11Steps} scheme from the \eqn{colorBlindness} R package. To obtain the color gradient, a PCA is carried out on the values from \eqn{env1} and \eqn{env2} and colors are calibrated on PC1. The obtained gradient is used to represent habitat conditions in the geographic and environmental space (left and right figure, respectively).
#'
#' @return A \eqn{data.frame} object reporting the environmental information for the species distribution. Each line reports a georeferenced cell with resolution derived from \eqn{env1} and \eqn{env2}. Furthermore, a figure is printed in the working directory reporting the species distribution on the left and the bi-dimensional envelope on the right.
#' @export
#'
#' @examples
envelope2d <- function(shape = NULL, env1 = NULL, env2 = NULL, samp = NULL,
                       figname = "envelopeR2d", env1_lab = NULL, env2_lab = NULL,
                       extent = NULL, a_lab = "(a)", b_lab = "(b)", cell_cex = 0.1,
                       col_rev = FALSE, env2_lab_line = 2.8, scale_pos = "bottomleft",
                       density = FALSE) {

  # Dependencies
  if(!require(terra, quietly = T)) install.packages("terra")
  if(!require(sp, quietly = T)) install.packages("sp")
  if(!require(prettymapr, quietly = T)) install.packages("prettymapr")
  if(!require(MASS, quietly = T)) install.packages("MASS")
  if(!require(colorBlindness, quietly = T)) install.packages("colorBlindness")

  # Loading the environmental rasters, the coastline, defining the colorRampPalette

  tryCatch({
    env1 <- terra::rast(env1)
  }, error = function(e){})

  crsM <- terra::crs(env1)

  tryCatch({
    env2 <- terra::rast(env2)
  }, error = function(e){})

  tryCatch({
    coast <- terra::vect("./../Data/ne_10m_coastline.shp")
  }, error = function(e){})

  rbPal <- colorRampPalette(colorBlindness::ModifiedSpectralScheme11Steps)

  # Environmental information at the sampled sites (if not null)
  if(!is.null(samp)) {
    samp <- read.csv(file = samp, header = TRUE, sep = ";")
    colnames(samp) <- c("x","y")
    samp$env1 <- terra::extract(x = env1, y = samp)[, 2]
    samp$env2 <- terra::extract(x = env2, y = samp[, 1:2])[, 2]
    samp <- na.omit(samp)
    samp$smp <- "Y"
  } else {
      samp <- data.frame("x" = NA,"y" = NA, env1 = "NA", env2 = "NA", "smp" = NA)
      }

  # Environmental information across the species distribution (if not null)
  if(!is.null(shape)) {
    shape <- terra::vect(x = shape)
    env1 <- terra::crop(x = env1, y = shape)
    env1 <- terra::mask(x = env1, mask = shape, inverse = FALSE, updatevalue = NA)
    env1_xy <- as.data.frame(terra::xyFromCell(object = env1, cell = terra::cells(x = env1)))
    env1_xy$env1 <- terra::extract(x = env1, y = env1_xy)[, 2]
    env2 <- terra::crop(x = env2, y = shape)
    env2 <- terra::mask(x = env2, mask = shape, inverse = FALSE, updatevalue = NA)
    env2_xy <- as.data.frame(terra::xyFromCell(object = env2, cell = terra::cells(x = env2)))
    env2_xy$env2 <- terra::extract(x = env2, y = env2_xy)[, 2]
    env12_xy <- cbind.data.frame(env1_xy, env2_xy)
    env12_xy <- na.omit(env12_xy)
    env12_xy <- env12_xy[, c("x", "y", "env1", "env2")]
    env12_xy$smp <- "N"
    env12_xy <- rbind.data.frame(samp, env12_xy)
    env12_xy <- na.omit(env12_xy)
    pca.df <- env12_xy[, c("env1","env2")]
    pca.df$env1 <- as.numeric(pca.df$env1)
    pca.df$env2 <- as.numeric(pca.df$env2)
    pca <- prcomp(pca.df, center = T, scale. = T)# PCA
    cat("\n"); print(summary(pca)); cat("\n")
    pca <- as.data.frame(pca$x)
    env12_xy <- cbind.data.frame(env12_xy, pca)
  } else if (!is.null(extent)) {
    # Environmental information across a specified extent (if not null)
    env1 <- terra::crop(x = env1, y = terra::ext(extent))
    env1_xy <- as.data.frame(terra::xyFromCell(object = env1, cell = terra::cells(x = env1)))
    env1_xy$env1 <- terra::extract(x = env1, y = env1_xy)[, 2]
    env2 <- terra::crop(x = env2, y = terra::ext(extent))
    env2_xy <- as.data.frame(terra::xyFromCell(object = env2, cell = terra::cells(x = env2)))
    env2_xy$env2 <- terra::extract(x = env2, y = env2_xy)[, 2]
    env12_xy <- cbind.data.frame(env1_xy, env2_xy)
    env12_xy <- na.omit(env12_xy)
    env12_xy <- env12_xy[, c("x", "y", "env1", "env2")]
    env12_xy$smp <- "N"
    env12_xy <- rbind.data.frame(samp, env12_xy)
    env12_xy <- na.omit(env12_xy)
    pca.df <- env12_xy[, c("env1","env2")]
    pca.df$env1 <- as.numeric(pca.df$env1)
    pca.df$env2 <- as.numeric(pca.df$env2)
    pca <- prcomp(pca.df, center = T, scale. = T) # PCA
    cat("\n"); print(summary(pca)); cat("\n")
    pca <- as.data.frame(pca$x)
    env12_xy <- cbind.data.frame(env12_xy, pca)
    coast <- terra::crop(x = coast, terra::ext(extent))
  } else {
    # Environmental information across the input raster data
    env1_xy <- as.data.frame(terra::xyFromCell(object = env1, cell = terra::cells(x = env1)))
    env1_xy$env1 <- terra::extract(x = env1, y = env1_xy)[, 2]
    env2_xy <- as.data.frame(terra::xyFromCell(object = env2, cell = terra::cells(x = env2)))
    env2_xy$env2 <- terra::extract(x = env2, y = env2_xy)[, 2]
    env12_xy <- cbind.data.frame(env1_xy, env2_xy)
    env12_xy <- na.omit(env12_xy)
    env12_xy <- env12_xy[, c("x", "y", "env1", "env2")]
    env12_xy$smp <- "N"
    env12_xy <- rbind.data.frame(samp, env12_xy)
    env12_xy <- na.omit(env12_xy)
    pca.df <- env12_xy[, c("env1","env2")]
    pca.df$env1 <- as.numeric(pca.df$env1)
    pca.df$env2 <- as.numeric(pca.df$env2)
    pca <- prcomp(pca.df, center = T, scale. = T) # PCA
    cat("\n"); print(summary(pca)); cat("\n")
    pca <- as.data.frame(pca$x)
    env12_xy <- cbind.data.frame(env12_xy, pca)
  }

  # Defining the extent of the map
  extent <- terra::ext(env1)

  # Color
  if (!isFALSE(col_rev)) {
    rbPal <- colorRampPalette(rev(colorBlindness::ModifiedSpectralScheme11Steps))
  }
  env12_xy$Col <- rbPal(100)[as.numeric(cut(env12_xy$PC1, breaks = 100))]

  # Subsetting sampled locations
  Y <- subset(env12_xy, env12_xy$smp == "Y")

  # Plotting
  jpeg(paste0(figname, ".jpeg"), width = 12, height = 6, units = 'in', res = 900)
  par(mar=c(1, 0, 1, 0), mfrow = c(1, 2), oma=c(4, 5, 4, 5), mgp=c(1.9, 0.9, 0))

  plot(env12_xy$x, env12_xy$y, t = "n", xlab = "", ylab = "", pch = 16, cex = cell_cex, las = 1, cex.axis = 0.8, xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]))
  points(env12_xy$x, env12_xy$y, pch = 16, cex = cell_cex, col = env12_xy$Col)
  plot(coast, add = TRUE, lwd = .5, col = "gray23")
  points(Y$x, Y$y, cex = 1, pch = 21, bg = Y$Col, lwd = 1.5)
  prettymapr::addscalebar(widthhint = 0.25, unitcategory = "metric", htin = 0.1, padin = c(0.25, 0.2), style = "bar", bar.cols = c("black", "black"), lwd = 1, linecol = "black", tick.cex = 0.7, labelpadin = 0.08, label.cex = 1, label.col = "black", pos = scale_pos)
  legend("topright", a_lab, cex = 1.8, bty = "n")
  box(lwd=2)

  plot(env12_xy$env1, env12_xy$env2, pch=16, cex=cell_cex, col = env12_xy$Col, yaxt = "n", xlab = "", ylab="", cex.axis = 0.8)

  # Density
  if (!isFALSE(density)) {
    dens <- MASS::kde2d(as.numeric(env12_xy$env1), as.numeric(env12_xy$env2))
    levels <- quantile(as.vector(dens$z))
    tryCatch({
      terra::contour(dens, levels = levels, labels = NULL, add = TRUE, labcex = 1, col = rev(gray.colors(length(levels))))
    }, error = function(e){})
  }

  # Quadrants and sampling points
  abline(v = median(as.numeric(env12_xy$env1)), lty=5, col="gray23")
  abline(h = median(as.numeric(env12_xy$env2)), lty=5, col="gray23")
  points(Y$env1, Y$env2, cex = 1, pch = 21, bg = Y$Col, lwd = 1.5)

  legend("topright", b_lab, cex = 1.8, bty = "n")
  axis(side = 4, las=1, cex.axis = 0.8)
  box(lwd=2)

  mtext(text = "Longitude", side = 1, line = 1.8, outer = T, at = .25, cex = 1.5)
  mtext(text = "Latitude", side = 2, line = 2.8, outer = T, at = .5, cex = 1.5)
  mtext(text = env1_lab, side = 1, line = 1.8, outer = T, at = .75, cex = 1.5)
  mtext(text = env2_lab, side = 4, line = env2_lab_line, outer = T, at = .5, cex = 1.5)

  dev.off()

  # Return
  return(env12_xy)

}

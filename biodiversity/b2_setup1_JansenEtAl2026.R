############################################################
# LOAD ENVIRONMENTAL CONTEXT RASTERS
############################################################
r.stack <- rast(file.path(env.derived,
                          paste0("Circumpolar_EnvData_", res, "_shelf_mask_unscaled_variables.tif")))
r2 <- r.stack$depth
# r2[1:350, ] <- NA
r.stack.subset <- r.stack[[c(1:3, 14, 19, 23, 28, 41, 45:47)]]
# rm(r.stack)


############################################################
# LOAD COASTLINE
############################################################
coast.unprojected <- vect(paste0(usr.dropbox.dir,
                                 "data_environmental/raw/antarctic_coastline_2023/add_coastline_high_res_polygon_v7_8.shp"))
coast.proj <- project(coast.unprojected, r2)

############################################################
# LOAD SAMPLE LOCATIONS
############################################################
sampling_dir   <- paste0(usr.dropbox.dir, "data_products/modelling_files/circum_antarctic")
dat <- readRDS(file.path(sampling_dir, paste0("cover_modelling_inputs_", res, ".rds")))
xy <- dat$cell_metrics[,6:7]


############################################################
# INSET MAP SETUP
############################################################

# Zoom extents
AP.ext <- ext(-2750000, -2250000, 1050000, 1900000)
PB.ext <- ext(2100000, 2600000, 300000, 1500000)
RS.ext <- ext(-450000, 500000, -1900000, -1250000)
MG.ext <- ext(1300000, 1800000, -2300000, -1900000)

# Shift positions for inset plotting
AP.dx <-  2700000
PB.dx <- -1100000
RS.dx <-  -275000
MG.dx <-  -300000

AP.dy <- -400000
PB.dy <- -500000
RS.dy <- 1100000
MG.dy <-  850000

# Legend positions
l.x  <- -2800000
l.y  <- -1750000
l.xt <- l.x + 480000
l.yt <- l.y + 50000

l.x2 <- -1600000
l.y2 <- -1950000
plg.c <- list(size = c(0.3, 1), ext = c(-5000000,1800000,-2200000,-2100000), loc = "bottom", horiz=TRUE)

# Title position
loc.main <- c(0, 2800000)

# Coastline insets
AP.coast <- rescale(crop(coast.proj, AP.ext), fx = 2, x0 = mean(AP.ext[1:2]), y0 = mean(AP.ext[3:4]))
PB.coast <- rescale(crop(coast.proj, PB.ext), fx = 2, x0 = mean(PB.ext[1:2]), y0 = mean(PB.ext[3:4]))
RS.coast <- rescale(crop(coast.proj, RS.ext), fx = 2, x0 = mean(RS.ext[1:2]), y0 = mean(RS.ext[3:4]))
MG.coast <- rescale(crop(coast.proj, MG.ext), fx = 2, x0 = mean(MG.ext[1:2]), y0 = mean(MG.ext[3:4]))

AP.coast.b <- shift(AP.coast, AP.dx, AP.dy)
PB.coast.b <- shift(PB.coast, PB.dx, PB.dy)
RS.coast.b <- shift(RS.coast, RS.dx, RS.dy)
MG.coast.b <- shift(MG.coast, MG.dx, MG.dy)

# Sample points for inset maps
sp.samples <- vect(xy, crs = crs(r2))

AP.samples <- rescale(crop(sp.samples, AP.ext), fx = 2)
PB.samples <- rescale(crop(sp.samples, PB.ext), fx = 2)
RS.samples <- rescale(crop(sp.samples, RS.ext), fx = 2)
MG.samples <- rescale(crop(sp.samples, MG.ext), fx = 2)

AP.b.samples <- shift(AP.samples, AP.dx, AP.dy)
PB.b.samples <- shift(PB.samples, PB.dx, PB.dy)
RS.b.samples <- shift(RS.samples, RS.dx, RS.dy)
MG.b.samples <- shift(MG.samples, MG.dx, MG.dy)

## crop, shift and zoom into inset maps
crop_and_shift <- function(ra, ra.polygons=NA){
  ## extract and zoom into inset maps
  AP <- rescale(crop(ra, AP.ext), fx=2, x0=mean(AP.ext[1:2]), y0=mean(AP.ext[3:4]))
  PB <- rescale(crop(ra, PB.ext), fx=2, x0=mean(PB.ext[1:2]), y0=mean(PB.ext[3:4]))
  RS <- rescale(crop(ra, RS.ext), fx=2, x0=mean(RS.ext[1:2]), y0=mean(RS.ext[3:4]))
  MG <- rescale(crop(ra, MG.ext), fx=2, x0=mean(MG.ext[1:2]), y0=mean(MG.ext[3:4]))
  
  ## shift location on where to plot the inset map
  AP.b <<- shift(AP, AP.dx, AP.dy)
  PB.b <<- shift(PB, PB.dx, PB.dy)
  RS.b <<- shift(RS, RS.dx, RS.dy)
  MG.b <<- shift(MG, MG.dx, MG.dy)
  
  if(!is.na(ra.polygons)){
    ## extract and zoom into inset maps
    AP <- rescale(crop(ra.polygons, AP.ext), fx=2, x0=mean(AP.ext[1:2]), y0=mean(AP.ext[3:4]))
    PB <- rescale(crop(ra.polygons, PB.ext), fx=2, x0=mean(PB.ext[1:2]), y0=mean(PB.ext[3:4]))
    RS <- rescale(crop(ra.polygons, RS.ext), fx=2, x0=mean(RS.ext[1:2]), y0=mean(RS.ext[3:4]))
    MG <- rescale(crop(ra.polygons, MG.ext), fx=2, x0=mean(MG.ext[1:2]), y0=mean(MG.ext[3:4]))
    ## shift location on where to plot the inset map
    AP.polygons <<- shift(AP, AP.dx, AP.dy)
    PB.polygons <<- shift(PB, PB.dx, PB.dy)
    RS.polygons <<- shift(RS, RS.dx, RS.dy)
    MG.polygons <<- shift(MG, MG.dx, MG.dy)
  }
}
## to help with plotting correctly
get_breaks_from_levels <- function(x) {
  vals <- levels(x)[[1]]$value
  vals <- sort(vals)
  c(min(vals) - 0.5, head(vals, -1) + 0.5, max(vals) + 0.5)
}
## plot inset maps and lines
add.plotsandlines <- function(inset.cols=cols, plot.points=FALSE, range=NULL, breaks=NULL, polygons=FALSE, border.col=NA, fill.col="red", coast=TRUE, alpha=NULL, fill_range=FALSE, maxcell=500000){
  ## add coastline
  coast.col="grey40"
  land.col="grey60"
  plot(coast.proj[values(coast.proj)=="land"], col=land.col, add=TRUE)
  lines(coast.proj, col=coast.col)
  
  ## add inset maps
  plot(AP.b, col=inset.cols, add=TRUE, colNA="white", legend=FALSE, range=range, breaks=breaks, fill_range=fill_range, maxcell=maxcell)
  if(coast){
    plot(AP.coast.b[values(AP.coast.b)=="land"], col=land.col, add=TRUE)
    plot(AP.coast.b, border=coast.col, add=TRUE)
  }
  if(polygons) plot(AP.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
  if(plot.points) {points(AP.b.samples, col="black", pch=19)}
  plot(PB.b, col=inset.cols, add=TRUE, colNA="white", legend=FALSE, range=range, breaks=breaks, fill_range=fill_range, maxcell=maxcell)
  if(coast){
    plot(PB.coast.b[values(PB.coast.b)=="land"], col=land.col, add=TRUE)
    plot(PB.coast.b, border=coast.col, add=TRUE)
  }
  if(polygons) plot(PB.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
  if(plot.points) {points(PB.b.samples, col="black", pch=19)}
  plot(RS.b, col=inset.cols, add=TRUE, colNA="white", legend=FALSE, range=range, breaks=breaks, fill_range=fill_range, maxcell=maxcell)
  if(coast){
    plot(RS.coast.b[values(RS.coast.b)=="land"], col=land.col, add=TRUE)
    plot(RS.coast.b, border=coast.col, add=TRUE)
  }
  if(polygons) plot(RS.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
  if(plot.points) {points(RS.b.samples, col="black", pch=19)}
  plot(MG.b, col=inset.cols, add=TRUE, colNA="white", legend=FALSE, range=range, breaks=breaks, fill_range=fill_range, maxcell=maxcell)
  if(coast){
    plot(MG.coast.b[values(MG.coast.b)=="land"], col=land.col, add=TRUE)
    plot(MG.coast.b, border=coast.col, add=TRUE)
  }
  if(polygons) plot(MG.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
  if(plot.points) {points(MG.b.samples, col="black", pch=19)}
  
  ## lines around those boxes and the zoom areas
  lines(AP.ext, lty=1)
  lines(PB.ext, lty=1)
  lines(RS.ext, lty=1)
  lines(MG.ext, lty=1)
  lines(ext(AP.b), lty=1)
  lines(ext(PB.b), lty=1)
  lines(ext(RS.b), lty=1)
  lines(ext(MG.b), lty=1)
  
  ## draw lines from areas to zooms
  lines(matrix(c(AP.ext[2], ext(AP.b)[1], mean(c(AP.ext[3],AP.ext[4])), mean(c(ext(AP.b)[3],ext(AP.b)[4]))), ncol=2), lty=1)
  lines(matrix(c(PB.ext[1], ext(PB.b)[2], mean(c(PB.ext[3],PB.ext[4])), mean(c(ext(PB.b)[3],ext(PB.b)[4]))), ncol=2), lty=1)
  lines(matrix(c(mean(c(RS.ext[1],RS.ext[2])), mean(c(ext(RS.b)[1],ext(RS.b)[2])), RS.ext[4], ext(RS.b)[3]), ncol=2), lty=1)
  lines(matrix(c(mean(c(MG.ext[1],MG.ext[2])), mean(c(ext(MG.b)[1],ext(MG.b)[2])), MG.ext[4], ext(MG.b)[3]), ncol=2), lty=1)
}

## plot inset maps and lines
add.polygons <- function(border.col=NA, fill.col="red", alpha=NULL){
  ## add inset polygons
  plot(AP.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
  plot(PB.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
  plot(RS.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
  plot(MG.polygons, add=TRUE, col=fill.col, border=border.col, alpha=alpha)
}

############################################################
# SETUP COLOURS
############################################################

better.ramp <- grDevices::colorRampPalette(rev(blues9))
better.ramp <- viridis
lcex <- 1

# Colours for discrete hotspot categories
perc.ramp <- grDevices::colorRampPalette(c("grey90", "grey65", rev(viridis(3))))
cols <- perc.ramp(5)
leg.steps <- c(1, 2, 2.75, 3)

# Colours for overlap categories
cols.discr <- c("#0072B2", "#009E73", "#E69F00")

# Continuous colour scales
cols.c <- viridis(99)
cols.grey <- grey.colors(99, end = 0.99)

# Continuous - logged
log_data <- log10(seq(0.01, 1, length.out = 100))
normalized_log_data <- (log_data - min(log_data)) / (max(log_data) - min(log_data))
cols.c.log <- viridis(99)[as.numeric(cut(normalized_log_data, breaks = 99))]
cols.grey.rev.log <- rev(grey.colors(99, end = 0.99))[as.numeric(cut(normalized_log_data, breaks = 99))]

# Continuous - power transformed
power_data <- seq(0, 1, length.out = 100)^0.5
normalized_power_data <- (power_data - min(power_data)) / (max(power_data) - min(power_data))
cols.c.power <- viridis(99)[as.numeric(cut(normalized_power_data, breaks = 99))]
cols.grey.rev.power <- rev(grey.colors(99, end = 0.99))[as.numeric(cut(normalized_power_data, breaks = 99))]

# Coast/land colours
coast.col <- "grey40"
land.col  <- "grey60"

# Violin plot colours
col.in  <- cols.discr[1]
col.out <- cols[1]

# Plotting ranges
range.pa <- c(0, 1)
range.ab <- c(0, 100)

# Common plotting parameters
plot.par <- par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))

xlim=c(-2800000,3000000)
ylim=c(-2700000,2580000)
maxcell=10000000

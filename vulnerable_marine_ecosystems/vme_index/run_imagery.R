################################################################################
#                                 INIT
################################################################################
# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index")

# Source method variables and load packages
source("config_imagery.R")
# Source utils functions
source("utils.R")
# Source raster_tot
#source("../data_preparation/csv_2_raster.R")

################################################################################
#                                 LOAD DATA
################################################################################
# Read raster
#raster_bio <- stack(path_bio_data)
bio_data <- read.csv("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210916_biigle_vme_cover.csv")
names(bio_data) <- gsub(x = names(bio_data), pattern = "\\.", replacement = "-")
bio_data = bio_data[bio_data$survey == "PS81", ]
# Curate data
#data_ccamlr <- curate_ccamlr_registry(ccamlr_registry)
# Read VME Taxa Scores
vulnerability_scores_taxa <- read.csv(path_taxa_scores)
#write.csv(vulnerability_scores_taxa,path_taxa_scores)
#vulnerability_scores_taxa$morpho_taxon <- paste(vulnerability_scores_taxa$ï..Parent, vulnerability_scores_taxa$Taxon, sep="_")
head(vulnerability_scores_taxa)

################################################################################
#                   COMPUTE VULNERABILITY SCORE OF EACH TAXA
################################################################################
cat("Computing vulnerability score of each taxa using",
    toupper(agg_vulnerability_score),
    "aggregation method...")
if (agg_vulnerability_score == "mean") {
  result <- apply(vulnerability_scores_taxa[, 2:ncol(vulnerability_scores_taxa)], 1,
                  function(x) mean(x, na.rm=TRUE))
} else if (agg_vulnerability_score == "quadratic_mean") {
  quadratic_mean <- function(x) {
    x_ <- x[!is.na(x)]
    return(sqrt(sum(x_^2)/length(x_)))
  }
  result <- apply(vulnerability_scores_taxa[, 2:ncol(vulnerability_scores_taxa)], 1,
                  function(x) quadratic_mean(x))
} else {
  cat("\tERROR: Unknown aggregation method:", agg_vulnerability_score,
      "\n\t\tPlease choose among:", c("mean", "quadratic_mean"))
  exit()
}
vulnerability_scores_taxa$score <- result
vulnerability_scores_taxa <- vulnerability_scores_taxa[order(vulnerability_scores_taxa$score,
                                                             decreasing = TRUE),]
print(vulnerability_scores_taxa[ , c("morpho_taxon", "score")])

################################################################################
#                           COMPUTE ABUNDANCE SCORES
################################################################################
cat("Computing abundance scores...")
morpho_taxa_list <- colnames(bio_data)[6:ncol(bio_data)]
cat("\tGrouping abundance data of each taxon into", n_abundance_categories,
    "categories using Jenks breaks method...")
# Get Abundance scores
df_abundance_score <- data.frame(bio_data)
for (taxon in morpho_taxa_list) {
  abundance_taxon <- bio_data[ , taxon]
  jenks_breaks <- get_jenks_breaks(abundance_taxon, n_abundance_categories)
  abundance_scores <- apply_jenks_breaks(abundance_taxon, jenks_breaks, n_abundance_categories)
  df_abundance_score[ , taxon] <- abundance_scores
}
print(head(df_abundance_score))

################################################################################
#                           COMPUTE VME INDEX
################################################################################
cat("Computing VME indexes...")
cat("\tModulating abundance scores with vulnerability scores...")
df_indicator_abundance_score <- data.frame(df_abundance_score)
for (taxon in morpho_taxa_list) {
  taxon_score <- vulnerability_scores_taxa[vulnerability_scores_taxa$morpho_taxon == taxon, "score"]
  df_indicator_abundance_score[ , taxon] <- df_abundance_score[ , taxon] * taxon_score
}
print(head(df_indicator_abundance_score))
cat("\tAggregate scores across taxa for each records using", toupper(vme_index_agg),"method...")
if (vme_index_agg == "max") {
  result <- apply(df_indicator_abundance_score[, morpho_taxa_list], 1,
                  function(x) max(x, na.rm=TRUE))
} else if (vme_index_agg == "mean") {
  result <- apply(df_indicator_abundance_score[, morpho_taxa_list], 1,
                  function(x) mean(x, na.rm=TRUE))
} else if (vme_index_agg == "median") {
  result <- apply(df_indicator_abundance_score[, morpho_taxa_list], 1,
                  function(x) median(x, na.rm=TRUE))
} else {
  cat("\tERROR: Unknown aggregation method:", agg_vulnerability_score,
      "\n\t\tPlease choose among:", c("mean", "median", "max"))
  exit()
}
df_indicator_abundance_score$vme_index <- result
print(head(df_indicator_abundance_score["vme_index"]))
# TODO
# Get Jenks breaks for vme indexes
vme_index_vals <- as.numeric(df_indicator_abundance_score$vme_index)
vme_index_breaks <- getJenksBreaks(vme_index_vals, n_index_categories+1)

# Get final VME index
df_vme_index <- data.frame(df_indicator_abundance_score)
df_vme_index <- df_vme_index[ , !(colnames(df_vme_index) %in% morpho_taxa_list)]
for (row in 1:nrow(df_vme_index)) {
  vme_index_area <- df_vme_index[row, "vme_index"]
  for (vme_index_break_idx in 1:(length(vme_index_breaks)-1)) {
    vme_index_break <- vme_index_breaks[vme_index_break_idx]
    if (vme_index_area >= vme_index_break) {
      df_vme_index[row, "vme_index"] <- vme_index_break_idx
    } else {
      break
    }
  }
}

################################################################################
#                           MAPPING
################################################################################
# Create new variable
pts <- df_vme_index
# Projection
coordinates(pts) <- c("longitude", "latitude")
projection(pts) <- "+proj=longlat +datum=WGS84"
pts <- SOproj(pts)

cat("Cropping and resampling (",
    resolution_raster, "x", resolution_raster, "m ) the reference raster...")
# Crop reference raster according to pts extent

mask <- buffer(pts[1, ], 10000)

raster_ref <- crop(SmallBathy, mask)
# Set raster resolution
if (resolution_raster < res(raster_ref)[1]) {
  resample_fact <- res(raster_ref)[1] / resolution_raster
  raster_ref <- disaggregate(raster_ref,
                             fact=resample_fact,
                             method='bilinear')
} else if (resolution_raster > res(raster_ref)[1]) {
  resample_fact <- resolution_raster / res(raster_ref)[1]
  raster_ref <- aggregate(raster_ref,
                          fact=resample_fact,
                          method='bilinear')
} else {
  raster_ref <- raster_ref
}

cat("Using the projection system:", proj4string(raster_ref), "...")
# Enforce same projection
proj4string(pts) <- proj4string(raster_ref)

# Rasterise data
# TODO: change fun
cat("Rasterizing data ...")
raster_vme <- rasterize(pts,
                        raster_ref,
                        "vme_index",
                        fun=mean)


blue.col <- colorRampPalette(c("darkblue", "lightblue"))
blue.br <- seq(from=-11000, to=0, by=1000)
yellow.col <- colorRampPalette(c("yellow", "red"))
yellow.br <- seq(from=1, to=max(n_index_categories), by=1)

par(mar=c(3,3,3,3))
plot(raster_ref, col=blue.col(22), breaks=blue.br, legend=FALSE)
plot(rasterToPolygons(raster_vme), add=TRUE, col=yellow.col(n_index_categories), breaks=yellow.br, legend=FALSE)
legend("left", legend=unique(pts$vme_index), col=yellow.col(n_index_categories), pch=15, border="black", title = "VME index")









# Stack rasters
raster_tot <- stack(raster_bio, raster_survey)
cat("Number of layers:", nlayers(raster_tot), "...")

# Save raster
if (save_raster) {
  cat("Saving raster:", path_out, "...")
  writeRaster(raster_tot,
              filename=path_out,
              suffix=names(raster_tot),
              bylayer=TRUE,
              format="GTiff",
              overwrite=TRUE)
}


















# Get CCAMLR coordinates
coords <- data.frame(x = as.numeric(df_vme_index$longitude), 
                     y= as.numeric(df_vme_index$latitude))
coordinates(coords) <- c("x","y")
projection(coords) <- "+proj=longlat +datum=WGS84"
coords <- SOproj(coords)
SOmap_auto(coords, input_points = FALSE, input_lines = FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)

# Create mask around one of the records, with some buffer
mask <- buffer(coords[20, ], 50000)
#mask <- buffer(coords[1, ], 50000)
#mask <- buffer(coords[2, ], 50000)
#mask <- buffer(coords[75, ], 50000)
#mask <- buffer(coords[48, ], 50000)
SOmap_auto(coords, input_points = FALSE, input_lines = FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)
plot(mask, add=TRUE)

# Crop raster around this mask
raster.crop <- crop(SmallBathy, mask)
# Aesthetic
blue.col <- colorRampPalette(c("darkblue", "lightblue"))
blue.br <- seq(from=-11000, to=0, by=1000)
# Plot
plot(raster.crop, col=blue.col(22), breaks=blue.br, legend=FALSE)
plot(coords, add = TRUE, pch = 19, col = 3, cex=1)

# Create Risk areas circle shapes
df_risk_areas <- data.frame(long = as.numeric(df_vme_index$LongitudeMid),
                            lat = as.numeric(df_vme_index$LatitudeMid))
sf_risk_areas <- st_as_sf(df_risk_areas, coords = c("long", "lat"), crs="+proj=longlat +datum=WGS84")
# Buffer circles by 1 nmile radius
sf_circles <- st_buffer(sf_risk_areas, dist = 1852)
sf_circles <- SOproj(sf_circles)
sf_circles$vme <- as.numeric(df_vme_index$vme_index)
yellow.col <- colorRampPalette(c("yellow", "orange"))
yellow.br <- seq(from=1, to=max(n_index_categories), by=1)
plot(sf_circles, add=TRUE, col=yellow.col(n_index_categories), legend=TRUE)
# Final resolution
if (resolution.raster < res(raster.crop)[1]) {
  resample_fact <- res(raster.crop)[1] / resolution.raster
  raster.crop.res <- disaggregate(raster.crop, fact=resample_fact, method='bilinear')
} else if (resolution.raster > res(raster.crop)[1]) {
  resample_fact <- resolution.raster / res(raster.crop)[1]
  raster.crop.res <- aggregate(raster.crop, fact=resample_fact, method='bilinear')
} else {
  raster.crop.res <- raster.crop
}
crs(raster.crop.res) <- "+proj=longlat +datum=WGS84"
# Rasterize
raster_circles <- rasterize(sf_circles, raster.crop.res, field="vme", fun="max") #, fun=mean)
unique(values(raster_circles))
unique(sf_circles$vme)
crs(raster_circles) <- "+proj=longlat +datum=WGS84"
# TODO: MOve between sf_circles and raster_circles
yellow.col <- colorRampPalette(c("yellow", "orange"))
yellow.br <- seq(from=1, to=max(n_index_categories)+1, by=1)
par(mar=c(2, 2, 2, 2))
plot(raster.crop.res, col=blue.col(22), breaks=blue.br, legend=FALSE)
plot(rasterToPolygons(raster_circles), add=TRUE, col=yellow.col(n_index_categories), breaks=yellow.br, legend=FALSE, border="black")
legend("topright", legend=unique(raster_circles), col=yellow.col(n_index_categories), pch=15, border="black", title = "VME index")
plot(sf_circles, add=TRUE, col=yellow.col(6), legend=TRUE)

par(mar=c(2, 2, 2, 2))
plot(raster.crop.res, col=blue.col(22), breaks=blue.br, legend=FALSE)
plot(rasterToPolygons(raster_circles), add=TRUE, col=yellow.col(n_index_categories), breaks=yellow.br, legend=FALSE, border="black")




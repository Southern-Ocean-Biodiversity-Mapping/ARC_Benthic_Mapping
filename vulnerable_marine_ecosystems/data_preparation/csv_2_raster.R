################################################################################
#                                 TODO
################################################################################
# 

################################################################################
#                                 INIT
################################################################################
# Packages
library(CCAMLRGIS)
library(raster)
library(SOmap)
library(sp)


# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/data_preparation")

# Path towards BIIGLE report
path_biigle <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210909_biigle_vme_cover.csv"

# Raster resolution
resolution_raster <- 500

# Output filename
path_out <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210914_raster_cover/raster_"

# Save raster as file
save_raster = FALSE

################################################################################
#                                 XX
################################################################################
# Read BIIGLE data
cat("Reading data:", path_biigle, "...")
df_biigle <- read.csv(path_biigle)
# Remove filename column
cat("Cleaning data ...")
df_biigle <- df_biigle[,!names(df_biigle) %in% c("filename")]
# Convert survey ID to numeric
df_biigle$survey <- as.numeric(df_biigle$survey)
# Create new variable
pts <- df_biigle
# Projection
coordinates(pts) <- c("longitude", "latitude")
projection(pts) <- "+proj=longlat +datum=WGS84"
pts <- SOproj(pts)

cat("Cropping and resampling (",
    resolution_raster, "x", resolution_raster, "m ) the reference raster...")
# Crop reference raster according to pts extent
raster_ref <- crop(SmallBathy, pts)
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
cat("Rasterizing data ...")
raster_bio <- rasterize(pts,
                        raster_ref,
                        colnames(df_biigle)[!names(df_biigle) %in% c("longitude", "latitude", "survey")],
                        fun=sum)
raster_survey <- rasterize(pts,
                           raster_ref,
                           c("survey"),
                           fun="first")
names(raster_survey) <- "survey"

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


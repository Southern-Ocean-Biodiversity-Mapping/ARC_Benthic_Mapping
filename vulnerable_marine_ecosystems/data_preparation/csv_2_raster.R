################################################################################
#                                 TODO
################################################################################
# - Account for BIIGLE 254
# - Add other covariates: date, image quality, acqusition method

################################################################################
#                                 INIT
################################################################################
# Packages
library(CCAMLRGIS)
library(raster)
library(SOmap)
library(sp)
library(proj4)
library(dplyr)

# Set working directory
setwd("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/data_preparation")

# Path towards data
path_data <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210920"
fname_bio <- "bio_data.csv"
fname_src <- "bio_data_source.csv"

# Raster resolution
resolution_raster <- 500

# Output filename
path_out <- "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20210930_raster_cover/raster_"

# Save raster as file
save_raster = FALSE
# Save df as file
save_df = TRUE

################################################################################
#                                 READ INPUT DATA
################################################################################
# Create output folder
dir.create(dirname(path_out), showWarnings = FALSE)

# Read Bio data
path_bio_data <- paste(path_data, fname_bio, sep="/")
cat("\nReading data:", path_bio_data, "...\n")
df <- read.csv(path_bio_data)
names(df) <- gsub(x = names(df), pattern = "\\.", replacement = "--")

# Read source data
path_src <- paste(path_data, fname_src, sep="/")
cat("\nReading data:", path_src, "...\n")
df_src <- read.csv(path_src)
df_src$morpho_taxon <- gsub(x = df_src$morpho_taxon, pattern = "-", replacement = "--")

# Get taxa names for each source
list_taxa_coralnet = df_src[df_src$source == "coralnet", "morpho_taxon"]
list_taxa_biigle839 = df_src[df_src$source == "biigle839", "morpho_taxon"]
list_taxa_biigle254 = df_src[df_src$source == "biigle254", "morpho_taxon"]

################################################################################
#                                 PREPARE RASTER
################################################################################

# Get reference
raster_ref <- SmallBathy

# Projection
cat("\nProjection:", proj4string(raster_ref), "...\n")
df$proj_coord_x <- project(df[,c("longitude", "latitude")],
                           proj=crs(raster_ref))$x
df$proj_coord_y <- project(df[,c("longitude", "latitude")],
                            proj=crs(raster_ref))$y
# Create new variable
pts <- df
coordinates(pts) <- c("proj_coord_x", "proj_coord_y")

# Crop reference raster according to pts extent
raster_ref <- crop(raster_ref, pts)
cat("\nCropping and resampling (",
    resolution_raster, "x", resolution_raster, "m ) the reference raster...\n")
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

# Get cell IDs
df$cellID <- raster::extract(raster_ref,
                     df[,c("proj_coord_x", "proj_coord_y")],
                     cellnumbers=TRUE)[,1]
list_cells <- unique(df$cellID)

# Compute percentage cover per cell
cat("\nComputing percentage cover per raster cell ...\n")
list_colnames_not_sum <- c("survey", "longitude", "latitude", "filename", "proj_coord_x", "proj_coord_y")
df_to_sum <- df[colnames(df)[!names(df) %in% list_colnames_not_sum]]
df_not_to_sum <- df[, c("cellID", list_colnames_not_sum)]
# Remove duplicates
df_not_to_sum <- df_not_to_sum[!duplicated(df_not_to_sum$cellID),]
# Sum Bio data within each raster cell
df_sum <- as.data.frame(df_to_sum
                        %>% group_by(cellID)
                        %>% summarise(across(everything(),
                                             function(x,...){if (!all(is.na(x))){sum(na.omit(x))} else{NA}})))

# Aggregate for CCAMLR taxa
list_taxonomic <- colnames(df_sum)[!(colnames(df_sum) %in% c("n_annotation", "area_pix", "area", "cellID"))]
list_taxonomic <- unique(list_taxonomic %>% strsplit( "--" ) %>% sapply(tail, 1))
df_sum_taxonomic <- df_sum[, c("n_annotation", "area_pix", "area", "cellID")]
for (taxa in list_taxonomic) {
  list_morpho_taxa_cur <- colnames(df_sum)[endsWith(colnames(df_sum), taxa)]
  if (length(list_morpho_taxa_cur) > 1) {
    df_sum_taxonomic[, taxa] <- rowSums(df_sum[, list_morpho_taxa_cur], na.rm=TRUE)
  } else {
    df_sum_taxonomic[, taxa] <- df_sum[, list_morpho_taxa_cur]
  }
}

################################################################################
#                                 COMPUTE PERCENT COVER
################################################################################

# Normalise to compute percentage cover
df_sum[, list_taxa_coralnet] <- df_sum[, list_taxa_coralnet] * 100. / df_sum[, "n_annotation"]
df_sum[, list_taxa_biigle839] <- df_sum[, list_taxa_biigle839] * 100. / df_sum[, "area_pix"]
#df_sum[, list_taxa_biigle254] <- df_sum[, list_taxa_biigle254] * 100. / df_sum[, "area_pix"]
df_sum_taxonomic[, c("bryozoans", "porifera", "stylasterids")] <- df_sum_taxonomic[, c("bryozoans", "porifera", "stylasterids")] * 100. / df_sum_taxonomic[, "n_annotation"]
list_taxonomic_biigle <- list_taxonomic[!(list_taxonomic %in% c("bryozoans", "porifera", "stylasterids"))]
df_sum_taxonomic[, list_taxonomic_biigle] <- df_sum_taxonomic[, list_taxonomic_biigle] * 100. / df_sum_taxonomic[, "area_pix"]
#df_sum_taxonomic[, list_taxa_biigle254] <- df_sum_taxonomic[, list_taxa_biigle254] * 100. / df_sum_taxonomic[, "area_pix"]

################################################################################
#                                CLEAN DATA FRAMES
################################################################################

# Join dataframes
df_to_rasterize <- left_join(df_sum, df_not_to_sum)
df_to_rasterize$survey <- as.numeric(as.factor(df_to_rasterize$survey))

rm_row_na <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Remove rows where BIIGLE and Coralnet are not both present
df_to_rasterize <- rm_row_na(df_to_rasterize, list_taxa_biigle839)
df_to_rasterize <- rm_row_na(df_to_rasterize, list_taxa_coralnet)
df_to_rasterize <- rm_row_na(df_to_rasterize, list_taxa_biigle254)

df_to_rasterize_ccamlr <- left_join(df_sum_taxonomic, df_not_to_sum)
df_to_rasterize_ccamlr$survey <- as.numeric(as.factor(df_to_rasterize_ccamlr$survey))

list_taxa_biigle839_ccamlr <- unique(sapply(strsplit(list_taxa_biigle839,"--"), `[`, 2))
list_taxa_coralnet_ccamlr <- unique(sapply(strsplit(list_taxa_coralnet,"--"), `[`, 2))
#list_taxa_biigle254_ccamlr <- unique(sapply(strsplit(list_taxa_biigle254,"--"), `[`, 2))

df_to_rasterize_ccamlr <- rm_row_na(df_to_rasterize_ccamlr, list_taxa_biigle839_ccamlr)
df_to_rasterize_ccamlr <- rm_row_na(df_to_rasterize_ccamlr, list_taxa_coralnet_ccamlr)
#df_to_rasterize_ccamlr <- rm_row_na(df_to_rasterize_ccamlr, list_taxa_biigle254_ccamlr)

################################################################################
#                                CLEAN DATA FRAMES
################################################################################

if (save_df) {
  cat("\nSaving dataframe:", path_out, "...\n")
  write.csv(df_to_rasterize,
            paste0(path_out, "data.csv"),
            row.names = FALSE)
}

if (save_df) {
  cat("\nSaving CCAMLR dataframe:", path_out, "...\n")
  write.csv(df_to_rasterize_ccamlr,
            paste0(path_out, "data_ccamlr.csv"),
            row.names = FALSE)
}

################################################################################
#                                RASTERIZE
################################################################################

pts_to_rasterize <- df_to_rasterize
coordinates(pts_to_rasterize) <- c("proj_coord_x", "proj_coord_y")

# Rasterise data
cat("\nRasterizing data ...\n")
list_colnames_to_rasterise <- c("survey",
                                "area",
                                list_taxa_biigle254,
                                list_taxa_biigle839,
                                list_taxa_coralnet)
raster_biodata <- rasterize(pts_to_rasterize,
                            raster_ref,
                            colnames(df_to_rasterize)[names(df_to_rasterize) %in% list_colnames_to_rasterise],
                            fun='first')

#cellnb.pos <- Which(raster_biodata$octocorals_fleshy_mushroom..alcyonacea > 0, cells = TRUE)
cat("\nNumber of layers:", nlayers(raster_biodata), "...\n")

# Save raster
if (save_raster) {
  cat("\nSaving raster:", path_out, "...\n")
  writeRaster(raster_biodata,
              filename=path_out,
              suffix=names(raster_biodata),
              bylayer=TRUE,
              format="GTiff",
              overwrite=TRUE)
}

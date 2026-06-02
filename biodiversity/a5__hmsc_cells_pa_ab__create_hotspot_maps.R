## a5__hmsc_cells_pa_ab__create_hotspot_maps.R

############################################################
# HMSC BIODIVERSITY HOTSPOT ANALYSIS PIPELINE
# ----------------------------------------------------------
# This script:
#
# 1) Loads predictive species distribution maps for multiple models
# 2) Calculates derived biodiversity metrics:
#    - Richness (sum of species probabilities)
#    - Total abundance
# 3) Converts these to percentile surfaces
# 4) Identifies hotspots (top X% of values)
# 5) Combines richness and abundance into biodiversity hotspots
# 6) Calculates multi-model consensus maps
#
# Designed for:
# - Circum-Antarctic seafloor biodiversity predictions
# - Outputs from HMSC modelling pipeline (assembled rasters)
#
############################################################


############################
# 1) SETUP
############################

library(terra)
library(dplyr)

# User environment setup
usr <- "VM"
source("0_SourceFile.R")

# load custom functions
source("0_Functions.R")

# Spatial resolution label (for consistency with upstream pipeline)
res <- "2km"

# Model identifiers (8 alternative environmental model formulations)
model_ids <- c(
  "npp_cafe","npp_cbpm","npp_eppl","npp_vpmg",
  "fam_cafe","fam_cbpm","fam_eppl","fam_vpmg"
)

# Source dropbox folder with original files
src_base <- file.path(usr.dropbox.dir, "data_products/predictive_maps/circum_antarctic")
# Local directory with copied files for faster processing
dst_base <- file.path(usr.main.dir, "4_model_prediction/copy_of_predictive_maps/circum_antarctic")

# ############################################################
# # COPY PREDICTIVE MAPS TO LOCAL DISK
# # ----------------------------------------------------------
# # Copy all predictive map files from the synced
# # storage location to a local directory, preserving the folder
# # structure for each model.
# ############################################################
# # Create destination root if needed
# dir.create(dst_base, recursive = TRUE, showWarnings = FALSE)
# # Find all files under source directory
# # Change pattern if only certain file types should be copied
# all_files <- list.files(src_base,
#   recursive = TRUE, full.names = TRUE, include.dirs = FALSE
# )
# # Keep only GeoTIFF files
# all_files <- all_files[grepl("\\.tif$", all_files, ignore.case = TRUE)]
# 
# # Build matching destination paths while preserving subfolder structure
# rel_paths <- substring(all_files, nchar(src_base) + 2)
# dst_files <- file.path(dst_base, rel_paths)
# 
# # Create all required destination subdirectories
# dst_dirs <- unique(dirname(dst_files))
# invisible(lapply(dst_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
# 
# # Copy files
# ok <- file.copy(from = all_files, to = dst_files, overwrite = FALSE, copy.mode = TRUE)
# # Report results
# message(sprintf("Copied %d of %d files", sum(ok), length(ok)))
# # Show any failures
# if (any(!ok)) {
#   message("The following files were not copied:")
#   print(all_files[!ok])
# }
# ## SYNC LOCAL FILE MTIMES TO MATCH SOURCE FILES
# # List tif files
# src_files <- list.files(src_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
# dst_files <- list.files(dst_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
# # Relative paths
# src_rel <- substring(src_files, nchar(src_base) + 2)
# dst_rel <- substring(dst_files, nchar(dst_base) + 2)
# # Match source files to local files
# m <- match(src_rel, dst_rel)
# # Keep only files that exist in both locations
# keep <- !is.na(m)
# # Apply source mtime to destination file
# src_mtime <- file.info(src_files[keep])$mtime
# for (i in seq_along(src_mtime)) {
#   Sys.setFileTime(dst_files[m[keep][i]], src_mtime[i])
# }
# message(sprintf("Updated timestamps for %d files", sum(keep)))
# ############################################################

############################################################
# CHECK WHETHER LOCAL COPIED FILES MATCH SOURCE FILES
# ----------------------------------------------------------
# Compares .tif files between source and local copy using:
# - relative file path
# - file size
# - modification time
#
# Returns a table of files that are:
# - missing locally
# - missing in source
# - different in size
# - different in modification time
############################################################
# List tif files in each location
src_files <- list.files(src_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
dst_files <- list.files(dst_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
# Convert to relative paths so files can be matched independent of root folder
src_rel <- substring(src_files, nchar(src_base) + 2)
dst_rel <- substring(dst_files, nchar(dst_base) + 2)
# Build metadata tables
src_info <- data.frame(
  rel_path = src_rel,
  src_file = src_files,
  src_size = file.info(src_files)$size,
  src_mtime = file.info(src_files)$mtime,
  stringsAsFactors = FALSE
)
dst_info <- data.frame(
  rel_path = dst_rel,
  dst_file = dst_files,
  dst_size = file.info(dst_files)$size,
  dst_mtime = file.info(dst_files)$mtime,
  stringsAsFactors = FALSE
)
# Merge source and destination tables by relative path
cmp <- merge(src_info, dst_info, by = "rel_path", all = TRUE)
# Identify problems
cmp$missing_in_local  <- is.na(cmp$dst_file)
cmp$missing_in_source <- is.na(cmp$src_file)
cmp$size_diff <- FALSE
cmp$mtime_diff <- FALSE
both_present <- !cmp$missing_in_local & !cmp$missing_in_source
cmp$size_diff[both_present] <- cmp$src_size[both_present] != cmp$dst_size[both_present]
cmp$mtime_diff[both_present] <- cmp$src_mtime[both_present] != cmp$dst_mtime[both_present]
# Files that differ in any way
cmp_issues <- cmp[
  cmp$missing_in_local |
    cmp$missing_in_source |
    cmp$size_diff |
    cmp$mtime_diff,
]
# Report result
if (nrow(cmp_issues) == 0) {
  message("Local predictive maps match source files (by path, size, and modification time).")
} else {
  message("Differences detected between source and local predictive maps:")
  print(cmp_issues[, c("rel_path", "missing_in_local", "missing_in_source", "size_diff", "mtime_diff")])
}

###################################################################
# Base directory containing predictive maps for all models
if (nrow(cmp_issues) == 0) {pred_base_dir <- dst_base}
# pred_base_dir <- file.path(
#   usr.dropbox.dir, "data_products/predictive_maps/circum_antarctic"
# )

############################
# 2) PROCESS EACH MODEL
############################
# Store outputs for all models in a list
model_outputs <- list()

for (nm in model_ids) {
  t0 <- Sys.time()
  
  message("======================================")
  message("Processing model: ", nm)
  
  # Directory containing predictions for this model
  model_dir <- file.path(pred_base_dir, paste0("hmsc_with_",nm))
  
  ############################################################
  # 2.1 LOAD SPECIES-LEVEL PREDICTIONS
  ############################################################
  message("loading files")
  # Each species has its own raster file:
  #   PA_<species>.tif with 5 layers (mean, median, SE, etc.)
  pa_files <- list.files(model_dir,
    pattern = "^PA_.*\\.tif$", full.names = TRUE
  )
  
  if (length(pa_files) == 0) {
    warning("No PA files found for model: ", nm)
    next
  }
  
  # Stack all species maps, extracting the median prediction layer for each species
  pa_median <- rast(pa_files, lyrs=2)
  
  ############################################################
  # 2.2 CALCULATE SPECIES RICHNESS
  ############################################################
  message("extracting richness")
  # Richness = sum of species occurrence probabilities
  richness <- sum(pa_median, na.rm = TRUE)
  
  ############################################################
  # 2.3 LOAD TOTAL ABUNDANCE
  ############################################################
  message("extracting abundance")
  
  # Abundance raster also contains multiple layers
  # Extract median layer (layer 2)
  ab_file <- file.path(model_dir, "total_abundance.tif")
  
  if (!file.exists(ab_file)) {
    warning("No abundance file found for model: ", nm)
    next
  }
  abundance <- rast(ab_file)[[2]]
  
  ############################################################
  # 2.4 CALCULATE PERCENTILE SURFACES
  ############################################################
  message("calculating percentiles")
  
  # Convert raw values into relative rank (0–1)
  richness_pct  <- calc_percentile_raster(richness)
  abundance_pct <- calc_percentile_raster(abundance)
  
  ############################################################
  # 2.5 IDENTIFY HOTSPOTS
  ############################################################
  
  # Define hotspots as top 10% of values
  richness_hot  <- topX.calc(richness,  0.10)
  abundance_hot <- topX.calc(abundance, 0.10)
  
  # Combined biodiversity hotspot:
  # cells that are hotspots for BOTH richness and abundance
  combined_hot <- (richness_hot + abundance_hot) == 2
  
  ############################################################
  # 2.6 STORE RESULTS
  ############################################################
  
  model_outputs[[nm]] <- list(
    richness       = richness,
    abundance      = abundance,
    richness_pct   = richness_pct,
    abundance_pct  = abundance_pct,
    richness_hot   = richness_hot,
    abundance_hot  = abundance_hot,
    combined_hot   = combined_hot
  )
  message(sprintf("Elapsed time: %.2f sec", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}


############################
# 3) WRITE OUTPUTS PER MODEL
############################

for (nm in names(model_outputs)) {
  
  message("Saving outputs for model: ", nm)
  
  # Create output directory
  out_dir <- file.path(src_base, paste0("hmsc_with_", nm), "hotspots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  obj <- model_outputs[[nm]]
  
  # Save continuous maps
  writeRaster(obj$richness,  file.path(out_dir, "richness.tif"), overwrite = TRUE)
  writeRaster(obj$abundance, file.path(out_dir, "abundance.tif"), overwrite = TRUE)
  
  # Save percentiles
  writeRaster(obj$richness_pct,  file.path(out_dir, "richness_percentiles.tif"), overwrite = TRUE)
  writeRaster(obj$abundance_pct, file.path(out_dir, "abundance_percentiles.tif"), overwrite = TRUE)
  
  # Save hotspot masks
  writeRaster(obj$richness_hot,  file.path(out_dir, "richness_hotspots.tif"), overwrite = TRUE)
  writeRaster(obj$abundance_hot, file.path(out_dir, "abundance_hotspots.tif"), overwrite = TRUE)
  writeRaster(obj$combined_hot,  file.path(out_dir, "combined_hotspots.tif"), overwrite = TRUE)
}


############################
# 4) MULTI-MODEL CONSENSUS
############################

# Combine hotspot maps across models to assess agreement

message("Calculating multi-model consensus...")

# --- richness consensus ---
rich_stack <- rast(lapply(model_outputs, function(x) x$richness_hot))
rich_consensus <- sum(rich_stack, na.rm = TRUE)

# --- abundance consensus ---
ab_stack <- rast(lapply(model_outputs, function(x) x$abundance_hot))
ab_consensus <- sum(ab_stack, na.rm = TRUE)

# --- combined biodiversity consensus ---
comb_stack <- rast(lapply(model_outputs, function(x) x$combined_hot))
comb_consensus <- sum(comb_stack, na.rm = TRUE)


############################
# 5) SAVE CONSENSUS MAPS
############################

cons_dir <- file.path(src_base, "hmsc_model_consensus")
dir.create(cons_dir, showWarnings = FALSE)

writeRaster(rich_consensus,
            file.path(cons_dir, "richness_hotspot_consensus.tif"),
            overwrite = TRUE)

writeRaster(ab_consensus,
            file.path(cons_dir, "abundance_hotspot_consensus.tif"),
            overwrite = TRUE)

writeRaster(comb_consensus,
            file.path(cons_dir, "combined_hotspot_consensus.tif"),
            overwrite = TRUE)


# ############################
# # 6) DEFINE ROBUST HOTSPOTS
# ############################
# 
# # Identify cells consistently identified across models
# # Example: at least 5 out of 8 models agree
# 
# robust_hotspots <- comb_consensus >= 5
# 
# writeRaster(
#   robust_hotspots,
#   file.path(cons_dir, "combined_hotspots_robust_5models.tif"),
#   overwrite = TRUE
# )




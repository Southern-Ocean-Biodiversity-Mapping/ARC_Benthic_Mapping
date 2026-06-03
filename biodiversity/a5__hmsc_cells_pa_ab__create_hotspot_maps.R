############################################################
# HMSC BIODIVERSITY HOTSPOT ANALYSIS PIPELINE
# ----------------------------------------------------------
# This script:
#
# 1) Loads predictive species distribution maps for multiple models
# 2) Derives biodiversity metrics and uncertainty layers
# 3) Calculates circumpolar percentile surfaces
# 4) Calculates planning-domain-specific percentile surfaces
# 5) Creates hotspot classes from raw-value hotspot thresholds
# 6) Creates richness-only and abundance-only contrast layers
# 7) Writes hotspot threshold rasters and class rasters
# 8) Creates polygons for the top 10% threshold on richness,
#    abundance and biodiversity
#
# Designed for:
# - Circum-Antarctic seafloor biodiversity predictions
# - Outputs assembled from the HMSC modelling pipeline
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
source("0_Functions.R")

# Spatial resolution label
res <- "2km"

# Model identifiers
model_ids <- c(
  "npp_cafe", "npp_cbpm", "npp_eppl", "npp_vpmg",
  "fam_cafe", "fam_cbpm", "fam_eppl", "fam_vpmg"
)

# Source dropbox folder with original files
src_base <- file.path(
  usr.dropbox.dir,
  "data_products/predictive_maps/circum_antarctic"
)

# Local directory with copied files for faster processing
dst_base <- file.path(
  usr.main.dir,
  "4_model_prediction/copy_of_predictive_maps/circum_antarctic"
)

# Planning-domain definitions
domain_ids <- c(1, 3, 4, 7, 8, 9)
domain_name_field <- "Name"

domain_shp <- file.path(
  usr.dropbox.dir,
  "data_environmental/raw/ccamlr_data_2025/geographical_data/mpapd/CCAMLR_MPAPD_EPSG6932.shp"
)

# # Beta-diversity directory
# # The script searches for files containing the model identifier.
# beta_dir <- file.path(
#   usr.dropbox.dir,
#   "data_products/predictive_maps/circum_antarctic_beta_diversity"
# )
#
# # Beta-diversity layer index to use if the file has multiple layers
# beta_layer_index <- 1


############################
# 2) HELPER FUNCTIONS
############################

# ----------------------------------------------------------
# Function: positive_difference
# ----------------------------------------------------------
# Returns the positive component of a raster difference
# Negative values are set to zero
# ----------------------------------------------------------
positive_difference <- function(ra1, ra2) {
  out <- ra1 - ra2
  out[out < 0] <- 0
  out[is.na(ra1)] <- NA
  out[is.na(ra2)] <- NA
  out
}

# ----------------------------------------------------------
# Function: safe_divide_raster
# ----------------------------------------------------------
# Divides one raster by another
# Zero or negative denominator values are set to NA
# ----------------------------------------------------------
safe_divide_raster <- function(num, den) {
  den2 <- den
  den2[den2 <= 0] <- NA
  out <- num / den2
  out[is.infinite(out)] <- NA
  out
}

# ----------------------------------------------------------
# Function: write_hotspot_stack
# ----------------------------------------------------------
# Writes a standard hotspot stack with named layers
# ----------------------------------------------------------
write_hotspot_stack <- function(hot_stack, filename, overwrite = TRUE) {
  names(hot_stack) <- c("top25", "top10", "top5", "top1", "intensity")
  writeRaster(hot_stack, filename = filename, overwrite = overwrite)
}

# ----------------------------------------------------------
# Function: hotspot_stack_to_class
# ----------------------------------------------------------
# Extracts the hotspot-class raster from a hotspot stack
# Classes correspond to the nested hotspot thresholds:
#   0 = outside top 25%
#   1 = top 25%
#   2 = top 10%
#   3 = top 5%
#   4 = top 1%
# ----------------------------------------------------------
hotspot_stack_to_class <- function(hot_stack) {
  out <- hot_stack[["intensity"]]
  names(out) <- "class"
  out
}

# ----------------------------------------------------------
# Function: empty_named_raster
# ----------------------------------------------------------
# Creates a single-layer NA raster with a specified name
# ----------------------------------------------------------
empty_named_raster <- function(template, layer_name) {
  out <- rast(template)
  out[] <- NA
  names(out) <- layer_name
  out
}

# # ----------------------------------------------------------
# # Function: find_beta_file
# # ----------------------------------------------------------
# # Finds an optional beta-diversity file for a model
# # ----------------------------------------------------------
# find_beta_file <- function(model_id, beta_dir) {
#   if (!dir.exists(beta_dir)) {
#     return(NA_character_)
#   }
#
#   beta_files <- list.files(beta_dir, pattern = "\\.tif$", full.names = TRUE)
#   if (length(beta_files) == 0) {
#     return(NA_character_)
#   }
#
#   hit <- grep(model_id, basename(beta_files), ignore.case = TRUE)
#
#   if (length(hit) == 1) {
#     return(beta_files[hit])
#   }
#
#   if (length(beta_files) == 1) {
#     return(beta_files[1])
#   }
#
#   NA_character_
# }


############################
# 3) CHECK LOCAL COPY OF PREDICTIVE MAPS
############################

# # COPY PREDICTIVE MAPS TO LOCAL DISK
# # ----------------------------------------------------------
# # Copy all predictive map files from the synced storage
# # location to a local directory, preserving the folder
# # structure for each model.
# ############################################################
# dir.create(dst_base, recursive = TRUE, showWarnings = FALSE)
#
# all_files <- list.files(
#   src_base,
#   recursive = TRUE,
#   full.names = TRUE,
#   include.dirs = FALSE
# )
#
# all_files <- all_files[grepl("\\.tif$", all_files, ignore.case = TRUE)]
#
# rel_paths <- substring(all_files, nchar(src_base) + 2)
# dst_files <- file.path(dst_base, rel_paths)
#
# dst_dirs <- unique(dirname(dst_files))
# invisible(lapply(dst_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
#
# ok <- file.copy(from = all_files, to = dst_files, overwrite = FALSE, copy.mode = TRUE)
# message(sprintf("Copied %d of %d files", sum(ok), length(ok)))
#
# if (any(!ok)) {
#   message("The following files were not copied:")
#   print(all_files[!ok])
# }
#
# src_files <- list.files(src_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
# dst_files <- list.files(dst_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
#
# src_rel <- substring(src_files, nchar(src_base) + 2)
# dst_rel <- substring(dst_files, nchar(dst_base) + 2)
#
# m <- match(src_rel, dst_rel)
# keep <- !is.na(m)
#
# src_mtime <- file.info(src_files[keep])$mtime
# for (i in seq_along(src_mtime)) {
#   Sys.setFileTime(dst_files[m[keep][i]], src_mtime[i])
# }
#
# message(sprintf("Updated timestamps for %d files", sum(keep)))

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

src_files <- list.files(
  src_base,
  pattern = "\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

dst_files <- list.files(
  dst_base,
  pattern = "\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

src_rel <- substring(src_files, nchar(src_base) + 2)
dst_rel <- substring(dst_files, nchar(dst_base) + 2)

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

cmp <- merge(src_info, dst_info, by = "rel_path", all = TRUE)

cmp$missing_in_local  <- is.na(cmp$dst_file)
cmp$missing_in_source <- is.na(cmp$src_file)
cmp$size_diff <- FALSE
cmp$mtime_diff <- FALSE

both_present <- !cmp$missing_in_local & !cmp$missing_in_source
cmp$size_diff[both_present]  <- cmp$src_size[both_present]  != cmp$dst_size[both_present]
cmp$mtime_diff[both_present] <- cmp$src_mtime[both_present] != cmp$dst_mtime[both_present]

cmp_issues <- cmp[
  cmp$missing_in_local |
    cmp$missing_in_source |
    cmp$size_diff |
    cmp$mtime_diff,
]

if (nrow(cmp_issues) == 0) {
  message("Local predictive maps match source files (by path, size, and modification time).")
  pred_base_dir <- dst_base
} else {
  message("Differences detected between source and local predictive maps.")
  message("Using source predictive maps.")
  pred_base_dir <- src_base
}


############################
# 4) LOAD PLANNING DOMAINS
############################
if (!file.exists(domain_shp)) {
  stop("Planning-domain shapefile not found: ", domain_shp)
}
domain_raw <- vect(domain_shp)

# Assign source CRS if missing
domain_crs <- crs(domain_raw)
if (is.na(domain_crs) || domain_crs == "") {
  crs(domain_raw) <- "EPSG:6932"
}

# Domain raster alignment is based on the first successfully loaded model
domain_template_raster <- NULL


############################
# 5) PROCESS EACH MODEL
############################
for (nm in model_ids) {
  t0 <- Sys.time()
  
  message("======================================")
  message("Processing model: ", nm)
  message("======================================")
  
  model_dir <- file.path(pred_base_dir, paste0("hmsc_with_", nm))
  if (!dir.exists(model_dir)) {
    warning("Model directory not found: ", model_dir)
    next
  }
  
  out_dir <- file.path(src_base, paste0("hmsc_with_", nm))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_dir_cont   <- file.path(out_dir, "continuous")
  out_dir_pct    <- file.path(out_dir, "percentiles")
  out_dir_hot    <- file.path(out_dir, "hotspots")
  out_dir_domain <- file.path(out_dir, "domains")
  out_dir_poly   <- file.path(out_dir, "polygons")
  
  dir.create(out_dir_cont,   recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_pct,    recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_hot,    recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_domain, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_poly,   recursive = TRUE, showWarnings = FALSE)
  
  ############################################################
  # 5.1 LOAD SPECIES-LEVEL PREDICTIONS
  ############################################################
  
  pa_files <- list.files(
    model_dir,
    pattern = "^PA_.*\\.tif$",
    full.names = TRUE
  )
  
  if (length(pa_files) == 0) {
    warning("No PA files found for model: ", nm)
    next
  }
  
  pa_mean   <- rast(pa_files, lyrs = 1)
  pa_median <- rast(pa_files, lyrs = 2)
  pa_se     <- rast(pa_files, lyrs = 3)
  pa_5      <- rast(pa_files, lyrs = 4)
  pa_95     <- rast(pa_files, lyrs = 5)
  # species_names <- sub("^PA_|\\.tif$", "", basename(pa_files))
  # names(pa_mean)   <- species_names
  # names(pa_median) <- species_names
  # names(pa_se)     <- species_names
  # names(pa_5)      <- species_names
  # names(pa_95)     <- species_names
  
  ab_file <- file.path(model_dir, "total_abundance.tif")
  if (!file.exists(ab_file)) {
    warning("No total abundance file found for model: ", nm)
    next
  }
  ab_stack <- rast(ab_file)
  
  # Convert abundance predictions from 540-point scale to 100-point scale
  abundance_mean   <- ab_stack[[1]] / 5.4
  abundance_median <- ab_stack[[2]] / 5.4
  abundance_se     <- ab_stack[[3]] / 5.4
  abundance_5      <- ab_stack[[4]] / 5.4
  abundance_95     <- ab_stack[[5]] / 5.4
  # Cap median abundance at 100 to retain interpretation as percent cover
  abundance_median_100 <- abundance_median
  abundance_median_100[abundance_median_100 > 100] <- 100
  names(abundance_median_100) <- "abundance_median_100"
  
  ############################################################
  # 5.2 DERIVE BIODIVERSITY METRICS
  ############################################################

  # Summed richness from species-level occurrence probabilities
  richness_mean   <- sum(pa_mean, na.rm = TRUE)
  richness_median <- sum(pa_median, na.rm = TRUE)
  
  # Bootstrap-based uncertainty for richness
  richness_boot <- bootstrap_richness_from_species_medians(
    pa_median = pa_median,
    n_boot = 1000,
    chunk_size = 10000,
    seed = 2
  )
  
  richness_boot_median <- richness_boot$median
  richness_boot_sd <- richness_boot$sd
  
  # Ratio of richness to abundance, and vice-versa
  richness_per_abundance <- safe_divide_raster(richness_median, abundance_median_100)
  abundance_per_richness <- safe_divide_raster(abundance_median_100, richness_median)
  
  ############################################################
  # 5.3 BETA-DIVERSITY INPUT
  ############################################################
  
  # beta_file <- find_beta_file(nm, beta_dir)
  # beta_median <- NULL
  #
  # if (!is.na(beta_file) && file.exists(beta_file)) {
  #   beta_stack <- rast(beta_file)
  #   if (nlyr(beta_stack) >= beta_layer_index) {
  #     beta_median <- beta_stack[[beta_layer_index]]
  #   }
  # }
  
  ############################################################
  # 5.4 PREPARE DOMAIN GEOMETRY
  ############################################################
  if (is.null(domain_template_raster)) {
    domain_template_raster <- richness_median
  }
  domain_dat <- project(domain_raw, crs(richness_median))
  
  ############################################################
  # 5.5 CIRCUMPOLAR PERCENTILE SURFACES
  ############################################################
  message("Calculating circumpolar percentile surfaces")
  
  richness_pct  <- calc_percentile_rank_raster(richness_median)
  abundance_pct <- calc_percentile_rank_raster(abundance_median_100)
  
  # Combined biodiversity score based on matched percentile scales
  biodiversity_raw <- richness_pct + abundance_pct
  biodiversity_pct <- calc_percentile_rank_raster(biodiversity_raw)
  
  # Fine contrasts based on raw percentile differences
  richness_only_pct  <- positive_difference(richness_pct, abundance_pct)
  abundance_only_pct <- positive_difference(abundance_pct, richness_pct)
  
  richness_per_abundance_pct <- calc_percentile_rank_raster(richness_per_abundance)
  abundance_per_richness_pct <- calc_percentile_rank_raster(abundance_per_richness)
  
  # if (!is.null(beta_median)) {
  #   beta_pct <- calc_percentile_rank_raster(beta_median)
  # }
  
  ############################################################
  # 5.6 HOTSPOT THRESHOLDS AND CLASSES
  ############################################################
  message("Creating hotspot thresholds and classes")
  
  # Hotspot thresholds are calculated on the raw analysis surfaces
  richness_hot   <- calc_hotspot_products(richness_median)
  abundance_hot  <- calc_hotspot_products(abundance_median_100)
  biodiversity_hot <- calc_hotspot_products(biodiversity_raw)

  richness_only_class  <- positive_difference(richness_hot$intensity, abundance_hot$intensity)
  abundance_only_class <- positive_difference(abundance_hot$intensity, richness_hot$intensity)
  
  # if (!is.null(beta_median)) {
  #   beta_hot <- calc_hotspot_products(beta_median)
  #   beta_class <- hotspot_stack_to_class(beta_hot)
  # }
  
  ############################################################
  # 5.10 DOMAIN-SPECIFIC PERCENTILES AND HOTSPOTS
  ############################################################
  message("Creating domain-specific percentile and hotspot products")
  
  domain_richness_pct_list <- vector("list", length(domain_ids))
  domain_abundance_pct_list <- vector("list", length(domain_ids))
  domain_biodiversity_pct_list <- vector("list", length(domain_ids))
  
  domain_richness_hot_top25_list <- vector("list", length(domain_ids))
  domain_richness_hot_top10_list <- vector("list", length(domain_ids))
  domain_richness_hot_top5_list  <- vector("list", length(domain_ids))
  domain_richness_hot_top1_list  <- vector("list", length(domain_ids))
  domain_richness_hot_class_list  <- vector("list", length(domain_ids))
  
  domain_abundance_hot_top25_list <- vector("list", length(domain_ids))
  domain_abundance_hot_top10_list <- vector("list", length(domain_ids))
  domain_abundance_hot_top5_list  <- vector("list", length(domain_ids))
  domain_abundance_hot_top1_list  <- vector("list", length(domain_ids))
  domain_abundance_hot_class_list  <- vector("list", length(domain_ids))
  
  domain_biodiversity_hot_top25_list <- vector("list", length(domain_ids))
  domain_biodiversity_hot_top10_list <- vector("list", length(domain_ids))
  domain_biodiversity_hot_top5_list  <- vector("list", length(domain_ids))
  domain_biodiversity_hot_top1_list  <- vector("list", length(domain_ids))
  domain_biodiversity_hot_class_list  <- vector("list", length(domain_ids))
  
  for (i in seq_along(domain_ids)) {
    dom_id <- domain_ids[i]
    message("  Domain ", dom_id)
    
    dom_poly <- domain_dat[
      domain_dat[[domain_name_field]] == dom_id
    ]
    
    if (nrow(dom_poly) == 0) {
      warning("No geometry found for domain: ", dom_id)
      # 
      # domain_richness_pct_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_abundance_pct_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_biodiversity_pct_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # 
      # domain_richness_hot_top25_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_richness_hot_top10_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_richness_hot_top5_list[[i]]  <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_richness_hot_top1_list[[i]]  <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # 
      # domain_abundance_hot_top25_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_abundance_hot_top10_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_abundance_hot_top5_list[[i]]  <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_abundance_hot_top1_list[[i]]  <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # 
      # domain_biodiversity_hot_top25_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_biodiversity_hot_top10_list[[i]] <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_biodiversity_hot_top5_list[[i]]  <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # domain_biodiversity_hot_top1_list[[i]]  <- empty_named_raster(domain_template_raster, paste0("domain_", dom_id))
      # 
      next
    }
    
    dom_richness <- mask(richness_median, dom_poly)
    dom_abundance <- mask(abundance_median_100, dom_poly)
    
    dom_richness_pct <- calc_percentile_rank_raster(dom_richness)
    dom_abundance_pct <- calc_percentile_rank_raster(dom_abundance)
    
    dom_biodiversity_raw <- dom_richness_pct + dom_abundance_pct
    dom_biodiversity_pct <- calc_percentile_rank_raster(dom_biodiversity_raw)
    
    dom_richness_hot <- calc_hotspot_products(dom_richness)
    dom_abundance_hot <- calc_hotspot_products(dom_abundance)
    dom_biodiversity_hot <- calc_hotspot_products(dom_biodiversity_raw)
    
    names(dom_richness_pct) <- paste0("domain_", dom_id)
    names(dom_abundance_pct) <- paste0("domain_", dom_id)
    names(dom_biodiversity_pct) <- paste0("domain_", dom_id)
    domain_richness_pct_list[[i]] <- dom_richness_pct
    domain_abundance_pct_list[[i]] <- dom_abundance_pct
    domain_biodiversity_pct_list[[i]] <- dom_biodiversity_pct
    
    domain_richness_hot_top25_list[[i]] <- dom_richness_hot[["top25"]]
    domain_richness_hot_top10_list[[i]] <- dom_richness_hot[["top10"]]
    domain_richness_hot_top5_list[[i]]  <- dom_richness_hot[["top5"]]
    domain_richness_hot_top1_list[[i]]  <- dom_richness_hot[["top1"]]
    domain_richness_hot_class_list[[i]] <- dom_richness_hot[["intensity"]]
    
    domain_abundance_hot_top25_list[[i]] <- dom_abundance_hot[["top25"]]
    domain_abundance_hot_top10_list[[i]] <- dom_abundance_hot[["top10"]]
    domain_abundance_hot_top5_list[[i]]  <- dom_abundance_hot[["top5"]]
    domain_abundance_hot_top1_list[[i]]  <- dom_abundance_hot[["top1"]]
    domain_abundance_hot_class_list[[i]] <- dom_abundance_hot[["intensity"]]
    
    domain_biodiversity_hot_top25_list[[i]] <- dom_biodiversity_hot[["top25"]]
    domain_biodiversity_hot_top10_list[[i]] <- dom_biodiversity_hot[["top10"]]
    domain_biodiversity_hot_top5_list[[i]]  <- dom_biodiversity_hot[["top5"]]
    domain_biodiversity_hot_top1_list[[i]]  <- dom_biodiversity_hot[["top1"]]
    domain_biodiversity_hot_class_list[[i]] <- dom_biodiversity_hot[["intensity"]]
    
    rm(dom_poly,
      dom_richness, dom_abundance,
      dom_richness_pct, dom_abundance_pct, dom_biodiversity_pct,
      dom_richness_hot, dom_abundance_hot, dom_biodiversity_hot
    )
    gc()
  }
  
  domain_richness_pct <- rast(domain_richness_pct_list)
  domain_abundance_pct <- rast(domain_abundance_pct_list)
  domain_biodiversity_pct <- rast(domain_biodiversity_pct_list)
  
  domain_richness_hot_top25 <- rast(domain_richness_hot_top25_list)
  domain_richness_hot_top10 <- rast(domain_richness_hot_top10_list)
  domain_richness_hot_top5  <- rast(domain_richness_hot_top5_list)
  domain_richness_hot_top1  <- rast(domain_richness_hot_top1_list)
  domain_richness_hot_class <- rast(domain_richness_hot_class_list)
  
  domain_abundance_hot_top25 <- rast(domain_abundance_hot_top25_list)
  domain_abundance_hot_top10 <- rast(domain_abundance_hot_top10_list)
  domain_abundance_hot_top5  <- rast(domain_abundance_hot_top5_list)
  domain_abundance_hot_top1  <- rast(domain_abundance_hot_top1_list)
  domain_abundance_hot_class <- rast(domain_abundance_hot_class_list)
  
  domain_biodiversity_hot_top25 <- rast(domain_biodiversity_hot_top25_list)
  domain_biodiversity_hot_top10 <- rast(domain_biodiversity_hot_top10_list)
  domain_biodiversity_hot_top5  <- rast(domain_biodiversity_hot_top5_list)
  domain_biodiversity_hot_top1  <- rast(domain_biodiversity_hot_top1_list)
  domain_biodiversity_hot_class <- rast(domain_biodiversity_hot_class_list)
  
  names(domain_richness_pct) <- paste0("domain_", domain_ids)
  names(domain_abundance_pct) <- paste0("domain_", domain_ids)
  names(domain_biodiversity_pct) <- paste0("domain_", domain_ids)
  
  names(domain_richness_hot_top25) <- paste0("domain_", domain_ids)
  names(domain_richness_hot_top10) <- paste0("domain_", domain_ids)
  names(domain_richness_hot_top5)  <- paste0("domain_", domain_ids)
  names(domain_richness_hot_top1)  <- paste0("domain_", domain_ids)
  names(domain_richness_hot_class) <- paste0("domain_", domain_ids)
  
  names(domain_abundance_hot_top25) <- paste0("domain_", domain_ids)
  names(domain_abundance_hot_top10) <- paste0("domain_", domain_ids)
  names(domain_abundance_hot_top5)  <- paste0("domain_", domain_ids)
  names(domain_abundance_hot_top1)  <- paste0("domain_", domain_ids)
  names(domain_abundance_hot_class) <- paste0("domain_", domain_ids)
  
  names(domain_biodiversity_hot_top25) <- paste0("domain_", domain_ids)
  names(domain_biodiversity_hot_top10) <- paste0("domain_", domain_ids)
  names(domain_biodiversity_hot_top5)  <- paste0("domain_", domain_ids)
  names(domain_biodiversity_hot_top1)  <- paste0("domain_", domain_ids)
  names(domain_biodiversity_hot_class) <- paste0("domain_", domain_ids)
  
  ############################################################
  # 5.11 POLYGONS FOR TOP 10% THRESHOLD
  ############################################################
  message("Creating polygons for top 10% threshold")
  
  richness_top10_poly <- hotspot_polygons_from_class(richness_hot[["top10"]], class_field = "top10")
  abundance_top10_poly <- hotspot_polygons_from_class(abundance_hot[["top10"]], class_field = "top10")
  biodiversity_top10_poly <- hotspot_polygons_from_class(biodiversity_hot[["top10"]], class_field = "top10")
  
  ############################################################
  # 5.12 WRITE CONTINUOUS OUTPUTS
  ############################################################
  message("Writing continuous outputs")
  writeRaster(richness_mean,
    filename = file.path(out_dir_cont, "richness_mean.tif"), overwrite = TRUE)
  
  writeRaster(richness_median,
    filename = file.path(out_dir_cont, "richness_median.tif"), overwrite = TRUE)
  
  writeRaster(richness_boot_median,
    filename = file.path(out_dir_cont, "richness_bootstrap_median.tif"), overwrite = TRUE)
  
  writeRaster(richness_boot_sd,
    filename = file.path(out_dir_cont, "richness_bootstrap_sd.tif"), overwrite = TRUE)
  
  
  writeRaster(abundance_mean,
    filename = file.path(out_dir_cont, "abundance_mean.tif"), overwrite = TRUE)
  
  writeRaster(abundance_median,
    filename = file.path(out_dir_cont, "abundance_median.tif"), overwrite = TRUE)
  
  writeRaster(abundance_median_100,
    filename = file.path(out_dir_cont, "abundance_median_100.tif"), overwrite = TRUE)
  
  writeRaster(abundance_se,
    filename = file.path(out_dir_cont, "abundance_se.tif"), overwrite = TRUE)
  
  
  writeRaster(richness_per_abundance,
    filename = file.path(out_dir_cont, "richness_per_abundance.tif"), overwrite = TRUE)
  
  writeRaster(abundance_per_richness,
    filename = file.path(out_dir_cont, "abundance_per_richness.tif"), overwrite = TRUE)
  
  # if (!is.null(beta_median)) {
  #   writeRaster(beta_median,
  #     filename = file.path(out_dir_cont, "beta_diversity_median.tif"), overwrite = TRUE)
  # }
  
  ############################################################
  # 5.13 WRITE PERCENTILE OUTPUTS
  ############################################################
  message("Writing percentile outputs")
  
  writeRaster(richness_pct,
    filename = file.path(out_dir_pct, "richness_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(abundance_pct,
    filename = file.path(out_dir_pct, "abundance_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(biodiversity_pct,
    filename = file.path(out_dir_pct, "biodiversity_percentiles.tif"), overwrite = TRUE)
  
  
  writeRaster(richness_only_pct,
    filename = file.path(out_dir_pct, "richness_only_percentile_difference.tif"), overwrite = TRUE)
  
  writeRaster(abundance_only_pct,
    filename = file.path(out_dir_pct, "abundance_only_percentile_difference.tif"), overwrite = TRUE)
  
  
  writeRaster(richness_per_abundance_pct,
    filename = file.path(out_dir_pct, "richness_per_abundance_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(abundance_per_richness_pct,
    filename = file.path(out_dir_pct, "abundance_per_richness_percentiles.tif"), overwrite = TRUE)
  
  # if (!is.null(beta_median)) {
  #   writeRaster(beta_pct,
  #     filename = file.path(out_dir_pct, "beta_diversity_percentiles.tif"), overwrite = TRUE)
  # }
  
  ############################################################
  # 5.14 WRITE HOTSPOT OUTPUTS
  ############################################################
  message("Writing hotspot outputs")
  
  write_hotspot_stack(richness_hot,
    filename = file.path(out_dir_hot, "richness_hotspots.tif"), overwrite = TRUE)
  
  write_hotspot_stack(abundance_hot,
    filename = file.path(out_dir_hot, "abundance_hotspots.tif"), overwrite = TRUE)
  
  write_hotspot_stack(biodiversity_hot,
    filename = file.path(out_dir_hot, "biodiversity_hotspots.tif"), overwrite = TRUE)
  
  writeRaster(richness_only_class,
    filename = file.path(out_dir_hot, "richness_only_hotspots.tif"), overwrite = TRUE)
  
  writeRaster(abundance_only_class,
    filename = file.path(out_dir_hot, "abundance_only_hotspots.tif"), overwrite = TRUE)
  
  # if (!is.null(beta_median)) {
  #   writeRaster(beta_class,
  #     filename = file.path(out_dir_hot, "beta_diversity_hotspot_classes.tif"), overwrite = TRUE)
  # }
  
  ############################################################
  # 5.16 WRITE DOMAIN OUTPUTS
  ############################################################
  message("Writing domain outputs")
  
  writeRaster(domain_richness_pct,
    filename = file.path(out_dir_domain, "domain_richness_percentiles.tif"), overwrite = TRUE)

  writeRaster(domain_abundance_pct,
    filename = file.path(out_dir_domain, "domain_abundance_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(domain_biodiversity_pct,
    filename = file.path(out_dir_domain, "domain_biodiversity_percentiles.tif"), overwrite = TRUE)
  
  
  writeRaster(domain_richness_hot_top25,
    filename = file.path(out_dir_domain, "domain_richness_hotspot_top25.tif"), overwrite = TRUE)
  
  writeRaster(domain_richness_hot_top10,
    filename = file.path(out_dir_domain, "domain_richness_hotspot_top10.tif"), overwrite = TRUE)
  
  writeRaster(domain_richness_hot_top5,
    filename = file.path(out_dir_domain, "domain_richness_hotspot_top5.tif"), overwrite = TRUE)
  
  writeRaster(domain_richness_hot_top1,
    filename = file.path(out_dir_domain, "domain_richness_hotspot_top1.tif"), overwrite = TRUE)

  writeRaster(domain_richness_hot_class,
    filename = file.path(out_dir_domain, "domain_richness_hotspot_class.tif"), overwrite = TRUE)
  

  writeRaster(domain_abundance_hot_top25,
    filename = file.path(out_dir_domain, "domain_abundance_hotspot_top25.tif"), overwrite = TRUE)
  
  writeRaster(domain_abundance_hot_top10,
    filename = file.path(out_dir_domain, "domain_abundance_hotspot_top10.tif"), overwrite = TRUE)
  
  writeRaster(domain_abundance_hot_top5,
    filename = file.path(out_dir_domain, "domain_abundance_hotspot_top5.tif"), overwrite = TRUE)
  
  writeRaster(domain_abundance_hot_top1,
    filename = file.path(out_dir_domain, "domain_abundance_hotspot_top1.tif"), overwrite = TRUE)

  writeRaster(domain_abundance_hot_class,
    filename = file.path(out_dir_domain, "domain_abundance_hotspot_class.tif"), overwrite = TRUE)
  

  writeRaster(domain_biodiversity_hot_top25,
    filename = file.path(out_dir_domain, "domain_biodiversity_hotspot_top25.tif"), overwrite = TRUE)
  
  writeRaster(domain_biodiversity_hot_top10,
    filename = file.path(out_dir_domain, "domain_biodiversity_hotspot_top10.tif"), overwrite = TRUE)
  
  writeRaster(domain_biodiversity_hot_top5,
    filename = file.path(out_dir_domain, "domain_biodiversity_hotspot_top5.tif"), overwrite = TRUE)
  
  writeRaster(domain_biodiversity_hot_top1,
    filename = file.path(out_dir_domain, "domain_biodiversity_hotspot_top1.tif"), overwrite = TRUE)

  writeRaster(domain_biodiversity_hot_class,
    filename = file.path(out_dir_domain, "domain_biodiversity_hotspot_class.tif"), overwrite = TRUE)
  

  ############################################################
  # 5.17 WRITE POLYGON OUTPUTS
  ############################################################
  message("Writing polygon outputs")
  
    writeVector(richness_top10_poly,
      filename = file.path(out_dir_poly, "richness_top10_hotspots.gpkg"), overwrite = TRUE)

    writeVector(abundance_top10_poly,
      filename = file.path(out_dir_poly, "abundance_top10_hotspots.gpkg"), overwrite = TRUE)

    writeVector(biodiversity_top10_poly,
      filename = file.path(out_dir_poly, "biodiversity_top10_hotspots.gpkg"), overwrite = TRUE)

  ############################################################
  # 5.19 CLEAN MEMORY
  ############################################################
  
  rm(
    pa_mean, pa_median, pa_se, pa_5, pa_95,
    ab_stack, abundance_mean, abundance_median, abundance_median_100, abundance_se, abundance_5, abundance_95,
    richness_mean, richness_median, richness_boot, richness_boot_median, richness_boot_sd,
    richness_per_abundance, abundance_per_richness,
    richness_pct, abundance_pct, biodiversity_pct, biodiversity_raw,
    richness_only_pct, abundance_only_pct,
    richness_per_abundance_pct, abundance_per_richness_pct,
    richness_hot, abundance_hot, biodiversity_hot,
    richness_only_class, abundance_only_class,
    domain_richness_pct, domain_abundance_pct, domain_biodiversity_pct,
    domain_richness_hot_top25, domain_richness_hot_top10, domain_richness_hot_top5, domain_richness_hot_top1, domain_richness_class,
    domain_abundance_hot_top25, domain_abundance_hot_top10, domain_abundance_hot_top5, domain_abundance_hot_top1, domain_abundance_class,
    domain_biodiversity_hot_top25, domain_biodiversity_hot_top10, domain_biodiversity_hot_top5, domain_biodiversity_hot_top1, domain_biodiversity_class,
    domain_richness_pct_list, domain_abundance_pct_list, domain_biodiversity_pct_list,
    domain_richness_hot_top25_list, domain_richness_hot_top10_list, domain_richness_hot_top5_list, domain_richness_hot_top1_list, domain_richness_class_list,
    domain_abundance_hot_top25_list, domain_abundance_hot_top10_list, domain_abundance_hot_top5_list, domain_abundance_hot_top1_list, domain_abundance_class_list,
    domain_biodiversity_hot_top25_list, domain_biodiversity_hot_top10_list, domain_biodiversity_hot_top5_list, domain_biodiversity_hot_top1_list, domain_biodiversity_class_list,
    richness_top10_poly, abundance_top10_poly, biodiversity_top10_poly
  )
  
  # if (!is.null(beta_median)) {rm(beta_median, beta_pct, beta_hot, beta_class)}
  
  gc()
  
  message(sprintf(
    "Finished %s in %.2f minutes",
    nm,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}


































































































































































############################################################
# HMSC BIODIVERSITY HOTSPOT ANALYSIS PIPELINE
# ----------------------------------------------------------
# This script:
#
# 1) Loads predictive species distribution maps for multiple models
# 2) Derives biodiversity metrics and uncertainty layers
# 3) Calculates circumpolar percentile surfaces
# 4) Calculates planning-domain-specific percentile surfaces
# 5) Creates simplified percentile class rasters
# 6) Creates multi-threshold hotspot intensity rasters
# 7) Creates richness-only and abundance-only hotspot layers
# 8) Creates uncertainty-based hotspot layers
# 9) Converts hotspot rasters to polygons
# 10) Calculates multi-model consensus products
#
# Designed for:
# - Circum-Antarctic seafloor biodiversity predictions
# - Outputs assembled from the HMSC modelling pipeline
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
source("0_Functions.R")

# Spatial resolution label
res <- "2km"

# Model identifiers
model_ids <- c(
  "npp_cafe", "npp_cbpm", "npp_eppl", "npp_vpmg",
  "fam_cafe", "fam_cbpm", "fam_eppl", "fam_vpmg"
)

# Source dropbox folder with original files
src_base <- file.path(
  usr.dropbox.dir,
  "data_products/predictive_maps/circum_antarctic"
)

# Local directory with copied files for faster processing
dst_base <- file.path(
  usr.main.dir,
  "4_model_prediction/copy_of_predictive_maps/circum_antarctic"
)

# Planning-domain definitions
domain_ids <- c(1, 3, 4, 7, 8, 9)
domain_name_field <- "Name"

domain_shp <- file.path(
  usr.dropbox.dir,
  "data_environmental/raw/ccamlr_data_2025/geographical_data/mpapd/CCAMLR_MPAPD_EPSG6932.shp"
)

# # Beta-diversity directory
# # The script searches for files containing the model identifier.
# beta_dir <- file.path(
#   usr.dropbox.dir,
#   "data_products/predictive_maps/circum_antarctic_beta_diversity"
# )
#
# # Beta-diversity layer index to use if the file has multiple layers
# beta_layer_index <- 1


############################
# 2) HELPER FUNCTIONS
############################

# ----------------------------------------------------------
# Function: positive_difference
# ----------------------------------------------------------
# Returns the positive component of a raster difference
# Negative values are set to zero
# ----------------------------------------------------------
positive_difference <- function(ra1, ra2) {
  out <- ra1 - ra2
  out[out < 0] <- 0
  out[is.na(ra1)] <- NA
  out[is.na(ra2)] <- NA
  out
}

# ----------------------------------------------------------
# Function: safe_divide_raster
# ----------------------------------------------------------
# Divides one raster by another
# Zero or negative denominator values are set to NA
# ----------------------------------------------------------
safe_divide_raster <- function(num, den) {
  den2 <- den
  den2[den2 <= 0] <- NA
  out <- num / den2
  out[is.infinite(out)] <- NA
  out
}

# ----------------------------------------------------------
# Function: write_hotspot_stack
# ----------------------------------------------------------
# Writes a standard hotspot stack with named layers
# ----------------------------------------------------------
write_hotspot_stack <- function(hot_stack, filename, overwrite = TRUE) {
  names(hot_stack) <- c("top25", "top10", "top5", "top1", "intensity")
  writeRaster(hot_stack, filename = filename, overwrite = overwrite)
}

# # ----------------------------------------------------------
# # Function: find_beta_file
# # ----------------------------------------------------------
# # Finds an optional beta-diversity file for a model
# # ----------------------------------------------------------
# find_beta_file <- function(model_id, beta_dir) {
#   if (!dir.exists(beta_dir)) {
#     return(NA_character_)
#   }
#
#   beta_files <- list.files(beta_dir, pattern = "\\.tif$", full.names = TRUE)
#   if (length(beta_files) == 0) {
#     return(NA_character_)
#   }
#
#   hit <- grep(model_id, basename(beta_files), ignore.case = TRUE)
#
#   if (length(hit) == 1) {
#     return(beta_files[hit])
#   }
#
#   if (length(beta_files) == 1) {
#     return(beta_files[1])
#   }
#
#   NA_character_
# }


############################
# 3) CHECK LOCAL COPY OF PREDICTIVE MAPS
############################

# # COPY PREDICTIVE MAPS TO LOCAL DISK
# # ----------------------------------------------------------
# # Copy all predictive map files from the synced storage
# # location to a local directory, preserving the folder
# # structure for each model.
# ############################################################
# dir.create(dst_base, recursive = TRUE, showWarnings = FALSE)
#
# all_files <- list.files(
#   src_base,
#   recursive = TRUE,
#   full.names = TRUE,
#   include.dirs = FALSE
# )
#
# all_files <- all_files[grepl("\\.tif$", all_files, ignore.case = TRUE)]
#
# rel_paths <- substring(all_files, nchar(src_base) + 2)
# dst_files <- file.path(dst_base, rel_paths)
#
# dst_dirs <- unique(dirname(dst_files))
# invisible(lapply(dst_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
#
# ok <- file.copy(from = all_files, to = dst_files, overwrite = FALSE, copy.mode = TRUE)
# message(sprintf("Copied %d of %d files", sum(ok), length(ok)))
#
# if (any(!ok)) {
#   message("The following files were not copied:")
#   print(all_files[!ok])
# }
#
# src_files <- list.files(src_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
# dst_files <- list.files(dst_base, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
#
# src_rel <- substring(src_files, nchar(src_base) + 2)
# dst_rel <- substring(dst_files, nchar(dst_base) + 2)
#
# m <- match(src_rel, dst_rel)
# keep <- !is.na(m)
#
# src_mtime <- file.info(src_files[keep])$mtime
# for (i in seq_along(src_mtime)) {
#   Sys.setFileTime(dst_files[m[keep][i]], src_mtime[i])
# }
#
# message(sprintf("Updated timestamps for %d files", sum(keep)))

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

src_files <- list.files(
  src_base,
  pattern = "\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

dst_files <- list.files(
  dst_base,
  pattern = "\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

src_rel <- substring(src_files, nchar(src_base) + 2)
dst_rel <- substring(dst_files, nchar(dst_base) + 2)

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

cmp <- merge(src_info, dst_info, by = "rel_path", all = TRUE)

cmp$missing_in_local  <- is.na(cmp$dst_file)
cmp$missing_in_source <- is.na(cmp$src_file)
cmp$size_diff <- FALSE
cmp$mtime_diff <- FALSE

both_present <- !cmp$missing_in_local & !cmp$missing_in_source
cmp$size_diff[both_present]  <- cmp$src_size[both_present]  != cmp$dst_size[both_present]
cmp$mtime_diff[both_present] <- cmp$src_mtime[both_present] != cmp$dst_mtime[both_present]

cmp_issues <- cmp[
  cmp$missing_in_local |
    cmp$missing_in_source |
    cmp$size_diff |
    cmp$mtime_diff,
]

if (nrow(cmp_issues) == 0) {
  message("Local predictive maps match source files (by path, size, and modification time).")
  pred_base_dir <- dst_base
} else {
  message("Differences detected between source and local predictive maps.")
  message("Using source predictive maps.")
  pred_base_dir <- src_base
}


############################
# 4) LOAD PLANNING DOMAINS
############################

if (!file.exists(domain_shp)) {
  stop("Planning-domain shapefile not found: ", domain_shp)
}

domain_raw <- vect(domain_shp)

# Assign source CRS if missing
domain_crs <- crs(domain_raw)
if (is.na(domain_crs) || domain_crs == "") {
  crs(domain_raw) <- "EPSG:6932"
}

# Domain raster alignment is based on the first successfully loaded model
domain_template_raster <- NULL

# Objects for multi-model consensus
consensus_files <- list()


############################
# 5) PROCESS EACH MODEL
############################

for (nm in model_ids) {
  t0 <- Sys.time()
  
  message("======================================")
  message("Processing model: ", nm)
  message("======================================")
  
  model_dir <- file.path(pred_base_dir, paste0("hmsc_with_", nm))
  if (!dir.exists(model_dir)) {
    warning("Model directory not found: ", model_dir)
    next
  }
  
  out_dir <- file.path(src_base, paste0("hmsc_with_", nm))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_dir_cont   <- file.path(out_dir, "continuous")
  out_dir_pct    <- file.path(out_dir, "percentiles")
  out_dir_class  <- file.path(out_dir, "classes")
  out_dir_hot    <- file.path(out_dir, "hotspots")
  out_dir_domain <- file.path(out_dir, "domains")
  out_dir_poly   <- file.path(out_dir, "polygons")
  
  dir.create(out_dir_cont,   recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_pct,    recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_class,  recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_hot,    recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_domain, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_poly,   recursive = TRUE, showWarnings = FALSE)
  
  ############################################################
  # 5.1 LOAD SPECIES-LEVEL PREDICTIONS
  ############################################################
  
  pa_files <- list.files(
    model_dir,
    pattern = "^PA_.*\\.tif$",
    full.names = TRUE
  )
  
  if (length(pa_files) == 0) {
    warning("No PA files found for model: ", nm)
    next
  }
  
  pa_mean   <- rast(pa_files, lyrs = 1)
  pa_median <- rast(pa_files, lyrs = 2)
  pa_se     <- rast(pa_files, lyrs = 3)
  pa_5      <- rast(pa_files, lyrs = 4)
  pa_95     <- rast(pa_files, lyrs = 5)
  
  # species_names <- sub("^PA_|\\.tif$", "", basename(pa_files))
  # names(pa_mean)   <- species_names
  # names(pa_median) <- species_names
  # names(pa_se)     <- species_names
  # names(pa_5)      <- species_names
  # names(pa_95)     <- species_names
  
  ab_file <- file.path(model_dir, "total_abundance.tif")
  if (!file.exists(ab_file)) {
    warning("No total abundance file found for model: ", nm)
    next
  }
  
  ab_stack <- rast(ab_file)
  # Convert abundance predictions from 540-point scale to 100-point scale
  abundance_mean   <- ab_stack[[1]] / 5.4
  abundance_median <- ab_stack[[2]] / 5.4
  abundance_se     <- ab_stack[[3]] / 5.4
  abundance_5      <- ab_stack[[4]] / 5.4
  abundance_95     <- ab_stack[[5]] / 5.4
  # Cap median abundance at 100 to retain interpretation as percent cover
  abundance_median_100 <- abundance_median
  abundance_median_100[abundance_median_100 > 100] <- 100
  names(abundance_median_100) <- "abundance_median_100"
  
  ############################################################
  # 5.2 DERIVE BIODIVERSITY METRICS
  ############################################################
  
  # Summed richness from species-level occurrence probabilities
  richness_mean   <- sum(pa_mean, na.rm = TRUE)
  richness_median <- sum(pa_median, na.rm = TRUE)
  
  # Bootstrap-based uncertainty for richness
  richness_boot <- bootstrap_richness_from_species_medians(
    pa_median = pa_median,
    n_boot = 1000,
    chunk_size = 10000,
    seed = 2
  )
  
  richness_boot_median <- richness_boot$median
  richness_boot_sd <- richness_boot$sd
  
  # Ratio of richness to abundance, and vice-versa
  richness_per_abundance <- safe_divide_raster(richness_median, abundance_median_100)
  abundance_per_richness <- safe_divide_raster(abundance_median_100, richness_median)
  
  ############################################################
  # 5.3 BETA-DIVERSITY INPUT
  ############################################################
  
  # beta_file <- find_beta_file(nm, beta_dir)
  # beta_median <- NULL
  #
  # if (!is.na(beta_file) && file.exists(beta_file)) {
  #   beta_stack <- rast(beta_file)
  #   if (nlyr(beta_stack) >= beta_layer_index) {
  #     beta_median <- beta_stack[[beta_layer_index]]
  #   }
  # }
  
  ############################################################
  # 5.4 PREPARE DOMAIN GEOMETRY
  ############################################################
  if (is.null(domain_template_raster)) {
    domain_template_raster <- richness_median
  }
  domain_dat <- project(domain_raw, crs(richness_median))
  
  ############################################################
  # 5.5 CIRCUMPOLAR PERCENTILE SURFACES
  ############################################################
  message("Calculating circumpolar percentile surfaces")
  
  richness_pct  <- calc_percentile_rank_raster(richness_median)
  abundance_pct <- calc_percentile_rank_raster(abundance_median_100)
  
  biodiversity_raw <- richness_pct + abundance_pct
  biodiversity_pct <- calc_percentile_rank_raster(biodiversity_raw)

  # Fine contrasts based on raw percentile differences
  richness_only_pct  <- positive_difference(richness_pct, abundance_pct)
  abundance_only_pct <- positive_difference(abundance_pct, richness_pct)
                                            
  richness_per_abundance_pct <- calc_percentile_rank_raster(richness_per_abundance)
  abundance_per_richness_pct <- calc_percentile_rank_raster(abundance_per_richness)
  
  # if (!is.null(beta_median)) {
  #   beta_pct <- calc_percentile_rank_raster(beta_median)
  # }
  
  ############################################################
  # 5.6 SIMPLIFIED PERCENTILE CLASSES
  ############################################################
  message("Creating simplified percentile classes")
  
  richness_class <- simplify_percentiles(richness_pct)
  abundance_class <- simplify_percentiles(abundance_pct)
  biodiversity_class <- simplify_percentiles(biodiversity_pct)
  richness_only_class  <- positive_difference(richness_class, abundance_class)
  abundance_only_class <- positive_difference(abundance_class, richness_class)
  
  # if (!is.null(beta_median)) {
  #   beta_class <- simplify_percentiles(beta_pct)
  # }

  ############################################################
  # 5.10 DOMAIN-SPECIFIC PERCENTILES AND HOTSPOTS
  ############################################################
  message("Creating domain-specific percentile and hotspot products")
  
  domain_richness_pct_list <- vector("list", length(domain_ids))
  domain_abundance_pct_list <- vector("list", length(domain_ids))
  domain_biodiversity_pct_list <- vector("list", length(domain_ids))
  
  domain_richness_class_list <- vector("list", length(domain_ids))
  domain_abundance_class_list <- vector("list", length(domain_ids))
  domain_biodiversity_class_list <- vector("list", length(domain_ids))
  
  domain_richness_hot_top10_list <- vector("list", length(domain_ids))
  domain_abundance_hot_top10_list <- vector("list", length(domain_ids))
  domain_biodiversity_hot_top10_list <- vector("list", length(domain_ids))
  
  domain_richness_poly_list <- list()
  domain_abundance_poly_list <- list()
  domain_biodiversity_poly_list <- list()
  
  for (i in seq_along(domain_ids)) {
    dom_id <- domain_ids[i]
    message("  Domain ", dom_id)
    
    dom_poly <- domain_dat[
      domain_dat[[domain_name_field]] == dom_id
    ]
    
    if (nrow(dom_poly) == 0) {
      warning("No geometry found for domain: ", dom_id)
      # dom_na <- mask(domain_template_raster, domain_template_raster)
      # names(dom_na) <- paste0("domain_", dom_id)
      # 
      # domain_richness_pct_list[[i]] <- dom_na
      # domain_abundance_pct_list[[i]] <- dom_na
      # domain_biodiversity_pct_list[[i]] <- dom_na
      # 
      # domain_richness_class_list[[i]] <- dom_na
      # domain_abundance_class_list[[i]] <- dom_na
      # domain_biodiversity_class_list[[i]] <- dom_na
      # 
      # domain_richness_hot_top10_list[[i]] <- dom_na
      # domain_abundance_hot_top10_list[[i]] <- dom_na
      # domain_biodiversity_hot_top10_list[[i]] <- dom_na
      # 
      # rm(dom_na)
      # next
    }
    
    dom_richness <- mask(richness_median, dom_poly)
    dom_abundance <- mask(abundance_median_100, dom_poly)
    
    dom_richness_pct <- calc_percentile_rank_raster(dom_richness)
    dom_abundance_pct <- calc_percentile_rank_raster(dom_abundance)
    
    dom_biodiversity_raw <- dom_richness_pct + dom_abundance_pct
    dom_biodiversity_pct <- calc_percentile_rank_raster(dom_biodiversity_raw)
    
    dom_richness_class <- simplify_percentiles(dom_richness_pct)
    dom_abundance_class <- simplify_percentiles(dom_abundance_pct)
    dom_biodiversity_class <- simplify_percentiles(dom_biodiversity_pct)
    
    names(dom_richness_pct) <- paste0("domain_", dom_id)
    names(dom_abundance_pct) <- paste0("domain_", dom_id)
    names(dom_biodiversity_pct) <- paste0("domain_", dom_id)
    
    names(dom_richness_class) <- paste0("domain_", dom_id)
    names(dom_abundance_class) <- paste0("domain_", dom_id)
    names(dom_biodiversity_class) <- paste0("domain_", dom_id)
    
    domain_richness_pct_list[[i]] <- dom_richness_pct
    domain_abundance_pct_list[[i]] <- dom_abundance_pct
    domain_biodiversity_pct_list[[i]] <- dom_biodiversity_pct
    
    domain_richness_class_list[[i]] <- dom_richness_class
    domain_abundance_class_list[[i]] <- dom_abundance_class
    domain_biodiversity_class_list[[i]] <- dom_biodiversity_class
    
    rm(
      dom_poly,
      dom_richness, dom_abundance,
      dom_richness_pct, dom_abundance_pct, dom_biodiversity_pct,
      dom_richness_class, dom_abundance_class, dom_biodiversity_class,
    )
    gc()
  }
  
  domain_richness_pct <- rast(domain_richness_pct_list)
  domain_abundance_pct <- rast(domain_abundance_pct_list)
  domain_biodiversity_pct <- rast(domain_biodiversity_pct_list)
  
  domain_richness_class <- rast(domain_richness_class_list)
  domain_abundance_class <- rast(domain_abundance_class_list)
  domain_biodiversity_class <- rast(domain_biodiversity_class_list)
  
  names(domain_richness_pct) <- paste0("domain_", domain_ids)
  names(domain_abundance_pct) <- paste0("domain_", domain_ids)
  names(domain_biodiversity_pct) <- paste0("domain_", domain_ids)
  
  names(domain_richness_class) <- paste0("domain_", domain_ids)
  names(domain_abundance_class) <- paste0("domain_", domain_ids)
  names(domain_biodiversity_class) <- paste0("domain_", domain_ids)
  
  ############################################################
  # 5.12 WRITE CONTINUOUS OUTPUTS
  ############################################################
  message("Writing continuous outputs")
  
  writeRaster(richness_mean, 
    filename = file.path(out_dir_cont, "richness_mean.tif"), overwrite = TRUE)
  
  writeRaster(richness_median,
    filename = file.path(out_dir_cont, "richness_median.tif"), overwrite = TRUE)
  
  writeRaster(richness_boot_median,
    filename = file.path(out_dir_cont, "richness_bootstrap_median.tif"), overwrite = TRUE)
  
  writeRaster(richness_boot_sd,
    filename = file.path(out_dir_cont, "richness_bootstrap_sd.tif"), overwrite = TRUE)
  
  
  writeRaster(abundance_mean,
    filename = file.path(out_dir_cont, "abundance_mean.tif"), overwrite = TRUE)
  
  writeRaster(abundance_median,
    filename = file.path(out_dir_cont, "abundance_median.tif"), overwrite = TRUE)
  
  writeRaster(abundance_median_100,
    filename = file.path(out_dir_cont, "abundance_median_100.tif"), overwrite = TRUE)
  
  writeRaster(abundance_se,
    filename = file.path(out_dir_cont, "abundance_se.tif"), overwrite = TRUE)
  
  
  writeRaster(richness_per_abundance,
    filename = file.path(out_dir_cont, "richness_per_abundance.tif"), overwrite = TRUE)
  
  writeRaster(abundance_per_richness,
    filename = file.path(out_dir_cont, "abundance_per_richness.tif"), overwrite = TRUE)

  # if (!is.null(beta_median)) {
  #   writeRaster( beta_median,
  #     filename = file.path(out_dir_cont, "beta_diversity_median.tif"), overwrite = TRUE)
  # }
  
  ############################################################
  # 5.13 WRITE PERCENTILE OUTPUTS
  ############################################################
  message("Writing percentile outputs")
  
  writeRaster(richness_pct,
    filename = file.path(out_dir_pct, "richness_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(abundance_pct,
    filename = file.path(out_dir_pct, "abundance_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(biodiversity_pct,
    filename = file.path(out_dir_pct, "biodiversity_percentiles.tif"), overwrite = TRUE)
  
  
  writeRaster(richness_only_pct,
    filename = file.path(out_dir_pct, "richness_only_percentile_difference.tif"), overwrite = TRUE)

  writeRaster(abundance_only_pct,
    filename = file.path(out_dir_pct, "abundance_only_percentile_difference.tif"), overwrite = TRUE)
  
  
  writeRaster(richness_per_abundance_pct,
    filename = file.path(out_dir_pct, "richness_per_abundance_percentiles.tif"), overwrite = TRUE)

  writeRaster(abundance_per_richness_pct,
    filename = file.path(out_dir_pct, "abundance_per_richness_percentiles.tif"), overwrite = TRUE)

  # if (!is.null(beta_median)) {
  #   writeRaster(beta_pct,
  #     filename = file.path(out_dir_pct, "beta_diversity_percentiles.tif"), overwrite = TRUE)
  # }
  
  ############################################################
  # 5.14 WRITE CLASS OUTPUTS
  ############################################################
  message("Writing class outputs")
  
  writeRaster(richness_class,
    filename = file.path(out_dir_class, "richness_percentile_classes.tif"), overwrite = TRUE)
  
  writeRaster(abundance_class,
    filename = file.path(out_dir_class, "abundance_percentile_classes.tif"), overwrite = TRUE)
  
  writeRaster(biodiversity_class,
    filename = file.path(out_dir_class, "biodiversity_percentile_classes.tif"), overwrite = TRUE)
  
  writeRaster(richness_only_class,
    filename = file.path(out_dir_class, "richness_only_percentile_classes.tif"), overwrite = TRUE)
  
  writeRaster(abundance_only_class,
    filename = file.path(out_dir_class, "abundance_only_percentile_classes.tif"), overwrite = TRUE)
  
  # if (!is.null(beta_median)) {
  #   writeRaster(beta_class,
  #     filename = file.path(out_dir_class, "beta_diversity_percentile_classes.tif"), overwrite = TRUE)
  # }
  
  ############################################################
  # 5.16 WRITE DOMAIN OUTPUTS
  ############################################################
  
  message("Writing domain outputs")
  
  writeRaster(domain_richness_pct,
    filename = file.path(out_dir_domain, "domain_richness_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(domain_abundance_pct,
    filename = file.path(out_dir_domain, "domain_abundance_percentiles.tif"), overwrite = TRUE)
  
  writeRaster(domain_biodiversity_pct,
    filename = file.path(out_dir_domain, "domain_biodiversity_percentiles.tif"), overwrite = TRUE)
  
  
  writeRaster(domain_richness_class,
    filename = file.path(out_dir_domain, "domain_richness_percentile_classes.tif"), overwrite = TRUE)
  
  writeRaster(domain_abundance_class,
    filename = file.path(out_dir_domain, "domain_abundance_percentile_classes.tif"), overwrite = TRUE)
  
  writeRaster(domain_biodiversity_class,
    filename = file.path(out_dir_domain, "domain_biodiversity_percentile_classes.tif"), overwrite = TRUE)

  ############################################################
  # 5.18 STORE FILE PATHS FOR CONSENSUS
  ############################################################
  
  consensus_files[[nm]] <- list(
    circ = list(
      richness_top10 = richness_hot_file,
      abundance_top10 = abundance_hot_file,
      biodiversity_top10 = biodiversity_hot_file,
      richness_intensity = richness_hot_file,
      abundance_intensity = abundance_hot_file,
      biodiversity_intensity = biodiversity_hot_file
    ),
    domains = list(
      richness_top10 = file.path(out_dir_domain, "domain_richness_hotspot_top10.tif"),
      abundance_top10 = file.path(out_dir_domain, "domain_abundance_hotspot_top10.tif"),
      biodiversity_top10 = file.path(out_dir_domain, "domain_biodiversity_hotspot_top10.tif"),
      richness_intensity = file.path(out_dir_domain, "domain_richness_hotspot_intensity.tif"),
      abundance_intensity = file.path(out_dir_domain, "domain_abundance_hotspot_intensity.tif"),
      biodiversity_intensity = file.path(out_dir_domain, "domain_biodiversity_hotspot_intensity.tif")
    )
  )
  
  ############################################################
  # 5.19 CLEAN MEMORY
  ############################################################
  
  rm(
    pa_mean, pa_median, pa_se, pa_5, pa_95,
    ab_stack, abundance_mean, abundance_median, abundance_median_100, abundance_se, abundance_5, abundance_95,
    richness_mean, richness_median, richness_boot, richness_boot_median, richness_boot_sd,
    richness_pct, abundance_pct, biodiversity_pct, biodiversity_raw,
    richness_class, abundance_class, biodiversity_class,
    richness_per_abundance, richness_per_abundance_pct,
    richness_only_class, abundance_only_class,
    domain_richness_pct, domain_abundance_pct, domain_biodiversity_pct,
    domain_richness_class, domain_abundance_class, domain_biodiversity_class,
  )
  
  # if (!is.null(beta_median)) {rm(beta_median, beta_pct, beta_class, beta_hot)}
  
  gc()
  
  message(sprintf(
    "Finished %s in %.2f minutes",
    nm,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}


############################
# 6) MULTI-MODEL CONSENSUS
############################

message("======================================")
message("Creating multi-model consensus products")
message("======================================")

processed_models <- names(consensus_files)

if (length(processed_models) == 0) {
  stop("No models were processed successfully. Consensus products were not created.")
}

cons_dir <- file.path(src_base, "consensus")
dir.create(cons_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------------------------------------
# 6.1 CIRCUMPOLAR CONSENSUS
#-----------------------------------------------------------

message("Creating circumpolar consensus rasters")

rich_top10_stack <- rast(lapply(
  processed_models,
  function(x) rast(consensus_files[[x]]$circ$richness_top10)[[2]]
))

abund_top10_stack <- rast(lapply(
  processed_models,
  function(x) rast(consensus_files[[x]]$circ$abundance_top10)[[2]]
))

biodiv_top10_stack <- rast(lapply(
  processed_models,
  function(x) rast(consensus_files[[x]]$circ$biodiversity_top10)[[2]]
))

rich_top10_consensus <- sum(rich_top10_stack, na.rm = TRUE)
abund_top10_consensus <- sum(abund_top10_stack, na.rm = TRUE)
biodiv_top10_consensus <- sum(biodiv_top10_stack, na.rm = TRUE)

rich_intensity_stack <- rast(lapply(
  processed_models,
  function(x) rast(consensus_files[[x]]$circ$richness_intensity)[[5]]
))

abund_intensity_stack <- rast(lapply(
  processed_models,
  function(x) rast(consensus_files[[x]]$circ$abundance_intensity)[[5]]
))

biodiv_intensity_stack <- rast(lapply(
  processed_models,
  function(x) rast(consensus_files[[x]]$circ$biodiversity_intensity)[[5]]
))

rich_intensity_consensus <- sum(rich_intensity_stack, na.rm = TRUE)
abund_intensity_consensus <- sum(abund_intensity_stack, na.rm = TRUE)
biodiv_intensity_consensus <- sum(biodiv_intensity_stack, na.rm = TRUE)

writeRaster(
  rich_top10_consensus,
  filename = file.path(cons_dir, "circumpolar_richness_hotspot_top10_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  abund_top10_consensus,
  filename = file.path(cons_dir, "circumpolar_abundance_hotspot_top10_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  biodiv_top10_consensus,
  filename = file.path(cons_dir, "circumpolar_biodiversity_hotspot_top10_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  rich_intensity_consensus,
  filename = file.path(cons_dir, "circumpolar_richness_hotspot_intensity_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  abund_intensity_consensus,
  filename = file.path(cons_dir, "circumpolar_abundance_hotspot_intensity_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  biodiv_intensity_consensus,
  filename = file.path(cons_dir, "circumpolar_biodiversity_hotspot_intensity_consensus.tif"),
  overwrite = TRUE
)


#-----------------------------------------------------------
# 6.2 DOMAIN-SPECIFIC CONSENSUS
#-----------------------------------------------------------

message("Creating domain-specific consensus rasters")

domain_rich_top10_consensus_list <- vector("list", length(domain_ids))
domain_abund_top10_consensus_list <- vector("list", length(domain_ids))
domain_biodiv_top10_consensus_list <- vector("list", length(domain_ids))

domain_rich_intensity_consensus_list <- vector("list", length(domain_ids))
domain_abund_intensity_consensus_list <- vector("list", length(domain_ids))
domain_biodiv_intensity_consensus_list <- vector("list", length(domain_ids))

for (i in seq_along(domain_ids)) {
  dom_id <- domain_ids[i]
  message("  Domain ", dom_id)
  
  rich_top10_dom_stack <- rast(lapply(
    processed_models,
    function(x) rast(consensus_files[[x]]$domains$richness_top10)[[i]]
  ))
  
  abund_top10_dom_stack <- rast(lapply(
    processed_models,
    function(x) rast(consensus_files[[x]]$domains$abundance_top10)[[i]]
  ))
  
  biodiv_top10_dom_stack <- rast(lapply(
    processed_models,
    function(x) rast(consensus_files[[x]]$domains$biodiversity_top10)[[i]]
  ))
  
  rich_int_dom_stack <- rast(lapply(
    processed_models,
    function(x) rast(consensus_files[[x]]$domains$richness_intensity)[[i]]
  ))
  
  abund_int_dom_stack <- rast(lapply(
    processed_models,
    function(x) rast(consensus_files[[x]]$domains$abundance_intensity)[[i]]
  ))
  
  biodiv_int_dom_stack <- rast(lapply(
    processed_models,
    function(x) rast(consensus_files[[x]]$domains$biodiversity_intensity)[[i]]
  ))
  
  domain_rich_top10_consensus_list[[i]] <- sum(rich_top10_dom_stack, na.rm = TRUE)
  domain_abund_top10_consensus_list[[i]] <- sum(abund_top10_dom_stack, na.rm = TRUE)
  domain_biodiv_top10_consensus_list[[i]] <- sum(biodiv_top10_dom_stack, na.rm = TRUE)
  
  domain_rich_intensity_consensus_list[[i]] <- sum(rich_int_dom_stack, na.rm = TRUE)
  domain_abund_intensity_consensus_list[[i]] <- sum(abund_int_dom_stack, na.rm = TRUE)
  domain_biodiv_intensity_consensus_list[[i]] <- sum(biodiv_int_dom_stack, na.rm = TRUE)
  
  names(domain_rich_top10_consensus_list[[i]]) <- paste0("domain_", dom_id)
  names(domain_abund_top10_consensus_list[[i]]) <- paste0("domain_", dom_id)
  names(domain_biodiv_top10_consensus_list[[i]]) <- paste0("domain_", dom_id)
  
  names(domain_rich_intensity_consensus_list[[i]]) <- paste0("domain_", dom_id)
  names(domain_abund_intensity_consensus_list[[i]]) <- paste0("domain_", dom_id)
  names(domain_biodiv_intensity_consensus_list[[i]]) <- paste0("domain_", dom_id)
}

domain_rich_top10_consensus <- rast(domain_rich_top10_consensus_list)
domain_abund_top10_consensus <- rast(domain_abund_top10_consensus_list)
domain_biodiv_top10_consensus <- rast(domain_biodiv_top10_consensus_list)

domain_rich_intensity_consensus <- rast(domain_rich_intensity_consensus_list)
domain_abund_intensity_consensus <- rast(domain_abund_intensity_consensus_list)
domain_biodiv_intensity_consensus <- rast(domain_biodiv_intensity_consensus_list)

names(domain_rich_top10_consensus) <- paste0("domain_", domain_ids)
names(domain_abund_top10_consensus) <- paste0("domain_", domain_ids)
names(domain_biodiv_top10_consensus) <- paste0("domain_", domain_ids)

names(domain_rich_intensity_consensus) <- paste0("domain_", domain_ids)
names(domain_abund_intensity_consensus) <- paste0("domain_", domain_ids)
names(domain_biodiv_intensity_consensus) <- paste0("domain_", domain_ids)

writeRaster(
  domain_rich_top10_consensus,
  filename = file.path(cons_dir, "domain_richness_hotspot_top10_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  domain_abund_top10_consensus,
  filename = file.path(cons_dir, "domain_abundance_hotspot_top10_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  domain_biodiv_top10_consensus,
  filename = file.path(cons_dir, "domain_biodiversity_hotspot_top10_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  domain_rich_intensity_consensus,
  filename = file.path(cons_dir, "domain_richness_hotspot_intensity_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  domain_abund_intensity_consensus,
  filename = file.path(cons_dir, "domain_abundance_hotspot_intensity_consensus.tif"),
  overwrite = TRUE
)

writeRaster(
  domain_biodiv_intensity_consensus,
  filename = file.path(cons_dir, "domain_biodiversity_hotspot_intensity_consensus.tif"),
  overwrite = TRUE
)

message("All hotspot products complete.")

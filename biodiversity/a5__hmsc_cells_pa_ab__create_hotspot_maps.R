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
  "fam_cafe", "fam_cbpm", "fam_eppl", "fam_vpmg",
  "npp_and_fam_cafe","npp_and_fam_cbpm","npp_and_fam_eppl","npp_and_fam_vpmg"
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

# COPY PREDICTIVE MAPS TO LOCAL DISK
# ----------------------------------------------------------
# Copy all predictive map files from the synced storage
# location to a local directory, preserving the folder
# structure for each model.
############################################################
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
# all_files <- all_files[which(grepl("npp_and_fam_cbpm/continuous", all_files, ignore.case = TRUE))]
# 
# rel_paths <- substring(all_files, nchar(src_base) + 2)
# dst_files <- file.path(dst_base, rel_paths)
# 
# dst_dirs <- unique(dirname(dst_files))
# invisible(lapply(dst_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
# 
# ok <- file.copy(from = all_files, to = dst_files, overwrite = TRUE, copy.mode = TRUE)
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
  
  richness_per_abundance_hot <- calc_hotspot_products(richness_per_abundance)
  abundance_per_richness_hot <- calc_hotspot_products(abundance_per_richness)
  
  ############################################################
  # 5.7 COMBINED HOTSPOT FILE FOR PLOTTING
  ############################################################
  message("Creating combined labelled hotspot file for plotting")
  top_percentiles <- build_all_hotspots_stack(
    richness_median = richness_median,
    abundance_median_100 = abundance_median_100,
    richness_per_abundance = richness_per_abundance,
    abundance_per_richness = abundance_per_richness,
    richness_hot_class = richness_hot[["intensity"]],
    abundance_hot_class = abundance_hot[["intensity"]],
    biodiversity_hot_class = biodiversity_hot[["intensity"]],
    richness_only_class = richness_only_class,
    abundance_only_class = abundance_only_class,
    richness_per_abundance_class = richness_per_abundance_hot[["intensity"]],
    abundance_per_richness_class = abundance_per_richness_hot[["intensity"]]
  )
  
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

  writeRaster(richness_per_abundance_hot,
              filename = file.path(out_dir_hot, "richness_per_abundance_hotspots.tif"), overwrite = TRUE)
  
  writeRaster(abundance_per_richness_hot,
              filename = file.path(out_dir_hot, "abundance_per_richness_hotspots.tif"), overwrite = TRUE)
  
  writeRaster(top_percentiles,
    filename = file.path(out_dir_hot, "all_hotspots_file.tif"), overwrite = TRUE)
  
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
    richness_per_abundance_hot, abundance_per_richness_hot, top_percentiles,
    domain_richness_pct, domain_abundance_pct, domain_biodiversity_pct,
    domain_richness_hot_top25, domain_richness_hot_top10, domain_richness_hot_top5, domain_richness_hot_top1,
    domain_abundance_hot_top25, domain_abundance_hot_top10, domain_abundance_hot_top5, domain_abundance_hot_top1,
    domain_biodiversity_hot_top25, domain_biodiversity_hot_top10, domain_biodiversity_hot_top5, domain_biodiversity_hot_top1,
    domain_richness_pct_list, domain_abundance_pct_list, domain_biodiversity_pct_list,
    domain_richness_hot_top25_list, domain_richness_hot_top10_list, domain_richness_hot_top5_list, domain_richness_hot_top1_list,
    domain_abundance_hot_top25_list, domain_abundance_hot_top10_list, domain_abundance_hot_top5_list, domain_abundance_hot_top1_list,
    domain_biodiversity_hot_top25_list, domain_biodiversity_hot_top10_list, domain_biodiversity_hot_top5_list, domain_biodiversity_hot_top1_list,
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
# 6) POST-PROCESSING: ADD HOTSPOT THRESHOLD LABELS
############################################################
# This section adds human-readable category labels to hotspot
# class rasters.
#
# Labels are based on the raw richness and abundance values
# associated with each hotspot class:
#
#   0 = outside top 25%
#   1 = top 25%
#   2 = top 10%
#   3 = top 5%
#   4 = top 1%
#
# For richness hotspot maps, labels are expressed as:
#   < / > X morphospecies
#
# For abundance hotspot maps, labels are expressed as:
#   < / > X %-cover
#
# For biodiversity hotspot maps, labels combine both:
#   < / > X %-cover & < / > Y morphospecies
#
# The labels are written into:
#   hotspots/all_hotspots_file.tif
#   domains/domain_richness_hotspot_class.tif
#   domains/domain_abundance_hotspot_class.tif
#   domains/domain_biodiversity_hotspot_class.tif
#
# This is a lightweight post-processing step. It does not
# recalculate hotspot classes.
############################################################


############################
# 6.1 HELPER FUNCTIONS
############################

# ----------------------------------------------------------
# Function: format_threshold_value
# ----------------------------------------------------------
# Formats numeric threshold values for legend labels.
# If the threshold is missing, returns "NA".
# ----------------------------------------------------------
format_threshold_value <- function(x, digits = 0) {
  if (is.na(x) || !is.finite(x)) {
    return("NA")
  }
  format(round(x, digits), trim = TRUE, nsmall = digits)
}


# ----------------------------------------------------------
# Function: min_surface_value_by_class
# ----------------------------------------------------------
# Returns the minimum value of a continuous raster surface
# within a specified hotspot class.
#
# Arguments:
#   surface      = continuous raster, e.g. richness or abundance
#   class_raster = hotspot class raster with integer classes
#   class_value  = hotspot class value to query
#
# This is used to derive the displayed threshold associated
# with each hotspot class.
# ----------------------------------------------------------
min_surface_value_by_class <- function(surface, class_raster, class_value) {
  vals <- values(c(surface, class_raster), mat = TRUE)
  
  surface_vals <- vals[, 1]
  class_vals   <- vals[, 2]
  
  sel <- which(
    !is.na(surface_vals) &
      is.finite(surface_vals) &
      !is.na(class_vals) &
      class_vals == class_value
  )
  
  if (length(sel) == 0) {
    return(NA_real_)
  }
  
  min(surface_vals[sel], na.rm = TRUE)
}


# ----------------------------------------------------------
# Function: hotspot_thresholds_from_classes
# ----------------------------------------------------------
# Calculates the raw threshold values associated with each
# hotspot class.
#
# For a normal hotspot raster:
#   class 1 = top 25%
#   class 2 = top 10%
#   class 3 = top 5%
#   class 4 = top 1%
#
# For contrast rasters such as richness_only or abundance_only,
# the maximum class may be less than 4. This function handles
# that automatically.
# ----------------------------------------------------------
hotspot_thresholds_from_classes <- function(surface, class_raster, max_class = NULL) {
  existing_classes <- sort(unique(values(class_raster, mat = FALSE, na.rm = TRUE)))
  existing_classes <- existing_classes[is.finite(existing_classes)]
  existing_classes <- existing_classes[existing_classes > 0]
  
  if (length(existing_classes) == 0) {
    return(numeric(0))
  }
  
  if (is.null(max_class)) {
    max_class <- max(existing_classes)
  }
  
  class_ids <- seq_len(max_class)
  
  thresholds <- sapply(class_ids, function(cl) {
    min_surface_value_by_class(
      surface = surface,
      class_raster = class_raster,
      class_value = cl
    )
  })
  
  names(thresholds) <- as.character(class_ids)
  thresholds
}


# ----------------------------------------------------------
# Function: make_single_metric_levels
# ----------------------------------------------------------
# Creates a raster attribute table for a single-metric hotspot
# raster, such as richness or abundance.
#
# Example output labels:
#   0 = < 18 morphospecies
#   1 = > 18 morphospecies
#   2 = > 24 morphospecies
#   3 = > 29 morphospecies
#   4 = > 37 morphospecies
# ----------------------------------------------------------
make_single_metric_levels <- function(
    class_raster,
    surface,
    unit_label,
    digits = 0,
    label_column = "label"
) {
  vals <- values(class_raster, mat = FALSE, na.rm = TRUE)
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) {
    return(data.frame(
      value = 0,
      label = paste0("No data"),
      stringsAsFactors = FALSE
    ))
  }
  
  max_class <- max(vals, na.rm = TRUE)
  thresholds <- hotspot_thresholds_from_classes(
    surface = surface,
    class_raster = class_raster,
    max_class = max_class
  )
  
  level_values <- 0:max_class
  level_labels <- character(length(level_values))
  
  # Class 0 is outside the top 25% hotspot threshold.
  level_labels[1] <- paste0(
    "< ",
    format_threshold_value(thresholds["1"], digits = digits),
    " ",
    unit_label
  )
  
  # Classes 1:max_class are increasingly intense hotspots.
  for (cl in seq_len(max_class)) {
    level_labels[cl + 1] <- paste0(
      "> ",
      format_threshold_value(thresholds[as.character(cl)], digits = digits),
      " ",
      unit_label
    )
  }
  
  out <- data.frame(
    value = level_values,
    label = level_labels,
    stringsAsFactors = FALSE
  )
  
  names(out)[2] <- label_column
  out
}


# ----------------------------------------------------------
# Function: make_biodiversity_levels
# ----------------------------------------------------------
# Creates raster attribute labels for biodiversity hotspot
# classes, using both the richness and abundance values within
# each biodiversity class.
#
# Example output labels:
#   0 = < 12 %-cover & < 18 morphospecies
#   1 = > 12 %-cover & > 18 morphospecies
#   2 = > 19 %-cover & > 24 morphospecies
#   3 = > 26 %-cover & > 29 morphospecies
#   4 = > 38 %-cover & > 37 morphospecies
# ----------------------------------------------------------
make_biodiversity_levels <- function(
    class_raster,
    richness_surface,
    abundance_surface,
    richness_digits = 0,
    abundance_digits = 0,
    label_column = "label"
) {
  vals <- values(class_raster, mat = FALSE, na.rm = TRUE)
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) {
    return(data.frame(
      value = 0,
      label = "No data",
      stringsAsFactors = FALSE
    ))
  }
  
  max_class <- max(vals, na.rm = TRUE)
  
  rich_thresholds <- hotspot_thresholds_from_classes(
    surface = richness_surface,
    class_raster = class_raster,
    max_class = max_class
  )
  
  abund_thresholds <- hotspot_thresholds_from_classes(
    surface = abundance_surface,
    class_raster = class_raster,
    max_class = max_class
  )
  
  level_values <- 0:max_class
  level_labels <- character(length(level_values))
  
  # Class 0 is outside the top 25% biodiversity hotspot threshold.
  level_labels[1] <- paste0(
    "< ",
    format_threshold_value(abund_thresholds["1"], digits = abundance_digits),
    " %-cover & < ",
    format_threshold_value(rich_thresholds["1"], digits = richness_digits),
    " morphospecies"
  )
  
  # Classes 1:max_class are increasingly intense biodiversity hotspots.
  for (cl in seq_len(max_class)) {
    level_labels[cl + 1] <- paste0(
      "> ",
      format_threshold_value(abund_thresholds[as.character(cl)], digits = abundance_digits),
      " %-cover & > ",
      format_threshold_value(rich_thresholds[as.character(cl)], digits = richness_digits),
      " morphospecies"
    )
  }
  
  out <- data.frame(
    value = level_values,
    label = level_labels,
    stringsAsFactors = FALSE
  )
  
  names(out)[2] <- label_column
  out
}


# ----------------------------------------------------------
# Function: set_layer_levels
# ----------------------------------------------------------
# Safely assigns category labels to a named layer within a
# multilayer SpatRaster.
# ----------------------------------------------------------
set_layer_levels <- function(x, layer_name, level_table) {
  if (!layer_name %in% names(x)) {
    warning("Layer not found: ", layer_name)
    return(x)
  }
  
  lyr <- x[[layer_name]]
  levels(lyr) <- level_table
  x[[layer_name]] <- lyr
  
  x
}


# ----------------------------------------------------------
# Function: assign_domain_levels
# ----------------------------------------------------------
# Assigns layer-specific category labels to a multilayer domain
# hotspot raster. Each layer corresponds to one planning domain.
# ----------------------------------------------------------
assign_domain_levels <- function(
    domain_class_raster,
    richness_surface,
    abundance_surface,
    metric = c("richness", "abundance", "biodiversity"),
    richness_digits = 0,
    abundance_digits = 0
) {
  metric <- match.arg(metric)
  
  level_list <- vector("list", nlyr(domain_class_raster))
  
  for (i in seq_len(nlyr(domain_class_raster))) {
    this_layer <- domain_class_raster[[i]]
    
    if (metric == "richness") {
      level_list[[i]] <- make_single_metric_levels(
        class_raster = this_layer,
        surface = richness_surface,
        unit_label = "morphospecies",
        digits = richness_digits,
        label_column = names(domain_class_raster)[i]
      )
    }
    
    if (metric == "abundance") {
      level_list[[i]] <- make_single_metric_levels(
        class_raster = this_layer,
        surface = abundance_surface,
        unit_label = "%-cover",
        digits = abundance_digits,
        label_column = names(domain_class_raster)[i]
      )
    }
    
    if (metric == "biodiversity") {
      level_list[[i]] <- make_biodiversity_levels(
        class_raster = this_layer,
        richness_surface = richness_surface,
        abundance_surface = abundance_surface,
        richness_digits = richness_digits,
        abundance_digits = abundance_digits,
        label_column = names(domain_class_raster)[i]
      )
    }
  }
  
  levels(domain_class_raster) <- level_list
  domain_class_raster
}


############################
# 6.2 APPLY LABELS TO EACH MODEL
############################

for (nm in model_ids) {
  
  message("======================================")
  message("Adding hotspot threshold labels for model: ", nm)
  message("======================================")
  
  model_out_dir <- file.path(src_base, paste0("hmsc_with_", nm))
  
  cont_dir   <- file.path(model_out_dir, "continuous")
  hot_dir    <- file.path(model_out_dir, "hotspots")
  domain_dir <- file.path(model_out_dir, "domains")
  
  required_files <- c(
    file.path(cont_dir, "richness_median.tif"),
    file.path(cont_dir, "abundance_median_100.tif"),
    file.path(hot_dir, "all_hotspots_file.tif"),
    file.path(domain_dir, "domain_richness_hotspot_class.tif"),
    file.path(domain_dir, "domain_abundance_hotspot_class.tif"),
    file.path(domain_dir, "domain_biodiversity_hotspot_class.tif")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    warning(
      "Skipping label assignment for ", nm, ". Missing files:\n",
      paste(missing_files, collapse = "\n")
    )
    next
  }
  
  ############################################################
  # Load continuous surfaces
  ############################################################
  
  richness_median <- rast(file.path(cont_dir, "richness_median.tif"))
  abundance_median_100 <- rast(file.path(cont_dir, "abundance_median_100.tif"))
  
  
  ############################################################
  # Label circumpolar hotspot display stack
  ############################################################
  
  message("  Labelling circumpolar hotspot display stack")
  
  top_percentiles <- rast(file.path(hot_dir, "all_hotspots_file.tif"))
  
  # Richness hotspots
  top_percentiles <- set_layer_levels(
    x = top_percentiles,
    layer_name = "richness",
    level_table = make_single_metric_levels(
      class_raster = top_percentiles[["richness"]],
      surface = richness_median,
      unit_label = "morphospecies",
      digits = 0,
      label_column = "richness"
    )
  )
  
  # Abundance hotspots
  top_percentiles <- set_layer_levels(
    x = top_percentiles,
    layer_name = "abundance",
    level_table = make_single_metric_levels(
      class_raster = top_percentiles[["abundance"]],
      surface = abundance_median_100,
      unit_label = "%-cover",
      digits = 0,
      label_column = "abundance"
    )
  )
  
  # Biodiversity hotspots
  top_percentiles <- set_layer_levels(
    x = top_percentiles,
    layer_name = "biodiversity",
    level_table = make_biodiversity_levels(
      class_raster = top_percentiles[["biodiversity"]],
      richness_surface = richness_median,
      abundance_surface = abundance_median_100,
      richness_digits = 0,
      abundance_digits = 0,
      label_column = "biodiversity"
    )
  )
  
  # Richness-only contrast hotspots
  if ("richness_only" %in% names(top_percentiles)) {
    top_percentiles <- set_layer_levels(
      x = top_percentiles,
      layer_name = "richness_only",
      level_table = make_biodiversity_levels(
        class_raster = top_percentiles[["richness_only"]],
        richness_surface = richness_median,
        abundance_surface = abundance_median_100,
        richness_digits = 0,
        abundance_digits = 0,
        label_column = "richness_only"
      )
    )
  }
  
  # Abundance-only contrast hotspots
  if ("abundance_only" %in% names(top_percentiles)) {
    top_percentiles <- set_layer_levels(
      x = top_percentiles,
      layer_name = "abundance_only",
      level_table = make_biodiversity_levels(
        class_raster = top_percentiles[["abundance_only"]],
        richness_surface = richness_median,
        abundance_surface = abundance_median_100,
        richness_digits = 0,
        abundance_digits = 0,
        label_column = "abundance_only"
      )
    )
  }
  
  # Richness per abundance hotspots
  if ("richness_per_abundance" %in% names(top_percentiles)) {
    rich_per_abund <- rast(file.path(cont_dir, "richness_per_abundance.tif"))
    
    top_percentiles <- set_layer_levels(
      x = top_percentiles,
      layer_name = "richness_per_abundance",
      level_table = make_single_metric_levels(
        class_raster = top_percentiles[["richness_per_abundance"]],
        surface = rich_per_abund,
        unit_label = "morphospecies for each 1% cover",
        digits = 1,
        label_column = "richness_per_abundance"
      )
    )
  }
  
  # Abundance per richness hotspots
  if ("abundance_per_richness" %in% names(top_percentiles)) {
    abund_per_rich <- rast(file.path(cont_dir, "abundance_per_richness.tif"))
    
    top_percentiles <- set_layer_levels(
      x = top_percentiles,
      layer_name = "abundance_per_richness",
      level_table = make_single_metric_levels(
        class_raster = top_percentiles[["abundance_per_richness"]],
        surface = abund_per_rich,
        unit_label = "%-cover for each morphospecies",
        digits = 1,
        label_column = "abundance_per_richness"
      )
    )
  }
  
  writeRaster(
    top_percentiles,
    filename = file.path(hot_dir, "all_hotspots_file_labelled.tif"),
    overwrite = TRUE
  )
  
  message("  Wrote labelled circumpolar file: ",
          file.path(hot_dir, "all_hotspots_file.tif"))
  
  
  ############################################################
  # Label domain-specific hotspot class rasters
  ############################################################
  
  message("  Labelling domain-specific hotspot rasters")
  
  domain_richness_hot_class <- rast(
    file.path(domain_dir, "domain_richness_hotspot_class.tif")
  )
  
  domain_abundance_hot_class <- rast(
    file.path(domain_dir, "domain_abundance_hotspot_class.tif")
  )
  
  domain_biodiversity_hot_class <- rast(
    file.path(domain_dir, "domain_biodiversity_hotspot_class.tif")
  )
  
  domain_richness_hot_class <- assign_domain_levels(
    domain_class_raster = domain_richness_hot_class,
    richness_surface = richness_median,
    abundance_surface = abundance_median_100,
    metric = "richness",
    richness_digits = 0,
    abundance_digits = 0
  )
  
  domain_abundance_hot_class <- assign_domain_levels(
    domain_class_raster = domain_abundance_hot_class,
    richness_surface = richness_median,
    abundance_surface = abundance_median_100,
    metric = "abundance",
    richness_digits = 0,
    abundance_digits = 0
  )
  
  domain_biodiversity_hot_class <- assign_domain_levels(
    domain_class_raster = domain_biodiversity_hot_class,
    richness_surface = richness_median,
    abundance_surface = abundance_median_100,
    metric = "biodiversity",
    richness_digits = 0,
    abundance_digits = 0
  )
  
  writeRaster(
    domain_richness_hot_class,
    filename = file.path(domain_dir, "domain_richness_hotspot_class_labelled.tif"),
    overwrite = TRUE
  )
  
  writeRaster(
    domain_abundance_hot_class,
    filename = file.path(domain_dir, "domain_abundance_hotspot_class_labelled.tif"),
    overwrite = TRUE
  )
  
  writeRaster(
    domain_biodiversity_hot_class,
    filename = file.path(domain_dir, "domain_biodiversity_hotspot_class_labelled.tif"),
    overwrite = TRUE
  )
  
  message("  Wrote labelled domain hotspot files")
  
  ############################################################
  # Clean memory
  ############################################################
  
  rm(
    richness_median,
    abundance_median_100,
    top_percentiles,
    domain_richness_hot_class,
    domain_abundance_hot_class,
    domain_biodiversity_hot_class
  )
  
  if (exists("rich_per_abund")) {
    rm(rich_per_abund)
  }
  
  if (exists("abund_per_rich")) {
    rm(abund_per_rich)
  }
  
  gc()
  
  message("Finished hotspot label assignment for model: ", nm)
}



# ############################################################
# # OLD TEMPORARY LOOP: CREATE all_hotspots_file.tif ONLY
# ############################################################
# 
# for (nm in model_ids) {
#   
#   message("======================================")
#   message("Creating all_hotspots_file for model: ", nm)
#   message("======================================")
#   
#   model_out_dir <- file.path(src_base, paste0("hmsc_with_", nm))
#   cont_dir <- file.path(model_out_dir, "continuous")
#   hot_dir  <- file.path(model_out_dir, "hotspots")
#   
#   required_files <- c(
#     file.path(cont_dir, "richness_median.tif"),
#     file.path(cont_dir, "abundance_median_100.tif"),
#     file.path(cont_dir, "richness_per_abundance.tif"),
#     file.path(cont_dir, "abundance_per_richness.tif"),
#     file.path(hot_dir, "richness_hotspots.tif"),
#     file.path(hot_dir, "abundance_hotspots.tif"),
#     file.path(hot_dir, "biodiversity_hotspots.tif"),
#     file.path(hot_dir, "richness_only_hotspots.tif"),
#     file.path(hot_dir, "abundance_only_hotspots.tif"),
#     file.path(hot_dir, "richness_per_abundance_hotspots.tif"),
#     file.path(hot_dir, "abundance_per_richness_hotspots.tif")
#   )
#   
#   missing_files <- required_files[!file.exists(required_files)]
#   
#   if (length(missing_files) > 0) {
#     warning(
#       "Skipping ", nm, ". Missing files:\n",
#       paste(missing_files, collapse = "\n")
#     )
#     next
#   }
#   
#   richness_median <- rast(file.path(cont_dir, "richness_median.tif"))
#   abundance_median_100 <- rast(file.path(cont_dir, "abundance_median_100.tif"))
#   richness_per_abundance <- rast(file.path(cont_dir, "richness_per_abundance.tif"))
#   abundance_per_richness <- rast(file.path(cont_dir, "abundance_per_richness.tif"))
#   
#   richness_hot <- rast(file.path(hot_dir, "richness_hotspots.tif"))
#   abundance_hot <- rast(file.path(hot_dir, "abundance_hotspots.tif"))
#   biodiversity_hot <- rast(file.path(hot_dir, "biodiversity_hotspots.tif"))
#   richness_only_class <- rast(file.path(hot_dir, "richness_only_hotspots.tif"))
#   abundance_only_class <- rast(file.path(hot_dir, "abundance_only_hotspots.tif"))
#   richness_per_abundance_hot <- rast(file.path(hot_dir, "richness_per_abundance_hotspots.tif"))
#   abundance_per_richness_hot <- rast(file.path(hot_dir, "abundance_per_richness_hotspots.tif"))
#   
#   top_percentiles <- build_all_hotspots_stack(
#     richness_median = richness_median,
#     abundance_median_100 = abundance_median_100,
#     richness_per_abundance = richness_per_abundance,
#     abundance_per_richness = abundance_per_richness,
#     richness_hot_class = richness_hot[["intensity"]],
#     abundance_hot_class = abundance_hot[["intensity"]],
#     biodiversity_hot_class = biodiversity_hot[["intensity"]],
#     richness_only_class = richness_only_class,
#     abundance_only_class = abundance_only_class,
#     richness_per_abundance_class = richness_per_abundance_hot[["intensity"]],
#     abundance_per_richness_class = abundance_per_richness_hot[["intensity"]]
#   )
#   
#   writeRaster(
#     top_percentiles,
#     filename = file.path(hot_dir, "all_hotspots_file.tif"),
#     overwrite = TRUE
#   )
#   
#   message("Wrote: ", file.path(hot_dir, "all_hotspots_file.tif"))
#   
#   rm(
#     richness_median, abundance_median_100,
#     richness_per_abundance, abundance_per_richness,
#     richness_hot, abundance_hot, biodiversity_hot,
#     richness_only_class, abundance_only_class,
#     richness_per_abundance_hot, abundance_per_richness_hot,
#     top_percentiles
#   )
#   gc()
# }


















##


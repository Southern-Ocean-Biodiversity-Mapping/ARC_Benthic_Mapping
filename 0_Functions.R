## 0_Functions.R

# ----------------------------------------------------------
# Function: topX.calc
# ----------------------------------------------------------
# Identifies the top (or bottom) X proportion of raster values
#
# Arguments:
#   ra  : input raster
#   x   : proportion threshold (e.g. 0.1 = top 10%)
#   low : if TRUE selects lowest X%, otherwise highest
#
# Returns:
#   Boolean raster (1 = hotspot, NA = masked cells)
# ----------------------------------------------------------
topX.calc <- function(ra, x = 0.1, low = FALSE) {
  # total number of values
  n_cells <- sum(!is.na(values(ra)))
  # Select exactly top X% cells
  topX <- selectHighest(ra, ceiling(x * n_cells), low = low)
  # Convert NA (outside selection) to 0
  topX[is.na(topX)] <- 0
  # Restore NA where original raster had NA
  topX[is.na(ra)] <- NA
  
  return(topX)
}

# # ----------------------------------------------------------
# # Function: calc_percentile_raster
# # ----------------------------------------------------------
# # Converts raster values into percentiles (0–1)
# #
# # Each cell value represents its rank relative to all cells
# #
# # Returns:
# #   Raster with values between 0 and 1
# # ----------------------------------------------------------
# calc_percentile_raster <- function(ra) {
#   vals <- values(ra)
#   # Rank values (ignoring NA)
#   ranks <- rank(vals, na.last = "keep")
#   # Convert ranks to percentiles
#   percentiles <- ranks / length(na.omit(vals))
#   # Assign back to raster
#   out <- setValues(ra, percentiles)
#   
#   return(out)
# }

# ----------------------------------------------------------
# Function: calc_percentile_rank_raster
# ----------------------------------------------------------
# Convert raster values to integer percentile ranks (1-100)
# Cells with NA remain NA
# -----------------------------------------------------------
calc_percentile_rank_raster <- function(ra) {
  vals <- values(ra)
  out_vals <- rep(NA_real_, length(vals))
  
  sel <- which(!is.na(vals))
  if (length(sel) == 0) {
    return(setValues(ra, out_vals))
  }
  
  ranks <- rank(vals[sel], ties.method = "average")
  pct <- ceiling(100 * ranks / length(sel))
  pct[pct < 1] <- 1
  pct[pct > 100] <- 100
  
  out_vals[sel] <- pct
  setValues(ra, out_vals)
}

# # ----------------------------------------------------------
# # Function: simplify_percentiles
# # ----------------------------------------------------------
# # Convert percentile ranks into simplified classes:
# # 0 = 0-75th percentile
# # 1 = >75th percentile
# # 2 = >90th percentile
# # 3 = >95th percentile
# # 4 = >99th percentile
# # ----------------------------------------------------------
# simplify_percentiles <- function(ra) {
#   vals <- values(ra)
#   out_vals <- rep(NA_real_, length(vals))
#   
#   sel <- which(!is.na(vals))
#   if (length(sel) == 0) {
#     return(setValues(ra, out_vals))
#   }
#   
#   out_vals[sel] <- 0
#   out_vals[sel][vals[sel] > 75] <- 1
#   out_vals[sel][vals[sel] > 90] <- 2
#   out_vals[sel][vals[sel] > 95] <- 3
#   out_vals[sel][vals[sel] > 99] <- 4
#   
#   setValues(ra, out_vals)
# }

# ----------------------------------------------------------
# Function: calc_hotspot_products
# ----------------------------------------------------------
# Create hotspot products for a raster using four thresholds:
# top 25%, top 10%, top 5%, top 1%
# The intensity layer is the sum of these four binary layers
# ----------------------------------------------------------
calc_hotspot_products <- function(ra, low = FALSE) {
  top25 <- topX.calc(ra, x = 0.25, low = low)
  top10 <- topX.calc(ra, x = 0.10, low = low)
  top5  <- topX.calc(ra, x = 0.05, low = low)
  top1  <- topX.calc(ra, x = 0.01, low = low)
  
  intensity <- sum(top1, top5, top10, top25, na.rm = TRUE)
  intensity[is.na(ra)] <- NA
  
  out <- c(top25, top10, top5, top1, intensity)
  names(out) <- c("top25", "top10", "top5", "top1", "intensity")
  out
}

# ----------------------------------------------------------
# Function: hotspot_polygons_from_class
# ----------------------------------------------------------
# Convert a classified hotspot raster to polygons
# Class 0 is excluded from polygon output
# ----------------------------------------------------------
hotspot_polygons_from_class <- function(ra, class_field = "class_id") {
  ra_poly <- ra
  ra_poly[ra_poly <= 0] <- NA
  
  vals <- values(ra_poly)
  if (all(is.na(vals))) {
    return(NULL)
  }
  
  poly <- as.polygons(ra_poly, dissolve = TRUE, values = TRUE)
  poly <- poly[!is.na(values(poly)), ]
  names(poly)[1] <- class_field
  poly
}

# ----------------------------------------------------------
# Function: bootstrap_richness_from_species_medians
# ----------------------------------------------------------
# Bootstrap uncertainty for species richness
#
# Richness is calculated as the sum of species-level median
# occurrence probabilities. Uncertainty is estimated from
# bootstrap resampling of morphospecies with replacement.
#
# For each bootstrap replicate:
# - morphospecies are resampled with replacement
# - resampled median occurrence probabilities are summed
#
# Returned rasters:
# - median richness across bootstrap replicates
# - standard deviation across bootstrap replicates
# ----------------------------------------------------------
bootstrap_richness_from_species_medians <- function(
    pa_median,
    n_boot = 1000,
    chunk_size = 10000,
    seed = 2
) {
  set.seed(seed)
  
  # Extract raster values as a matrix:
  # rows = cells, columns = species
  vals <- values(pa_median, mat = TRUE)
  
  if (is.null(dim(vals))) {
    stop("Species median raster stack must contain multiple layers.")
  }
  
  n_cells <- nrow(vals)
  n_species <- ncol(vals)
  
  # Record cells with complete species predictions
  complete_rows <- complete.cases(vals)
  
  # Prepare output vectors
  boot_median_vals <- rep(NA_real_, n_cells)
  boot_sd_vals <- rep(NA_real_, n_cells)
  
  if (sum(complete_rows) == 0) {
    boot_median <- setValues(pa_median[[1]], boot_median_vals)
    boot_sd <- setValues(pa_median[[1]], boot_sd_vals)
    return(list(median = boot_median, sd = boot_sd))
  }
  
  vals_complete <- vals[complete_rows, , drop = FALSE]
  
  # Create bootstrap count matrix
  # Each column gives the number of times each species is sampled
  # in one bootstrap replicate
  boot_counts <- replicate(
    n_boot,
    tabulate(sample.int(n_species, size = n_species, replace = TRUE), nbins = n_species)
  )
  
  # Process cells in chunks
  idx_complete <- which(complete_rows)
  n_complete <- length(idx_complete)
  
  chunk_starts <- seq(1, n_complete, by = chunk_size)
  
  for (s in chunk_starts) {
    e <- min(s + chunk_size - 1, n_complete)
    chunk_idx <- s:e
    
    # Cell x species matrix for this chunk
    vals_chunk <- vals_complete[chunk_idx, , drop = FALSE]
    
    # Cell x bootstrap matrix of richness replicates
    # Each bootstrap richness value is the sum of resampled species medians
    rich_boot_chunk <- vals_chunk %*% boot_counts
    
    # Per-cell bootstrap summaries
    boot_median_vals[idx_complete[chunk_idx]] <- apply(rich_boot_chunk, 1, median)
    boot_sd_vals[idx_complete[chunk_idx]] <- apply(rich_boot_chunk, 1, sd)
  }
  
  boot_median <- setValues(pa_median[[1]], boot_median_vals)
  boot_sd <- setValues(pa_median[[1]], boot_sd_vals)
  
  names(boot_median) <- "richness_bootstrap_median"
  names(boot_sd) <- "richness_bootstrap_sd"
  
  list(median = boot_median, sd = boot_sd)
}

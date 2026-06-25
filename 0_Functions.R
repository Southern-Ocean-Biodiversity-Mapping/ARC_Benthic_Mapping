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


# ----------------------------------------------------------
# Function: set_hotspot_plot_levels
# ----------------------------------------------------------
# Sets categorical labels for hotspot-class rasters.
#
# Input hotspot raster should contain integer classes:
#   0 = outside top hotspot class
#   1 = top 25%
#   2 = top 10%
#   3 = top 5%
#   4 = top 1%
#
# It automatically handles rasters that only have 0:2 or 0:3.
# ----------------------------------------------------------
set_hotspot_plot_levels <- function(hotspot,
                                    metrics,
                                    units,
                                    digits = 0,
                                    positive_sign = ">",
                                    zero_sign = "<",
                                    connector = " & ") {
  
  old_name <- names(hotspot)
  
  # Get classes present in the hotspot raster
  ff <- terra::freq(hotspot)
  vals <- sort(unique(ff$value[!is.na(ff$value)]))
  
  positive_classes <- sort(vals[vals > 0], decreasing = TRUE)
  
  if (length(positive_classes) == 0) {
    warning("No positive hotspot classes found in: ", old_name)
    return(hotspot)
  }
  
  # Allow one sign for all metrics, or one sign per metric
  if (length(positive_sign) == 1) {
    positive_sign <- rep(positive_sign, length(metrics))
  }
  if (length(zero_sign) == 1) {
    zero_sign <- rep(zero_sign, length(metrics))
  }
  
  if (length(metrics) != length(units)) {
    stop("metrics and units must have the same length.")
  }
  if (length(metrics) != length(positive_sign)) {
    stop("positive_sign must have length 1 or the same length as metrics.")
  }
  if (length(metrics) != length(zero_sign)) {
    stop("zero_sign must have length 1 or the same length as metrics.")
  }
  
  # Calculate minimum metric value within each positive hotspot class
  get_thresholds <- function(metric) {
    zz <- terra::zonal(metric, hotspot, fun = "min", na.rm = TRUE)
    threshold_vals <- setNames(zz[[2]], zz[[1]])
    round(threshold_vals[as.character(positive_classes)], digits)
  }
  
  thresholds <- lapply(metrics, get_thresholds)
  
  # Labels for positive hotspot classes, in descending class order
  positive_labels_desc <- vapply(seq_along(positive_classes), function(i) {
    paste(
      mapply(
        function(th, unit, sign) paste0(sign, " ", th[i], " ", unit),
        thresholds,
        units,
        positive_sign
      ),
      collapse = connector
    )
  }, character(1))
  
  # Label for class 0 uses the class-1 threshold
  zero_label <- paste(
    mapply(
      function(th, unit, sign) paste0(sign, " ", th[length(th)], " ", unit),
      thresholds,
      units,
      zero_sign
    ),
    collapse = connector
  )
  
  # terra levels need to be supplied in ascending raster-value order
  level_values <- sort(c(0, positive_classes))
  level_labels <- character(length(level_values))
  
  level_labels[level_values == 0] <- zero_label
  
  for (i in seq_along(positive_classes)) {
    this_class <- positive_classes[i]
    level_labels[level_values == this_class] <- positive_labels_desc[i]
  }
  
  hotspot <- as.factor(hotspot)
  levels(hotspot) <- data.frame(
    value = level_values,
    label = level_labels
  )
  
  names(hotspot) <- old_name
  hotspot
}


# ----------------------------------------------------------
# Function: build_all_hotspots_stack
# ----------------------------------------------------------
# Builds the seven-layer all_hotspots_file object with plot labels.
# ----------------------------------------------------------
build_all_hotspots_stack <- function(richness_median,
                                     abundance_median_100,
                                     richness_per_abundance,
                                     abundance_per_richness,
                                     richness_hot_class,
                                     abundance_hot_class,
                                     biodiversity_hot_class,
                                     richness_only_class,
                                     abundance_only_class,
                                     richness_per_abundance_class,
                                     abundance_per_richness_class) {
  
  top_percentiles <- c(
    richness_hot_class,
    abundance_hot_class,
    biodiversity_hot_class,
    richness_only_class,
    abundance_only_class,
    richness_per_abundance_class,
    abundance_per_richness_class
  )
  
  names(top_percentiles) <- c(
    "richness",
    "abundance",
    "biodiversity",
    "richness_only",
    "abundance_only",
    "richness_per_abundance",
    "abundance_per_richness"
  )
  
  top_percentiles[["biodiversity"]] <- set_hotspot_plot_levels(
    hotspot = top_percentiles[["biodiversity"]],
    metrics = list(
      abundance = abundance_median_100,
      richness = richness_median
    ),
    units = c("%-cover", "morphospecies"),
    digits = 0,
    positive_sign = ">",
    zero_sign = "<"
  )
  
  top_percentiles[["richness"]] <- set_hotspot_plot_levels(
    hotspot = top_percentiles[["richness"]],
    metrics = list(
      richness = richness_median
    ),
    units = c("morphospecies"),
    digits = 0,
    positive_sign = ">",
    zero_sign = "<"
  )
  
  top_percentiles[["abundance"]] <- set_hotspot_plot_levels(
    hotspot = top_percentiles[["abundance"]],
    metrics = list(
      abundance = abundance_median_100
    ),
    units = c("%-cover"),
    digits = 0,
    positive_sign = ">",
    zero_sign = "<"
  )
  
  top_percentiles[["richness_only"]] <- set_hotspot_plot_levels(
    hotspot = top_percentiles[["richness_only"]],
    metrics = list(
      abundance = abundance_median_100,
      richness = richness_median
    ),
    units = c("%-cover", "morphospecies"),
    digits = 0,
    positive_sign = ">",
    zero_sign = "<"
  )
  
  top_percentiles[["abundance_only"]] <- set_hotspot_plot_levels(
    hotspot = top_percentiles[["abundance_only"]],
    metrics = list(
      abundance = abundance_median_100,
      richness = richness_median
    ),
    units = c("%-cover", "morphospecies"),
    digits = 0,
    positive_sign = ">",
    zero_sign = "<"
  )
  
  top_percentiles[["richness_per_abundance"]] <- set_hotspot_plot_levels(
    hotspot = top_percentiles[["richness_per_abundance"]],
    metrics = list(
      richness_per_abundance = richness_per_abundance
    ),
    units = c("morphospecies for each 1% cover"),
    digits = 1,
    positive_sign = ">",
    zero_sign = "<"
  )
  
  top_percentiles[["abundance_per_richness"]] <- set_hotspot_plot_levels(
    hotspot = top_percentiles[["abundance_per_richness"]],
    metrics = list(
      abundance_per_richness = abundance_per_richness
    ),
    units = c("%-cover per morphospecies"),
    digits = 1,
    positive_sign = ">",
    zero_sign = "<"
  )
  
  top_percentiles
}



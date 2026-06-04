############################################################
# Environmental Gaps Analysis
# Setup and Data Loading
#
# Purpose
# -------
# This section defines project paths, loads required input
# data, standardises metadata fields used by downstream code,
# and prepares raster-derived environmental values needed for
# gap analyses and hypervolume workflows.
#
# Inputs
# ------
# - 0_SourceFile.R
# - cover_modelling_inputs_<res>.rds
# - Circumpolar_EnvData_<res>_shelf_mask_scaled.tif
#
# Outputs created here (if absent)
# --------------------------------
# - Circumpolar_EnvData_<res>_env_values_scaled.RData
#   containing:
#     * env_values : data.frame of all raster cell values
#     * na.sel     : integer vector of rows with missing values
#
############################################################


############################
# 0) USER SETTINGS
############################

# Select user profile defined in 0_SourceFile.R
# Valid options depend on that source file.
# Examples:
# usr <- "VM"
# usr <- "SJ"
usr <- "JJ"

# Load user-specific root directories and projection settings
source("0_SourceFile.R")

# Resolution label used in file names
res <- "2km"


############################
# 1) PACKAGES
############################

needed_pkgs <- c("terra", "hypervolume")

# Install any missing packages if required
# install.packages(setdiff(needed_pkgs, rownames(installed.packages())))

invisible(lapply(needed_pkgs, require, character.only = TRUE))


############################
# 2) DIRECTORY DEFINITIONS
############################

# Directory containing modelling-ready sampling inputs
sampling_dir <- file.path(
  usr.dropbox.dir,
  "data_products", "modelling_files", "circum_antarctic"
)

# Directory containing environmental rasters and derived raster products
env.derived <- file.path(
  usr.dropbox.dir,
  "data_environmental", "derived"
)

# Directory for gap-analysis outputs
gap_output_dir <- file.path(
  usr.dropbox.dir,
  "data_products", "environmental_gaps", "circum_antarctic"
)

# Create output directory if required
if (!dir.exists(gap_output_dir)) {
  dir.create(gap_output_dir, recursive = TRUE)
}


############################
# 3) DEFINE INPUT FILES
############################

# Sampling metadata and modelling inputs
sampling_file <- file.path(
  sampling_dir,
  paste0("cover_modelling_inputs_", res, ".rds")
)

# Scaled environmental raster stack used for environmental comparisons
scaled_raster_file <- file.path(
  env.derived,
  paste0("Circumpolar_EnvData_", res, "_shelf_mask_scaled.tif")
)

# Cached raster-value table used by hypervolume workflows
env_values_file <- file.path(
  env.derived,
  paste0("Circumpolar_EnvData_", res, "_env_values_scaled.RData")
)


############################
# 4) CHECK INPUT FILES
############################

if (!file.exists(sampling_file)) {
  stop("Sampling input file not found:\n  ", sampling_file)
}

if (!file.exists(scaled_raster_file)) {
  stop("Scaled environmental raster not found:\n  ", scaled_raster_file)
}


############################
# 5) LOAD SAMPLING METADATA
############################

# The modelling input object is expected to contain a component
# named `cell_metrics`, which stores one row per sampled cell and
# includes survey identifiers and projected coordinates.
dat <- readRDS(sampling_file)

if (!("cell_metrics" %in% names(dat))) {
  stop(
    "The sampling object does not contain a `cell_metrics` element.\n",
    "Expected object structure: dat$cell_metrics"
  )
}

img.metadata <- dat$cell_metrics

# Check that required fields are present
required_metadata_cols <- c("cell_id", "surveyID", "proj_coord_x", "proj_coord_y")
missing_metadata_cols <- setdiff(required_metadata_cols, names(img.metadata))

if (length(missing_metadata_cols) > 0) {
  stop(
    "The sampling metadata is missing required columns:\n  ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

# Standardise survey field for compatibility with legacy downstream code.
# A factor is used because later steps may rely on level-based indexing.
img.metadata$survey <- factor(img.metadata$surveyID)

# Ensure coordinates are numeric
img.metadata$proj_coord_x <- as.numeric(img.metadata$proj_coord_x)
img.metadata$proj_coord_y <- as.numeric(img.metadata$proj_coord_y)


############################
# 6) LOAD ENVIRONMENTAL RASTER STACK
############################

# The scaled environmental raster stack is used here because
# variables are already harmonised for multivariate analysis.
r.stack <- terra::rast(scaled_raster_file)

# Basic validation
if (terra::nlyr(r.stack) == 0) {
  stop("The environmental raster stack contains zero layers:\n  ", scaled_raster_file)
}

if (!("depth" %in% names(r.stack))) {
  stop(
    "Expected a layer named `depth` in the scaled environmental raster stack.\n",
    "Please verify the raster layer names in:\n  ", scaled_raster_file
  )
}


############################
# 7) PREPARE OR LOAD RASTER VALUE TABLE
############################

# `env_values` stores all raster cell values as a data.frame and
# `na.sel` stores row indices where at least one predictor is missing.
# These objects are reused throughout the Environmental Gaps script.

if (file.exists(env_values_file)) {
  
  load(env_values_file)
  
  # Validate loaded objects
  if (!exists("env_values") || !exists("na.sel")) {
    stop(
      "The cached environment-values file exists but does not contain ",
      "`env_values` and `na.sel`:\n  ", env_values_file
    )
  }
  
} else {
  
  # Extract raster values for all cells
  env_values.raw <- as.data.frame(terra::values(r.stack))
  
  # Identify rows where one or more layers are missing
  na.sel <- which(is.na(rowSums(env_values.raw)))
  
  # Store the table used by downstream scripts
  env_values <- env_values.raw
  
  # Save for reuse
  save(env_values, na.sel, file = env_values_file)
}

# Final consistency check
if (!is.data.frame(env_values)) {
  stop("`env_values` is not a data.frame after loading/creation.")
}

if (!is.numeric(na.sel)) {
  stop("`na.sel` is not numeric after loading/creation.")
}


############################
# 8) PREPARE COORDINATE TABLES
############################

# Full coordinate table used to identify raster cells sampled by all surveys
all.coords <- img.metadata[, c("proj_coord_x", "proj_coord_y")]

# Optional convenience object used in several parts of the script.
# This is retained to ease compatibility with existing code that
# expects x/y columns in a two-column data.frame.
colnames(all.coords) <- c("x", "y")


############################
# 9) OPTIONAL: REPORT BASIC INPUT SUMMARY
############################

message("Sampling metadata loaded: ", nrow(img.metadata), " sampled cells.")
message("Unique surveys: ", length(unique(img.metadata$survey)))
message("Environmental raster layers loaded: ", terra::nlyr(r.stack))
message("Raster-value table available with ", nrow(env_values), " rows.")


############################################################
# 10) DEFINE SURVEY GROUPS AND IDENTIFY SAMPLED RASTER CELLS
#
# Purpose
# -------
# This section:
#   1) Defines optional survey groupings (e.g., East vs West Antarctica)
#   2) Extracts projected coordinates for sampled cells
#   3) Converts coordinates to raster cell indices
#   4) Creates:
#        - cells           : all unique sampled raster cells
#        - cells.E         : all unique sampled raster cells in East group
#        - cells.W         : all unique sampled raster cells in West group
#        - cells.individual: list of unique sampled raster cells for each survey
#
# Notes
# -----
# - Survey grouping is defined using explicit survey names rather than
#   factor-level positions. This is more robust to changes in factor
#   ordering and metadata updates.
# - If East/West groupings are not required for a particular analysis,
#   the vectors `east_surveys` and `west_surveys` can be left empty.
# - Projected coordinates must be in the same CRS as the environmental
#   raster stack.
#
############################################################

############################
# 10.1) OPTIONAL: DEFINE SURVEY GROUPS
############################
# Define survey groups explicitly using survey identifiers present in `img.metadata$surveyID`.
west_surveys <- c("Antarctica_Peninsula_2015_JR15005", "Antarctica_Peninsula_2011_JR262", "Antarctica_Peninsula_2002_PS61",
                  "Antarctica_WeddellSea_1991_PS18", "Antarctica_Peninsula_2019_PS118", "Antarctica_WeddellSea_2016_PS96",
                  "Antarctica_WeddellSea_1989_PS14", "Antarctica_Peninsula_2013_PS81", "Antarctica_Peninsula_2017_JR17003",
                  "Antarctica_Peninsula_2017_JR17001", "Antarctica_Peninsula_2010_NBP1001", "Antarctica_WeddellSea_1985_PS06",
                  "Antarctica_Peninsula_2013_LMG1311", "Antarctica_Peninsula_2008_NBP0808", "Antarctica_Peninsula_2009_CRS")

east_surveys <- c("Antarctica_East_2014_NBP1402", "Antarctica_RossSea_2008_TAN0802", "Antarctica_RossSea_2015_NBP1502",  
                  "Antarctica_RossSea_2019_TAN1901", "Antarctica_RossSea_2018_TAN1802", "Antarctica_East_2011_AA2011")

# Basic checks on survey-group definitions
all_survey_ids <- unique(img.metadata$surveyID)

############################
# 10.2) EXTRACT PROJECTED COORDINATES
############################
# Coordinate table for all sampled records
all.coords <- as.data.frame(img.metadata[, c("proj_coord_x", "proj_coord_y")])
colnames(all.coords) <- c("x", "y")

# Coordinate tables for East and West survey groups
sel.E <- which(img.metadata$surveyID %in% east_surveys)
sel.W <- which(img.metadata$surveyID %in% west_surveys)

E.coords <- as.data.frame(img.metadata[sel.E, c("proj_coord_x", "proj_coord_y"), drop = FALSE])
W.coords <- as.data.frame(img.metadata[sel.W, c("proj_coord_x", "proj_coord_y"), drop = FALSE])
colnames(E.coords) <- c("x", "y")
colnames(W.coords) <- c("x", "y")


############################
# 10.4) IDENTIFY UNIQUE SAMPLED CELLS FOR EACH SURVEY
############################
cells <- img.metadata$cell_id

# Preserve survey order as it appears in the metadata
survey_ids <- unique(img.metadata$surveyID)

# Initialise list to store one integer vector of raster cells per survey
cells.individual <- vector("list", length(survey_ids))
names(cells.individual) <- survey_ids

for (i in seq_along(survey_ids)) {
  loop.sel <- which(img.metadata$surveyID == survey_ids[i])
  cells.individual[[i]] <- img.metadata[loop.sel, c("cell_id"), drop = FALSE]
}


############################
# calculate hypervolumes for all surveys
############################
model_vars_fixed <- c("depth","depth2","logslope","tpi","distance2canyons","distance2canyons2",
                "seafloortemperature","seafloorcurrents_mean","seafloorcurrents_residual","seafloorsalinity")
model_vars_swap_npp_mean <- c("cafe_mean", "cbpm_mean", "eppl_mean", "vpmg_mean")
model_vars_swap_npp_sd   <- c("cafe_sd", "cbpm_sd", "eppl_sd", "vpmg_sd")
model_vars_swap_fam_flx <- c("log.flux.mean.cafe", "log.flux.mean.cbpm", "log.flux.mean.eppl", "log.flux.mean.vpmg")
model_vars_swap_fam_sed <- c("sed.mean.cafe", "sed.mean.cbpm", "sed.mean.eppl", "sed.mean.vpmg")
## build pairs
npp_pairs <- data.frame(mean = model_vars_swap_npp_mean, sd   = model_vars_swap_npp_sd, stringsAsFactors = FALSE)
fam_pairs <- data.frame(flux = model_vars_swap_fam_flx, sed  = model_vars_swap_fam_sed, stringsAsFactors = FALSE)
##
env_stack <- list()
env_values_raw <- list()
env_stack_fixed <- r.stack[[names(r.stack)%in%model_vars_fixed]]
env_values_fixed <- data.frame(values(env_stack_fixed))
for(i in 1:4){
  sel <- c(model_vars_swap_npp_mean[i], model_vars_swap_npp_sd[i], model_vars_swap_fam_flx[i], model_vars_swap_fam_sed[i])
  env_stack_loop <- r.stack[[names(r.stack)%in%sel]]
  env_values_raw[[i]] <- cbind(env_values_fixed, data.frame(values(env_stack_loop)))
  env_stack[[i]] <- c(env_stack_fixed, env_stack_loop)
}
# 
# na.sel <- which(is.na(rowSums(env_values_raw[[1]])))
# # na.sel <- unique(c(which(is.na(r.stack$seafloortemperature[])),which(is.na(r.stack$depth[]))))
# for(i in 1:nlyr(env.stack)){
#   print(i)
#   env.stack[[i]][na.sel] <- NA
# }
# ## prepare environmental data layers
# env_values <- data.frame(values(env.stack))
# #save(env_values, na.sel, file=paste0(env.dir,"Circumpolar_EnvData_2km_env_values.Rdata"))
# load(paste0(env.dir,"Circumpolar_EnvData_2km_env_values.Rdata"))


## all environmental variables (takes ~5min)
Ant_dat <- data.frame(extract(env_stack[[1]],cells))
na.sel <- which(is.na(rowSums(Ant_dat)))
if(any(is.na(rowSums(Ant_dat)))){
  Ant_dat <- Ant_dat[-na.sel,]
}
hv_comb.all <- hypervolume(Ant_dat, name='comb', verbose=FALSE, quantile.requested = 0.95)
comb_inout.all <- hypervolume_inclusion_test(hv_comb.all, env_values_raw[[1]][-na.sel,], verbose=FALSE)
#save(hv_comb.all, comb_inout.all, file=paste0(gap_output_dir, "Circumpolar_Analysis_GapHypervolume_2km_AllSurveys_AllVariables_CafeNPPandFAM.Rdata"))
load(paste0(gap_output_dir,"Circumpolar_Analysis_GapHypervolume_2km_AllSurveys_AllVariables_CafeNPPandFAM.Rdata"))
r.gaps.all <- env.stack$depth
values(r.gaps.all)[-na.sel] <- 0
values(r.gaps.all)[-na.sel][which(comb_inout.all)] <- 1
plot(r.gaps.all)


##################################
#### ASSESS OUTPUTS

r.gaps.all090   <- rast(paste0(DP.dir,"Circumpolar_Analysis_GapHypervolume_2km_All21SurveysCombined_AllVariablesExceptDepth2Canyons_Quantile090.tif")           )
r.gaps.all099 <- rast(paste0(DP.dir,"Circumpolar_Analysis_GapHypervolume_2km_All21SurveysCombined_AllVariablesExceptDepth2Canyons_Quantile099.tif"))
r.gaps.all100 <- rast(paste0(DP.dir,"Circumpolar_Analysis_GapHypervolume_2km_All21SurveysCombined_AllVariablesExceptDepth2Canyons_Quantile100.tif"))

## environment inside vs outside
r.stack.subset <- subset(r.stack, c(1,29,31,3,7,30,25,22,23,24,15,27,28))
env.inside <-  mask(r.stack.subset, r.gaps.all099, maskvalues = 0)
env.outside <- mask(r.stack.subset, r.gaps.all099, maskvalues = 1)

vals.inside <- as.data.frame(values(env.inside))
vals.outside <- as.data.frame(values(env.outside))

vals.inside$group <- "inside"
vals.outside$group <- "outside"
vals.all <- rbind(vals.inside, vals.outside)

## density plots
# library(ggplot2)
# library(reshape2)
# melted <- melt(vals.all, id.vars = "group")
# ggplot(melted, aes(x = value, fill = group)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~variable, scales = "free") +
#   theme_minimal()

##
par(mfrow=c(2,3))
plot(r.gaps.all090, xlim=c(-500000,500000), ylim=c(-2000000,-1300000), main="gaps090")
plot(r.gaps.all099, xlim=c(-500000,500000), ylim=c(-2000000,-1300000), main="gaps099")
plot(r.gaps.all100, xlim=c(-500000,500000), ylim=c(-2000000,-1300000), main="gaps100")
plot(r.gaps.all090, xlim=c(-200000,100000), ylim=c(-1800000,-1500000), main="gaps090")
plot(r.gaps.all099, xlim=c(-200000,100000), ylim=c(-1800000,-1500000), main="gaps099")
plot(r.gaps.all100, xlim=c(-200000,100000), ylim=c(-1800000,-1500000), main="gaps100")
points(-50000,-1600000, col="red", pch=16)
points(-130000,-1700000, col="blue", pch=16)

## checking for two locations:
location_coords <- data.frame(x=c(-50000,-130000), y=c(-1600000,-1700000))
location_env <- terra::extract(r.stack.subset, location_coords)
# Calculate z-scores or percentiles
z_scores <- sapply(names(location_env)[-1], function(var) {
  (location_env[[var]] - mean(vals.inside[[var]], na.rm = TRUE)) / sd(vals.inside[[var]], na.rm = TRUE)
})
z_scores


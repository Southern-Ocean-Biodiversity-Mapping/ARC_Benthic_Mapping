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
usr <- "VM"
# usr <- "SJ"
#usr <- "JJ"

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
sampling_file <- file.path(sampling_dir,
  paste0("cover_modelling_inputs_", res, ".rds"))

# Scaled environmental raster stack used for environmental comparisons
scaled_raster_file <- file.path(env.derived,
  paste0("Circumpolar_EnvData_", res, "_shelf_mask_scaled.tif"))

# Cached raster-value table used by hypervolume workflows
env_values_file <- file.path(env.derived,
  paste0("Circumpolar_EnvData_", res, "_env_values_scaled.RData"))

############################
# 4) LOAD SAMPLING METADATA
############################

# The modelling input object is expected to contain a component
# named `cell_metrics`, which stores one row per sampled cell and
# includes survey identifiers and projected coordinates.
dat <- readRDS(sampling_file)
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
# 5) LOAD ENVIRONMENTAL RASTER STACK AND VALUE TABLE
############################
# The scaled environmental raster stack is used here because
# variables are already harmonised for multivariate analysis.
r.stack <- terra::rast(scaled_raster_file)

# `env_values` stores all raster cell values as a data.frame and
# `na.sel` stores row indices where at least one predictor is missing.
# These objects are reused throughout the Environmental Gaps script.
load(env_values_file)

############################
# 6) PREPARE COORDINATE TABLES
############################
# Full coordinate table used to identify raster cells sampled by all surveys
all.coords <- img.metadata[, c("proj_coord_x", "proj_coord_y")]

# Optional convenience object used in several parts of the script.
# This is retained to ease compatibility with existing code that
# expects x/y columns in a two-column data.frame.
colnames(all.coords) <- c("x", "y")


############################################################
# 7) DEFINE SURVEY GROUPS (IF APPLICABLE) AND IDENTIFY SAMPLED RASTER CELLS
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
# 7.1) OPTIONAL: DEFINE SURVEY GROUPS
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
# 7.2) EXTRACT PROJECTED COORDINATES
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
# 7.4) IDENTIFY UNIQUE SAMPLED CELLS FOR EACH SURVEY
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
# 8) calculate hypervolumes for all surveys
############################
model_vars_fixed <- c("depth","slope","tpi","distance2canyons",
                "seafloortemperature","seafloorcurrents_mean","seafloorcurrents_residual","seafloorsalinity")
model_vars_swap_npp_mean <- c("cafe_mean", "cbpm_mean", "eppl_mean", "vpmg_mean")
model_vars_swap_fam_flx <- c("log.flux.mean.cafe", "log.flux.mean.cbpm", "log.flux.mean.eppl", "log.flux.mean.vpmg")
model_vars_swap_fam_sed <- c("sed.mean.cafe", "sed.mean.cbpm", "sed.mean.eppl", "sed.mean.vpmg")
## build pairs
npp_fam <- data.frame(npp = model_vars_swap_npp_mean, flux = model_vars_swap_fam_flx, sed  = model_vars_swap_fam_sed, stringsAsFactors = FALSE)
##
model.names <- c("CafeNPPandFAM","CbpmNPPandFAM","EpplNPPandFAM","VpmgNPPandFAM")
##
env_stack <- list()
env_values_raw <- list()
env_stack_fixed <- r.stack[[names(r.stack)%in%model_vars_fixed]]
env_values_fixed <- data.frame(values(env_stack_fixed))
for(i in 1:4){
  ## select raster layers
  env_stack_loop <- r.stack[[names(r.stack)%in%npp_fam[i,]]]
  ## extract values
  env_values_raw[[i]] <- cbind(env_values_fixed, data.frame(values(env_stack_loop)))
  ## prepare raster stack
  env_stack[[i]] <- c(env_stack_fixed, env_stack_loop)
}

#### prep and save environmental variables
for(i in 1:4){
  print(i)
  Ant_dat <- data.frame(extract(env_stack[[i]],cells))
  na.sel <- which(is.na(rowSums(Ant_dat)))
  if(any(is.na(rowSums(Ant_dat)))){
    Ant_dat <- Ant_dat[-na.sel,]
  }
  env.na.sel <- which(is.na(rowSums(env_values_raw[[i]])))
  env_values_red <- env_values_raw[[i]][-env.na.sel,]
  print("env-prep done")
  save(env_values_red, env.na.sel, file=file.path(env.derived,paste0("Circumpolar_EnvData_2km_env_values_scaled_forGaps_",model.names[i],".Rdata")))
}

## ~1h per loop
for(i in 1:4){
  print(i)
  ## load environment
  load(file.path(env.derived,paste0("Circumpolar_EnvData_2km_env_values_scaled_forGaps_",model.names[i],".Rdata")))

  ####  calculate hypervolume, takes ~1h
  ## option SVM is not good, way to narrow volume
  file.basename <- file.path(gap_output_dir, paste0("Circumpolar_Analysis_GapHypervolume_2km_AllSurveys_AllVariables_ExceptDepth2_",model.names[i]))
  hv_comb.all <- hypervolume(Ant_dat, name='comb', verbose=FALSE, quantile.requested = 0.95)
  print("hypervolume calculated")
  ## 
  ptm <- proc.time()
  comb_inout.all <- hypervolume_inclusion_test(hv_comb.all, env_values_red, verbose=FALSE)#, reduction.factor = 0.1)
  print(runtime <- proc.time() - ptm)
  save(hv_comb.all, comb_inout.all, env.na.sel, file=paste0(file.basename,".Rdata"))
  print("inclusion test finished")
  
  #### save tif files
  r.gaps.all <- rast(r.stack$depth)
  r.gaps.all[-env.na.sel] <- 0
  values(r.gaps.all)[-env.na.sel][which(comb_inout.all)] <- 1
  writeRaster(r.gaps.all, filename=paste0(file.basename,".tif"), overwrite=TRUE)
}


## INDIVIDUAL SURVEYS ? per loop
for(i in 1:4){
  message(i)
  for(j in 1:length(cells.individual)){
    s.name <- names(cells.individual)[j]
    cells.ind <- as.vector(cells.individual[[j]])$cell_id
    print(j)
  #### prep and save environmental variables
  Ant_dat <- data.frame(extract(env_stack[[i]],cells.ind))
  na.sel <- which(is.na(rowSums(Ant_dat)))
  if(any(is.na(rowSums(Ant_dat)))){
    Ant_dat <- Ant_dat[-na.sel,]
  }
  env.na.sel <- which(is.na(rowSums(env_values_raw[[i]])))
  env_values_red <- env_values_raw[[i]][-env.na.sel,]
  print("env-prep done")
  save(env_values_red, env.na.sel, file=file.path(env.derived,paste0("Circumpolar_EnvData_2km_env_values_scaled_forGaps_",model.names[i],"_",s.name,".Rdata")))
  
  ####  calculate hypervolume, takes ~1h
  ## option SVM is not good, way to narrow volume
  file.basename <- file.path(gap_output_dir, paste0("Circumpolar_Analysis_GapHypervolume_2km_IndividualSurveys_AllVariables_ExceptDepth2_",model.names[i],"_",s.name))
  hv_comb.all <- hypervolume(Ant_dat, name='comb', verbose=FALSE, quantile.requested = 0.95)
  print("hypervolume calculated")
  ## 
  ptm <- proc.time()
  comb_inout.all <- hypervolume_inclusion_test(hv_comb.all, env_values_red, verbose=FALSE)#, reduction.factor = 0.1)
  print(runtime <- proc.time() - ptm)
  save(hv_comb.all, comb_inout.all, env.na.sel, file=paste0(file.basename,".Rdata"))
  print("inclusion test finished")
  
  #### save tif files
  r.gaps.all <- rast(r.stack$depth)
  r.gaps.all[-env.na.sel] <- 0
  values(r.gaps.all)[-env.na.sel][which(comb_inout.all)] <- 1
  writeRaster(r.gaps.all, filename=paste0(file.basename,".tif"), overwrite=TRUE)
  }
}



# ## ~1h per loop
# for(i in 1:4){
#   print(i)
#   load(file.path(env.derived,paste0("Circumpolar_EnvData_2km_env_values_scaled_forGaps_",model.names[i],".Rdata")))
#   ####  calculate hypervolume, takes ~1h
#   file.basename <- file.path(gap_output_dir, paste0("Circumpolar_Analysis_GapHypervolume_2km_AllSurveys_AllVariables_",model.names[i]))
#   load(paste0(file.basename,".Rdata"))
#   #### save tif files
#   r.gaps.all <- rast(r.stack$depth)
#   r.gaps.all[-env.na.sel] <- 0
#   values(r.gaps.all)[-env.na.sel][which(comb_inout.all)] <- 1
#   writeRaster(r.gaps.all, filename=paste0(file.basename,".tif"), overwrite=TRUE)
# }
############################
# environmental conditions inside/outside
############################




# 
# ##################################
# #### ASSESS OUTPUTS
# 
# r.gaps.all090   <- rast(paste0(DP.dir,"Circumpolar_Analysis_GapHypervolume_2km_All21SurveysCombined_AllVariablesExceptDepth2Canyons_Quantile090.tif")           )
# r.gaps.all099 <- rast(paste0(DP.dir,"Circumpolar_Analysis_GapHypervolume_2km_All21SurveysCombined_AllVariablesExceptDepth2Canyons_Quantile099.tif"))
# r.gaps.all100 <- rast(paste0(DP.dir,"Circumpolar_Analysis_GapHypervolume_2km_All21SurveysCombined_AllVariablesExceptDepth2Canyons_Quantile100.tif"))
# 
# ## environment inside vs outside
# r.stack.subset <- subset(r.stack, c(1,29,31,3,7,30,25,22,23,24,15,27,28))
# env.inside <-  mask(r.stack.subset, r.gaps.all099, maskvalues = 0)
# env.outside <- mask(r.stack.subset, r.gaps.all099, maskvalues = 1)
# 
# vals.inside <- as.data.frame(values(env.inside))
# vals.outside <- as.data.frame(values(env.outside))
# 
# vals.inside$group <- "inside"
# vals.outside$group <- "outside"
# vals.all <- rbind(vals.inside, vals.outside)
# 
# ## density plots
# # library(ggplot2)
# # library(reshape2)
# # melted <- melt(vals.all, id.vars = "group")
# # ggplot(melted, aes(x = value, fill = group)) +
# #   geom_density(alpha = 0.5) +
# #   facet_wrap(~variable, scales = "free") +
# #   theme_minimal()
# 
# ##
# par(mfrow=c(2,3))
# plot(r.gaps.all090, xlim=c(-500000,500000), ylim=c(-2000000,-1300000), main="gaps090")
# plot(r.gaps.all099, xlim=c(-500000,500000), ylim=c(-2000000,-1300000), main="gaps099")
# plot(r.gaps.all100, xlim=c(-500000,500000), ylim=c(-2000000,-1300000), main="gaps100")
# plot(r.gaps.all090, xlim=c(-200000,100000), ylim=c(-1800000,-1500000), main="gaps090")
# plot(r.gaps.all099, xlim=c(-200000,100000), ylim=c(-1800000,-1500000), main="gaps099")
# plot(r.gaps.all100, xlim=c(-200000,100000), ylim=c(-1800000,-1500000), main="gaps100")
# points(-50000,-1600000, col="red", pch=16)
# points(-130000,-1700000, col="blue", pch=16)
# 
# ## checking for two locations:
# location_coords <- data.frame(x=c(-50000,-130000), y=c(-1600000,-1700000))
# location_env <- terra::extract(r.stack.subset, location_coords)
# # Calculate z-scores or percentiles
# z_scores <- sapply(names(location_env)[-1], function(var) {
#   (location_env[[var]] - mean(vals.inside[[var]], na.rm = TRUE)) / sd(vals.inside[[var]], na.rm = TRUE)
# })
# z_scores


library(dplyr)
library(raster)
library(usdm)
library(blockCV)
library(sf)

# Load bio data
df_bio = read.csv("biodata_step4.csv")

# Load metadata
load('/perm_storage/shared_space/BioMAS/ARC_Data/annotation/Circumpolar_Annotation_Env_Data_500m.RData')
# Get metadata
df_metadata = cell_metadata_env[c("cellID", "lon", "lat", "proj_coord_x",
                                  "proj_coord_y", "image_quality_score", "year",
                                  "gear")]
# Load env data
fname_env<-'/perm_storage/shared_space/BioMAS/environmental_data/Circumpolar_EnvData_500m_shelf_mask_scaled.tif' 
raster_env=stack(fname_env)
df_env = raster::extract(x=raster_env,
                         y=cbind(df_metadata$proj_coord_x, df_metadata$proj_coord_y),
                         cellnumbers=TRUE)
df_env_ = cbind(cbind(df_metadata$proj_coord_x, df_metadata$proj_coord_y),
                data.frame(df_env))
names(df_env_)[1:2] <- c('proj_coord_x', 'proj_coord_y')
df_env_$geomorphology <- as.factor(as.character(df_env_$geomorphology))

# Get cells of interest
cell_metadata = df_metadata$cellID
cell_bio = df_bio$cellID
cell_of_interest = intersect(cell_metadata, cell_bio)
df_metadata_cOi = df_metadata[df_metadata$cellID %in% cell_of_interest, ]
df_env_cOi = df_env_[df_env_$cells %in% cell_of_interest, ]

# Get covariates of interest
cov_lst = c("depth",
            "depth2",
            "tpi11",
            "logslope",
            "tpi",
            "seafloorcurrents_mean",
            "seafloorcurrents_residual",
            "seafloorcurrents_absolute",
            "seafloortemperature",
            "seafloorsalinity",
            "test_settle08",
            "test_susp08",
            "distance2canyons",
            "distance2canyons2")
df_env_cOi_ = df_env_cOi[, c(colnames(df_env_cOi)[1:3], cov_lst)]
df_env_cOi_[df_env_cOi_=='NaN'] <- NA
df_env_cOi_ <- dplyr::mutate_all(df_env_cOi_, function(x) as.numeric(as.character(x)))

# VIF
vif_ = usdm::vifstep(df_env_cOi_[cov_lst], th=10)
vif_
df_env_vif = usdm::exclude(df_env_cOi_[cov_lst], vif_)
df_env_clean = cbind(df_env_cOi_[1:3], df_env_vif)

# Check correlations
suppressWarnings(chart.Correlation(df_env_clean[, 4:ncol(df_env_clean)], histogram=TRUE, pch=19))

# Remove NAs
df_env_clean = df_env_clean[complete.cases(df_env_clean), ]
names(df_env_clean)[names(df_env_clean) == 'cells'] <- 'cellID'
df_metadata_clean = df_metadata_cOi[df_metadata_cOi$cellID %in% df_env_clean$cellID, ]
df_bio_clean = df_bio[df_bio$cellID %in% df_env_clean$cellID, ]

# BlockCV
sf_dat <- st_as_sf(df_metadata_clean,
                   coords = c("proj_coord_x", "proj_coord_y"),
                   crs = crs(raster_env))
spatial_block <- spatialBlock(speciesData = sf_dat,
                              species = NULL,
                              theRange = 2000000,
                              k = 5,
                              selection = "random",
                              iteration = 100,
                              seed=7109)
df_tmp <- as.data.frame(sf_dat)[c('cellID')]
df_tmp$fold <- spatial_block$foldID
df_metadata_clean <- merge(df_tmp, df_metadata_clean, by="cellID")

# Order by cellID
df_bio_clean = df_bio_clean[order(df_bio_clean$cellID), ]
df_env_clean = df_env_clean[order(df_env_clean$cellID), ]
df_metadata_clean = df_metadata_clean[order(df_metadata_clean$cellID), ]

# Save
save(df_bio_clean, df_metadata_clean, df_env_clean,
     file ="modelling_data.RData")

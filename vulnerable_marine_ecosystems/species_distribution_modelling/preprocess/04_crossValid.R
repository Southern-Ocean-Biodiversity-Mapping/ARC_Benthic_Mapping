# Todo:
#   

#https://cran.r-project.org/web/packages/blockCV/vignettes/BlockCV_for_SDM.html
# install.packages("blockCV", dependencies = TRUE)
library(blockCV)
library(raster)
library(sf)

# import raster data
raster_env <- raster::brick("Circumpolar_EnvData_bathy500m_shelf_gebco2020_depth.grd")
# load data
load("03_data.RData")

df_pa_ <- merge(df_pa, df_metadata, by="cellID")

sf_bio <- st_as_sf(df_pa_,
                   coords = c("proj_coord_x", "proj_coord_y"),
                   crs = crs(raster_env))

sb <- blockCV::spatialBlock(speciesData = sf_bio,
                   species = NULL,
                   #species = "pa_data",
                   #rasterLayer = raster_env,
                   theRange = 100000,
                   k = 5,
                   selection = "random",
                   iteration = 100,
                   seed=7109)

df_ <- as.data.frame(sf_bio)[c('cellID')]
df_$fold <- sb$foldID

df_metadata_ <- merge(df_, df_metadata, by="cellID")
#sf_result <- st_as_sf(df_metadata_,
#                      coords = c("proj_coord_x", "proj_coord_y"),
#                      crs = crs(raster_env))
#plot(sf_result["fold"])

# Order by cellID
df_abd = df_abd[order(df_abd$cellID), ]
df_pa = df_pa[order(df_pa$cellID), ]
df_metadata_ = df_metadata_[order(df_metadata_$cellID), ]
df_env_vif = df_env_vif[order(df_env_vif$cellID), ]
df_env_scaled = df_env_scaled[order(df_env_scaled$cellID), ]

save(df_abd, df_pa, df_metadata_, df_env_vif, df_env_scaled, scale_means, scale_vars, file ="04_data.RData")

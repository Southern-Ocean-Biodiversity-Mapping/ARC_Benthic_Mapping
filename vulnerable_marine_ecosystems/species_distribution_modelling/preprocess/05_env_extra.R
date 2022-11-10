# Todo:
#   Todo on Desktop R Studio
#   Correct waom4k_settle6test

library(tidyverse)
library(raster)
library(rgdal)
library(rasterVis)
library(stringr)

# Load data
load("04_data.RData")

# Params
#cell_centroid = 65476765
#cell_centroid = "random"
#if (cell_centroid == "random") {
#  cell_centroid = sample(df_metadata_$cellID, 1)
#}
n_cell = 20
env.derived <-  "C:/Users/cgros/data/env_layers/derived"
n_area = 10

# Get cells of interest
size_half = 500 * n_cell / 2
lst_area <- c()
for (i_area in 1:n_area) {
  cell_centroid = sample(df_metadata_$cellID, 1)
  print(cell_centroid)
  x_centroid = df_metadata_[which(df_metadata_$cellID == cell_centroid), "proj_coord_x"]
  y_centroid = df_metadata_[which(df_metadata_$cellID == cell_centroid), "proj_coord_y"]
  proj_x_00 = c(x_centroid - size_half)
  proj_y_00 = c(y_centroid - size_half)
  lst_proj_x = c()
  lst_proj_y = c()
  for (x_cell in 0:(n_cell-1)) {
    x_cur = proj_x_00 + x_cell * 500
    for (y_cell in 0:(n_cell-1)) {
      y_cur = proj_y_00 + y_cell * 500
      lst_proj_x = append(lst_proj_x, x_cur)
      lst_proj_y = append(lst_proj_y, y_cur)
    }
  }
  if (i_area == 1) {
    df_extrapol = data.frame(proj_coord_x=lst_proj_x, proj_coord_y=lst_proj_y)
  }
  else {
    df_extrapol = rbind(df_extrapol,
                        data.frame(proj_coord_x=lst_proj_x, proj_coord_y=lst_proj_y))
  }
}
head(df_extrapol)

# Load raster data
# Get file names of all environmental rasters and bricks and load into one big stack----
# All files with "gri" extension
env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
# Subset to  "shelf" files
env_list<-env_list[grep("shelf", env_list)]
# For the single rasters layer names are missing. Extract from file name.
env_names<-gsub(".*_|\\..*","",env_list)
# Stack all environmental layers and make sure they have appropriate names
env_stack<-stack(env_list)
names(env_stack)
names(env_stack)[1:5]<-env_names[1:5]
names(env_stack)[14:22]<-paste(rep(c("CARS_NO3", "CARS_O2", "CARS_PO4"),each=3),c("mean", "seas_range", "std_dev"), sep="_")
names(env_stack)[23] <-"distance2canyons"
names(env_stack)[34]<-"NPP_su_mean"
names(env_stack)[35:40]<-c("ssh_mean","ssh_sd","ssh_sp_mean","ssh_sp_sd","ssh_su_mean","ssh_su_sd")
names(env_stack)[41:46]<-c("sst_mean","sst_sd","sst_sp_mean","sst_sp_sd","sst_su_mean","sst_su_sd")
names(env_stack)[47:56]<-c("waom2k_seafloorcurrents", "waom2k_seafloortemperature", "waom4k_seafloorcurrents_absolute", "waom4k_seafloorcurrents_mean", 
                           "waom4k_seafloorcurrents_residual", "waom4k_seafloorsalinity", "waom4k_seafloortemperature",
                           "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08")

#add environmental data with non-conformant names- 
#### remember to update column index if changes!!!
env_stack<-stack( env_stack,
                  raster(paste0(env.derived, "/Circumpolar_EnvData_geomorphology")))
names(env_stack)[57]<-"geomorph"

geomorph_cat<-levels(env_stack[[57]])[[1]]
names(env_stack)

#extract environmental data
df_extrapol <- cbind(df_extrapol, raster::extract(env_stack, df_extrapol, cellnumbers=TRUE))
names(df_extrapol)

save(df_abd, df_pa, df_metadata_, df_env_vif, df_env_scaled, scale_means, scale_vars, df_extrapol,
     file ="05_data.RData")

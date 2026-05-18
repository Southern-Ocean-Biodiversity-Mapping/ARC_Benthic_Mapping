##########
## WHAT THIS SCRIPT DOES:
## - reading in fitted models (from "biodiversity_models_b_fit...") and environmental data
## - predicting responses and CIs for each model
## - saving the predictions and plotting the results 
##########


##### Setting up----
library(raster)
# library(readxl)
# library(readr)
# library(dplyr)
# library(data.table)
# library(proj4)
# library(stringr)
# library(RColorBrewer)
# library(SOmap)

user = "Jan"
#user = "charley"
#user="nicole"
if (user == "Jan") {
  sci.dir <-      "C:/Users/jjansen/Desktop/science/"
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  #bio.dir <-      paste0(sci.dir,"data_biological/")
  ## remote repository (DOESN'T WORK YET):
  # env.dir <- "https://data.imas.utas.edu.au/data_transfer/admin/files/EnvironmentalData/"
  ## common paths (after "sci.dir")
  tools.dir <-    paste0(sci.dir,"SouthernOceanBiodiversityMapping/Useful_Functions_Tools/")
  ARC_Data.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Data/")
} 
if (user == "charley") {
  sci.dir <- "C:/Users/cgros/code/IMAS/"
  ARC_Data.dir <- paste0(sci.dir,"ARC_Data/")
  env.derived <-  "C:/Users/cgros/data/SO_env_layers/derived/"
  tools.dir <-    paste0(sci.dir,"Useful_Functions_Tools/")
}
if (user == "nicole") {
  sci.dir <-    "C:/Users/hillna/OneDrive - University of Tasmania/UTAS_work/Projects/Benthic Diversity ARC/"
  ARC_Data.dir <- paste0(sci.dir,"Analysis/ARC_Data/")
  env.derived <-  paste0(sci.dir,"data_environmental/derived/")
  tools.dir <-    paste0(sci.dir,"Analysis/Useful_Functions_Tools/")
}
##############################################################################################################
##############################################################################################################
biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")
load(file=paste0(biodiv.dir,"/biodiversity_env_dat.Rdata"))
load(file=paste0(biodiv.dir,"/biodiversity_bio_dat.Rdata"))

## load hmsc output
load(paste0(biodiv.dir,"pred_files/model_1_pa_chains_2_thin_10_samples_100_pred_x1_y5.Rdata"))

# ## functions
# source(paste0(tools.dir,"SOmap_functions_JJ.R"))
# 
# ## projection
# stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# 
# ## bathymetry
# ## from "ReadIn_Circumpolar_Environmental_Data.Rmd"
r2 <- raster(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
# 
## load coastline
load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))


hmsc_preds_bryo <- rast(paste0(biodiv.dir,list.files(biodiv.dir, pattern=".tif")))
names(hmsc_preds_bryo) <- dimnames(predY.mean)[[2]]
par(mfrow=c(2,2))
plot(hmsc_preds_bryo, range=c(0,1), xlim=c(-2800000,-2000000), ylim=c(1400000,2400000))
plot(coast.proj, add=TRUE)





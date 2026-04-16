
library(terra)
library(viridis)
library(adespatial)

##########################
sci.dir <-      "C:/Users/jjansen/OneDrive - University of Tasmania/science/"
env.derived <-  paste0(sci.dir,"data_environmental/derived/")
bio.dir <-      paste0(sci.dir,"data_biological/")
pred.dir <- paste0(bio.dir,"circumpolar_prediction_outputs/")
tools.dir <-    paste0(sci.dir,"SouthernOceanBiodiversityMapping/Useful_Functions_Tools/")
DP.dir <- paste0(sci.dir,"DP190101858_MappingAntarcticSeafloorBiodiversity/")
ARC_Data.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Data/")

#########################

##########################
#### read in and compare hmsc species predictions, prepared in "display_predictive_maps_setup.R"
#ra.list <- list.files(pred.dir, pattern=".tif")
ra.list.pa <- list.files(paste0(pred.dir,"2km_model_cells_1_pa/"), pattern=".tif", full.names=TRUE)
ra.list.ab <- list.files(paste0(pred.dir,"2km_model_cells_1_abundcondpres/"), pattern=".tif", full.names=TRUE)

#ra.list.pa.all <- ra.list.pa[grep("1_pa_",ra.list.pa)]
ra.list.pa.env <- ra.list.pa[grep("envonly",basename(ra.list.pa))]

## P/A maps
#hmsc.maps.pa <- rast(paste0(pred.dir,"old/",ra.list.pa))
hmsc.maps.pa.env <- rast(ra.list.pa.env)

#hmsc.maps.pa.median     <- subset(hmsc.maps.pa, seq(2,nlyr(hmsc.maps.pa),by=4))
hmsc.maps.pa.env.median <- subset(hmsc.maps.pa.env, seq(2,nlyr(hmsc.maps.pa.env),by=5))

##############################
# #### testing on a small region first
# test.ra <- crop(hmsc.maps.pa.env.median, ext(-200000,100000,-1700000,-1400000))
# plot(test.ra[[1]])
# 
# dat <- values(test.ra)
# na.sel <- which(is.na(rowSums(dat)))
# dat.clean <- dat[-na.sel,]
# ## ~ 1.5min for 22k rows
# ## ~ 4min for 35k rows
# ## 50k rows crashes, 110k cells isn't possible (90GB vector)
# s.time <- Sys.time()
# dat.beta <- adespatial::beta.div(dat.clean)
# Sys.time()-s.time
# 
# test.ra.beta <- rast(test.ra[[1]])
# values(test.ra.beta)[-na.sel] <- dat.beta$LCBD
# plot(test.ra.beta)
# 
# #### sampling randomly and calculating betadiv works, maps looks the same as above
# ## sample cells regularly to allow computation
# ## 20 runs with every 20th cell
# dat.beta.lcbd <- rep(NA, nrow(dat.clean))
# step_size <- 20
# sample.size <- nrow(dat.clean)/step_size
# for(i in 1:step_size){
#   print(i)
#   if(i == 1){
#     selected_indices <- sample(1:nrow(dat.clean), sample.size, replace = FALSE)
#   }else{
#     # Generate a random subset of indices for the current loop, excluding already selected indices
#     available_indices <- setdiff(1:nrow(dat.clean), selected_indices)
#     indices <- sample(available_indices, sample.size, replace = FALSE)
#     #indices <- seq(i, nrow(dat.clean), by = step_size) ## regular sampling introduces artefacts
#     # Update the selected indices
#     selected_indices <- c(selected_indices, indices)
#   }
#   ##
#   dat.clean.loop <- dat.clean[selected_indices,]
#   dat.beta.loop <- adespatial::beta.div(dat.clean.loop)
#   dat.beta.lcbd[selected_indices] <- dat.beta.loop$LCBD
# }
# 
# looped.ra.beta <- rast(test.ra[[1]])
# values(looped.ra.beta)[-na.sel] <- dat.beta.lcbd
# plot(looped.ra.beta)
########################################
#### testing on a large region
test.ra <- crop(hmsc.maps.pa.env.median, ext(-500000,500000,-2100000,-1200000))
plot(test.ra[[1]])

dat <- values(test.ra)
na.sel <- which(is.na(rowSums(dat)))
dat.clean <- dat[-na.sel,]

## 20 runs with every 20th cell
#dat.beta.lcbd <- rep(NA, nrow(dat.clean))
step_size <- 20
sample.size <- floor(nrow(dat.clean)/step_size)
sample.size1 <- nrow(dat.clean)-(step_size-1)*sample.size
for(i in 1:step_size){
  print(i)
  if(i == 1){
    selected_indices <- sample(1:nrow(dat.clean), sample.size1, replace = FALSE)
  }else{
    # Generate a random subset of indices for the current loop, excluding already selected indices
    available_indices <- setdiff(1:nrow(dat.clean), selected_indices)
    indices <- sample(available_indices, sample.size, replace = FALSE)
    #indices <- seq(i, nrow(dat.clean), by = step_size) ## regular sampling introduces artefacts
    # Update the selected indices
    selected_indices <- c(selected_indices, indices)
  }
  ##
  dat.clean.loop <- dat.clean[selected_indices,]
  dat.beta.loop <- adespatial::beta.div(dat.clean.loop)
  save(dat.beta.loop, file=paste0(pred.dir,"2km_model_cells_1_pa/","BetaDiversity_TestSample",sprintf("%03d", i),".Rdata"))
}

dat.beta.lcbd[selected_indices] <- dat.beta.loop$LCBD


ross.ra.beta <- rast(test.ra[[1]])
values(ross.ra.beta)[-na.sel] <- dat.beta.lcbd
plot(ross.ra.beta)



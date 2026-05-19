##############################################################################################################
##### SELECT MODEL SETUP
##############################################################################################################
library(terra)
library(Hmsc)
usr <- "VM"
source("0_SourceFile.R")

env_dir   <- paste0(usr.dropbox.dir, "data_environmental/derived")
model_dir   <- paste0(usr.dropbox.dir, "data_products/modelling_files/circum_antarctic")
pred_dir   <- paste0(usr.main.dir, "4_model_prediction/pred_files/")

res <- "2km"

thin     <- 10
samples  <- 800
nChains  <- 4

modelspec <- "_envonly"

se <- function(x) {
  sd(x) / sqrt(length(x))
}
p5 <- function(x, z=1.96) {
  mean(x) - z * sd(x) / sqrt(length(x))
}
p95 <- function(x, z=1.96) {
  mean(x) + z * sd(x) / sqrt(length(x))
}
# '%!in%' <- function(x,y)!('%in%'(x,y))

# nm <- "fam_cafe"   # <<<<<< CHANGE THIS PER RUN
# model_file <- file.path(
#   "/pvol/2_fitting_and_running_models/",
#   paste0(res, "_model_cells_", nm,
#     "_chains_", nChains,
#     "_thin_", thin,
#     "_samples_", samples,
#     ".Rdata"
#   )
# )
# load(model_file)  # loads: models, XFormula, runtime
# # Environment-only model (PA)
# pa <- models$mENV
# # Abundance model
# ab <- models$mAB
# rm(models)

model_ids <- c(
  "npp_cafe","npp_cbpm","npp_eppl","npp_vpmg",
  "fam_cafe","fam_cbpm","fam_eppl","fam_vpmg"
)

#############################################################
###### prepare data (RUN ONCE)
# # Load prediction dataframe (already scaled)
# pred_df <- readRDS(file.path(env_dir, paste0("Circumpolar_EnvData_", res, "_shelf_mask_scaled_dataframe.rds")))
# 
# ## Identify valid shelf cells
# r.stack <- rast(file.path(env_dir, paste0("Circumpolar_EnvData_", res, "_shelf_mask_unscaled_variables.tif")))
# r2 <- r.stack$depth
# sel.not.na <- which(!is.na(r2[]))
# 
# ## add sampling effort (Set constant value)
# pred_df$cover_points_scorable <- 540
# 
# ##
# grid <- pred_df
# sel <- which(complete.cases(grid))
# 
# ## Extract coordinates
# xy.grid.raw <- pred_df[, c("proj_coord_x", "proj_coord_y")]
# 
# ## prepare data to predict on
# XData.grid <- grid[sel, ]
# xy.grid    <- xy.grid.raw[sel, ]

# ## Save for reuse
# save(sel, sel.not.na, file = file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_sel.Rdata")))
# save(XData.grid, xy.grid, file = file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_grid.Rdata")))

load(file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_sel.Rdata")))
load(file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_grid.Rdata")))


#############################################################
## TILE SETUP (KEEP 100 km STRUCTURE)
# Tile boundaries (100 km)
xmin <- seq(-3000000, 2900000, by = 100000)
xmax <- seq(-2900000, 3000000, by = 100000)
ymin <- seq(-3000000, 2900000, by = 100000)
ymax <- seq(-2900000, 3000000, by = 100000)

#############################################################
## LOAD TILE LOOKUP
## create a look-up table to check which cells we need to predict
## keep in mind that the raster starts at the bottom left and the matrix start filling in values from the top left!
# cells_with_data <- matrix(NA, nrow=length(ymin), ncol=length(xmin))
# for(i in 1:length(xmin)){
#   print(i)
#   for(k in 1:length(ymin)){
#     sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
#                         xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
#     ## fill the matrix from the bottom up
#     #x.sel <- (length(xmin):1)[i]
#     #y.sel <- (length(ymin):1)[k]
#     if(length(sel.loop>0)){
#       cells_with_data[k,i] <- 1
#     }
#   }}
# save(cells_with_data, file=file.path(model_dir, paste0(res,"_model_100km_cells_with_data.Rdata")))
load(file = file.path(model_dir, paste0(res, "_model_100km_cells_with_data.Rdata")))

cell.sel.v <- which(!is.na(cells_with_data))
cell.sel.df <- which(!is.na(cells_with_data), arr.ind = TRUE)

#############################################################
## RUN ALL MODELS
#############################################################
#############################################################
for (nm in model_ids) {
  ## parallel processing: PER CELL that contains values
  library(doParallel)
  library(foreach)
  parallel::detectCores()
  #UseCores = parallel::detectCores() - 1
  UseCores = 4
  c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
  registerDoParallel(c1)
  getDoParWorkers()
  
  message("====================================")
  message("Running prediction for model: ", nm)
  message("====================================")
  
  model_file <- file.path(model.dir, paste0("2_fitting_and_running_models/", res, "_model_cells_", nm, "_chains_", nChains, "_thin_", thin, "_samples_", samples, ".Rdata")
  )
  load(model_file)
  pa <- models$mENV
  ab <- models$mAB
  rm(models)
  
  # predictor_names <- colnames(pa$XData)
  # grid <- pred_df[, predictor_names, drop = FALSE]
  # sel <- which(complete.cases(grid))
  # XData.grid <- grid[sel, ]
  # xy.grid    <- xy.grid.raw[sel, ]

  #iterations <- 60
  
  message("setup complete, parallel loops over prediction areas now")
  
  ptm = proc.time()
  foreach(j=1:length(cell.sel.v), .packages = c("Hmsc")) %dopar%{ #3:length(xmin)
  # for(j in 1:length(cell.sel.v)){
    i <- cell.sel.df[j,2]
    k <- cell.sel.df[j,1]
    ## select cells in tile
    sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
                        xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
    if (length(sel.loop) == 0) {
      return(NULL)
    }
    
    if (j %% 50 == 0) message("Tile ", j)
    ## subset predictors + coordinates
    XData.grid.loop <- XData.grid[sel.loop,]
    xy.grid.loop <- xy.grid[sel.loop,]
    
    ## setup prediction - pa
    Gradient.pa = prepareGradient(pa, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
    predY.loop.pa <- predict(pa, Gradient=Gradient.pa, expected=TRUE) ## this gives probabilities instead of integer outcomes
    rm(Gradient.pa)
    # mat.names <- dimnames(predY.loop.pa[[1]])
    # predY.loop.pa <- array(unlist(predY.loop.pa), c(nrow(xy.grid.loop), ncol(pa$Y), samples*nChains), dimnames(predY.loop.pa[[1]]))
    
    ## setup prediction - abund
    Gradient.ab = prepareGradient(ab, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
    predY.loop.ab <- predict(ab, Gradient=Gradient.ab, expected=TRUE) ## this gives probabilities instead of integer outcomes
    rm(Gradient.ab)  
    # mat.names <- dimnames(predY.loop.ab[[1]])
    # predY.loop.ab <- array(unlist(predY.loop.ab), c(nrow(xy.grid.loop), ncol(ab$Y), samples*nChains), dimnames(predY.loop.ab[[1]]))
    
    # Convert matrices
    predY.loop.pa.array = simplify2array(predY.loop.pa)
    predY.loop.ab.array = simplify2array(predY.loop.ab)
    
    # Get posterior median
    predY.pa.mean   = as.data.frame(apply(predY.loop.pa.array, 1:2, mean))
    predY.pa.median = as.data.frame(apply(predY.loop.pa.array, 1:2, median))
    predY.ab.mean   = as.data.frame(apply(predY.loop.ab.array, 1:2, mean))
    predY.ab.median = as.data.frame(apply(predY.loop.ab.array, 1:2, median))
    
    # Posterior width
    predY.pa.se = as.data.frame(apply(predY.loop.pa.array, 1:2, se))
    predY.pa.5  = as.data.frame(apply(predY.loop.pa.array, 1:2, p5))
    predY.pa.95 = as.data.frame(apply(predY.loop.pa.array, 1:2, p95))
    predY.ab.se = as.data.frame(apply(predY.loop.ab.array, 1:2, se))
    predY.ab.5  = as.data.frame(apply(predY.loop.ab.array, 1:2, p5))
    predY.ab.95 = as.data.frame(apply(predY.loop.ab.array, 1:2, p95))
    
    ## save-string for 100km cell tiles
    dat.name <- paste0(pred_dir, res,"_model_cells_",nm, "_chains_",as.character(nChains),"_thin_", as.character(thin),"_samples_", as.character(samples),"_pred_")
    run.name <- sprintf("%06d",cell.sel.v[j])
    
    save(predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.5, predY.ab.95,
         predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95,
         sel.loop, XData.grid.loop, xy.grid.loop,
         file=paste0(dat.name,modelspec,run.name,".Rdata"))
    rm(predY.loop.pa, predY.loop.ab, predY.loop.pa.array, predY.loop.ab.array,
       predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.5, predY.ab.95,
       predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95)
  }
  computational.time = proc.time() - ptm
  message("Finished model: ", nm)
  parallel::stopCluster(cl = c1)
}




























































## Predicting a p/a combined with abundance hmsc onto the entire Antarctic continental shelf
## prediction is splitted into chunks of small areas to speed up the process

library(Hmsc)
library(terra)
'%!in%' <- function(x,y)!('%in%'(x,y))

se <- function(x) {
  sd(x) / sqrt(length(x))
}
p5 <- function(x, z=1.96) {
  mean(x) - z * sd(x) / sqrt(length(x))
}
p95 <- function(x, z=1.96) {
  mean(x) + z * sd(x) / sqrt(length(x))
}

##############################################################################################################
## select resolution
#res <- "500m"
res <- "2km"

## select which model to run
#model.sel <- 1 ## full model
model.sel <- 2 ## environment only model

if(model.sel==2){
  modelspec <- "_envonly"
} else modelspec <- ""
##############################################################################################################
env.derived <- "/pvol3TB/data_environmental/"

r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth
## create an empty raster to fill for mapping
empty.ra <- rast(r2)
empty.ra[] <- NA

## specify model to load
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 800 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
nChains = 4

## presence absence model
modeltype = 1
model = 1
load(paste0("/pvol3TB/2_fitting_and_running_models/",res,"_model_cells_",model,"_pa_chains_4_thin_",thin,"_samples_",samples,".Rdata"))
pa <- models[[model.sel]]

## abundance model
model = 1
load(paste0("/pvol3TB/2_fitting_and_running_models/",res,"_model_cells_",model,"_abundcondpres_chains_4_thin_",thin,"_samples_",samples,".Rdata"))
ab <- models[[model.sel]]

rm(models)

## predicted vs observed
#pred.pa <- computePredictedValues(pa)


#############################################################
###### prepare data (RUN ONCE)
## load per cell environmental data
load(paste0("/pvol3TB/data_environmental/Circumpolar_EnvData_",res,"_shelf_mask_scaled_dataframe.Rdata"))
## identify which cells we have ignored before
sel.not.na <- which(!is.na(r2[]))
#depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+seafloortemperature+seafloorcurrents_mean+seafloorcurrents_residual+seafloorsalinity+npp_mean+log.flux.mean+sed.mean+log(cover_points_scorable)
## relevant env data
grid <- pred_stack.dat[,c(1,29,31,3,7,30,25,22,23,24,15,27,28)]
## spatial information
#xy <- crds(r2)
xy.grid.raw <- crds(r2) #xy #[sel.not.na,]
## remove NAs
sel <- which(complete.cases(grid))
XData.grid <- grid[sel,]
xy.grid <- data.frame(xy.grid.raw[sel,])

## save data
save(sel, sel.not.na, file=paste0("/pvol3TB/4_model_prediction/hmsc_",res,"_model_cell_sel.Rdata"))
save(XData.grid, xy.grid, file=paste0("/pvol3TB/4_model_prediction/hmsc_",res,"_model_cell_grid.Rdata"))
rm(xy.grid.raw, grid, pred_stack.dat)

#############################################################
## load data
load(paste0("/pvol3TB/4_model_prediction/hmsc_",res,"_model_cell_sel.Rdata"))
load(paste0("/pvol3TB/4_model_prediction/hmsc_",res,"_model_cell_grid.Rdata"))

## predict to 5 images per cell
XData.grid$cover_points_scorable <- 540

# ## size of prediction boxes (100km most time efficient)
xmin <- seq(-3000000,2900000, by=100000)
xmax <- seq(-2900000,3000000, by=100000)
ymin <- seq(-3000000,2900000, by=100000)
ymax <- seq(-2900000,3000000, by=100000)

# plot(r2)
# for(i in 2:length(xmin)){
#   abline(h=ymin[i])
#   abline(v=xmin[i])
# }

#############################################################
##### best/fastest resolution to run on the VM is 100km cells
## at 150km cells, speed is ? tiles per hour (3 cores not working!), with 374 tiles in total ->  ?days
## at 100km cells, speed is 4 tiles per hour (4 cores), with 537 tiles in total ->  5.5days
## at 50km cells, speed is 10 tiles per hour (10 cores), with 2.3k tiles in total ->  9.5days
## at 2km cells, speed is ? tiles per hour (32 cores), with 1m tiles in total -> ? days
## at 10km cells, speed is 32 tiles per hour (32 cores), with 43k tiles in total -> 55 days
## at 20km cells, speed is 32 tiles per hour (32 cores), with 11.8k tiles in total -> 15 days
## at 30km cells, speed is 20 tiles per hour (20 cores, error on 24 cores), with 5.5k tiles in total -> 11.5 days days

## create a look-up table to check which cells we need to predict
## keep in mind that the raster starts at the bottom left and the matrix start filling in values from the top left!
# cells_with_data <- matrix(NA, nrow=length(ymin), ncol=length(xmin))
# for(i in 1:length(xmin)){
#   print(i)
#   for(k in 1:length(ymin)){
#     sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
#                         xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
#     ## fill the matrix from the bottom up
#     #x.sel <- (length(xmin):1)[i]
#     #y.sel <- (length(ymin):1)[k]
#     if(length(sel.loop>0)){
#       cells_with_data[k,i] <- 1
#     }
#   }}
# save(cells_with_data, file=paste0("/pvol/4_model_prediction/",res,"_model_150km_cells_with_data.Rdata"))
# load(file=paste0("/pvol/4_model_prediction/",res,"_model_150km_cells_with_data.Rdata"))
# save(cells_with_data, file=paste0("/pvol/4_model_prediction/",res,"_model_50km_cells_with_data.Rdata"))
# load(file=paste0("/pvol/4_model_prediction/",res,"_model_50km_cells_with_data.Rdata"))
# save(cells_with_data, file=paste0("/pvol/4_model_prediction/",res,"_model_30km_cells_with_data.Rdata"))
# load(file=paste0("/pvol/4_model_prediction/",res,"_model_30km_cells_with_data.Rdata"))
# save(cells_with_data, file=paste0("/pvol/4_model_prediction/",res,"_model_20km_cells_with_data.Rdata"))
# load(file=paste0("/pvol/4_model_prediction/",res,"_model_20km_cells_with_data.Rdata"))
# save(cells_with_data, file=paste0("/pvol/4_model_prediction/",res,"_model_10km_cells_with_data.Rdata"))
# load(file=paste0("/pvol/4_model_prediction/",res,"_model_10km_cells_with_data.Rdata"))

# save(cells_with_data, file=paste0("/pvol3TB/4_model_prediction/",res,"_model_100km_cells_with_data.Rdata"))
load(file=paste0("/pvol3TB/4_model_prediction/",res,"_model_100km_cells_with_data.Rdata"))

#################
## hmsc seems to predict onto a new level of a non-spatial random effect
## the prepareGradient function creates an object with 20 levels instead of the original 19 in the hmsc model object
## I don't understand why it's displaying a warning message though
# library(dplyr)
# XData.grid.loop$surveyID <- factor("PS96", levels=m$ranLevels$surveyID$pi)

##### Preds for only Ross Sea and Denman:
# plot(r2)
# for(i in 2:length(xmin)){
#   abline(h=ymin[i])
#   abline(v=xmin[i])
# }
# abline(v=xmin[26], col="blue")
# abline(v=xmin[36], col="blue")
# abline(v=xmin[55], col="blue")
# abline(v=xmin[60], col="blue")
# 
# abline(h=ymin[22], col="red")
# abline(h=ymin[30], col="red")
# abline(h=ymin[20], col="red")
# abline(h=ymin[8], col="red")
# 
# ## Ross Sea
# cells_with_data[,c(1:25,37:60)] <- NA
# cells_with_data[c(1:7,21:60),] <- NA
# 
# ## Denman
# cells_with_data[,c(1:54)] <- NA
# cells_with_data[c(1:21,31:60),] <- NA
# 

#############################################################
## parallel processing: PER CELL that contains values
library(doParallel)
library(foreach)
parallel::detectCores()
#UseCores = parallel::detectCores() - 1
UseCores = 2
c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
registerDoParallel(c1)
getDoParWorkers()

cell.sel.v <- which(!is.na(cells_with_data))
cell.sel.df <- which(!is.na(cells_with_data), arr.ind = TRUE)

#iterations <- 60

ptm = proc.time()
#foreach(j=1:iterations) %dopar%{ #3:length(xmin)
foreach(j=1:length(cell.sel.v)) %dopar%{ #3:length(xmin)
#foreach(j=1:96) %dopar%{ #3:length(xmin)
  #library(Hmsc)
  i <- cell.sel.df[j,2]
  k <- cell.sel.df[j,1]
  sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
                      xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
  print(i)
  ##
  XData.grid.loop <- XData.grid[sel.loop,]
  ## moving the survey-ID into the fixed component doesn't seem to work
  #XData.grid.loop$surveyID <- factor("PS96", levels=m$ranLevels$surveyID$pi)
  xy.grid.loop <- xy.grid[sel.loop,]
 
  ## setup prediction - pa
  Gradient.pa = prepareGradient(pa, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
  predY.loop.pa <- predict(pa, Gradient=Gradient.pa, expected=TRUE) ## this gives probabilities instead of integer outcomes
  rm(Gradient.pa)
  #predY.loop.pa <- predict(pa, XData=XData.grid.loop, expected=TRUE)
  
  # mat.names <- dimnames(predY.loop.pa[[1]])
  # predY.loop.pa <- array(unlist(predY.loop.pa), c(nrow(xy.grid.loop), ncol(pa$Y), samples*nChains), dimnames(predY.loop.pa[[1]]))

  ## setup prediction - abund
  Gradient.ab = prepareGradient(ab, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
  predY.loop.ab <- predict(ab, Gradient=Gradient.ab, expected=TRUE) ## this gives probabilities instead of integer outcomes
  rm(Gradient.ab)  
  #predY.loop.ab <- predict(ab, XData=XData.grid.loop, expected=TRUE)
  
  # mat.names <- dimnames(predY.loop.ab[[1]])
  # predY.loop.ab <- array(unlist(predY.loop.ab), c(nrow(xy.grid.loop), ncol(ab$Y), samples*nChains), dimnames(predY.loop.ab[[1]]))
  
  # Multiply matrices
  predY.loop = simplify2array(predY.loop.pa) * simplify2array(predY.loop.ab)
  predY.loop.pa.array = simplify2array(predY.loop.pa)
  
  # Get posterior median
  predY.pa.mean = as.data.frame(apply(predY.loop.pa.array, 1:2, mean))
  predY.pa.median = as.data.frame(apply(predY.loop.pa.array, 1:2, median))
  predY.mean = as.data.frame(apply(predY.loop, 1:2, mean))
  predY.median = as.data.frame(apply(predY.loop, 1:2, median))
  
  # Posterior width
  predY.pa.se = as.data.frame(apply(predY.loop.pa.array, 1:2, se))
  predY.pa.5 = as.data.frame(apply(predY.loop.pa.array, 1:2, p5))
  predY.pa.95 = as.data.frame(apply(predY.loop.pa.array, 1:2, p95))
  predY.se = as.data.frame(apply(predY.loop, 1:2, se))
  predY.5 = as.data.frame(apply(predY.loop, 1:2, p5))
  predY.95 = as.data.frame(apply(predY.loop, 1:2, p95))
 
  # ## predict
  # ## derived values
  # predY.mean <- apply(predY.loop, 1:2, mean)
  # predY.sd <- apply(predY.loop, 1:2, sd)
  # dimnames(predY.mean) <- dimnames(predY.sd) <- mat.names
  
  ## save-string for 10km cell tiles
  dat.name <- paste0("/pvol3TB/4_model_prediction/pred_files",modelspec,"/",res,"/abundance/",
                     res,"_model_cells_",as.character(model), "_abundcondpres_chains_",as.character(nChains),"_thin_", ... = as.character(thin),"_samples_", as.character(samples),"_pred_")
  run.name <- sprintf("%06d",cell.sel.v[j])
  ## save-string for 50km cell tiles
  # dat.name <- paste0("/pvol/4_model_prediction/pred_files/",res,"/abundance/",
  #                    res,"_model_cells_",as.character(model), "_",c("pa","abundance")[modeltype],
  #                    "_chains_",as.character(nChains),"_thin_", ... = as.character(thin),"_samples_", as.character(samples),"_pred_")
  # run.name <- sprintf("%05d",cell.sel.v[j])
  ## save and then remove objects
  # save(predY.loop, file=paste0(dat.name,"fulldat_",run.name,".Rdata"))
  save(predY.mean, predY.median, predY.se, predY.5, predY.95,
       predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95,
       sel.loop, XData.grid.loop, xy.grid.loop,
       file=paste0(dat.name,modelspec,run.name,".Rdata"))
  rm(predY.loop.pa, predY.loop.ab, predY.loop, 
     predY.mean, predY.median, predY.se, predY.5, predY.95,
     predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95)
}
computational.time = proc.time() - ptm
parallel::stopCluster(cl = c1)


###############################
## If not all areas/cell have run successfully, we can use the code below to only run the unsuccessful ones again

## identifying which cells/regions are already predicted and saved to file:
pred.list <- list.files(paste0("/pvol3TB/4_model_prediction/pred_files",modelspec,"/2km/abundance/"))
pred.list <- pred.list[grep("abundcondpres",pred.list)]
#pred.list <- pred.list[grep("abundance",pred.list)]
# pred.list2 <- pred.list[grep("pa_only",pred.list)]
# pred.list3 <- pred.list[-grep("pa_only",pred.list)]
#pred.list <- pred.list[-grep("fulldat",pred.list)]
# pred.list.numbers <- c(substr(pred.list3,63,67),substr(pred.list2,71,75))
pred.list.numbers <- c(substr(pred.list,75,80))
pred.list.numbers <- as.numeric(sub("^0+", "", pred.list.numbers) )
# pred.list.numbers <- substr(pred.list2,71,75)
# pred.list.numbers <- as.numeric(sub("^0+", "", pred.list.numbers) )

## LAST FILE WITHOUT PA SAVED SEPARATELY IS "..._02331.Rdata"

try.again.v <- which(cell.sel.v%!in%pred.list.numbers)

ptm = proc.time()
#for(j in try.again.v[1:300]){
foreach(j=try.again.v) %dopar% {
  i <- cell.sel.df[j,2]
  k <- cell.sel.df[j,1]
  sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
                      xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
  message(j)
  ##
  XData.grid.loop <- XData.grid[sel.loop,]
  xy.grid.loop <- xy.grid[sel.loop,]

  ## name for saving:
  dat.name <- paste0("/pvol3TB/4_model_prediction/pred_files",modelspec,"/",res,"/abundance/",res,"_model_cells_", as.character(model), 
                     "_abundcondpres_chains_",as.character(nChains),
                     "_thin_", ... = as.character(thin),
                     "_samples_", as.character(samples),
                     "_pred_")
  run.name <- sprintf("%05d",cell.sel.v[j])

  ## setup prediction - pa
  Gradient.pa = prepareGradient(pa, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
  predY.loop.pa <- predict(pa, Gradient=Gradient.pa, expected=TRUE) ## this gives probabilities instead of integer outcomes
  rm(Gradient.pa)
  ## mat.names <- dimnames(predY.loop.pa[[1]])
  ## predY.loop.pa <- array(unlist(predY.loop.pa), c(nrow(xy.grid.loop), ncol(pa$Y), samples*nChains), dimnames(predY.loop.pa[[1]]))
  
  ## setup prediction - abund
  Gradient.ab = prepareGradient(ab, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
  predY.loop.ab <- predict(ab, Gradient=Gradient.ab, expected=TRUE) ## this gives probabilities instead of integer outcomes
  rm(Gradient.ab)  
  # mat.names <- dimnames(predY.loop.ab[[1]])
  # predY.loop.ab <- array(unlist(predY.loop.ab), c(nrow(xy.grid.loop), ncol(ab$Y), samples*nChains), dimnames(predY.loop.ab[[1]]))
  
  # Multiply matrices
  predY.loop = simplify2array(predY.loop.pa) * simplify2array(predY.loop.ab)
  predY.loop.pa.array = simplify2array(predY.loop.pa)
  
  # Get posterior median
  predY.pa.mean = as.data.frame(apply(predY.loop.pa.array, 1:2, mean))
  predY.pa.median = as.data.frame(apply(predY.loop.pa.array, 1:2, median))
  predY.mean = as.data.frame(apply(predY.loop, 1:2, mean))
  predY.median = as.data.frame(apply(predY.loop, 1:2, median))
  
  # Posterior width
  predY.pa.5 = as.data.frame(apply(predY.loop.pa.array, 1:2, p5))
  predY.pa.95 = as.data.frame(apply(predY.loop.pa.array, 1:2, p95))
  predY.5 = as.data.frame(apply(predY.loop, 1:2, p5))
  predY.95 = as.data.frame(apply(predY.loop, 1:2, p95))

  ## save
  save(predY.mean, predY.median, predY.5, predY.95,
       predY.pa.mean, predY.pa.median, predY.pa.5, predY.pa.95,
       sel.loop, XData.grid.loop, xy.grid.loop,
       file=paste0(dat.name,modelspec,run.name,".Rdata"))
  rm(predY.loop.pa, predY.loop.ab, predY.loop, 
     predY.mean, predY.median, predY.5, predY.95,
     predY.pa.mean, predY.pa.median, predY.pa.5, predY.pa.95)
}
computational.time = proc.time() - ptm
parallel::stopCluster(cl = c1)

#############################################################
#############################################################
#############################################################

288
306
Error in { : task 1 failed - "non-numeric argument to binary operator"
  In addition: Warning messages:
    1: executing %dopar% sequentially: no parallel backend registered 
  2: In sprintf("new_unit", 1:(nyNew)) :
    one argument not used by format 'new_unit'
  3: In sprintf("new_unit", 1:(nyNew)) :
    one argument not used by format 'new_unit'
  4: In sprintf("new_unit", 1:(nyNew)) :
    one argument not used by format 'new_unit'
  5: In sprintf("new_unit", 1:(nyNew)) :
    one argument not used by format 'new_unit'
  Session restored from your saved work on 2025-Jan-15 15:13:37 UTC (5 hours ago)
























#############################################################
#############################################################
#############################################################
##### OLD CODE

## 10h on the laptop for 15 species
## parallel processing:
library(doParallel)
library(foreach)
UseCores = parallel::detectCores() - 1
UseCores=8
c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
registerDoParallel(c1)
getDoParWorkers()

## 1500s for 8 cells with 8 cores and 50km cells
## 1640s for 12 cells with 12 cores  and 50km cells

ptm = proc.time()
foreach(i=1:60) %dopar%{ #1:length(xmin)) %dopar%{
  #library(Hmsc)
  for(k in 1:length(ymin)){
    print(paste0("i = ",i,"; y = ",k))
    sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
                        xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
    XData.grid.loop <- XData.grid[sel.loop,]
    xy.grid.loop <- xy.grid[sel.loop,]
    
    ## setup prediction
    Gradient = prepareGradient(m, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
    ## predict
    predY.loop <- predict(m, Gradient=Gradient)
    mat.names <- dimnames(predY.loop[[1]])
    predY.loop <- array(unlist(predY.loop), c(nrow(xy.grid.loop),ncol(m$Y),samples*nChains), dimnames(predY.loop[[1]]))
    predY.mean <- apply(predY.loop, 1:2, mean)
    predY.sd <- apply(predY.loop, 1:2, sd)
    dimnames(predY.mean) <- dimnames(predY.sd) <- mat.names
    
    dat.name <- paste0("/pvol3TB/4_model_prediction/pred_files/",res,"/",
                       res,"_model_", as.character(model), "_",
                       c("pa","abundance")[modeltype],
                       "_chains_",as.character(nChains),
                       "_thin_", ... = as.character(thin),
                       "_samples_", as.character(samples),
                       "_pred_")
    run.name <- paste0("x",i,"_y",k)
    save(predY.loop, file=paste0(dat.name,"fulldat_",run.name,".Rdata"))
    save(predY.mean, predY.sd, sel.loop, XData.grid.loop, xy.grid.loop,
         file=paste0(dat.name,run.name,".Rdata"))
    rm(predY.loop, predY.mean, predY.sd)
  }
}
computational.time = proc.time() - ptm

parallel::stopCluster(cl = c1)

# ptm = proc.time()
# for(i in 1:length(xmin)){
#   message(paste0("x = ",i))
#   for(k in 1:length(ymin)){
#     print(paste0("y = ",k))
#     sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
#                         xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
#     XData.grid.loop <- XData.grid[sel.loop,]
#     xy.grid.loop <- xy.grid[sel.loop,]
# 
#     ## setup prediction
#     Gradient = prepareGradient(m, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
#     ## predict
#     predY.loop <- predict(m, Gradient=Gradient)
#     mat.names <- dimnames(predY.loop[[1]])
#     predY.loop <- array(unlist(predY.loop), c(nrow(xy.grid.loop),ncol(m$Y),samples*nChains), dimnames(predY.loop[[1]]))
#     predY.mean <- apply(predY.loop, 1:2, mean)
#     predY.sd <- apply(predY.loop, 1:2, sd)
#     dimnames(predY.mean) <- dimnames(predY.sd) <- mat.names
#     
#     dat.name <- paste0("/pvol/4_model_prediction/pred_files/","model_", as.character(model), "_",
#                        c("pa","abundance")[modeltype],
#                        "_chains_",as.character(nChains),
#                        "_thin_", ... = as.character(thin),
#                        "_samples_", as.character(samples),
#                        "_pred_")
#     run.name <- paste0("x",i,"_y",k)
#     save(predY.loop, file=paste0(dat.name,"fulldat_",run.name,".Rdata"))
#     save(predY.mean, predY.sd, sel.loop, XData.grid.loop, xy.grid.loop,
#          file=paste0(dat.name,run.name,".Rdata"))
#     rm(predY.loop, predY.mean, predY.sd)
#   }
# }
# computational.time = proc.time() - ptm





















# ## setup prediction
# Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(cellID = xy.grid))
# ## predict
# ptm = proc.time()
# predY <- predict(m, Gradient=Gradient)
# computational.time = proc.time() - ptm
# 
# 
# ## size of prediction boxes
# xmin <- seq(-3000000,2000000, by=1000000)
# xmax <- seq(-2000000,3000000, by=1000000)
# ymin <- seq(-3000000,2000000, by=1000000)
# ymax <- seq(-2000000,3000000, by=1000000)
# ## spatial data
# xy.grid.raw <- coordinates(r2)[which(!is.na(r2[])),]
# ## remove NAs
# sel <- which(!complete.cases(pred_stack.dat))
# XData.grid <- grid[-sel,]
# xy.grid <- xy.grid.raw[-sel,]
# 
# rm(xy.grid.raw, grid, r2, pred_stack.dat)
# 
# ptm = proc.time()
# for(i in 1:length(xmin)){
#   print(i)
#   for(k in 1:length(ymin)){
#     sel2 <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
#                   xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
#     XData.grid.loop <- XData.grid[sel2,]
#     xy.grid.loop <- xy.grid[sel2,]
#     ## setup prediction
#     Gradient = prepareGradient(m, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
#     ## predict
#     predY <- predict(m, Gradient=Gradient)
#    }
# }
# computational.time = proc.time() - ptm





















# 
# 
# ## expected count
# EpredY = apply(abind(predY, along=3, c(1,2), mean))
# ## expected probability of occurrence
# EpredY = apply(abind(predY, along=3, c(1,2), FUN=function(a){mean(a>0)}))
# 



























# 
# 
# thin = 1
# samples = 1000
# nChains = 2
# comp.time = matrix(nrow=2, ncol=3)
# for (modeltype in 1:2){
#   for (model in 1:3){
#     filename = file.path(biodiv.dir, paste("model_",as.character(model),"_",
#                                            c("pa","abundance")[modeltype],
#                                            "_chains_",as.character(nChains),
#                                            "_thin_", as.character(thin),"_samples_",
#                                            as.character(samples),
#                                            ".Rdata",sep = ""))
#     # filename = file.path(biodiv.dir, paste("model_", "pa", "_thin_", ... = as.character(thin),
#     #                                        "_samples_", as.character(samples), ".Rdata", sep = ""))
#     load(filename)
#     comp.time[modeltype,model] = computational.time[1]
#     mpost = convertToCodaObject(m)
#     es.beta = effectiveSize(mpost$Beta)
#     ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
#     es.gamma = effectiveSize(mpost$Gamma)
#     ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
#     # es.rho = effectiveSize(mpost$Rho)
#     # ge.rho = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
#     es.V = effectiveSize(mpost$V)
#     ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
#     if (model==2){
#       es.omega = NA
#       ge.omega = NA
#     } else {
#       es.omega = effectiveSize(mpost$Omega[[1]])
#       ge.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
#     }
#     mixing = list(es.beta=es.beta, ge.beta=ge.beta,
#                   es.gamma=es.gamma, ge.gamma=ge.gamma,
#                   # es.rho=es.rho, ge.rho=ge.rho,
#                   es.V=es.V, ge.V=ge.V,
#                   es.omega=es.omega, ge.omega=ge.omega)
#     filename = file.path(biodiv.dir, paste("mixing_",as.character(model),"_",
#                                            c("pa","abundance")[modeltype],
#                                            "_chains_",as.character(nChains),
#                                            "_thin_", as.character(thin),"_samples_",
#                                            as.character(samples),
#                                            ".Rdata",sep = ""))
#     save(file=filename, mixing)
#   }}
# 
# 
# 
# 
# 
# #setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
# library(Hmsc)
# library(colorspace)
# library(vioplot)
# 
# #include in samples_list and thin_list only those models that you have actually fitted!
# samples_list = 1000 #c(5,250,250,250)
# thin_list = 1 #c(1,1,10,100)
# nst = length(thin_list)
# nChains = 2
# 
# ma = NULL
# na = NULL
# for (Lst in 1:nst) {
#   thin = thin_list[Lst]
#   samples = samples_list[Lst]
#   
#   filename = file.path(biodiv.dir, paste("model_",as.character(model),"_",
#                                          c("pa","abundance")[modeltype],
#                                          "_chains_",as.character(nChains),
#                                          "_thin_", as.character(thin),"_samples_",
#                                          as.character(samples),
#                                          ".Rdata",sep = ""))
#   load(filename)
#   mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
#   psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
#   tmp = summary(psrf.beta)
#   if(is.null(ma)){
#     ma=psrf.beta[,1]
#     na = paste0(as.character(thin),",",as.character(samples))
#   } else {
#     ma = cbind(ma,psrf.beta[,1])
#     if(j==1){
#       na = c(na,paste0(as.character(thin),",",as.character(samples)))
#     } else {
#       na = c(na,"")
#     }
#   }
# }
# 
# pdf(file=paste("MCMC_convergence.pdf"))
# par(mfrow=c(2,1))
# vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
# vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
# dev.off()
# 
# 
# 
# 
# 
# 
# ################################
# ##### evaluating model fit #####
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## script S3 - fitting a model
# #setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
# 
# load(file = "models/unfitted_models") #models, modelnames
# 
# samples_list = c(5,250,250,250,250,250)
# thin_list = c(1,1,10,100,1000,10000)
# nChains = 4
# for(Lst in 1:length(samples_list)){
#   thin = thin_list[Lst]
#   samples = samples_list[Lst]
#   print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
#   nm = length(models)
#   for (model in 1:nm) {
#     print(paste0("model = ",modelnames[model]))
#     m = models[[model]]
#     m = sampleMcmc(m, samples = samples, thin=thin,
#                    adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
#                    transient = ceiling(0.5*samples*thin),
#                    nChains = nChains) 
#     models[[model]] = m
#   }
#   filename = paste("models/models_thin_", as.character(thin),
#                    "_samples_", as.character(samples),
#                    "_chains_",as.character(nChains),
#                    ".Rdata",sep = "")
#   save(models,modelnames,file=filename)
# }
# 
# 
# 







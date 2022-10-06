# ##### Setting up----
# library(raster)
# library(readxl)
# library(readr)
# library(dplyr)
# library(data.table)
# library(proj4)
# library(stringr)
# library(RColorBrewer)
# library(SOmap)

library(Hmsc)
library(terra)

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

biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity")

##############################################################################################################
##############################################################################################################
## get file names of all environmental rasters and bricks and load into one big stack----
#all files with "gri" extension
env_list<-list.files(path = env.derived, pattern="gri$",  full.names=TRUE) 
#subset to  "shelf" files
env_list<-env_list[grep(".500m_shelf_scaled", env_list)]
#for the single rasters layer names are missing. Extract from file name.
env_names<-gsub(".*_|\\..*","",env_list)
#stack all environmental layers and make sure they have appropriate names (currently manual and a bit messy!)
env_stack_scaled<-stack(env_list)
names(env_stack_scaled) <- env_names
names(env_stack_scaled)[10:17]<-c("waom4k_seafloorcurrents_absolute", "waom4k_seafloorcurrents_mean", 
                                  "waom4k_seafloorcurrents_residual", "waom4k_seafloorsalinity", "waom4k_seafloortemperature",
                                  "waom4k_test_flux08","waom4k_test_settle08","waom4k_test_susp08")



## load data (generated in "biodiversity_models_a_prep.R" in folder ARC_Benthic_Mapping/biodiversity)
load(file=paste0(biodiv.dir,"/biodiversity_bio_dat.Rdata"))
load(file=paste0(biodiv.dir,"/biodiversity_env_dat.Rdata"))

# dat.hmsc <- cbind(dat_pa_clean, cell_metadata_env_clean[,c(21,22,43,69:ncol(cell_metadata_env_clean))])

metadat <- cell_metadata_env_clean_scaled

##############################################################################################################

## fitting an hmsc using Otsos course scripts and : https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13345&file=mee313345-sup-0002-AppendixS2.pdf


###########################
##### set up the data #####

## for simplicity, start by analysing only 10 species at 100 sites
# s <- sample(1:nrow(cell_metadata_env_clean_scaled),100)
# s2 <- s[order(s)]
# metadat <- metadat[s2,]
# metadat$cellID <- factor(metadat$cellID)
# Y <- dat_cov_pa[s2,c(1,4,6,7:11,13,15)]  ## species data

## or go with the full dataset here (minus species that are super rare)
# Y <- dat_cov_pa  ## species data
# Y <- dat_cov_pa[,-c(1, 25, 34, 40, 52, 68, 70, 74, 80, 85, 86, 90,100,102,103,105,109,110,114,115,118,120:130)] ## minus species at less than 10 sites
#Y <- dat_cov_pa[,-c(1,24,25,34,40,52,54,62,65,67:74,76,80,85,86,88:90,97,99:115,118:130)] ## minus species at less than 20 sites
## only Bryozoans
dat_cov_species[dat_cov_species==0] <- NA
Y <- dat_cov_species[,4] ##
# Y <- dat_cov_species[,c(4,5,13,16,17,20,21,30,31,74,61,97,104,105,107)] ##


## with NPP:
# XData <- metadat[,c(21:23,34,45,60:67,70:77)] ## environmental data, 78 is geomorph
## without NPP:
XData <- metadat[,c(21:23,34,60:67,70:77)] ## environmental data, 78 is geomorph
############################
##### set up the model #####

## study design - a random effect at the sample level (raster-cell)
studyDesign <- data.frame(cellID=metadat$cellID)#,
                          # imageQuality=metadat$image_quality_score,
                          # transectID=metadat$cover_cells_transect1)

## spatial random effect
xy <- metadat[,4:5]
colnames(xy) = c("x","y")
sRL = xy
rownames(sRL) = metadat$cellID
rL = HmscRandomLevel(sData=sRL)
rL$nfMin = 5
rL$nfMax = 10

# ## simple random effect
# rL = HmscRandomLevel(units=studyDesign$transectID_full)

XFormula = ~ depth + depth2 + slope + tpi + distance2canyons + waom4k_seafloortemperature + waom4k_seafloorcurrents_mean + waom4k_test_settle08#+ NPP_su_mean
## a few NPP values are NA, set to 0 (or the scaled equivalent of 0)
#XData$NPP_su_mean[which(is.na(XData$NPP_su_mean))] <- min(XData$NPP_su_mean, na.rm=T)

mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
             studyDesign = studyDesign, ranLevels = list(cellID=rL)) #ranLevels = list(transectID_full=rL))
mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
            studyDesign = studyDesign)
mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "lognormal poisson",
              studyDesign = studyDesign, ranLevels = list(cellID=rL)) #ranLevels = list(transectID_full=rL))

#######################################
##### run MCMC and save the model #####
## full model: 1000 iterations take about 45minutes
models <- list(mFULL, mENV, mSPACE)
modeltype = 2
model = 1
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 100 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
adaptNf = rep(ceiling(0.4*samples*thin),1)
nChains = 2
set.seed(1)
ptm = proc.time()
for(i in 2){
  print(i)
  models[[i]] <- sampleMcmc(models[[i]], samples = samples, thin = thin,
               adaptNf = adaptNf, transient = transient,
               nChains = nChains, nParallel = nChains,
               initPar = "fixed effects")
}
computational.time = proc.time() - ptm
filename = file.path(biodiv.dir, paste("model_", as.character(model), "_",
                                     c("pa","abundance")[modeltype], 
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples), ".Rdata", sep = ""))
# filename = file.path(biodiv.dir, paste("model_", "pa", "_thin_", ... = as.character(thin),
#                                        "_samples_", as.character(samples), ".Rdata", sep = ""))
save(models, file=filename, computational.time)
m <- models[[2]]

#######################################
##### examine parameter estimates #####

## extracting the posterior distribution from the model object
mpost = convertToCodaObject(m)

## viewing parameter estimates
summary(mpost$Beta)

## assess explanatory power
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY = preds)
MF

#######################################
##### examine MCMC convergence #####

### Graphical overview
plot(mpost$Beta)
### If good, the following is true:
### - different chain yielding the same results
### - chains rise and fall rapidly without apparent autocorrelation
### - the first half looks essentially identical to the second half

### quantitative overview:
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta)$psrf
### If good, the following is true:
### - effective sample size not too far away from sample size
### - potential scale reduction factors close to 1, indicating the multiple chains give consistent results

### if early iteration look different from late ones: discard more transient iterations
### if sampled chains show autocorrelation: increase number of samples or increase thinning interval (better)
### 1000 samples is usually sufficient

## wAIC

#########################################
##### predicting new sampling units #####
Gradient = constructGradient(m, focalVariable = "depth")
predY.1 = predict(m, Gradient = Gradient, expected = TRUE)
Gradient = constructGradient(m, focalVariable = "depth2")
predY.2 = predict(m, Gradient = Gradient, expected = TRUE)
Gradient = constructGradient(m, focalVariable = "slope")
predY.3 = predict(m, Gradient = Gradient, expected = TRUE)
Gradient = constructGradient(m, focalVariable = "tpi")
predY.4 = predict(m, Gradient = Gradient, expected = TRUE)
Gradient = constructGradient(m, focalVariable = "waom4k_seafloorcurrents_mean")
predY.5 = predict(m, Gradient = Gradient, expected = TRUE)
Gradient = constructGradient(m, focalVariable = "waom4k_seafloortemperature")
predY.6 = predict(m, Gradient = Gradient, expected = TRUE)
Gradient = constructGradient(m, focalVariable = "waom4k_test_settle08")
predY.7 = predict(m, Gradient = Gradient, expected = TRUE)

par(mfrow=c(4,4))
for(i in 1:15){
  plotGradient(m, Gradient, pred=predY.1, measure="Y", index=i, showData=TRUE, main="depth")
  }
par(mfrow=c(4,4))
for(i in 1:15){
  plotGradient(m, Gradient, pred=predY.3, measure="Y", index=i, showData=TRUE, main="slope")
}
par(mfrow=c(4,4))
for(i in 1:15){
  plotGradient(m, Gradient, pred=predY.5, measure="Y", index=i, showData=TRUE, main="seafloor currents")
}

par(mfrow=c(3,3))
plotGradient(m, Gradient, pred=predY.1, measure="Y", index=1, showData=TRUE)
plotGradient(m, Gradient, pred=predY.2, measure="Y", index=1, showData=TRUE)
plotGradient(m, Gradient, pred=predY.3, measure="Y", index=1, showData=TRUE)
plotGradient(m, Gradient, pred=predY.4, measure="Y", index=1, showData=TRUE)
plotGradient(m, Gradient, pred=predY.5, measure="Y", index=1, showData=TRUE)
plotGradient(m, Gradient, pred=predY.6, measure="Y", index=1, showData=TRUE)
plotGradient(m, Gradient, pred=predY.7, measure="Y", index=1, showData=TRUE)

#######################################################
##### evaluating model fit using cross validation #####
partition = createPartition(m, nfolds=2) ## use column to partition according to different hierarchies (e.g. leave an entire region out)
preds.cv = computePredictedValues(m, partition = partition)
MF.cv = evaluateModelFit(hM=m, predY = preds.cv)
MF.cv

####################################################################
load(paste0(biodiv.dir,"/model_1_abundance_chains_2_thin_10_samples_100.Rdata"))

load(paste0(biodiv.dir,"/biodiversity_pred_stack_scaled_dat.Rdata")) #(~5GB large)

# ###################### MAP THE HMSC !!!
# library(raster)
# r2 <- raster(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
# ## mapping the predictions
# # READING THE GRID DATA
# # grid = read.csv(file.path(dataDir, "grid_10000.csv"))
# # grid = grid[!(grid$Habitat=="Ma"),]
# 
# # depth+depth2+slope+tpi+distance2canyons+
# # waom4k_seafloortemperature+waom4k_seafloorcurrents_mean+waom4k_test_settle08
# ## env data
# grid <- pred_stack.dat[,c(1:3,6,7,17,14,11)]
# ## spatial data
# xy.grid.raw <- coordinates(r2)[which(!is.na(r2[])),] #~5GB large
# ## remove NAs
# sel <- which(!complete.cases(pred_stack.dat))
# XData.grid <- grid[-sel,]
# xy.grid <- xy.grid.raw[-sel,]
# rm(xy.grid.raw, grid, r2, pred_stack.dat)
# 
# ## focus on small areas
# ext <- c(-3000000,-2000000,2000000,3000000)
# ext.name <- "ext1"
# sel2 <- which(xy.grid[,1]>ext[1] & xy.grid[,1]<ext[2] & xy.grid[,2]>ext[3] & xy.grid[,2]<ext[4])
# XData.grid <- XData.grid[sel2,]
# xy.grid <- xy.grid[sel2,]
# 
# ## setup prediction
# Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(cellID = xy.grid))
# ## predict
# predY <- predict(m, Gradient=Gradient)
# 
# ## expected count
# EpredY = apply(abind(predY, along=3, c(1,2), mean))
# ## expected probability of occurrence
# EpredY = apply(abind(predY, along=3, c(1,2), FUN=function(a){mean(a>0)}))

## prepare data for prediction and mapping:
r2 <- rast(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
# plot(r2)
# abline(v=-3000000)
# abline(v=-2500000)
# abline(v=-2000000)
# abline(h=3000000)
# abline(h=2500000)
# abline(h=2000000)
# abline(h=1500000)
# abline(h=1000000)
# abline(h=500000)
# abline(h=0)

grid <- pred_stack.dat[,c(1:3,6,7,17,14,11)]
## spatial data
xy.grid.raw <- coordinates(r2)[which(!is.na(r2[])),]
## remove NAs
sel <- which(complete.cases(pred_stack.dat))
XData.grid <- grid[sel,]
xy.grid <- xy.grid.raw[sel,]

## first find which cells we have ignored before
sel.not.na <- which(!is.na(r2[]))
## create an empty raster to fill for mapping
empty.ra <- rast(r2)
empty.ra[] <- NA

## save data
save(sel, sel.not.na, file="biodiversity/hmsc_model_cell_sel.Rdata")
save(XData.grid, xy.grid, file="biodiversity/hmsc_model_cell_grid.Rdata")
rm(xy.grid.raw, grid, r2, pred_stack.dat)



#############################################################
## load data
load("biodiversity/model_1_pa_chains_2_thin_10_samples_100.Rdata")
m <- models[[2]]
load("biodiversity/hmsc_model_cell_sel.Rdata")
load("biodiversity/hmsc_model_cell_grid.Rdata")

r2 <- rast(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))
empty.ra <- rast(r2)
empty.ra[] <- NA


## size of prediction boxes
xmin <- seq(-3000000,2500000, by=500000)
xmax <- seq(-2500000,3000000, by=500000)
ymin <- seq(-3000000,2500000, by=500000)
ymax <- seq(-2500000,3000000, by=500000)

## 10h on the laptop for 15 species
ptm = proc.time()
for(i in 1:length(xmin)){
  message(paste0("x = ",i))
  for(k in 1:length(ymin)){
    print(paste0("y = ",k))
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

    dat.name <- paste0("biodiversity/pred_files/","model_", as.character(model), "_",
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




##########################################################################

##### puzzle the predictions back together, one species at a time
## loop over species and boxes individually to fill the raster and save
ptm = proc.time()
for(s in 1:ncol(m$Y)){
  sp.name <- dimnames(m$Y)[[2]][s]
  message(sp.name)
  hmsc.pred.ra <- empty.ra
  ## loop over x-axis
  for(i in 1:length(xmin)){
    message(paste0("x = ",i))
    ## loop over y-axis
    for(k in 1:length(ymin)){
      print(paste0("y = ",k))
      ## find datafile to load
      dat.name <- paste0("biodiversity/pred_files/","model_", as.character(model), "_",
                         c("pa","abundance")[modeltype], 
                         "_chains_",as.character(nChains),
                         "_thin_", ... = as.character(thin),
                         "_samples_", as.character(samples),
                         "_pred_")
      run.name <- paste0("x",i,"_y",k)
      load(file=paste0(dat.name,run.name,".Rdata"))
      
      ## which raster cells to fill with data
      sel.fill <- sel.not.na[sel[sel.loop]]
      ## fill data
      hmsc.pred.ra[sel.fill] <- predY.mean[,s]
    }
  }
  ## save raster file before moving on
  writeRaster(hmsc.pred.ra, filename=paste0(biodiv.dir,"/hmsc_pred_abund_",sp.name,".tif"), overwrite=TRUE)
}
computational.time = proc.time() - ptm

test <- rast(paste0(biodiv.dir,"/hmsc_pred_",sp.name,".tif"))

tifs <- list.files(biodiv.dir, pattern=".tif")
abund_sel <- grep("abund", tifs)
tifs_abund <- tifs[abund_sel]

maps_abund <- rast(paste0(biodiv.dir,"/", tifs_abund))
maps_pa <- rast(paste0(biodiv.dir,"/", tifs)[-abund_sel])

maps_abund_pa_BHBr <- maps_pa[[6]] * maps_abund[[2]]
maps_abund_pa_BHLettuce <- maps_pa[[9]] * maps_abund[[3]]
maps_abund_pa_UBS_B <- maps_pa[[15]] * maps_abund[[4]]

xlim <- c(-2800000,-2000000)
ylim <- c(1000000,2500000)
plot(maps_abund_pa_UBS_B, xlim=xlim, ylim=ylim, range=c(0,1), main="UBS_B")
plot(maps_abund_pa_BHBr, xlim=xlim, ylim=ylim, range=c(0,1), main="BH_Br")
plot(maps_abund_pa_BHLettuce, xlim=xlim, ylim=ylim, range=c(0,1), main="BH_Lettuce")



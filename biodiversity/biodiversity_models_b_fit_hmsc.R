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
library(Hmsc)

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
Y <- dat_cov_pa[,c(4,5,13,16,17,20,21,30,31,74,61,97,104,105,107)] ##
#Y <- dat_cov_pa[s2,c(4,5,13,16,17,20,21,30,31,74,61,97,104,105,107)] ##


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

mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
             studyDesign = studyDesign, ranLevels = list(cellID=rL)) #ranLevels = list(transectID_full=rL))
mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
            studyDesign = studyDesign)
mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "probit",
              studyDesign = studyDesign, ranLevels = list(cellID=rL)) #ranLevels = list(transectID_full=rL))

#######################################
##### run MCMC and save the model #####
## space only is very slow: 5min for 2 samples, 9min for 4 samples (-> 10k iteration will take a week!)
models <- list(mFULL, mENV, mSPACE)
modeltype = 1
model = 1
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 100 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
adaptNf = rep(ceiling(0.4*samples*thin),1)
nChains = 2
set.seed(1)
ptm = proc.time()
for(i in 1){
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
load(paste0(biodiv.dir,"/model_1_pa_chains_2_thin_10_samples_1000.Rdata"))

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
library(terra)
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
for(i in 3:length(xmin)){
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
for(s in 4){#ncol(m$Y)){
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
  writeRaster(hmsc.pred.ra, filename=paste0(biodiv.dir,"/hmsc_pred_",sp.name,".tif"), overwrite=TRUE)
}
computational.time = proc.time() - ptm

test <- rast(paste0(biodiv.dir,"/hmsc_pred_",sp.name,".tif"))


## map the predictions, e.g. South Orkneys at x=2, y=11
load(paste0(biodiv.dir,"/pred_files/model_1_pa_chains_2_thin_10_samples_100_pred_fulldat_x2_y11.Rdata"))
r2 <- raster(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))

sel.not.na <- which(!is.na(r2[]))
empty.ra <- rast(r2)
empty.ra[] <- NA
hmsc.pred.ra <- rast(list(empty.ra,empty.ra))



sel.fill <- sel.not.na[sel[sel.loop]]
hmsc.pred.ra[[1]][sel.fill] <- predY.mean[,1]
hmsc.pred.ra[[2]][sel.fill] <- predY.mean[,2]
for(i in 10:15){
  print(i)
hmsc.pred.ra.loop <- empty.ra
hmsc.pred.ra.loop[sel.fill] <- predY.mean[,i]
hmsc.pred.ra <- c(hmsc.pred.ra,hmsc.pred.ra.loop)
}
names(hmsc.pred.ra) <- dimnames(predY.mean)[[2]]
hmsc.pred.ra.crop <- crop(hmsc.pred.ra, ext(xmin[2],ymin[11],xmax[2],ymax[11]))

hmsc.pred.ra <- c(hmsc.pred.ra,empty.ra)
hmsc.pred.ra[[i]][sel.fill] <- predY.mean[,i]

hmsc.pred.ra$bryo <- empty.ra
hmsc.pred.ra$ 










# DEFINING THE NEW STUDY DESIGN
nyNew = nrow(grid)
StudyDesignNew = matrix(NA,nyNew,1, dimnames=list(NULL,names(m$studyDesign)))
StudyDesignNew[,1] = sprintf('new_Route_%.3d',1:nyNew)
StudyDesignNew[,2] = sprintf('new_Year_%.3d',1:nyNew)
StudyDesignNew = as.data.frame(StudyDesignNew)
StudyDesignAll=rbind(m$studyDesign,StudyDesignNew)
# DEFINING RANDOM EFFECTS THAT INCLUDE BOTH OLD (USED FOR MODEL FITTING)
# AND NEW (THE GRID DATA) UNITS
rL1 = m$ranLevels[[1]]
xyold = rL1$s
xy = grid[,1:2]
rownames(xy) = StudyDesignNew[,1]
colnames(xy) = colnames(xyold)
xyall = rbind(xyold,xy)
rL1$pi = StudyDesignAll[,1]
rL1$s = xyall
if(FALSE){
  predYR1 = predict(m, studyDesign=StudyDesignNew, XData=grid, ranLevels=list(Route=rL1),
                    expected=TRUE, predictEtaMean=TRUE)
  predYR = apply(abind(predYR1,along=3),c(1,2),mean)
  save(predYR,file="panels/predictions/predYR_thin_100_samples_1000_grid_10000.Rdata")
} else {
  load("panels/predictions/predYR_thin_100_samples_1000_grid_10000.Rdata")
  ta = 1:dim(predYR)[1]
  predYR = predYR[ta,]
  xy = grid[,1:2]
  xy = xy[ta,]
}
# COMPUTE SPECIES RICHNESS (S), COMMUNITY WEIGHTED MEANS (predT),
# REGIONS OF COMMON PROFILE (RCP)
S=rowSums(predYR)
predT = (predYR%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)
RCP = kmeans(predYR, 7)
RCP$cluster = as.factor(RCP$cluster)
# EXTRACT THE OCCURRENCE PROBABILITIES OF ONE EXAMPLE SPECIES
pred_Cm = predYR[,50]
# MAKE A DATAFRAME OF THE DATA TO BE PLOTTED
mapData=data.frame(xy,S,predT,pred_Cm,RCP$cluster)


sp <- ggplot(data = mapData, aes(x=x, y=y, color=pred_Cm))+geom_point(size=1)
sp + ggtitle("Predicted Corvus monedula occurrence") +
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  scale_color_gradient(low="blue", high="red", name ="Occurrence probability")
sp <- ggplot(data = mapData, aes(x=x, y=y, color=S))+geom_point(size=1)
sp + ggtitle("Predicted species richness") +
  xlab("East coordinate (km)") + ylab("North coordinate (km)") +
  scale_color_gradient(low="blue", high="red", name ="Species richness")
####################################################################





test <- predict(m)


computeAUC = function(Y, predY){
  ns = dim(Y)[2]
  AUC = rep(NaN,ns)
  ## take care that Y has only levels {0,1} as specified in auc() below
  Y <- ifelse(Y > 0, 1, 0)
  for (i in 1:ns){
    sel = !is.na(Y[,i])
    if(length(unique(Y[sel,i]))==2)
      AUC[i] = pROC::auc(Y[sel,i],predY[sel,i], levels=c(0,1),direction="<")
  }
  return(AUC)
}

computeTjurR2 = function(Y, predY) {
  ns = dim(Y)[2]
  R2 = rep(NaN, ns)
  for (i in 1:ns) {
    R2[i] = mean(predY[which(Y[, i] == 1), i]) - mean(predY[which(Y[,i] == 0), i])
  }
  return(R2)
}

test <- computeAUC(Y, preds.cv)








thin = 1
samples = 1000
nChains = 2
comp.time = matrix(nrow=2, ncol=3)
for (modeltype in 1:2){
  for (model in 1:3){
    filename = file.path(biodiv.dir, paste("model_",as.character(model),"_",
                                         c("pa","abundance")[modeltype],
                                         "_chains_",as.character(nChains),
                                         "_thin_", as.character(thin),"_samples_",
                                         as.character(samples),
                                         ".Rdata",sep = ""))
    # filename = file.path(biodiv.dir, paste("model_", "pa", "_thin_", ... = as.character(thin),
    #                                        "_samples_", as.character(samples), ".Rdata", sep = ""))
    load(filename)
    comp.time[modeltype,model] = computational.time[1]
    mpost = convertToCodaObject(m)
    es.beta = effectiveSize(mpost$Beta)
    ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    es.gamma = effectiveSize(mpost$Gamma)
    ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
    # es.rho = effectiveSize(mpost$Rho)
    # ge.rho = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
    es.V = effectiveSize(mpost$V)
    ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
    if (model==2){
      es.omega = NA
      ge.omega = NA
    } else {
      es.omega = effectiveSize(mpost$Omega[[1]])
      ge.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
    }
    mixing = list(es.beta=es.beta, ge.beta=ge.beta,
                  es.gamma=es.gamma, ge.gamma=ge.gamma,
                  # es.rho=es.rho, ge.rho=ge.rho,
                  es.V=es.V, ge.V=ge.V,
                  es.omega=es.omega, ge.omega=ge.omega)
    filename = file.path(biodiv.dir, paste("mixing_",as.character(model),"_",
                                          c("pa","abundance")[modeltype],
                                          "_chains_",as.character(nChains),
                                          "_thin_", as.character(thin),"_samples_",
                                          as.character(samples),
                                          ".Rdata",sep = ""))
    save(file=filename, mixing)
  }}





#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)
library(colorspace)
library(vioplot)

#include in samples_list and thin_list only those models that you have actually fitted!
samples_list = 1000 #c(5,250,250,250)
thin_list = 1 #c(1,1,10,100)
nst = length(thin_list)
nChains = 2

ma = NULL
na = NULL
for (Lst in 1:nst) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  filename = file.path(biodiv.dir, paste("model_",as.character(model),"_",
                                         c("pa","abundance")[modeltype],
                                         "_chains_",as.character(nChains),
                                         "_thin_", as.character(thin),"_samples_",
                                         as.character(samples),
                                         ".Rdata",sep = ""))
  load(filename)
  mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    if(is.null(ma)){
      ma=psrf.beta[,1]
      na = paste0(as.character(thin),",",as.character(samples))
    } else {
      ma = cbind(ma,psrf.beta[,1])
      if(j==1){
        na = c(na,paste0(as.character(thin),",",as.character(samples)))
      } else {
        na = c(na,"")
      }
    }
}

pdf(file=paste("MCMC_convergence.pdf"))
par(mfrow=c(2,1))
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
dev.off()






################################
##### evaluating model fit #####






















## script S3 - fitting a model
#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are

load(file = "models/unfitted_models") #models, modelnames

samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nChains = 4
for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  nm = length(models)
  for (model in 1:nm) {
    print(paste0("model = ",modelnames[model]))
    m = models[[model]]
    m = sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains) 
    models[[model]] = m
  }
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  save(models,modelnames,file=filename)
}










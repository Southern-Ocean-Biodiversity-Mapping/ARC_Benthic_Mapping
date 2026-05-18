# ##### Setting up----
#library(raster)
library(terra)
library(Hmsc)
# library(readr)
library(dplyr)
# library(data.table)
# library(proj4)
# library(stringr)
# library(RColorBrewer)
# library(SOmap)

'%!in%' <- function(x,y)!('%in%'(x,y))

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
distr.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/species_distributions")

##############################################################################################################

#res <- "500m"
res <- "2km"

##############################################################################################################
## load scaled environmental rasters:
env_stack_scaled <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled.tif"))

## load data
load(file=paste0(ARC_Data.dir,"Cell_level_env_",res,".Rdata"))
load(file=paste0(ARC_Data.dir,"Cell_level_bio_2pc_",res,".Rdata"))
cover_mod <- cover_mod.2km
count_mod <- count_mod.2km

# dat.hmsc <- cbind(dat_pa_clean, cell_metadata_env_clean[,c(21,22,43,69:ncol(cell_metadata_env_clean))])

r <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r$depth

#### check for NAs in the data:
## find NAs in waom data
waom.na.sel <- which(is.na(cell_metadata_env_scaled$seafloorcurrents_absolute))
npp.na.sel <- which(is.na(cell_metadata_env_scaled$npp_mean))
na.sel <- c(waom.na.sel,npp.na.sel)
metadat <- cell_metadata_env_scaled[-na.sel,]
metadat$cellID <- factor(metadat$cellID)
cover_cells <- cover_mod.2km[-na.sel,]

## presence absence data:
cov_pa.raw <- cover_cells[,-1]
cov_pa.raw[cov_pa.raw>0] <- 1

## remove rare species
cov_pa.raw2 <- cov_pa.raw[,-which(colSums(cov_pa.raw)<18)]

## combine UBS_B with Bryozoan_Hard_Branching_Antler
cov_pa.raw2$Bryozoan_Hard_Branching_Antler <- cov_pa.raw2$Bryozoan_Hard_Branching_Antler+cov_pa.raw2$UBS_B
cov_pa.raw2$Bryozoan_Hard_Branching_Antler[cov_pa.raw2$Bryozoan_Hard_Branching_Antler>1] <- 1
cov_pa.raw3 <- cov_pa.raw2[,-grep("UBS_B",names(cov_pa.raw2))]

## remove substrates, noid and unscorable
if(res=="2km") cov_pa <- cov_pa.raw3[,-c(grep("Sub",names(cov_pa.raw3)),
                                         grep("Unsco",names(cov_pa.raw3)),
                                         grep("NoID",names(cov_pa.raw3))[1:2])]

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

## or only PS81 data
# s <- which(metadat$counts_cells_survey=="PS81")
# Y <- cov_pa[s,]
# metadat <- metadat[s,]
# metadat$cellID <- factor(metadat$cellID)
# metadat$cover_cells_transect1 <- factor(metadat$cover_cells_transect1)

## or go with the full dataset here (minus species that are super rare)
## COMPARE THIS LIST TO THE EXCEL FILE WITH THE 2% CUTOFF!
Y <- cov_pa#[,-which((colSums(dat_cov_pa)/nrow(dat_cov_pa))<0.018)]
#Y <- Y[,1:15]
## species data
# Y <- dat_cov_pa[,-c(1, 25, 34, 40, 52, 68, 70, 74, 80, 85, 86, 90,100,102,103,105,109,110,114,115,118,120:130)] ## minus species at less than 10 sites
#Y <- dat_cov_pa[,-c(1,24,25,34,40,52,54,62,65,67:74,76,80,85,86,88:90,97,99:115,118:130)] ## minus species at less than 20 sites
## only Bryozoans
# Y <- c(grep("Bryo",names(dat_cov_pa)),grep("UBS_B",names(dat_cov_pa)))
#Y <- dat_cov_pa[s2,c(4,5,13,16,17,20,21,30,31,74,61,97,104,105,107)] ##


## with NPP:
# XData <- metadat[,c(21:23,34,45,60:67,70:77)] ## environmental data, 78 is geomorph
## without NPP:
#XData <- metadat[,c(23:36,38:72)] ## environmental data, 37 is geomorph

## XData only the variables we choose in XFormula below
model_vars <- c("depth","depth2","logslope","tpi","distance2canyons","distance2canyons2","seafloortemperature","seafloorcurrents_mean","seafloorcurrents_sd","npp_mean")
#XData <- dplyr::select(metadat, model_vars)
XData <- metadat[,which(names(metadat)%in%model_vars)]

############################
##### set up the model #####

## study design - a random effect at the sample level (raster-cell)
studyDesign <- data.frame(cellID=metadat$cellID, surveyID=metadat$cover_cells_survey, transectID=metadat$cover_cells_transect1,
                          gear=metadat$gear, year=as.factor(metadat$year))

## survey effect:
rL.s = HmscRandomLevel(units=levels(metadat$cover_cells_survey))
## transect effect:
rL.t = HmscRandomLevel(units=levels(metadat$cover_cells_transect1))
## gear effect:
rL.g = HmscRandomLevel(units=levels(metadat$gear))
## year effect:
rL.y = HmscRandomLevel(units=levels(as.factor(metadat$year)))

## spatial random effect
xy <- metadat[,4:5]
colnames(xy) = c("x","y")
sRL = xy
rownames(sRL) = metadat$cellID

# rL = HmscRandomLevel(sData=sRL)
# rL$nfMin = 5
# rL$nfMax = 10

## first specifying knots on a grid, keeping them to a minimum


########## knots at 200km distance, 250km min distance:
## add points between AP and Ross Sea
xy.knots <- rbind(xy,
                  c(-2143647,498436),
                  c(-1916289,355285),
                  c(-2000000,100000),
                  c(-1983655,-368890),
                  c(-1739456,-638350),
                  c(-1621567,-975176),
                  c(-1360527,-1160431),
                  c(-900000,-1300000),
                  c(-627930,-1345685))
## add points to East Antarctica
xy.knots <- rbind(xy.knots,
                  c(710952,-2154067),
                  c(1115144,-2288798),
                  c(2100359,-1724614),
                  c(2529812,-1017280),
                  c(2757170,-629930),
                  c(2824535,-318366),
                  c(2672963,-57326),
                  c(2656122,254237),
                  c(2487709,498436),
                  c(2386661,742635),
                  c(2268772,1256295),
                  c(2125621,1685748),
                  c(1645644,1812057),
                  c(820421,2056256),
                  c(1200000,2000000))
## add a point to Weddell Sea
xy.knots <- rbind(xy.knots,
                  c(-1503678,1054199),
                  c(-1436312,1300000),
                  c(-1958393,1298398))
Knots = constructKnots(xy.knots, knotDist = 200000, minKnotDist = 250000)
#Knots = constructKnots(xy, knotDist = 50000, minKnotDist = 2000000)
# plot(r2)
# points(xy.knots[,1],xy.knots[,2],pch=18)
# points(Knots[,1],Knots[,2],col='red',pch=18)
rL = HmscRandomLevel(sData=sRL, sMethod='GPP', sKnot=Knots)
rL$nfMax=10
rL.s$nfMax=10
rL.t$nfMax=10
rL.g$nfMax=10
rL.y$nfMax=10
#rL = setPriors(rL,nfMin=1,nfMax=1)

# ## simple random effect
# rL = HmscRandomLevel(units=studyDesign$transectID_full)

###
XFormula = ~depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+seafloortemperature+seafloorcurrents_mean+seafloorcurrents_sd+npp_mean

mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
             studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
            studyDesign = studyDesign)
mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "probit",
              studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))

#######################################
##### run MCMC and save the model #####
## 28h for full model with 10k iterations on 74 species, 455 samples, 9 covariates, 1 spatial variable using GPP
## space only is very slow: 5min for 2 samples, 9min for 4 samples (-> 10k iteration will take a week!)
models <- list(mFULL, mENV, mSPACE)
modeltype = 1
model = 1
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 500 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
adaptNf = rep(ceiling(0.4*samples*thin),1)
nChains = 2
set.seed(1)
ptm = proc.time()
for(i in 1){
  print(i)
  models[[i]] <- sampleMcmc(models[[i]], samples = samples, thin = thin,
               adaptNf = rep(adaptNf,length(models[[i]]$rLNames)), transient = transient,
               nChains = nChains, nParallel = nChains,initPar = "fixed effects")
}
computational.time = proc.time() - ptm
filename = file.path(biodiv.dir, paste(res,"_model_", as.character(model), "_",
                                     c("pa","abundance")[modeltype], 
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples), ".Rdata", sep = ""))
save(models, file=filename, computational.time)
m <- models[[1]]

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
# load(paste0(biodiv.dir,"/",res,"_model_1_pa_chains_4_thin_1_samples_1000.Rdata"))
load(paste0(biodiv.dir,"/",res,"_model_1_pa_chains_4_thin_10_samples_1000.Rdata"))

load(paste0(env.derived,"/Circumpolar_EnvData_",res,"_shelf_mask_scaled_dataframe.Rdata")) #(~5GB large)

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
r <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
# plot(r$depth)
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

grid <- pred_stack.dat[,c(1,48,50,3,14,44,41,46)]#[,c(1:3,6,7,17,14,11)]
## spatial data
xy.grid.raw <- crds(r$depth)
## remove NAs
sel <- which(complete.cases(grid))
XData.grid <- grid[sel,]
xy.grid <- xy.grid.raw[sel,]

## first find which cells we have ignored before
sel.not.na <- which(!is.na(r$depth[]))
## create an empty raster to fill for mapping
empty.ra <- rast(r$depth)
empty.ra[] <- NA

## save data
save(sel, sel.not.na, file=paste0("biodiversity/hmsc_",res,"_model_cell_sel.Rdata"))
save(XData.grid, xy.grid, file=paste0("biodiversity/hmsc_",res,"_model_cell_grid.Rdata"))
rm(xy.grid.raw, grid, r, pred_stack.dat)



#############################################################
## load data
load(paste0("biodiversity/",res,"_model_1_pa_chains_4_thin_10_samples_1000.Rdata"))
m <- models[[1]]
load(paste0(biodiv.dir,"/hmsc_",res,"_model_cell_sel.Rdata"))
load(paste0(biodiv.dir,"/hmsc_",res,"_model_cell_grid.Rdata"))

r <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
empty.ra <- rast(r$depth)
empty.ra[] <- NA


## size of prediction boxes
xmin <- seq(-3000000,2500000, by=500000)
xmax <- seq(-2500000,3000000, by=500000)
ymin <- seq(-3000000,2500000, by=500000)
ymax <- seq(-2500000,3000000, by=500000)

## size of prediction boxes
xmin <- seq(-3000000,2750000, by=250000)
xmax <- seq(-2750000,3000000, by=250000)
ymin <- seq(-3000000,2750000, by=250000)
ymax <- seq(-2750000,3000000, by=250000)

## size of prediction boxes
xmin <- seq(-3000000,2900000, by=100000)
xmax <- seq(-2900000,3000000, by=100000)
ymin <- seq(-3000000,2900000, by=100000)
ymax <- seq(-2900000,3000000, by=100000)

## size of prediction boxes
xmin <- seq(-3000000,2950000, by=50000)
xmax <- seq(-2950000,3000000, by=50000)
ymin <- seq(-3000000,2950000, by=50000)
ymax <- seq(-2950000,3000000, by=50000)

plot(r$depth)
for(i in 2:length(xmin)){
  abline(h=ymin[i])
  abline(v=xmin[i])
}

## we can reduce the runtime by 12h (for 14400 cells) if we skip over the empty cells
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
# save(cells_with_data, file=paste0("biodiversity/",res,"_model_50km_cells_with_data.Rdata"))
load(file=paste0("biodiversity/",res,"_model_50km_cells_with_data.Rdata"))

##################

## 10h on the laptop for 15 species
ptm = proc.time()
for(i in 30){#3:length(xmin)){
  message(paste0("x = ",i))
  for(k in 14){#1:length(ymin)){
    print(paste0("y = ",k))
    print(proc.time())
    sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
                        xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
    XData.grid.loop <- XData.grid[sel.loop,]
    xy.grid.loop <- xy.grid[sel.loop,]
    print(paste0(nrow(xy.grid.loop)," cells to predict on"))
    ## setup prediction
    Gradient = prepareGradient(m, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
    ## predict
    predY.loop <- predict(m, Gradient=Gradient)
    mat.names <- dimnames(predY.loop[[1]])
    predY.loop <- array(unlist(predY.loop), c(nrow(xy.grid.loop),ncol(m$Y),samples*nChains), dimnames(predY.loop[[1]]))

    predY.mean <- apply(predY.loop, 1:2, mean)
    predY.sd <- apply(predY.loop, 1:2, sd)
    dimnames(predY.mean) <- dimnames(predY.sd) <- mat.names

    dat.name <- paste0("biodiversity/pred_files/",res,"_model_", as.character(model), "_",
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






## parallel processing: PER CELL that contains values
library(doParallel)
library(foreach)
parallel::detectCores()
UseCores = parallel::detectCores() - 1
c1<-makeCluster(UseCores, outfile="")
registerDoParallel(c1)
getDoParWorkers()

cell.sel.v <- which(!is.na(cells_with_data))
cell.sel.df <- which(!is.na(cells_with_data), arr.ind = TRUE)

iterations <- 100

ptm = proc.time()
foreach(j=1:iterations) %dopar%{ #3:length(xmin)
  library(Hmsc)
  i <- cell.sel.df[j,2]
  k <- cell.sel.df[j,1]
  sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
                      xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
  print(i)
  ##
  XData.grid.loop <- XData.grid[sel.loop,]
  xy.grid.loop <- xy.grid[sel.loop,]
  ## setup prediction
  Gradient = prepareGradient(m, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
  ## predict
  predY.loop <- predict(m, Gradient=Gradient)
  mat.names <- dimnames(predY.loop[[1]])
  predY.loop <- array(unlist(predY.loop), c(nrow(xy.grid.loop),ncol(m$Y),samples*nChains), dimnames(predY.loop[[1]]))
  ## derived values
  predY.mean <- apply(predY.loop, 1:2, mean)
  predY.sd <- apply(predY.loop, 1:2, sd)
  dimnames(predY.mean) <- dimnames(predY.sd) <- mat.names
  ## save
  dat.name <- paste0("biodiversity/pred_files/",res,"/",res,"_model_", as.character(model), "_",
                     c("pa","abundance")[modeltype],
                     "_chains_",as.character(nChains),
                     "_thin_", ... = as.character(thin),
                     "_samples_", as.character(samples),
                     "_pred_")
  run.name <- sprintf("%05d",cell.sel.v[j])
  save(predY.loop, file=paste0(dat.name,"fulldat_",run.name,".Rdata"))
  save(predY.mean, predY.sd, sel.loop, XData.grid.loop, xy.grid.loop,
       file=paste0(dat.name,run.name,".Rdata"))
  rm(predY.loop, predY.mean, predY.sd)
}
computational.time = proc.time() - ptm
parallel::stopCluster(cl = c1)

pred.files <- list.files(paste0("biodiversity/pred_files/",res,"/"))
pred.files.longindices <- substr(pred.files, nchar(pred.files[1])-10, nchar(pred.files[1])-6)
pred.files.indices <- as.numeric(sub("^0+", "", pred.files.longindices))

cell.sel.v[which(cell.sel.v%!in%pred.files.indices)]


# ## parallel processing: PER COLUMN
# library(doParallel)
# library(foreach)
# parallel::detectCores()
# UseCores = parallel::detectCores() - 1
# c1<-makeCluster(UseCores, outfile="")
# registerDoParallel(c1)
# getDoParWorkers()
# 
# col.sel <- which(colSums(cells_with_data, na.rm=T)>0)
# iterations <- 8
# 
# ptm = proc.time()
# foreach(j=1:iterations) %dopar%{ #3:length(xmin)
#   i <- col.sel[j]
#   library(Hmsc)
#   print(i)
#   ys <- which(!is.na(cells_with_data[,i]))
#   for(k in ys){#1:length(ymin)){
#     sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
#                         xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
#     if(length(sel.loop)==0){
#       print(paste0("skipping: x = ",i,"; y = ",k)) 
#     }
#     #print(paste0("x = ",i,"; y = ",k))
#     #print(proc.time())
#     XData.grid.loop <- XData.grid[sel.loop,]
#     xy.grid.loop <- xy.grid[sel.loop,]
#     #print(paste0(nrow(xy.grid.loop)," cells to predict on"))
#     ## setup prediction
#     Gradient = prepareGradient(m, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
#     ## predict
#     predY.loop <- predict(m, Gradient=Gradient)
#     mat.names <- dimnames(predY.loop[[1]])
#     predY.loop <- array(unlist(predY.loop), c(nrow(xy.grid.loop),ncol(m$Y),samples*nChains), dimnames(predY.loop[[1]]))
# 
#     predY.mean <- apply(predY.loop, 1:2, mean)
#     predY.sd <- apply(predY.loop, 1:2, sd)
#     dimnames(predY.mean) <- dimnames(predY.sd) <- mat.names
# 
#     dat.name <- paste0("biodiversity/pred_files/",res,"/",res,"_model_", as.character(model), "_",
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
# 
# parallel::stopCluster(cl = c1)

























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










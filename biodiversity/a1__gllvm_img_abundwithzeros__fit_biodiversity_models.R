## fitting an hmsc using Otsos book, course scripts and : https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13345&file=mee313345-sup-0002-AppendixS2.pdf

##############################################################################################################

library(gllvm)
library(corrplot)
library(terra)
'%!in%' <- function(x,y)!('%in%'(x,y))



!!! GLLVM package has a bug that hasnt been fixed (July 2023), it cant predict when a row-effect is specified !!!

## from https://github.com/JenniNiku/gllvm/issues/86
data(antTraits)
y <- as.matrix(antTraits$abund)
X <- scale(antTraits$env[, 1:3])
fake_site <- data.frame(site = rep(1:3, each = 10))

# Fit gllvm model
fit <- gllvm(y = y, X, 
             family = poisson(), 
             studyDesign = fake_site, 
             row.eff = ~ (1|site))

# predict with new data
xnew <- cbind(rnorm(10), rnorm(10), rnorm(10))
colnames(xnew) <- colnames(X)
predfit <- predict(fit, newX = xnew, type = "response", level = 0)
  
## Error in object$TMBfn$env$data$dr0 %*% object$params$row.params : 
##  non-conformable arguments 

!!! Maybe glmmTMB could work: https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html


##############################################################################################################

#res <- "500m"
res <- "2km"

##############################################################################################################
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

## load scaled environmental rasters:
env_stack_scaled <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled.tif"))

## load data
load(file=paste0(ARC_Data.dir,"Cell_level_env_",res,".Rdata"))
load(file=paste0(ARC_Data.dir,"Image_level_bio.Rdata"))

## check for NAs in the data:
## find NAs in waom and npp data
waom.na.sel <- cell_metadata_env_scaled$cellID[which(is.na(cell_metadata_env_scaled$seafloorcurrents_absolute))]
npp.na.sel <- cell_metadata_env_scaled$cellID[which(is.na(cell_metadata_env_scaled$npp_mean))]
## remove seamount transects (tan1802 & tan1901)
seamount.transects <- c("TAN1802_160","TAN1802_170","TAN1802_179","TAN1802_180","TAN1802_184","TAN1802_185","TAN1802_191","TAN1802_193",
                        "TAN1802_195","TAN1802_196","TAN1802_197","TAN1802_207","TAN1802_208","TAN1802_209","TAN1802_213","tan1901_209")
seamount.transects.sel <- cell_metadata_env_scaled$cellID[which(cell_metadata_env_scaled$cover_cells_transect1%in%seamount.transects)]
## combine to one vector
na.sel.chr <- unique(c(waom.na.sel,npp.na.sel,seamount.transects.sel))
na.sel <- which(img.metadata$cellID_2km%in%na.sel.chr)

## remove nas from data
cover_imgs <- cover_mod[-na.sel,]
metadat <- img.metadata[-na.sel,]
metadat$cellID_2km <- factor(metadat$cellID_2km)

## presence absence data:
cov_pa.raw <- cover_imgs[,-1]
cov_pa.raw[cov_pa.raw>0] <- 1

## abundance
cov_ab.raw <- cover_imgs[,-1]

## remove rare species
sel.rm <- c(which(colSums(cov_pa.raw)<=9))
cov_ab.raw2 <- cov_ab.raw[,-sel.rm]

## combine UBS_B with Bryozoan_Hard_Branching_Antler
cov_ab.raw2$'Bryozoa - Hard - Branching - Morphotype 1 - Antler' <- cov_ab.raw2$'Bryozoa - Hard - Branching - Morphotype 1 - Antler' + cov_ab.raw2$'Unidentified Biological Matrix - Bryozoan associated'
cov_ab.raw3 <- cov_ab.raw2[,-grep("Unidentified Biological Matrix - Bryozoan associated",names(cov_ab.raw2))]

## remove substrates, noid and unscorable
substrate <- c("Biologenic Rubble", "Boulders", "Cobbles", "Pebble / Gravel", "Rock", "Sand / Mud")

if(res=="2km") cov_ab <- cov_ab.raw3[,-c(which(names(cov_ab.raw3)%in%substrate),
                                         grep("Unsco",names(cov_ab.raw3)),
                                         grep("Unidentifiable",names(cov_ab.raw3)))]

#cov_offset <- metadat$cover_points_scorable

## subsample to fewer sites and species
sel.subset <- c(1:71, seq(72, nrow(cov_ab), by=5))
cov_ab <- cov_ab[sel.subset,1:20]
metadat <- metadat[sel.subset,]
# sel.subset.subset <- which(rowSums(cov_ab)!=0)
# cov_ab <- cov_ab[sel.subset.subset,]
# metadat <- metadat[sel.subset.subset,]
message("USING A SUBSET ONLY")

###########################
##### set up the data #####
Y <- cov_ab
## XData only the variables we choose in XFormula below
model_vars <- c("depth","depth2","logslope","tpi","distance2canyons","distance2canyons2","seafloortemperature","seafloorcurrents_mean","npp_mean","seafloorsalinity")
#XData <- dplyr::select(metadat, model_vars)
XData.cells <- cell_metadata_env_scaled[match(metadat$cellID_2km, cell_metadata_env_scaled$cellID),]
XData <- XData.cells[,which(names(XData.cells)%in%model_vars)]
XData$transectID <- factor(metadat$transectID_full)

XFormula = ~depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+seafloortemperature+seafloorcurrents_mean+seafloorsalinity+npp_mean

metadat$Filename.standardised <- factor(metadat$Filename.standardised)
metadat$survey <- factor(metadat$survey)
metadat$transectID_full <- factor(metadat$transectID_full)

## study design
studyDesign <- data.frame(transectID=metadat$transectID_full)
#studyDesign <- data.frame(filename=metadat$Filename.standardised, surveyID=metadat$survey, transectID=metadat$transectID_full)
#studyDesign <- data.frame(filename=metadat$Filename.standardised, surveyID=metadat$survey, transectID=metadat$transectID_full, gear=metadat$gear)
# XData$surveyID=metadat$survey

# ## fit zero-inflated poisson, AIC 12898.89 for 20 species subset
# fit.gllvm.zip.env_only <- gllvm(y=Y, X=XData, formula=XFormula, family="ZIP")
# ## neg.binomial is better: , AIC 12553.66 for 20 species subset
# fit.gllvm.nb.env_only <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial")
# 
# ## adding random effects on survey doesn't improve AIC: 12559.51
# fit.gllvm.nb.re_on_survey <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|surveyID))
# ## adding random effects on site doesn't improve AIC: 12555.17
# fit.gllvm.nb.re_on_site <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|filename))
# ## adding random effects on gear doesn't improve AIC: 12552.39
# fit.gllvm.nb.re_on_gear <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|gear))
# ## adding random effects on transect improves AIC: 12538.53
# fit.gllvm.nb.re_on_transect <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|transectID))
# 
# ## adding 1LV improves AIC: AIC=12500.53
# fit.gllvm.nb.re_on_transect_plus_1lv <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|transectID), num.lv = 1)
# ## adding a second LV makes things worse: AIC=12538.53
# fit.gllvm.nb.re_on_transect_plus_2lv <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|transectID), num.lv = 2)

## 1 LV is better than none, but we can't use that for prediction
#fit.gllvm.nb.re_on_transect_plus_1lv <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|transectID), num.lv = 1, sd.errors = TRUE) #

fit.gllvm.nb.re_on_transect <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", row.eff=~(1|transectID), sd.errors=TRUE, num.lv=0)
fit <- fit.gllvm.nb.re_on_transect

## full model takes around 2h to fit
fit.gllvm.nb.re_on_transect <- gllvm(y=Y, X=XData, formula=XFormula, family="negative.binomial", studyDesign=studyDesign, row.eff=~(1|transectID), sd.errors=TRUE)
fit <- fit.gllvm.nb.re_on_transect
save(fit, file="biodiversity/image_model_gllvm_re_on_transectID.Rdata")

load("biodiversity/image_model_gllvm_re_on_transectID.Rdata")

par(mfrow=c(3,2))
plot(fit)

coefplot(fit, which.Xcoef = 1:2)
coefplot(fit, which.Xcoef = 3:4)
coefplot(fit, which.Xcoef = 5:6)
coefplot(fit, which.Xcoef = 7:8)
coefplot(fit, which.Xcoef = 9:10)

## Correlation matrix for model with predictors shows correlation patterns between species when the effect of the predictors are taken into account
crx <- getResidualCor(fit)
par(mfrow=c(1,1))
corrplot(crx[1:40,1:40], diag = FALSE, type = "lower", method = "square", tl.srt = 25)


























############################
##### gllvm #####

## study design - a random effect at the sample level (raster-cell)
# studyDesign <- data.frame(cellID=metadat$cellID)#,
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

### using the standard algorithm:
## 5min per iteration to fit the spatial model, which is way too long (~1 month for 10k iterations)
#rL = HmscRandomLevel(sData=sRL)

### trying NNGP:
## doesn't work, the error message is: "Failed updaters and their counts in chain 1  ( 15  attempts)"
# rL = HmscRandomLevel(sData=sRL, sMethod="NNPG")
# rL = setPriors(rL,nfMin=1,nfMax=1)

### trying GPP:
## ~ 40s for 10 iterations; 60s for 50 iterations, 84s for 100 iterations; using knots at 500km distance -> ~1h20min for 10k iterations
## ~ 3.5min for 10 iterations; 5min for 50 iterations; 6.6min for 100 iterations; using knots at 250km distance -> ~2h40min for 10k iterations
## ~ 7min for 2 iterations; 15min for 100 iterations using knots at 200km distance -> ~13h20min for 10k iterations
## BUT, on a 250km grid, with the full dataset, the predictions will take 41 days!!!
## first specifying knots on a grid
# xy.knots <- rbind(xy,c(2900000,0)) ## add a point to the right to allow mapping of East Antarctica
# Knots = constructKnots(xy.knots, knotDist = 250000, minKnotDist = 2500000)
# #Knots = constructKnots(xy, knotDist = 50000, minKnotDist = 2000000)
# plot(xy.knots[,1],xy.knots[,2],pch=18, asp=1)
# points(Knots[,1],Knots[,2],col='red',pch=18)

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
rL = HmscRandomLevel(sData=sRL, sMethod='GPP', sKnot=Knots)
rL$nfMax=10
rL.s$nfMax=10
rL.t$nfMax=10
rL.g$nfMax=10
rL.y$nfMax=10


# mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
#              studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
# mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
#             studyDesign = studyDesign)
# mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "lognormal poisson",
#               studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
# mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
#              studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
# mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
#             studyDesign = studyDesign)
# mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "probit",
#               studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
             studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
            studyDesign = studyDesign)
mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "lognormal poisson",
              studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))

#setPriors(mFULL, offset = Normal(1, 1e-6)))
#######################################
##### run MCMC and save the model #####
#######################################
## 5 days on the VM to fit all three 2km models with 8000 iterations and 75 species

models <- list(mFULL, mENV, mSPACE)
modeltype = 2
model = 3
thin = 1  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 50 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
adaptNf = rep(ceiling(0.4*samples*thin),1)
nChains = 4
filename.string <- paste(res,"_model_", as.character(model), "_",
                         c("pa","abundance")[modeltype], 
                         "_chains_",as.character(nChains),
                         "_thin_", ... = as.character(thin),
                         "_samples_", as.character(samples), sep = "")

set.seed(2)
ptm = proc.time()
for(i in 1){
  print(i)
  print(proc.time())
  models[[i]] <- sampleMcmc(models[[i]], samples = samples, thin = thin, transient = transient,
                            nChains = nChains, nParallel = nChains,
                            #initPar = "fixed effects",
                            updater = list(GammaEta = FALSE))
}
computational.time = proc.time() - ptm
filename = file.path(paste(filename.string, ".Rdata", sep = ""))
save(models, file=filename, computational.time)

################################
##### evaluating model fit #####
################################
MF <- list()
for(i in 1:3){
  preds = computePredictedValues(models[[i]])
  MF[[i]] = evaluateModelFit(hM=models[[i]], predY = preds)
}
MF
filename2 = file.path(paste(filename.string, "_MF.Rdata", sep = ""))
save(MF, file=filename2)


####################################################### TAKES X-TIMES (FOLDS) LONGER THAN THE MODEL FITTING!!!
##### evaluating model fit using cross validation #####
partition = createPartition(models[[1]], nfolds=2) ## use column to partition according to different hierarchies (e.g. leave an entire region out)
# preds.cv = computePredictedValues(m, partition = partition)
# MF.cv = evaluateModelFit(hM=m, predY = preds.cv)
# MF.cv

## This takes a long time, 5 days for the full model & 1h for the environment only model
ptm = proc.time()
MF.cv = list()
for(i in 1){
  preds = computePredictedValues(models[[i]], partition=partition)
  MF.cv[[i]] = evaluateModelFit(hM=models[[i]], predY = preds)
}
computational.time = proc.time() - ptm
filename2 = file.path(paste(res,"_model_", as.character(model), "_",
                            c("pa","abundance")[modeltype], 
                            "_chains_",as.character(nChains),
                            "_thin_", ... = as.character(thin),
                            "_samples_", as.character(samples), "_2foldcv.Rdata", sep = ""))
save(MF.cv, file=filename2, computational.time)
# for(i in 2){
#   preds = computePredictedValues(models[[i]], partition=partition)
#   MF.cv[[i]] = evaluateModelFit(hM=models[[i]], predY = preds)
# }
# save(MF.cv, file=filename2, computational.time)















# ##### Setting up----
library(terra)
library(PerformanceAnalytics) ## plotting correlations
library(MASS)    ## glm.nb
library(glmmTMB) ## model fitting for glmm
library(buildmer)## stepwise model selection for glmmTMB
library(DHARMa)  ## model diagnostics https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(mgcv)    ## gams

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

###############################################################################

#res <- "500m"
res <- "2km"

###############################################################################
## load scaled environmental rasters:
env_stack_scaled <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled.tif"))

## load data
load(file=paste0(ARC_Data.dir,"Cell_level_env_",res,".Rdata"))
load(file=paste0(ARC_Data.dir,"Cell_level_bio_2pc_",res,".Rdata"))

r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth

## check for NAs in the data:
## find NAs in waom and npp data
waom.na.sel <- which(is.na(cell_metadata_env_scaled$seafloorcurrents_absolute))
npp.na.sel <- which(is.na(cell_metadata_env_scaled$npp_mean))
## remove seamount transects (tan1802 & tan1901)
seamount.transects <- c("TAN1802_160","TAN1802_170","TAN1802_179","TAN1802_180","TAN1802_184","TAN1802_185","TAN1802_191","TAN1802_193",
                        "TAN1802_195","TAN1802_196","TAN1802_197","TAN1802_207","TAN1802_208","TAN1802_209","TAN1802_213","tan1901_209")
seamount.transects.sel <- which(cell_metadata_env_scaled$cover_cells_transect1%in%seamount.transects)
## combine to one vector
na.sel <- unique(c(waom.na.sel,npp.na.sel,seamount.transects.sel))

## remove nas from data
cover_cells <- cover_mod.2km[-na.sel,]
metadat <- cell_metadata_env_scaled[-na.sel,]
metadat$cellID <- factor(metadat$cellID)
cover_groups <- cover_groupings[-na.sel,]
rm(cover_mod.2km, cover_groupings, cell_metadata_env_scaled)

## presence absence data:
cov_pa.raw <- cover_cells[,-1]
cov_pa.raw[cov_pa.raw>0] <- 1

## abundance conditional on presence
cov_ab.raw <- cover_cells[,-1]
cov_ab.raw[cov_pa.raw==0] <- NA

## remove rare species
cov_pa.raw2 <- cov_pa.raw[,-which(colSums(cov_pa.raw)<=9)]
cov_ab.raw2 <- cov_ab.raw[,-which(colSums(cov_pa.raw)<=9)]

## combine UBS_B with Bryozoan_Hard_Branching_Antler
cov_pa.raw2$Bryozoan_Hard_Branching_Antler <- cov_pa.raw2$Bryozoan_Hard_Branching_Antler+cov_pa.raw2$UBS_B
cov_pa.raw3 <- cov_pa.raw2[,-grep("UBS_B",names(cov_pa.raw2))]
cov_ab.raw2$Bryozoan_Hard_Branching_Antler <- cov_ab.raw2$Bryozoan_Hard_Branching_Antler+cov_ab.raw2$UBS_B
cov_ab.raw3 <- cov_ab.raw2[,-grep("UBS_B",names(cov_ab.raw2))]

## remove substrates, noid and unscorable
if(res=="2km") cov_ab <- cov_ab.raw3[,-c(grep("Sub",names(cov_ab.raw3)),
                                         grep("Unsco",names(cov_ab.raw3)),
                                         grep("NoID",names(cov_ab.raw3))[1:2])]
if(res=="2km") cov_pa <- cov_pa.raw3[,-c(grep("Sub",names(cov_pa.raw3)),
                                         grep("Unsco",names(cov_pa.raw3)),
                                         grep("NoID",names(cov_pa.raw3))[1:2])]
cov_pa[cov_pa>0] <- 1

## divide by number of scorable points per cell
cov_ab_perc <- (cov_ab/metadat$cover_points_scorable)*100

cov_ab_perc_scaledlog <- scale(log(cov_ab_perc))
hist(cov_ab_perc_scaledlog)
# par(mfrow=c(4,4))
# for(i in 1:ncol(cov_ab_perc_scaledlog)){
#   hist(cov_ab_perc_scaledlog[,i], main=colnames(cov_ab_perc_scaledlog)[i])
# }



###############################################################################
## set up dataframe
dat.glm <- cbind(cover_groups, metadat)
dat.glm$image_quality_score <- dat.glm$image_quality_score-min(dat.glm$image_quality_score)+1
dat.glm$annotated_area <- dat.glm$cover_area * (dat.glm$cover_points_scorable/dat.glm$cover_points_total)

covars <- "depth + depth2 + logslope + tpi + distance2canyons + distance2canyons2 + seafloortemperature + seafloorcurrents_mean + npp_mean + seafloorsalinity"
dat.for.fit <- dat.glm[,-83]#[c(4,32,78,80,34,45,79,77,73,57,76,82,83)]

## we need to remove PS06 from CV, because it cannot be predicted from most folds when only two transects exist
dat.for.cv <- dat.for.fit[-which(dat.for.fit$cover_cells_survey=="PS06"),]

###############################################################################
## partition for 5-fold cross validation (same as hmsc models on VM):
partition <- c(3,2,3,4,2,3,5,5,1,3,1,1,4,4,4,4,2,1,1,3,5,4,4,4,4,5,3,3,4,4,5,4,3,4,4,3,4,2,1,1,1,4,5,5,2,4,1,3,1,5,5,2,4,3,4,1,1,3,3,1,2,5,1,
5,2,3,1,3,4,2,3,4,1,1,1,3,3,3,3,2,2,3,3,3,2,2,2,2,5,5,5,5,2,2,2,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,1,1,2,2,2,
1,1,1,1,2,2,2,2,4,4,4,4,4,4,1,1,1,1,4,4,4,5,5,5,2,2,2,2,1,1,1,1,3,3,3,1,2,3,3,4,4,3,3,3,1,1,1,2,1,2,2,5,4,2,5,5,5,3,5,4,2,2,2,
2,2,4,1,5,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,2,2,2,5,5,1,3,3,4,4,4,2,2,5,5,3,4,2,1,1,5,5,3,5,1,1,3,3,2,2,2,2,2,4,4,2,2,5,3,4,4,1,4,
4,4,1,5,3,2,2,2,4,3,3,1,5,5,4,3,3,4,1,1,2,2,1,1,3,3,5,4,4,5,5,5,3,3,3,1,4,4,3,1,2,2,5,5,5,3,3,3,3,2,3,2,5,5,5,2,4,4,5,5,1,2,3,
4,1,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,3,3,5,5,5,5,1,2,4,1,2,4,1,1,5,4,3,4,3,4,3,2,3,1,3,5,5,4,5,2,2,1,1,2,5,2,1,2,4,5,4,2,3,1,4,
5,2,5,4,2,2,5,4,2,4,2,1,5,3,5,2,1,5,2,4,3,1,1,5,4,3,1,4,5,5,2,5,1,2,2,4,3,3,3,1,5,5,5,2,1,1,5,5,5,3,4,4,2,2,4,4,1,1,3,3,4,3,2,
4,3,4,5,5,3,4,4,1,1,3,3,3)
## again, remove PS06
partition.cv <- partition[-which(dat.for.fit$cover_cells_survey=="PS06")]

## partition for 10-fold cross validation (same as hmsc models on VM):
#set.seed(2)
#partition10 <- 
## again, remove PS06
partition.cv <- partition[-which(dat.for.fit$cover_cells_survey=="PS06")]


###############################################################################
##### modelling richness
#### GLMs
summary(fit.glm.r <- stepAIC(glm.nb(formula(paste0("richness~",covars,"+log(cover_points_scorable)")), data=dat.for.fit)))
summary(fit.glm.r.o <- stepAIC(glm.nb(formula(paste0("richness~",covars,"+offset(log(cover_points_scorable))")), data=dat.for.fit)))

## diagnostics
fit.sim <- simulateResiduals(fit.glm.r)
plot(fit.sim)
fit.sim.o <- simulateResiduals(fit.glm.r.o)
plot(fit.sim.o)

#### GLMMs
# glmm.formula.r <- paste0("richness~",covars," + (1|cover_cells_survey) +offset(log(cover_points_scorable))") # + (1|gear) + (1|cover_cells_transect1) 
# summary(fit.glmm.r <- glmmTMB(formula(glmm.formula.r), data=dat.for.fit, family=nbinom1))
# summary(fit.glmm.r <- update(fit.glmm.r, ".~.-distance2canyons2"))
# summary(fit.glmm.r <- update(fit.glmm.r, ".~.-distance2canyons"))
#glmm.formula.r <- paste0("richness~",covars,"+log(cover_points_scorable) + (1|cover_cells_survey)") # + (1|gear) + (1|cover_cells_transect1)
glmm.formula.r <- paste0("richness~",covars,"+(1|cover_cells_survey)+offset(log(cover_points_scorable))") # + (1|gear) + (1|cover_cells_transect1)
summary(fit.glmm.r <- glmmTMB(formula(glmm.formula.r), data=dat.for.fit, family=nbinom1))
summary(fit.glmm.r <- update(fit.glmm.r, ".~.-distance2canyons2"))
summary(fit.glmm.r <- update(fit.glmm.r, ".~.-distance2canyons"))

## diagnostics
fit.sim <- simulateResiduals(fit.glmm.r)
plot(fit.sim)
par(mfrow=c(3,3))
plotResiduals(fit.sim, dat.for.fit$depth)#[-c(70,71)])
plotResiduals(fit.sim, dat.for.fit$logslope)#[-c(70,71)])
plotResiduals(fit.sim, dat.for.fit$tpi11)#[-c(70,71)])

#### calculate CV
preds.glm.r.cv <- list()
preds.glmm.r.cv <- list()
for(i in 1:5){
  message(i)
  dat.glm.train <- dat.for.cv[partition.cv!=i,]
  dat.glm.test <- dat.for.cv[partition.cv==i,]
  ## glm
  fit.glm.r.cv <- glm.nb(fit.glm.r.o$call$formula, data=dat.glm.train)
  preds.glm.r.cv[[i]] <-  predict(fit.glm.r.cv, dat.glm.test, type="response")
  print("glm done")
  ## glmm
  fit.glmm.r.cv <- glmmTMB(fit.glmm.r$call$formula, data=dat.glm.train, family=nbinom1)
  preds.glmm.r.cv[[i]] <-  predict(fit.glmm.r.cv, dat.glm.test, type="response")
}

#### plot cv
## glm
plot(fitted(fit.glm.r), fit.glm.r$y, col="red", xlab="predicted", ylab="observed", xlim=c(0,150), ylim=c(0,150))
abline(0,1,col="grey", lty=2)
for(i in 1:5){
  points(preds.glm.r.cv[[i]], dat.for.cv$richness[partition.cv==i])
}
legend("topleft", col=c("black","red"), legend=c("CV", "full"), pch=1)
## glmm
plot(fitted(fit.glmm.r), fit.glm.r$y, col="red", xlab="predicted", ylab="observed", xlim=c(0,50), ylim=c(0,50))
abline(0,1,col="grey", lty=2)
for(i in 1:5){
  points(preds.glmm.r.cv[[i]], dat.for.cv$richness[partition.cv==i])
}
legend("topleft", col=c("black","red"), legend=c("CV", "full"), pch=1)


###############################################################################
##### modelling total cover
glm.formula.c <-  paste0("cover_all~",covars,"+offset(log(cover_points_scorable))")
glmm.formula.c <- paste0("cover_all~",covars,"+offset(log(cover_points_scorable)) + (1|cover_cells_survey)") # + (1|gear) + (1|cover_cells_transect1)
glmm.formula.ct <- paste0("cover_all~",covars,"+offset(log(cover_points_scorable)) + (1|cover_cells_survey) + (1|cover_cells_transect1)")

#### glm
summary(fit.glm.c <- stepAIC(glm.nb(formula(glm.formula.c), data=dat.for.fit)))

#### glmm
summary(fit.glmm.c <- glmmTMB(formula(glmm.formula.c), data=dat.for.fit, family=nbinom1))
summary(fit.glmm.ct <- glmmTMB(formula(glmm.formula.ct), data=dat.for.fit, family=nbinom1))
summary(fit.glmm.ct <- update(fit.glmm.ct, ".~.-distance2canyons2"))
summary(fit.glmm.ct <- update(fit.glmm.ct, ".~.-seafloorsalinity"))
summary(fit.glmm.ct <- update(fit.glmm.ct, ".~.-npp_mean"))
## diagnostics
fit.sim <- simulateResiduals(fit.glmm.c)
plot(fit.sim)

#### calculate CV
preds.glm.c.cv <- list()
preds.glmm.c.cv <- list()
preds.glmm.ct.cv <- list()
for(i in 1:5){
  message(i)
  dat.glm.train <- dat.for.cv[partition.cv!=i,]
  dat.glm.test <- dat.for.cv[partition.cv==i,]
  ## glm
  fit.glm.c.cv <- glm.nb(fit.glm.c$call$formula, data=dat.glm.train)
  preds.glm.c.cv[[i]] <- predict(fit.glm.c.cv, dat.glm.test, type="response")
  print("glm done")
  ## glmm
  fit.glmm.c.cv <- glmmTMB(fit.glmm.c$call$formula, data=dat.glm.train, family=nbinom1)
  preds.glmm.c.cv[[i]] <- predict(fit.glmm.c.cv, dat.glm.test, type="response")
  ## glmm with transect random effect
  fit.glmm.ct.cv <- glmmTMB(fit.glmm.ct$call$formula, data=dat.glm.train, family=nbinom1)
  preds.glmm.ct.cv[[i]] <- predict(fit.glmm.ct.cv, dat.glm.test, type="response", allow.new.levels=TRUE)
}

#### plot cv
## glm
plot(fitted(fit.glm.c), fit.glm.c$y, col="red", xlab="predicted", ylab="observed", xlim=c(0,1700), ylim=c(0,1700), main="total cover GLM")
abline(0,1,col="grey", lty=2)
for(i in 1:5){
  points(preds.glm.c.cv[[i]], dat.for.cv$cover_all[partition.cv==i])
}
legend("bottomright", col=c("black","red"), legend=c("CV", "full"), pch=1)
## glmm
plot(fitted(fit.glmm.c), fit.glm.c$y, col="red", xlab="predicted", ylab="observed", xlim=c(0,1700), ylim=c(0,1700), main="total cover GLMM")
abline(0,1,col="grey", lty=2)
for(i in 1:5){
  points(preds.glmm.c.cv[[i]], dat.for.cv$cover_all[partition.cv==i])
}
legend("bottomright", col=c("black","red"), legend=c("CV", "full"), pch=1)
## glmm transects
plot(fitted(fit.glmm.ct), fit.glm.c$y, col="red", xlab="predicted", ylab="observed", xlim=c(0,1700), ylim=c(0,1700), main="total cover GLMM")
abline(0,1,col="grey", lty=2)
for(i in 1:5){
  points(preds.glmm.ct.cv[[i]], dat.for.cv$cover_all[partition.cv==i])
}
legend("bottomright", col=c("black","red"), legend=c("CV", "full"), pch=1)

## glmm transects vs without
plot(preds.glmm.c.cv[[1]], dat.for.cv$cover_all[partition.cv==1], col="red", xlab="predicted", ylab="observed", xlim=c(0,1700), ylim=c(0,1700), main="GLMM transects vs without")
abline(0,1,col="grey", lty=2)
for(i in 2:5){
  points(preds.glmm.c.cv[[i]], dat.for.cv$cover_all[partition.cv==i], col="red")
}
for(i in 1:5){
  points(preds.glmm.ct.cv[[i]], dat.for.cv$cover_all[partition.cv==i])
}
legend("bottomright", col=c("black","red"), legend=c("transects", "no transects"), pch=1)

###############################################################################
##### modelling SF cover
glm.formula.sf <- paste0("cover_SF~",covars,"+offset(log(cover_points_scorable))")
glmm.formula.sf <- paste0("cover_SF~",covars,"+offset(log(cover_points_scorable)) + (1|cover_cells_survey)") # + (1|gear) + (1|cover_cells_transect1)

#### glm
summary(fit.glm.sf <- glm.nb(formula(glm.formula.sf), data=dat.for.fit))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-tpi))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-seafloorsalinity))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-npp_mean))

#### glmm
summary(fit.glmm.sf <- glmmTMB(formula(glmm.formula.sf), data=dat.for.fit, family=nbinom1))
summary(fit.glmm.sf <- update(fit.glmm.sf, ".~.-distance2canyons2"))
summary(fit.glmm.sf <- update(fit.glmm.sf, ".~.-tpi"))
## diagnostics
fit.sim <- simulateResiduals(fit.glmm.sf)
plot(fit.sim)

#### calculate CV
preds.glm.sf.cv <- list()
preds.glmm.sf.cv <- list()
for(i in 1:5){
  message(i)
  dat.glm.train <- dat.for.cv[partition.cv!=i,]
  dat.glm.test <- dat.for.cv[partition.cv==i,]
  ## glm
  fit.glm.sf.cv <- glm.nb(fit.glm.sf$call$formula, data=dat.glm.train)
  preds.glm.sf.cv[[i]] <- predict(fit.glm.sf.cv, dat.glm.test, type="response")
  print("glm done")
  ## glmm
  fit.glmm.sf.cv <- glmmTMB(fit.glmm.sf$call$formula, data=dat.glm.train, family=nbinom1)
  preds.glmm.sf.cv[[i]] <- predict(fit.glmm.sf.cv, dat.glm.test, type="response")
}

#### plot cv
## glm
plot(fitted(fit.glm.sf), fit.glm.sf$y, col="red", xlab="predicted", ylab="observed", xlim=c(0,1700), ylim=c(0,1700), main="sf cover GLM")
abline(0,1,col="grey", lty=2)
for(i in 1:5){
  points(preds.glm.sf.cv[[i]], dat.for.cv$cover_SF[partition.cv==i])
}
legend("bottomright", col=c("black","red"), legend=c("CV", "full"), pch=1)
## glmm
plot(fitted(fit.glmm.sf), fit.glm.sf$y, col="red", xlab="predicted", ylab="observed", xlim=c(0,1700), ylim=c(0,1700), main="sf cover GLMM")
abline(0,1,col="grey", lty=2)
for(i in 1:5){
  points(preds.glmm.sf.cv[[i]], dat.for.cv$cover_SF[partition.cv==i])
}
legend("bottomright", col=c("black","red"), legend=c("CV", "full"), pch=1)


###############################################################################
##### modelling individual species - abundance
dat.full <- cbind(cov_ab,dat.glm[,c(32,78,80,34,45,79,77,73,57,76,82,19)])
dat.full.cv <- dat.full[-which(dat.for.fit$cover_cells_survey=="PS06"),]
glm.formula.sp <- paste0(names(cov_ab),"~",covars,"+offset(log(cover_points_scorable))")
glmm.formula.sp <- paste0(names(cov_ab),"~",covars,"+offset(log(cover_points_scorable)) + (1|cover_cells_survey)") # + (1|gear) + (1|cover_cells_transect1)

sp.not.na.sel.list <- list()
for(i in 1:83){
  sp.not.na.sel.list[[i]] <- which(!is.na(dat.full[,i]))
}

fit.glm.list <- list()
for(i in c(1:12,14:15,17:76,78:80,82)){
  print(i)
  summary(fit.glm.list[[i]] <- glm.nb(formula(glm.formula.sp[i]), data=dat.full))
}

fit.glmm.list <- list()
for(i in 1:83){
  print(i)
  summary(fit.glmm.list[[i]] <- glmmTMB(formula(glmm.formula.sp[i]), data=dat.full, family=nbinom1))
}

## CV (PROBLEMS WITH THE INDIVIDUAL GLMMS!!!)
preds.glm.ab.cv <- data.frame(matrix(NA, ncol=83, nrow=nrow(dat.full.cv)))
preds.glmm.ab.cv <- data.frame(matrix(NA, ncol=83, nrow=nrow(dat.full.cv)))
for(i in 1:5){
  message(i)
  dat.glm.train <- dat.full.cv[partition.cv!=i,]
  dat.glm.test <- dat.full.cv[partition.cv==i,]
  ## below is because of an error in predict.glmmTMB
  dat.glm.test$cover_cells_survey <- as.character(dat.glm.test$cover_cells_survey)
  ##
  ## glm
  for(j in 1:83){#c(1:12,14:15,17:20,22:58,61:70,73:74,76,78,80,82)){
    tryCatch({
    print(j)
    fit.glm.ab.cv.loop <- glm.nb(formula(glm.formula.sp[j]), data=dat.glm.train)
    preds.glm.ab.cv[partition.cv==i,j] <- predict(fit.glm.ab.cv.loop, dat.glm.test, type="response")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  print("glm iteration done")
  ## glmm
  for(k in 1:83){
    tryCatch({
    print(k)
    fit.glmm.ab.cv.loop <- glmmTMB(formula(glmm.formula.sp[k]), data=dat.glm.train, family=nbinom1)
    preds.glmm.ab.cv[partition.cv==i,k] <- predict(fit.glmm.ab.cv.loop, newdata=dat.glm.test, type="response")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  print("glmm iteration done")
}
## GLM "ERROR : missing value where TRUE/FALSE needed" in:
## i=1: 21,59,60,71,72,77,79,81,83
## i=2: 42,60,79:81,83
## i=3: 21,80,81,83
## i=4: 57,60,80,81,83
## i=5: 21,34,42,60,80,81
## GLM "ERROR : NA/NaN/Inf in 'x'" in:
## i=1: 75; i=2: 77; i=3: 13; i=4: 13,78
## GLM "ERROR : no valid set of coefficients has been found: please supply starting values" in:
## i=2: 13,16; i=3: 46; i=4: 16,77; i=5: 16
## GLMM errors in: (i=2,k=36), (i=4,k=40), (i=4,k=83)


###############################################################################
##### modelling individual species - p/a
glm.formula.sp.pa <- paste0(names(cov_pa),"~",covars,"+cover_points_scorable")
glmm.formula.sp.pa <- paste0(names(cov_pa),"~",covars,"+cover_points_scorable + (1|cover_cells_survey)") # + (1|gear) + (1|cover_cells_transect1)

dat.pa <- cbind(cov_pa,dat.glm)
dat.pa.cv <- dat.pa[-which(dat.for.fit$cover_cells_survey=="PS06"),]

fit.pa.glm.list <- list()
for(i in 1:83){
  print(i)
  summary(fit.pa.glm.list[[i]] <- glm(formula(glm.formula.sp.pa[i]), data=dat.pa, family="binomial"))
}

fit.pa.glmm.list <- list()
for(i in 1:83){
  print(i)
  summary(fit.pa.glmm.list[[i]] <- glmmTMB(formula(glmm.formula.sp.pa[i]), data=dat.pa, family="binomial"))
}

# ## CV
# preds.glm.pa.cv <- list()
# preds.glmm.pa.cv <- list()
# for(i in 1:5){
#   message(i)
#   dat.glm.train <- dat.pa.cv[partition.cv!=i,]
#   dat.glm.test <- dat.pa.cv[partition.cv==i,]
#   preds.glm.pa.cv[[i]] <- matrix(NA,nrow=nrow(dat.pa.cv), ncol=83)
#   preds.glmm.pa.cv[[i]] <- matrix(NA,nrow=nrow(dat.pa.cv), ncol=83)
#   ## glm
#   for(j in c(1:71,73:83)){
#     print(j)
#     fit.glm.pa.cv.loop <- glm.nb(formula(glm.formula.sp.pa[j]), data=dat.glm.train)
#     preds.glm.pa.cv[[i]][partition.cv==i,j] <- predict(fit.glm.pa.cv.loop, dat.glm.test, type="response")
#   }
#   print("glm iteration done")
#   ## glmm
#   for(k in c(1:71,73:83)){
#     print(k)
#     fit.glmm.pa.cv.loop <- glmmTMB(formula(glmm.formula.sp.pa[k]), data=dat.glm.train, family="binomial")
#     preds.glmm.pa.cv[[i]][partition.cv==i,k] <- predict(fit.glmm.pa.cv.loop, dat.glm.test, type="response")
#   }
#   print("glmm iteration done")
# }

## CV
preds.glm.pa.cv <- matrix(NA,nrow=nrow(dat.pa.cv), ncol=83)
preds.glmm.pa.cv <- matrix(NA,nrow=nrow(dat.pa.cv), ncol=83)
for(i in 1:5){
  message(i)
  dat.glm.train <- dat.pa.cv[partition.cv!=i,]
  dat.glm.test <- dat.pa.cv[partition.cv==i,]
  ## glm
  for(j in c(1:64, 66:71,73:83)){
    print(j)
    fit.glm.pa.cv.loop <- glm.nb(formula(glm.formula.sp.pa[j]), data=dat.glm.train)
    preds.glm.pa.cv[partition.cv==i,j] <- predict(fit.glm.pa.cv.loop, dat.glm.test, type="response")
  }
  print("glm iteration done")
  ## glmm
  for(k in c(1:71,73:83)){
    print(k)
    fit.glmm.pa.cv.loop <- glmmTMB(formula(glmm.formula.sp.pa[k]), data=dat.glm.train, family="binomial")
    preds.glmm.pa.cv[partition.cv==i,k] <- predict(fit.glmm.pa.cv.loop, dat.glm.test, type="response")
  }
  print("glmm iteration done")
}

###############################################################################
###############################################################################
#### save outputs
save(partition, dat.for.fit, partition.cv, dat.for.cv, dat.full, dat.full.cv, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_dat.Rdata"))
save(fit.glm.r, fit.glmm.r, preds.glm.r.cv, preds.glmm.r.cv, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_richness.Rdata"))
save(fit.glm.c, fit.glmm.c, fit.glmm.ct, preds.glm.c.cv, preds.glmm.c.cv, preds.glmm.ct.cv, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_totalcover.Rdata"))
save(fit.glm.sf, fit.glmm.sf, preds.glm.sf.cv, preds.glmm.sf.cv, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_sfcover.Rdata"))
save(fit.glm.list, fit.glmm.list, preds.glm.ab.cv, preds.glmm.ab.cv, sp.not.na.sel.list, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_speciescover.Rdata"))
save(fit.pa.glm.list, fit.pa.glmm.list, preds.glm.pa.cv, preds.glmm.pa.cv, sp.not.na.sel.list, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_speciespa.Rdata"))






























# ## mobile counts
# summary(fit.m <- glm.nb(formula(paste0("count_mobile~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
# summary(fit.m <- update(fit.m, .~.-waom4k_seafloortemperature))
# summary(fit.m <- update(fit.m, .~.-distance2canyons))
# summary(fit.m <- update(fit.m, .~.-waom4k_seafloorcurrents_residual))
# summary(fit.m <- update(fit.m, .~.-waom4k_seafloorsalinity))
# 
# ## echinoderms
# summary(fit.e <- glm.nb(formula(paste0("abund_echino~",covars,"+offset(log(n_not_na))")), data=dat.glm))
# summary(fit.e <- update(fit.e, .~.-waom4k_seafloortemperature))
# summary(fit.e <- update(fit.e, .~.-distance2canyons))
# summary(fit.e <- update(fit.e, .~.-waom4k_seafloorcurrents_residual))
# summary(fit.e <- update(fit.e, .~.-waom4k_seafloorsalinity))
# 
# ## crustacea
# summary(fit.cr <- glm.nb(formula(paste0("abund_crust~",covars,"+offset(log(n_not_na))")), data=dat.glm))
# summary(fit.cr <- update(fit.cr, .~.-distance2canyons))
# summary(fit.cr <- update(fit.cr, .~.-waom4k_seafloorsalinity))
# summary(fit.cr <- update(fit.cr, .~.-waom4k_seafloortemperature))
# summary(fit.cr <- update(fit.cr, .~.-waom4k_seafloorcurrents_residual))

##### modelling bryozoan cover

## B-abundance
summary(fit.glm.b <- glm.nb(formula(paste0("cover_B~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.glm.b <- update(fit.glm.b, .~.-waom4k_seafloorcurrents_absolute))
summary(fit.glm.b <- update(fit.glm.b, .~.-waom4k_seafloorsalinity))
summary(fit.glm.b <- update(fit.glm.b, .~.-tpi))

summary(fit.glmm.b <- glmmTMB(cover_B ~ depth + depth2 + logslope +
                                waom4k_seafloorcurrents_mean +
                                waom4k_test_susp08 +
                                waom4k_test_settle08 +
                                offset(log(cover_points_scorable)) +
                                (1|cover_cells_transect1) + (1|cover_cells_survey),
                              data=dat.glm, family=nbinom1)) #dat.glm[-c(70,71),]
fit.sim <- simulateResiduals(fit.glmm.b)
plot(fit.sim)

save(fit.glm.b, fit.glmm.b, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_bcover.Rdata"))

##### modelling sponge cover

## S-abundance
summary(fit.glm.s <- glm.nb(formula(paste0("cover_S~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.glm.s <- update(fit.glm.s, .~.-distance2canyons2))
summary(fit.glm.s <- update(fit.glm.s, .~.-waom4k_test_settle08))
summary(fit.glm.s <- update(fit.glm.s, .~.-waom4k_seafloorcurrents_absolute))

summary(fit.glmm.s <- glmmTMB(cover_S ~ depth + depth2 + logslope +
                                waom4k_seafloorcurrents_mean +
                                waom4k_seafloorcurrents_residual +
                                waom4k_seafloorsalinity +
                                waom4k_test_susp08 +
                                waom4k_test_settle08 +
                                distance2canyons +
                                offset(log(cover_points_scorable)) +
                                (1|cover_cells_transect1) + (1|cover_cells_survey),
                              data=dat.glm, family=nbinom1)) #dat.glm[-c(70,71),]
fit.sim <- simulateResiduals(fit.glmm.s)
plot(fit.sim)

save(fit.glm.s, fit.glmm.s, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_scover.Rdata"))

#######################################################
cor(cell_metadata_env[,c(21,22,43,70:76)])











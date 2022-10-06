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

##############################################
# load(file=paste0(biodiv.dir,"biodiversity_pred_stack_dat.Rdata"))
load(file=paste0(biodiv.dir,"biodiversity_pred_stack_scaled_dat.Rdata"))
pred_stack.dat$cover_points_scorable <- 108

#sel2 <- which(!is.na(rowSums(pred_stack.dat[,1:10])))

##### richness
load(file=paste0(biodiv.dir,"biodiversity_fit_richness.Rdata"))

## slope, tpi and waom have some NAs
sel.nas <- which(is.na(pred_stack.dat[,5])|is.na(pred_stack.dat[,7])|is.na(pred_stack.dat[,10]))
ptm = proc.time() ## 5-10min
#pred1 <- predict(fit, pred_stack.dat[sel2[1:1000],], type="response")#, se.fit=TRUE)
pred1 <- predict(fit.glmm.r, pred_stack.dat[-sel.nas,][1:6000000,], type="response")
pred2 <- predict(fit.glmm.r, pred_stack.dat[-sel.nas,][6000001:12000000,], type="response")
pred3 <- predict(fit.glmm.r, pred_stack.dat[-sel.nas,][12000001:18000000,], type="response")
pred4 <- predict(fit.glmm.r, pred_stack.dat[-sel.nas,][18000001:nrow(pred_stack.dat),], type="response")
computational.time = proc.time() - ptm

pred.glmm.r <- c(pred1, pred2, pred3, pred4)
save(pred.glmm.r, file=paste0(biodiv.dir,"biodiversity_fit_richness_glmm_pred.Rdata"))

ptm = proc.time() #20sec
pred1 <- predict(fit.glm.r, pred_stack.dat[1:6000000,], type="response", se.fit=TRUE)
pred2 <- predict(fit.glm.r, pred_stack.dat[6000001:12000000,], type="response", se.fit=TRUE)
pred3 <- predict(fit.glm.r, pred_stack.dat[12000001:18000000,], type="response", se.fit=TRUE)
pred4 <- predict(fit.glm.r, pred_stack.dat[18000001:nrow(pred_stack.dat),], type="response", se.fit=TRUE)
computational.time = proc.time() - ptm
pred.glm.r <- c(pred1$fit, pred2$fit, pred3$fit, pred4$fit)
pred.glm.r.se <- c(pred1$se.fit, pred2$se.fit, pred3$se.fit, pred4$se.fit)
save(pred.glm.r, pred.glm.r.se, file=paste0(biodiv.dir,"biodiversity_fit_richness_glm_pred.Rdata"))

##### total cover
load(file=paste0(biodiv.dir,"biodiversity_fit_totalcover.Rdata"))
pred_stack.dat$cover_points_scorable <- 100

ptm = proc.time() #20sec
pred1 <- predict(fit.glm.c, pred_stack.dat[1:6000000,], type="response")#, se.fit=TRUE)
pred2 <- predict(fit.glm.c, pred_stack.dat[6000001:12000000,], type="response")#, se.fit=TRUE)
pred3 <- predict(fit.glm.c, pred_stack.dat[12000001:18000000,], type="response")#, se.fit=TRUE)
pred4 <- predict(fit.glm.c, pred_stack.dat[18000001:nrow(pred_stack.dat),], type="response")#, se.fit=TRUE)
computational.time = proc.time() - ptm
pred.glm.c <- c(pred1, pred2, pred3, pred4)
# pred.glm.c <- c(pred1$fit, pred2$fit, pred3$fit, pred4$fit)
# pred.glm.c.se <- c(pred1$se.fit, pred2$se.fit, pred3$se.fit, pred4$se.fit)
# save(pred.glm.c, pred.glm.c.se, file=paste0(biodiv.dir,"biodiversity_fit_totalcover_glm_pred.Rdata"))
save(pred.glm.c, file=paste0(biodiv.dir,"biodiversity_fit_totalcover_glm_pred.Rdata"))

##### sf cover
load(file=paste0(biodiv.dir,"biodiversity_fit_sfcover.Rdata"))

ptm = proc.time() #20sec
pred1 <- predict(fit.glm.sf, pred_stack.dat[1:6000000,], type="response")#, se.fit=TRUE)
pred2 <- predict(fit.glm.sf, pred_stack.dat[6000001:12000000,], type="response")#, se.fit=TRUE)
pred3 <- predict(fit.glm.sf, pred_stack.dat[12000001:18000000,], type="response")#, se.fit=TRUE)
pred4 <- predict(fit.glm.sf, pred_stack.dat[18000001:nrow(pred_stack.dat),], type="response")#, se.fit=TRUE)
computational.time = proc.time() - ptm
pred.glm.sf <- c(pred1, pred2, pred3, pred4)
# pred.glm.sf <- c(pred1$fit, pred2$fit, pred3$fit, pred4$fit)
# pred.glm.sf.se <- c(pred1$se.fit, pred2$se.fit, pred3$se.fit, pred4$se.fit)
# save(pred.glm.sf, pred.glm.sf.se, file=paste0(biodiv.dir,"biodiversity_fit_totalcover_glm_pred.Rdata"))
save(pred.glm.sf, file=paste0(biodiv.dir,"biodiversity_fit_sfcover_glm_pred.Rdata"))

##### b cover
load(file=paste0(biodiv.dir,"biodiversity_fit_bcover.Rdata"))

ptm = proc.time() #20sec
pred1 <- predict(fit.glm.b, pred_stack.dat[1:6000000,], type="response")#, se.fit=TRUE)
pred2 <- predict(fit.glm.b, pred_stack.dat[6000001:12000000,], type="response")#, se.fit=TRUE)
pred3 <- predict(fit.glm.b, pred_stack.dat[12000001:18000000,], type="response")#, se.fit=TRUE)
pred4 <- predict(fit.glm.b, pred_stack.dat[18000001:nrow(pred_stack.dat),], type="response")#, se.fit=TRUE)
computational.time = proc.time() - ptm
pred.glm.b <- c(pred1, pred2, pred3, pred4)
# pred.glm.b <- c(pred1$fit, pred2$fit, pred3$fit, pred4$fit)
# pred.glm.b.se <- c(pred1$se.fit, pred2$se.fit, pred3$se.fit, pred4$se.fit)
# save(pred.glm.b, pred.glm.b.se, file=paste0(biodiv.dir,"biodiversity_fit_bryozoancover_glm_pred.Rdata"))
save(pred.glm.b, file=paste0(biodiv.dir,"biodiversity_fit_bcover_glm_pred.Rdata"))

# sel.nas <- which(is.na(pred_stack.dat[,5])|is.na(pred_stack.dat[,7])|is.na(pred_stack.dat[,10]))
ptm = proc.time() ## 5-10min
pred1 <- predict(fit.glmm.b, pred_stack.dat[1:6000000,], type="response") #[-sel.nas,]
pred2 <- predict(fit.glmm.b, pred_stack.dat[6000001:12000000,], type="response") #[-sel.nas,]
pred3 <- predict(fit.glmm.b, pred_stack.dat[12000001:18000000,], type="response") #[-sel.nas,]
pred4 <- predict(fit.glmm.b, pred_stack.dat[18000001:nrow(pred_stack.dat),], type="response") #[-sel.nas,]
computational.time = proc.time() - ptm

pred.glmm.b <- c(pred1, pred2, pred3, pred4)
save(pred.glmm.b, sel.nas, file=paste0(biodiv.dir,"biodiversity_fit_bcover_glmm_pred.Rdata"))

##### s cover
load(file=paste0(biodiv.dir,"biodiversity_fit_scover.Rdata"))

ptm = proc.time() #20sec
pred1 <- predict(fit.glm.s, pred_stack.dat[1:6000000,], type="response")#, se.fit=TRUE)
pred2 <- predict(fit.glm.s, pred_stack.dat[6000001:12000000,], type="response")#, se.fit=TRUE)
pred3 <- predict(fit.glm.s, pred_stack.dat[12000001:18000000,], type="response")#, se.fit=TRUE)
pred4 <- predict(fit.glm.s, pred_stack.dat[18000001:nrow(pred_stack.dat),], type="response")#, se.fit=TRUE)
computational.time = proc.time() - ptm
pred.glm.s <- c(pred1, pred2, pred3, pred4)
# pred.glm.s <- c(pred1$fit, pred2$fit, pred3$fit, pred4$fit)
# pred.glm.s.se <- c(pred1$se.fit, pred2$se.fit, pred3$se.fit, pred4$se.fit)
save(pred.glm.s, file=paste0(biodiv.dir,"biodiversity_fit_scover_glm_pred.Rdata"))

# sel.nas <- which(is.na(pred_stack.dat[,5])|is.na(pred_stack.dat[,7])|is.na(pred_stack.dat[,10]))
ptm = proc.time() ## 5-10min
pred1 <- predict(fit.glmm.s, pred_stack.dat[1:6000000,], type="response") #[-sel.nas,]
pred2 <- predict(fit.glmm.s, pred_stack.dat[6000001:12000000,], type="response") #[-sel.nas,]
pred3 <- predict(fit.glmm.s, pred_stack.dat[12000001:18000000,], type="response") #[-sel.nas,]
pred4 <- predict(fit.glmm.s, pred_stack.dat[18000001:nrow(pred_stack.dat),], type="response") #[-sel.nas,]
computational.time = proc.time() - ptm

pred.glmm.s <- c(pred1, pred2, pred3, pred4)
save(pred.glmm.s, file=paste0(biodiv.dir,"biodiversity_fit_scover_glmm_pred.Rdata"))


################################################################
load(paste0(biodiv.dir,"biodiversity_fit_richness_glm_pred.Rdata"))
load(paste0(biodiv.dir,"biodiversity_fit_totalcover_glm_pred.Rdata"))
load(paste0(biodiv.dir,"biodiversity_fit_sfcover_glm_pred.Rdata"))
load(paste0(biodiv.dir,"biodiversity_fit_bcover_glm_pred.Rdata"))
load(paste0(biodiv.dir,"biodiversity_fit_scover_glm_pred.Rdata"))
load(paste0(biodiv.dir,"biodiversity_fit_bcover_glmm_pred.Rdata"))
load(paste0(biodiv.dir,"biodiversity_fit_scover_glmm_pred.Rdata"))

library(viridis)
sel <- which(!is.na(env_stack_scaled$depth[]))

pred.ra.r.glm <- raster(env_stack_scaled$depth)
pred.ra.r.glm[sel] <- pred.glm.r
pred.ra.r.glm.se <- raster(env_stack_scaled$depth)
pred.ra.r.glm.se[sel] <- pred.glm.r.se
pred.ra.r.glm.uci <- raster(env_stack_scaled$depth)
pred.ra.r.glm.lci <- raster(env_stack_scaled$depth)
pred.ra.r.glm.uci[sel] <- pred.glm.r+1.96*pred.glm.r.se
pred.ra.r.glm.lci[sel] <- pred.glm.r-1.96*pred.glm.r.se
# pred.ra.r.glmm <- raster(env_stack_scaled$depth)
# pred.ra.r.glmm[sel] <- pred.glmm.r

pred.ra.c.glm <- raster(env_stack_scaled$depth)
pred.ra.c.glm[sel] <- pred.glm.c

pred.ra.sf.glm <- raster(env_stack_scaled$depth)
pred.ra.sf.glm[sel] <- pred.glm.sf

pred.ra.b.glm <- raster(env_stack_scaled$depth)
pred.ra.b.glm[sel] <- pred.glm.b
pred.ra.b.glmm <- raster(env_stack_scaled$depth)
pred.ra.b.glmm[sel] <- pred.glmm.b

pred.ra.s.glm <- raster(env_stack_scaled$depth)
pred.ra.s.glm[sel] <- pred.glm.s
pred.ra.s.glmm <- raster(env_stack_scaled$depth)
pred.ra.s.glmm[sel] <- pred.glmm.s


## LIMITING THE SCALE!!!
pred.ra.r.glm[pred.ra.r.glm[]>50] <- 50

pred.ra.r.glm.uci[pred.ra.r.glm.uci[]>50] <- 50
pred.ra.r.glm.uci@data@min <- 0
pred.ra.r.glm.uci@data@max <- 50

pred.ra.r.glm.lci[pred.ra.r.glm.lci[]>50] <- 50
pred.ra.r.glm.lci[pred.ra.r.glm.lci[]<0] <- 0
pred.ra.r.glm.lci@data@min <- 0
pred.ra.r.glm.lci@data@max <- 50

pred.ra.c.glm[pred.ra.c.glm[]>100] <- 100

pred.ra.b.glm[pred.ra.b.glm[]>100] <- 100
pred.ra.b.glmm[pred.ra.b.glmm[]>100] <- 100
pred.ra.s.glm[pred.ra.s.glm[]>100] <- 100

######
par(mfrow=c(1,1))
plot(pred.ra.r.glm, col=viridis(99))

## Bryozoans
cex.b <- dat_cov_sum$cover_B/cell_metadata_env_clean_scaled$cover_points_scorable
plot(pred.ra.b.glm, xlim=c(-2800000,-2200000),ylim=c(1200000,1900000), main="Bryozoan cover - GLM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.b*10, pch=1)
plot(pred.ra.b.glmm, xlim=c(-2800000,-2200000),ylim=c(1200000,1900000), main="Bryozoan cover - GLMM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.b*10, pch=1)

plot(pred.ra.b.glm, xlim=c(-400000,400000),ylim=c(-2000000,-1200000), main="Bryozoan cover - GLM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.b*10, pch=1)
plot(pred.ra.b.glmm, xlim=c(-400000,400000),ylim=c(-2000000,-1200000), main="Bryozoan cover - GLMM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.b*10, pch=1)


## Sponges
cex.s <- dat_cov_sum$cover_S/cell_metadata_env_clean_scaled$cover_points_scorable
plot(pred.ra.s.glm, xlim=c(-2800000,-2200000),ylim=c(1200000,1900000), main="Sponge cover - GLM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.s*10, pch=1)
plot(pred.ra.s.glmm, xlim=c(-2800000,-2200000),ylim=c(1200000,1900000), main="Sponge cover - GLMM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.s*10, pch=1)

plot(pred.ra.s.glm, xlim=c(-400000,400000),ylim=c(-2000000,-1200000), main="Sponge cover - GLM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.s*10, pch=1)
plot(pred.ra.s.glmm, xlim=c(-400000,400000),ylim=c(-2000000,-1200000), main="Sponge cover - GLMM")
points(cell_metadata_env_clean_scaled[,4:5], col="blue", cex=cex.s*10, pch=1)


## AP
xlim=c(-2750000,-2050000)
ylim=c(800000,1900000)
par(mfrow=c(1,3), mar=c(4,3,2,2), oma=c(0,0,0,2))
plot(pred.ra.r.glm, xlim=xlim, ylim=ylim, col=viridis(99), main="species richness per 10 sqm", axes=FALSE)
# plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)
plot(pred.ra.c.glm, xlim=xlim, ylim=ylim, col=viridis(99), main="total cover", axes=FALSE)
# plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)
plot(r2, xlim=xlim, ylim=ylim, col=viridis(99), main="depth in m", axes=FALSE)
# plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)

## Weddell Sea
xlim=c(-1200000,-500000)
ylim=c(850000,1950000)
par(mfrow=c(1,3), mar=c(4,3,2,2), oma=c(0,0,0,2))
plot(pred.ra.r.glm, xlim=xlim, ylim=ylim, col=viridis(99), main="species richness per 10 sqm", axes=FALSE)
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.7, pch=16)
plot(pred.ra.c.glm, xlim=xlim, ylim=ylim, col=viridis(99), main="total cover", axes=FALSE)
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.7, pch=16)
plot(r2, xlim=xlim, ylim=ylim, col=viridis(99), main="depth in m", axes=FALSE)
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.7, pch=16)

## Ross Sea
xlim=c(-700000,500000)
ylim=c(-2100000,-1200000)
par(mfrow=c(2,2), mar=c(4,3,2,2), oma=c(0,0,0,2))
plot(pred.ra.r.glm, xlim=xlim, ylim=ylim, col=viridis(99), main="species richness per 10 sqm", axes=FALSE)
plot(coast.proj, add=TRUE)
plot(pred.ra.c.glm, xlim=xlim, ylim=ylim, col=viridis(99), main="total cover", axes=FALSE)
plot(coast.proj, add=TRUE)
plot(r2, xlim=xlim, ylim=ylim, col=viridis(99), main="depth in m", axes=FALSE)
plot(coast.proj, add=TRUE)


#######################################


## AP - richness
xlim <- c(-2750000,-2050000)
ylim <- c(800000,1900000)
breaks <- seq(0,50, by=0.5)
col <- viridis(100)
at <- seq(0,50, by=10)
labels <- seq(0,50,by=10)

par(mfrow=c(1,3), mar=c(4,1,2,2), oma=c(0,0,0,2))
plot(pred.ra.r.glm.lci, xlim=xlim, ylim=ylim, col=col, main="lower confidence interval", axes=FALSE,
     breaks=breaks, legend=FALSE)
plot(pred.ra.r.glm.lci, legend.only=TRUE, breaks=breaks, col=col,
     axis.args=list(at=at, labels=labels))
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)

plot(pred.ra.r.glm, xlim=xlim, ylim=ylim, col=col, main="species richness per 10 sqm", axes=FALSE,
     breaks=breaks, legend=FALSE)
plot(pred.ra.r.glm, legend.only=TRUE, breaks=breaks, col=col,
     axis.args=list(at=at, labels=labels))
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)

plot(pred.ra.r.glm.uci, xlim=xlim, ylim=ylim, col=col, main="upper confidence interval", axes=FALSE,
     breaks=breaks, legend=FALSE)
plot(pred.ra.r.glm.uci, legend.only=TRUE, breaks=breaks, col=col,
     axis.args=list(at=at, labels=labels))
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)


## Weddell Sea
xlim=c(-1200000,-500000)
ylim=c(850000,1950000)
breaks <- seq(0,50, by=0.5)
col <- viridis(100)
at <- seq(0,50, by=10)
labels <- seq(0,50,by=10)

par(mfrow=c(1,3), mar=c(4,1,2,2), oma=c(0,0,0,2))
plot(pred.ra.r.glm.lci, xlim=xlim, ylim=ylim, col=col, main="lower confidence interval", axes=FALSE,
     breaks=breaks, legend=FALSE)
plot(pred.ra.r.glm.lci, legend.only=TRUE, breaks=breaks, col=col,
     axis.args=list(at=at, labels=labels))
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)
plot(pred.ra.r.glm, xlim=xlim, ylim=ylim, col=col, main="species richness per 10 sqm", axes=FALSE,
     breaks=breaks, legend=FALSE)
plot(pred.ra.r.glm, legend.only=TRUE, breaks=breaks, col=col,
     axis.args=list(at=at, labels=labels))
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)
plot(pred.ra.r.glm.uci, xlim=xlim, ylim=ylim, col=col, main="upper confidence interval", axes=FALSE,
     breaks=breaks, legend=FALSE)
plot(pred.ra.r.glm.uci, legend.only=TRUE, breaks=breaks, col=col,
     axis.args=list(at=at, labels=labels))
plot(coast.proj, add=TRUE)
points(cell_metadata_env_clean_scaled[,4:5], col="red", cex=0.5, pch=16)




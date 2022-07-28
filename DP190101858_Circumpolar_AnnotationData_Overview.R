### 1) ### Setting up----
library(raster)
library(readxl)
library(readr)
library(dplyr)
library(data.table)
library(proj4)
library(stringr)
library(RColorBrewer)
library(SOmap)

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

## functions
source(paste0(tools.dir,"SOmap_functions_JJ.R"))

## projection
stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

## bathymetry
## from "ReadIn_Circumpolar_Environmental_Data.Rmd"
r2 <- raster(paste0(env.derived,"Circumpolar_EnvData_500m_shelf_bathy_gebco_depth.grd"))

## load coastline
load(paste0(env.derived,"Circumpolar_Coastline.Rdata"))


##### load biological and environmental data
load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Data.Rdata"))
load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Env_Data.RData"))

##### we need spatial data to plot onto SO_map
spatial.dat <- data.frame(cell_metadata[,2:3])
names(spatial.dat) <- c("Longitude","Latitude")
coordinates(spatial.dat) <- c("Longitude","Latitude")
proj4string(spatial.dat) <- CRS("+proj=longlat +datum=WGS84")
polar.dat.cells <- spTransform(spatial.dat, CRS(stereo))


##### set colour parameters
## disrete colours:
leg.cols <- brewer.pal(12,"Paired")
ross_sea_cols <- leg.cols[c(4,8,11)]
wap_cols <- leg.cols[3:11]#leg.cols[c(4,8,11,6)]

## set continuous colors for raster:
# depth.pal <- colorRampPalette(c('black','grey10','grey20','navy','blue','dodgerblue','lightcyan'))
# range.pal <- colorRampPalette(c('black','grey10','grey20','navy','blue','dodgerblue','skyblue1','cadetblue1','aquamarine1','palegreen2','green2','OliveDrab1','yellow1','wheat1'))
#better.ramp2 <- grDevices::colorRampPalette(c("#54A3D1", "#54A3D1", "#54A3D1", "#54A3D1", "#54A3D1", "#54A3D1", "#60B3EB", "#60B3EB", "#78C8F0", "#98D1F5", "#B5DCFF", "#BDE1F0", "#CDEBFA", "#D6EFFF", "#EBFAFF")) #, "grey99", "grey90", "grey92", "grey94", "grey96", "white"
#blues <- c(blues9[9],blues9[9],rev(blues9))
blues <- rev(blues9)

## choose depending on the bathymetry file used (whether deeper than -3000 is NA or not)
better.ramp2 <- grDevices::colorRampPalette(blues)
cols <- better.ramp2(99)
breaks <- seq(-3000,0,length.out=100)
leg.args <- seq(-3000,0, by=500)

#full.ramp <- grDevices::colorRampPalette(c("black","grey15","grey25",blues))
#cols <- full.ramp(132)
#breaks <- c(seq(-7500,-3000,length.out=34),seq(-2500,0,length.out=100))
#leg.args <- seq(-7500,0, by=500)


West.cols.AP <- brewer.pal(9,"Set1")
West.cols.WS <- brewer.pal(5,"Set1")[-2]
East.cols <- brewer.pal(6,"Set1")[-2]

West.cols.leg.AP <- West.cols.AP
West.cols.leg.WS <- c(West.cols.WS,rep("white",5))
East.cols.leg <- c(rep("white",4),East.cols)

West.cruise.AP <- c("JR17001","CRS","LMG1311","PS81","PS61","JR17003","PS118","JR262","JR15005")
West.cruise.WS <- c("PS96","PS06","PS14","PS18","","","","","")
East.cruise <- c("","","","","NBP1402","AA2011","TAN0802","TAN1802","TAN1901")

## points plotted in the sequence of objects in names(polar.dat.list.clean)
col.v <- c(West.cols.AP[2],West.cols.WS[2:4],West.cols.AP[c(5,4,4)],West.cols.WS[1],West.cols.AP[7],East.cols[c(3:5,1:2)],West.cols.AP[c(8:9,1,6,3)])

pch <- 20
cex <- 2

text.plot <- function(transects,dat,polar.dat,col="black", adj=1.5, cex=0.7){
  for(i in 1:length(transects)){
    pts <- polar.dat[dat$transectID==transects[i]]
    text(pts[1],labels=transects[i], adj=adj, cex=cex,col=col)
  }}
leg.steps <- c(1,2,2.75,3)

#### Overview of sampling locations

#### Image survey sites on bathymetry
#pdf(file="Circumpolar_overview_images.pdf", width=10, height=10)
par(oma=c(0,0,0,0))
JJ_SOmap(ramp.col = better.ramp2, border=FALSE, graticules = TRUE, label = "depth\n(m)", leg.steps=leg.steps)
JJ_SOleg(col=West.cols.leg.AP, ticks=length(West.cols.leg.AP), border_width=0, label="Survey IDs\nAntarctic Pensinsula", position = "topleft",tlabs=West.cruise.AP, leg.steps=leg.steps)
JJ_SOleg(col=West.cols.leg.WS, ticks=length(West.cols.WS), border_width=0, label="Survey IDs\nWeddell and Lazarev Sea", position = "topright",tlabs=West.cruise.WS, leg.steps=leg.steps)
JJ_SOleg(col=East.cols.leg, ticks=length(East.cols), border_width=0, label="Survey IDs\nEast Antarctic", position = "bottomright",tlabs=East.cruise, leg.steps=leg.steps)
plot(polar.dat.cells, add=TRUE, pch=pch, cex=cex, col="black")
# for(i in 1:length(polar.dat.list.clean)){
#   plot(polar.dat.list.clean[[i]], add=TRUE, pch=pch, cex=cex, col=col.v[i])
# }
scalebar(1000000, type="bar", label=c("0","500","1000"), below="km")
#dev.off()

par(mar=c(4,3,3,1))
plot(r2, xlim=c(-500000,500000), ylim=c(-2250000,-1250000),col=cols, main="Ross Sea", axes=FALSE, breaks=breaks, legend=FALSE)
plot(coast.proj, add=TRUE)
scalebar(200000, type="bar", label=c("0","100","200"), below="km")
points(polar.dat.cells, pch=pch, col=ross_sea_cols[1], cex=cex)
plot(r2,legend.only=TRUE,col=cols, breaks=breaks, axis.args=list(at=leg.args,labels=leg.args))

colSums(cover_cells)


## points and unscorables per cell
n_total <- rowSums(cover_cells)
n_na <- cover_cells$Unscorable
n_not_na <- n_total - n_na

## images per cell
cell_metadata$cover_N

## names of faunal groups for cover_cells:
dataset.names <- names(cover_cells)

## selector for each faunal class
sel_S <- grep("S_",substr(dataset.names,1,2))
sel_O <- grep("O_",substr(dataset.names,1,2))
sel_B <- c(grep("B_",substr(dataset.names,1,2)),grep("BH_",substr(dataset.names,1,3)),grep("BS_",substr(dataset.names,1,3)))
sel_M <- grep("M_",substr(dataset.names,1,2))
sel_E <- grep("E_",substr(dataset.names,1,2))
sel_Asc <- grep("Asc_",substr(dataset.names,1,4))
sel_TW <- grep("WP_TubeSF",dataset.names)
sel_Hy <- grep("Hyd",dataset.names)

sel_SF <- c(sel_S,sel_O,sel_B,sel_Asc,sel_TW,sel_Hy)

sel_sed_soft <- c(grep("Fine",dataset.names),grep("PbGrv",dataset.names))
sel_sed_loose <- c(grep("Cbble",dataset.names),grep("BioRu",dataset.names),grep("BioShl",dataset.names),grep("BioOth",dataset.names)) 
sel_sed_hard <- c(grep("Bould",dataset.names),grep("Rock",dataset.names))
sel_sed <- grep("Sub_",dataset.names)

cover_cells_pa <- cover_cells
cover_cells_pa[cover_cells_pa>0] <- 1

abund_SF.prop <- rowSums(cover_cells[,sel_SF])/n_not_na
abund_SF <- rowSums(cover_cells[,sel_SF])
pa_SF <- abund_SF
pa_SF[abund_SF>0] <- 1
richness <- rowSums(cover_cells_pa[,-sel_sed])#/n_total
richness.l <- rowSums(cover_cells_pa[,-sel_sed])/log(n_total)
abund_all.prop <- rowSums(cover_cells[,-sel_sed])/n_not_na
abund_all <- rowSums(cover_cells[,-sel_sed])
## translate abundance into circle size, remember to use the area rather than radius or diameter
## area = pi*r^2
## r = sqrt(area/pi)
cex_abund_SF <- sqrt(50*(abund_SF.prop+0.001)/pi)
cex_abund_all <- sqrt(50*(abund_all.prop+0.001)/pi)
cex_richness <- sqrt(richness/pi)

## plot some regions
par(mar=c(4,3,3,1))
par(mfrow=c(1,1))
plot(r2, xlim=c(-500000,500000), ylim=c(-2250000,-1250000),col=cols, main="Ross Sea", axes=FALSE, breaks=breaks, legend=FALSE)
plot(coast.proj, add=TRUE)
scalebar(200000, type="bar", label=c("0","100","200"), below="km")
points(polar.dat.cells, pch=1, col="black", cex=cex_abund_SF)
plot(r2,legend.only=TRUE,col=cols, breaks=breaks, axis.args=list(at=leg.args,labels=leg.args))

xlim <- c(1250000,1800000)
ylim <- c(-2300000,-1900000)
plot(r2, xlim=xlim, ylim=ylim,col=cols, main="George V shelf", axes=FALSE, breaks=breaks, legend=FALSE)
plot(coast.proj, add=TRUE)
scalebar(200000, type="bar", label=c("0","100","200"), below="km")
points(polar.dat.cells, pch=1, col="black", cex=cex_abund_SF)
plot(r2,legend.only=TRUE,col=cols, breaks=breaks, axis.args=list(at=leg.args,labels=leg.args))

plot(r2, xlim=xlim, ylim=ylim,col=cols, main="George V shelf", axes=FALSE, breaks=breaks, legend=FALSE)
plot(coast.proj, add=TRUE)
scalebar(200000, type="bar", label=c("0","100","200"), below="km")
points(polar.dat.cells, pch=1, col="black", cex=cex_richness)
plot(r2,legend.only=TRUE,col=cols, breaks=breaks, axis.args=list(at=leg.args,labels=leg.args))

plot(r2, xlim=xlim, ylim=ylim,col=cols, main="George V shelf", axes=FALSE, breaks=breaks, legend=FALSE)
plot(coast.proj, add=TRUE)
scalebar(200000, type="bar", label=c("0","100","200"), below="km")
points(polar.dat.cells, pch=1, col="black", cex=cex_abund_all)
plot(r2,legend.only=TRUE,col=cols, breaks=breaks, axis.args=list(at=leg.args,labels=leg.args))


## plot some responses to the environment:
par(mfrow=c(2,2))
plot(abund_SF.prop~ cell_metadata_env$depth)
plot(abund_SF.prop~ cell_metadata_env$slope)
plot(abund_SF.prop~ cell_metadata_env$waom4k_seafloorcurrents_mean)
plot(abund_SF.prop~ cell_metadata_env$waom4k_test_settle08)

dat <- cbind(abund_SF,n_not_na,cell_metadata_env[,20:75])
dat$depth2 <- poly(dat$depth,2)[,2]
summary(fit <- glm(abund_SF.prop~depth+I(depth^2)+log(slope)+
                     waom4k_seafloorcurrents_mean+
                     waom4k_seafloorcurrents_residual+
                     waom4k_seafloortemperature+
                     waom4k_seafloorsalinity+
                     distance2canyons
                   ,data=dat, family="binomial"))

summary(fit <- glm(abund_SF.prop~depth+I(depth^2)+log(slope)+
                     waom4k_seafloorcurrents_mean+
                     distance2canyons
                   ,data=dat, family="binomial"))


summary(fit <- glm.nb(abund_SF~depth+I(depth^2)+log(slope)+
                     waom4k_seafloorcurrents_mean+
                     distance2canyons, data=dat))





################################################################
##### fitting a hurdle model on suspension feeder abundances
dat <- cbind(abund_SF,pa_SF, n_not_na, cell_metadata_env[,20:75])
## pa:
summary(fit <- glm(pa_SF~depth+I(depth^2)+log(slope)+
                     waom4k_seafloorcurrents_mean+
                     waom4k_seafloorcurrents_residual+
                     waom4k_seafloortemperature+
                     waom4k_seafloorsalinity+
                     waom4k_test_flux08+
                     waom4k_test_settle08+
                     distance2canyons
                   ,data=dat, family="binomial"))
summary(fit <- glm(pa_SF~depth+I(depth^2)+log(slope)+
                     waom4k_seafloorcurrents_mean+
                     waom4k_seafloorcurrents_residual+
                     distance2canyons
                   ,data=dat, family="binomial"))
## abund conditional on presence
dat.abund <- dat[dat$pa_SF==1,]
summary(fit <- glm.nb(abund_SF~depth+I(depth^2)+log(slope)+
                     waom4k_seafloorcurrents_mean+
                     waom4k_seafloorcurrents_residual+
                     waom4k_seafloortemperature+
                     waom4k_seafloorsalinity+
                     waom4k_test_flux08+
                     waom4k_test_settle08+
                     distance2canyons+ offset(log(n_not_na)),
                   data=dat.abund))
summary(fit <- glm.nb(abund_SF~depth+I(depth^2)+log(slope)+
                        waom4k_seafloorcurrents_mean+
                        waom4k_seafloortemperature+
                        waom4k_seafloorsalinity+
                        waom4k_test_flux08+
                        waom4k_test_settle08+
                        distance2canyons+ offset(log(n_not_na)),
                      data=dat.abund))
##check fit:
rootogram(fit)
nas <- which(is.na(rowSums(dat.abund%>%dplyr::select(depth,slope,
                                               waom4k_seafloorcurrents_mean,
                                               waom4k_test_flux08,
                                               waom4k_test_settle08))))

plot(fit$fitted.values, dat.abund$abund_SF[-nas])
















library(countreg)
library(pscl)
## first, fit the full model
summary(fit <- hurdle(abund_SF~
                        depth+depth2+log(slope)+
                        waom4k_seafloorcurrents_mean+
                        waom4k_seafloortemperature+
                        waom4k_seafloorsalinity+
                        waom4k_test_flux08+
                        waom4k_test_susp08+
                        # distance2canyons+
                        ice_prop+
                        ice_mean+
                        tpi,
                      data=dat, dist="negbin"))
## write separate model formulas and remove for each model the least relevant term 
summary(fit <- hurdle(abund_SF~
                        depth+depth2+log(slope)+
                        waom4k_seafloorcurrents_mean+
                        waom4k_seafloortemperature+
                        waom4k_seafloorsalinity+
                        #waom4k_test_settle08+
                        waom4k_test_susp08+
                        waom4k_test_flux08+
                        #distance2canyons+
                        ice_prop+
                        ice_mean+
                        tpi,
                      data=dat, dist="negbin"))
## backwards select terms
summary(fit<-hurdle(abund_SF~
                      depth+depth2+log(slope)+
                      waom4k_seafloorcurrents_mean+
                      waom4k_test_flux08+
                      ice_mean+
                      tpi
                    |
                      depth+depth2+log(slope)+
                      waom4k_seafloorcurrents_mean+
                      waom4k_test_susp08+
                      ice_mean,
                    data=dat, dist="negbin"))
## check fit:
rootogram(fit)

nas <- which(is.na(rowSums(dat%>%dplyr::select(depth,depth2,slope,
                                               waom4k_seafloorcurrents_mean,
                                               waom4k_test_flux08,
                                               waom4k_test_susp08,
                                               ice_mean,
                                               tpi))))

plot(fit$fitted.values, dat$abund_SF[-nas])


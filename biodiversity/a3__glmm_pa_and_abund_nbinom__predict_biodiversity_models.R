## Predicting a p/a hmsc onto the entire Antarctic continental shelf
## prediction is splitted into chunks of small areas to speed up the process

library(glmmTMB)
library(terra)
'%!in%' <- function(x,y)!('%in%'(x,y))

##############################################################################################################
#res <- "500m"
res <- "2km"
##############################################################################################################
sci.dir <- "C:/Users/jjansen/Desktop/science/"
DP.dir <-       paste0(sci.dir,"DP190101858_MappingAntarcticSeafloorBiodiversity/")
env.derived <- paste0(sci.dir,"data_environmental/derived/")
biodiv.dir <- paste0(sci.dir, "SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")

r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth
## create an empty raster to fill for mapping
empty.ra <- rast(r2)
empty.ra[] <- NA

## specify model to load
load(file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_dat.Rdata"))
load(file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_richness.Rdata"))
load(file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_totalcover.Rdata"))
load(file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_sfcover.Rdata"))
load(file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_speciescover.Rdata"))
load(file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_glmms_speciespa.Rdata"))

#############################################################
###### prepare data (RUN ONCE)
# ## load per cell environmental data
# load(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled_dataframe.Rdata"))
# ## identify which cells we have ignored before
# sel.not.na <- which(!is.na(r2[]))
# ## depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+
# ## seafloortemperature+seafloorcurrents_mean+seafloorsalinity+npp_mean
# ## relevant env data
# grid <- pred_stack.dat[,c(1,28,30,3,8,29,27,23,26,16)]
# ## spatial information
# #xy <- crds(r2)
# xy.grid.raw <- crds(r2) #xy #[sel.not.na,]
# ## remove NAs
# sel <- which(complete.cases(grid))
# XData.grid <- grid[sel,]
# xy.grid <- data.frame(xy.grid.raw[sel,])
# 
# ## save data
# save(sel, sel.not.na, file=paste0(biodiv.dir,"pred_",res,"_model_cell_sel.Rdata"))
# save(XData.grid, xy.grid, file=paste0(biodiv.dir,"pred_",res,"_model_cell_grid.Rdata"))
# rm(xy.grid.raw, grid, pred_stack.dat)

#############################################################
## load data
load(paste0(biodiv.dir,"pred_",res,"_model_cell_sel.Rdata"))
load(paste0(biodiv.dir,"pred_",res,"_model_cell_grid.Rdata"))

X <- XData.grid
X$cover_cells_transect1 <- as.factor("PS96_001")
X$transectID <- as.factor("PS96_001")
X$cover_cells_survey <- as.factor("PS96")
X$cover_points_scorable <- 540

#############################################################
## make predictions
## logslope has some negative infinite values: change to super low values instead
inf.sel <- which(is.infinite(X$logslope))
X$logslope[inf.sel] <- min(X$logslope[-inf.sel])

## takes about 10 seconds (30min per 100k points with se.fit)
ptm = proc.time()
preds.r <- predict(fit.glmm.r, newdata=X, type="response")
computational.time = proc.time() - ptm

## takes about 10 seconds (30min per 100k points with se.fit)
ptm = proc.time()
preds.c <- predict(fit.glmm.c, newdata=X, type="response")#, se.fit=TRUE)
computational.time = proc.time() - ptm

## takes about 10 seconds (30min per 100k points with se.fit)
ptm = proc.time()
preds.ct <- predict(fit.glmm.ct, newdata=X, type="response")#, se.fit=TRUE)
computational.time = proc.time() - ptm

## for se.fit, breaking it up into a loop to be able to check progress (takes 6h in total)
start.vec <- seq(1,900001, by=100000) #c(1,10001,20001)
stop.vec <- c(seq(100000, 900000, by=100000),991094) #c(10000,20000,30000)
ptm = proc.time()
preds.c <- predict(fit.glmm.c, newdata=X[start.vec[1]:stop.vec[1],], type="response", se.fit=TRUE)
for(i in 2:10){
  print(proc.time() - ptm)
  message(i)
  pred.loop <- predict(fit.glmm.c, newdata=X[start.vec[i]:stop.vec[i],], type="response", se.fit=TRUE)
  preds.c$fit <- c(preds.c$fit,pred.loop$fit)
  preds.c$se.fit <- c(preds.c$se.fit,pred.loop$se.fit)
}
computational.time = proc.time() - ptm
computational.time

# preds.r.se <- NA
# preds.r.se[1:100000] <- predict(fit.glmm.r, newdata=X[1:100000,], type="response", se.fit=TRUE)
# preds.c <- predict(fit.glmm.r, newdata=X, type="response", se.fit=TRUE)

#############################################################
## which cells do we need to fill with data
row.numbers <- as.numeric(rownames(xy.grid))
all.sel.ra <- sel.not.na[sel[row.numbers]]

preds.r.adj <- preds.r
ra.r <- r2
ra.r[] <- NA
ra.r[all.sel.ra] <- preds.r.adj

preds.c.adj <- preds.c$fit
preds.c.adj.se <- preds.c$se.fit
ra.c <- r2
ra.c[] <- NA
ra.c.se <- ra.c
ra.c[all.sel.ra] <- preds.c.adj
ra.c.se[all.sel.ra] <- preds.c.adj.se

preds.ct.adj <- preds.ct#$fit
#preds.ct.adj.se <- preds.ct$se.fit
ra.ct <- r2
ra.ct[] <- NA
#ra.ct.se <- ra.c
ra.ct[all.sel.ra] <- preds.ct.adj
#ra.ct.se[all.sel.ra] <- preds.ct.adj.se

plot(ra.r)
plot(ra.c)
plot(ra.ct)

## some crazy high richness and abundance numbers:
## need to adjust the numbers based on the environmental space
## same approach as in script "display_biodiversity_maps_setup.R"
gaps.hv.a <- rast(paste0(DP.dir,"Circumpolar_Analysis_GapHypervolume_2km_All21SurveysCombined_AllVariablesExceptCanyons.tif"))
gaps.hv.a.notzero.sel <- which(gaps.hv.a[]>0)

ra.r.corrected <- ra.r
# identify the maximum value predicted within the observed environmental space
max.range.within <- max(ra.r.corrected[gaps.hv.a.notzero.sel], na.rm=TRUE)
# identify cells in the raster that exceed that value
limit.sel <- which(ra.r.corrected[]>max.range.within)
# set values outside of environmental space to the maximum identified within
ra.r.corrected[limit.sel] <- max.range.within
# print max value
print(max.range.within)

ra.c.corrected <- ra.c
# identify the maximum value predicted within the observed environmental space
max.range.within <- max(ra.c.corrected[gaps.hv.a.notzero.sel], na.rm=TRUE)
# identify cells in the raster that exceed that value
limit.sel <- which(ra.c.corrected[]>max.range.within)
# set values outside of environmental space to the maximum identified within
ra.c.corrected[limit.sel] <- max.range.within
# print max value
print(max.range.within)
plot(ra.c.corrected)

ra.ct.corrected <- ra.ct
# identify the maximum value predicted within the observed environmental space
max.range.within <- max(ra.ct.corrected[gaps.hv.a.notzero.sel], na.rm=TRUE)
# identify cells in the raster that exceed that value
limit.sel <- which(ra.ct.corrected[]>max.range.within)
# set values outside of environmental space to the maximum identified within
ra.ct.corrected[limit.sel] <- max.range.within
# print max value
print(max.range.within)
plot(ra.ct.corrected)

#############################################################
## save output
base.str <- paste0(sci.dir,"data_biological/circumpolar_prediction_outputs/")
writeRaster(ra.r.corrected, paste0(base.str,res,"_model_cells_GLMM_richness.tif"))
writeRaster(ra.c.corrected, paste0(base.str,res,"_model_cells_GLMM_totalcover.tif"), overwrite=TRUE)
writeRaster(ra.c.se, paste0(base.str,res,"_model_cells_GLMM_totalcover_se.tif"))
writeRaster(ra.ct.corrected, paste0(base.str,res,"_model_cells_GLMM_totalcover_withtransects.tif"), overwrite=TRUE)


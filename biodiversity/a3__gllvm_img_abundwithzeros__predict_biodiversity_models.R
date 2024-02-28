## Predicting a p/a hmsc onto the entire Antarctic continental shelf
## prediction is splitted into chunks of small areas to speed up the process

library(gllvm)
library(terra)
'%!in%' <- function(x,y)!('%in%'(x,y))

##############################################################################################################
#res <- "500m"
res <- "2km"
##############################################################################################################
sci.dir <- "C:/Users/jjansen/Desktop/science/"
env.derived <- paste0(sci.dir,"data_environmental/derived/")
biodiv.dir <- paste0(sci.dir, "SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")

r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth
## create an empty raster to fill for mapping
empty.ra <- rast(r2)
empty.ra[] <- NA

## specify model to load
load("biodiversity/image_model_gllvm_re_on_transectID.Rdata")

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

# ## size of prediction boxes
# xmin <- seq(-3000000,2750000, by=250000)
# xmax <- seq(-2750000,3000000, by=250000)
# ymin <- seq(-3000000,2750000, by=250000)
# ymax <- seq(-2750000,3000000, by=250000)
xmin <- seq(-3000000,2950000, by=50000)
xmax <- seq(-2950000,3000000, by=50000)
ymin <- seq(-3000000,2950000, by=50000)
ymax <- seq(-2950000,3000000, by=50000)

plot(r2)
for(i in 2:length(xmin)){
  abline(h=ymin[i])
  abline(v=xmin[i])
}

X <- XData.grid
X$transectID <- as.factor("PS96_001")
preds <- predict.gllvm(fit, newX=X[1:10000,], type="response", level=0)




#############################################################
## RUN ONLY ONCE: 
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
# save(cells_with_data, file=paste0("/pvol/biodiversity_prediction/",res,"_model_50km_cells_with_data.Rdata"))
load(file=paste0("/pvol/biodiversity_prediction/",res,"_model_50km_cells_with_data.Rdata"))


#############################################################

## parallel processing: PER CELL that contains values
library(doParallel)
library(foreach)
parallel::detectCores()
#UseCores = parallel::detectCores() - 1
UseCores = 16
c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
registerDoParallel(c1)
getDoParWorkers()

cell.sel.v <- which(!is.na(cells_with_data))
cell.sel.df <- which(!is.na(cells_with_data), arr.ind = TRUE)

#iterations <- 60
ptm = proc.time()
#foreach(j=1:iterations) %dopar%{ #3:length(xmin)
foreach(j=1:length(cell.sel.v)) %dopar%{ #3:length(xmin)
  #library(Hmsc)
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
  dat.name <- paste0("/pvol/biodiversity_prediction/pred_files/",res,"/pa/",res,"_model_", as.character(model), "_",
                     c("pa","abundance")[modeltype],
                     "_chains_",as.character(nChains),
                     "_thin_", ... = as.character(thin),
                     "_samples_", as.character(samples),
                     "_pred_")
  run.name <- sprintf("%05d",cell.sel.v[j])
  ## not enough space to save the ~1GB size predY.loop files for each iteration
  #save(predY.loop, file=paste0(dat.name,"fulldat_",run.name,".Rdata"))
  ## saving only the derivatives
  save(predY.mean, predY.sd, sel.loop, XData.grid.loop, xy.grid.loop,
       file=paste0(dat.name,run.name,".Rdata"))
  rm(predY.loop, predY.mean, predY.sd)
}
computational.time = proc.time() - ptm
parallel::stopCluster(cl = c1)


## If not all areas/cell have run successfully, we can use the code below to only run the unsuccessful ones again

## identifying which cells/regions are already predicted and saved to file:
pred.list <- list.files("/pvol/biodiversity_prediction/pred_files/2km/pa/")
pred.list <- pred.list[grep("800_pred",pred.list)]
#pred.list <- pred.list[-grep("fulldat",pred.list)]
pred.list.numbers <- substr(pred.list,50,54)
pred.list.numbers <- as.numeric(sub("^0+", "", pred.list.numbers) )

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
  dat.name <- paste0("/pvol/biodiversity_prediction/pred_files/",res,"/",res,"_model_", as.character(model), "_",
                     c("pa","abundance")[modeltype],
                     "_chains_",as.character(nChains),
                     "_thin_", ... = as.character(thin),
                     "_samples_", as.character(samples),
                     "_pred_")
  run.name <- sprintf("%05d",cell.sel.v[j])
  #save(predY.loop, file=paste0(dat.name,"fulldat_",run.name,".Rdata"))
  save(predY.mean, predY.sd, sel.loop, XData.grid.loop, xy.grid.loop,
       file=paste0(dat.name,run.name,".Rdata"))
  rm(predY.loop, predY.mean, predY.sd)
}
computational.time = proc.time() - ptm
parallel::stopCluster(cl = c1)

#############################################################
#############################################################
#############################################################



























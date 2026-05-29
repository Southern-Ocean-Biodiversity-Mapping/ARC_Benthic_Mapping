##############################################################################################################
###### MODEL PREDICTIONS ONTO ANTARCTIC SHELF
###### All data accessed and written from/to the dropbox
##############################################################################################################

##############################################################################################################
##### SELECT MODEL SETUP
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
# cell_id <- grid$cell_id
# sel <- which(complete.cases(grid))
# 
# ## Extract coordinates
# xy.grid.raw <- pred_df[, c("proj_coord_x", "proj_coord_y")]
# 
# ## prepare data to predict on
# XData.grid <- grid[sel, ]
# xy.grid    <- xy.grid.raw[sel, ]
# 
# ## Save for reuse
# save(sel, sel.not.na, file = file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_sel.Rdata")))
# save(XData.grid, xy.grid, file = file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_grid.Rdata")))
# save(cell_id, file = file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_cell_ids.Rdata")))

load(file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_sel.Rdata")))
load(file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_model_cell_grid.Rdata")))
load(file.path(model_dir, paste0("4_model_prediction/hmsc_", res, "_cell_ids.Rdata")))

pred_df <- readRDS(file.path(usr.dropbox.dir, "data_products/modelling_files/circum_antarctic/prediction_grid_lookup.rds"))

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
for(nm in model_ids) {
  ## parallel processing: PER CELL that contains values
  library(doParallel)
  library(foreach)
  parallel::detectCores()
  #UseCores = parallel::detectCores() - 1
  UseCores = 12
  c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
  registerDoParallel(c1)
  getDoParWorkers()
  
  message("====================================")
  message("Running prediction for model: ", nm)
  message("====================================")
  
  model_file <- file.path(model_dir, "2_fitting_and_running_models/",
                          paste0(res, "_model_cells_", nm, "_chains_", nChains, "_thin_", thin, "_samples_", samples, ".Rdata"))
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
    if (length(sel.loop) == 0) {return(NULL)}
    if (j %% 50 == 0) message("Tile ", j)
    global_row <- sel[sel.loop]               # index in pred_df
    cell_id.loop <- pred_df$cell_id[global_row]
    
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
    
    #############################################################
    ## FORCE CONSISTENT DIMENSIONS
    #############################################################
    ## AB: handle all cases
    if (is.null(dim(predY.loop.ab.array))) 
      # vector → single cell, single variable
      predY.loop.ab.array <- array(predY.loop.ab.array, dim = c(1, 1, length(predY.loop.ab.array)))
    
    ## Get posterior mean/median/uncertainty
    predY.pa.mean   <- apply(predY.loop.pa.array, 1:2, mean)
    predY.pa.median <- apply(predY.loop.pa.array, 1:2, median)
    predY.ab.mean   <- apply(predY.loop.ab.array, 1:2, mean)
    predY.ab.median <- apply(predY.loop.ab.array, 1:2, median)
    
    predY.pa.se <- apply(predY.loop.pa.array, 1:2, se)
    predY.pa.5  <- apply(predY.loop.pa.array, 1:2, p5)
    predY.pa.95 <- apply(predY.loop.pa.array, 1:2, p95)
    predY.ab.se <- apply(predY.loop.ab.array, 1:2, se)
    predY.ab.5  <- apply(predY.loop.ab.array, 1:2, p5)
    predY.ab.95 <- apply(predY.loop.ab.array, 1:2, p95)
    
    ## save-string for 100km cell tiles
    dat.name <- paste0(pred_dir, res,"_model_cells_",nm, "_chains_",as.character(nChains),"_thin_", as.character(thin),"_samples_", as.character(samples),"_pred_")
    run.name <- sprintf("%06d",cell.sel.v[j])
    
    save(predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.5, predY.ab.95,
         predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95,
         global_row, cell_id.loop, sel.loop, XData.grid.loop, xy.grid.loop,
         file=paste0(dat.name,modelspec,run.name,".Rdata"))
    rm(predY.loop.pa, predY.loop.ab, predY.loop.pa.array, predY.loop.ab.array,
       predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.5, predY.ab.95,
       predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95)
  }
  computational.time = proc.time() - ptm
  message("Finished model: ", nm)
  parallel::stopCluster(cl = c1)
}


# #############################################################
# ## RUN UNSUCCESSUL ONES AGAIN
# done <- list.files("/pvol/4_model_prediction/pred_files")
# #extract the 6-digit numbers at end of filenames
# done_ids <- as.integer(sub(".*envonly(\\d{6})\\.Rdata$", "\\1", done))
# 
# # which requested cells are NOT present in files
# missing <- setdiff(cell.sel.v, done_ids)
# 
# missing.ids <- which(cell.sel.v%in%missing)
# 
# nm <- model_ids[1]
# for (nm in model_ids) {
#   message("====================================")
#   message("Running prediction for model: ", nm)
#   message("====================================")
#   
#   model_file <- file.path(model_dir, "2_fitting_and_running_models/",
#                           paste0(res, "_model_cells_", nm, "_chains_", nChains, "_thin_", thin, "_samples_", samples, ".Rdata")
#   )
#   load(model_file)
#   pa <- models$mENV
#   ab <- models$mAB
#   rm(models)
# 
#   message("setup complete, parallel loops over prediction areas now")
#   
#   ptm = proc.time()
#   for(j in missing.ids){
#     i <- cell.sel.df[j,2]
#     k <- cell.sel.df[j,1]
#     ## select cells in tile
#     sel.loop <- which(xy.grid[,1]>xmin[i] & xy.grid[,1]<xmax[i] &
#                         xy.grid[,2]>ymin[k] & xy.grid[,2]<ymax[k])
#     if (length(sel.loop) == 0) {
#       return(NULL)
#     }
#     
#     if (j %% 50 == 0) message("Tile ", j)
#     ## subset predictors + coordinates
#     XData.grid.loop <- XData.grid[sel.loop,]
#     xy.grid.loop <- xy.grid[sel.loop,]
#     
#     ## setup prediction - pa
#     Gradient.pa = prepareGradient(pa, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
#     predY.loop.pa <- predict(pa, Gradient=Gradient.pa, expected=TRUE) ## this gives probabilities instead of integer outcomes
#     rm(Gradient.pa)
#     # mat.names <- dimnames(predY.loop.pa[[1]])
#     # predY.loop.pa <- array(unlist(predY.loop.pa), c(nrow(xy.grid.loop), ncol(pa$Y), samples*nChains), dimnames(predY.loop.pa[[1]]))
#     
#     ## setup prediction - abund
#     Gradient.ab = prepareGradient(ab, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
#     predY.loop.ab <- predict(ab, Gradient=Gradient.ab, expected=TRUE) ## this gives probabilities instead of integer outcomes
#     rm(Gradient.ab)  
#     # mat.names <- dimnames(predY.loop.ab[[1]])
#     # predY.loop.ab <- array(unlist(predY.loop.ab), c(nrow(xy.grid.loop), ncol(ab$Y), samples*nChains), dimnames(predY.loop.ab[[1]]))
#     
#     # Convert matrices
#     predY.loop.pa.array = simplify2array(predY.loop.pa)
#     predY.loop.ab.array = simplify2array(predY.loop.ab)
#     
#     if(nrow(predY.loop.pa.array)==1){
#       # Get posterior median
#       predY.pa.mean   = as.data.frame(apply(predY.loop.pa.array[1,,], 1:2, mean))
#       predY.pa.median = as.data.frame(apply(predY.loop.pa.array[1,,], 1:2, median))
#       predY.ab.mean   = mean(predY.loop.ab.array)
#       predY.ab.median = median(predY.loop.ab.array)
#      # Posterior width
#       predY.pa.se = as.data.frame(apply(predY.loop.pa.array[1,,], 1:2, se))
#       predY.pa.5  = as.data.frame(apply(predY.loop.pa.array[1,,], 1:2, p5))
#       predY.pa.95 = as.data.frame(apply(predY.loop.pa.array[1,,], 1:2, p95))
#       predY.ab.se = se(predY.loop.ab.array)
#       predY.ab.5  = p5(predY.loop.ab.array)
#       predY.ab.95 = p95(predY.loop.ab.array)
#     }else{
#       # Get posterior median
#       predY.pa.mean   = as.data.frame(apply(predY.loop.pa.array, 1:2, mean))
#       predY.pa.median = as.data.frame(apply(predY.loop.pa.array, 1:2, median))
#       predY.ab.mean   = as.data.frame(apply(predY.loop.ab.array, 1:2, mean))
#       predY.ab.median = as.data.frame(apply(predY.loop.ab.array, 1:2, median))
#       # Posterior width
#       predY.pa.se = as.data.frame(apply(predY.loop.pa.array, 1:2, se))
#       predY.pa.5  = as.data.frame(apply(predY.loop.pa.array, 1:2, p5))
#       predY.pa.95 = as.data.frame(apply(predY.loop.pa.array, 1:2, p95))
#       predY.ab.se = as.data.frame(apply(predY.loop.ab.array, 1:2, se))
#       predY.ab.5  = as.data.frame(apply(predY.loop.ab.array, 1:2, p5))
#       predY.ab.95 = as.data.frame(apply(predY.loop.ab.array, 1:2, p95))
#     }
# 
#     ## save-string for 100km cell tiles
#     dat.name <- paste0(pred_dir, res,"_model_cells_",nm, "_chains_",as.character(nChains),"_thin_", as.character(thin),"_samples_", as.character(samples),"_pred_")
#     run.name <- sprintf("%06d",cell.sel.v[j])
#     
#     save(predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.5, predY.ab.95,
#          predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95,
#          sel.loop, XData.grid.loop, xy.grid.loop,
#          file=paste0(dat.name,modelspec,run.name,".Rdata"))
#     rm(predY.loop.pa, predY.loop.ab, predY.loop.pa.array, predY.loop.ab.array,
#        predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.5, predY.ab.95,
#        predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95)
#   }
#   computational.time = proc.time() - ptm
#   message("Finished model: ", nm)
# }




## a3__hmsc_cells_pa_ab__predict_biodiversity_models.R

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
model_dir <- paste0(usr.dropbox.dir, "data_products/modelling_files/circum_antarctic")
pred_dir  <- paste0(usr.main.dir, "4_model_prediction/pred_files/")
model_dir_local <- paste0(usr.main.dir, "2_fitting_and_running_models")

res <- "2km"

thin     <- 10
samples  <- 800
nChains  <- 4

modelspec <- "_envonly"

se <- function(x) {
  sd(x) / sqrt(length(x))
}
lwr95 <- function(x, z=1.96) {
  mean(x) - z * sd(x) / sqrt(length(x))
}
upr95 <- function(x, z=1.96) {
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
  "fam_cafe","fam_cbpm","fam_eppl","fam_vpmg",
  "npp_and_fam_cafe","npp_and_fam_cbpm","npp_and_fam_eppl","npp_and_fam_vpmg"
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
# grid <- pred_df[,!names(pred_df) %in% c("arag_sd","o2_mean","o2_sd","IBCSO_v2_2km_geomorph")]
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

pred_df <- readRDS(file.path(model_dir, "prediction_grid_lookup.rds"))

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
## RUN ALL MODELS FOR PA AND AB
#############################################################
#############################################################
for(nm in model_ids[6:8]) { #c(12,1:8)
  ## parallel processing: PER CELL that contains values
  library(doParallel)
  library(foreach)
  library(Hmsc)
  parallel::detectCores()
  #UseCores = parallel::detectCores() - 1
  UseCores = 12
  c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
  registerDoParallel(c1)
  getDoParWorkers()
  
  message("====================================")
  message("Running prediction for model: ", nm)
  message("====================================")
  
  model_file <- file.path(model_dir_local, paste0(res, "_model_cells_", nm, "_chains_", nChains, "_thin_", thin, "_samples_", samples, ".Rdata"))
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
  
  #############################################################
  ## IDENTIFY TILES STILL NEEDING PREDICTION
  #############################################################
  dat.name <- paste0(
    pred_dir, res, "_model_cells_", nm,
    "_chains_", as.character(nChains),
    "_thin_", as.character(thin),
    "_samples_", as.character(samples),
    "_pred_"
  )
  
  expected.files <- paste0(
    dat.name,
    modelspec,
    sprintf("%06d", cell.sel.v),
    ".Rdata"
  )
  
  tiles.to.run <- which(!file.exists(expected.files))
  
  message("Total tiles: ", length(cell.sel.v))
  message("Completed tiles: ", length(cell.sel.v) - length(tiles.to.run))
  message("Tiles still to run: ", length(tiles.to.run))
  
  if (length(tiles.to.run) == 0) {
    message("All tiles already completed for model: ", nm)
    parallel::stopCluster(cl = c1)
    next
  }
  
  ptm = proc.time()
  foreach(j = tiles.to.run, .packages = c("Hmsc"), .errorhandling = "pass") %dopar% {#3:length(xmin)
    
    run.name <- sprintf("%06d", cell.sel.v[j])
    out.file <- paste0(dat.name, modelspec, run.name, ".Rdata")
    err.file <- paste0(dat.name, modelspec, run.name, "_ERROR.txt")
    
    if (file.exists(out.file)) {
      message("Skipping existing tile: ", run.name)
      return(NULL)
    }

    tryCatch({
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
      predY.loop.pa <- predict(pa, Gradient=Gradient.pa, expected=TRUE) ## TRUE gives probabilities instead of integer outcomes
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
      predY.pa.ci.lwr95 <- apply(predY.loop.pa.array, 1:2, lwr95)
      predY.pa.ci.upr95 <- apply(predY.loop.pa.array, 1:2, upr95)
      predY.ab.se       <- apply(predY.loop.ab.array, 1:2, se)
      predY.ab.ci.lwr95 <- apply(predY.loop.ab.array, 1:2, lwr95)
      predY.ab.ci.upr95 <- apply(predY.loop.ab.array, 1:2, upr95)
      
      save(predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.ci.lwr95, predY.ab.ci.upr95,
           predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.ci.lwr95, predY.pa.ci.upr95,
           global_row, cell_id.loop, sel.loop, XData.grid.loop, xy.grid.loop,
           file=out.file)

      rm(predY.loop.pa, predY.loop.ab, predY.loop.pa.array, predY.loop.ab.array,
         predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.ci.lwr95, predY.ab.ci.upr95,
         predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.ci.lwr95, predY.pa.ci.upr95)
      
      gc()
      
    }, error = function(e) {
      writeLines(c(paste0("model: ", nm), paste0("tile_index_j: ", j), paste0("tile_file_id: ", cell.sel.v[j]), paste0("error: ", conditionMessage(e))), con = err.file)
      return(NULL)
    })
  }
  computational.time = proc.time() - ptm
  message("Finished model: ", nm)
  parallel::stopCluster(cl = c1)
}

# rclone copy ~/local_scratch/4_model_prediction/pred_files \
# dropbox:Data/4_model_prediction/pred_files \
# --progress



#############################################################
## RUN NPP_and_FAM model for richness and broad taxa predictions
#############################################################
#############################################################
for(nm in model_ids[9]) { 
  ## parallel processing: PER CELL that contains values
  library(doParallel)
  library(foreach)
  library(Hmsc)
  parallel::detectCores()
  #UseCores = parallel::detectCores() - 1
  UseCores = 12
  c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
  registerDoParallel(c1)
  getDoParWorkers()
  
  message("====================================")
  message("Running prediction for model: ", nm)
  message("====================================")
  
  model_file <- file.path(model_dir_local, paste0(res, "_model_cells_", nm, "_chains_", nChains, "_thin_", thin, "_samples_", samples, ".Rdata"))
  load(model_file)
  pa <- models$mENV
  rm(models)
  
  message("setup complete, parallel loops over prediction areas now")
  
  #############################################################
  ## IDENTIFY TILES STILL NEEDING PREDICTION
  #############################################################
  dat.name <- paste0(
    pred_dir, res, "_model_cells_", nm,
    "_chains_", as.character(nChains),
    "_thin_", as.character(thin),
    "_samples_", as.character(samples),
    "_biodiv_pred_"
  )
  
  expected.files <- paste0(
    dat.name,
    modelspec,
    sprintf("%06d", cell.sel.v),
    ".Rdata"
  )
  
  tiles.to.run <- which(!file.exists(expected.files))
  
  message("Total tiles: ", length(cell.sel.v))
  message("Completed tiles: ", length(cell.sel.v) - length(tiles.to.run))
  message("Tiles still to run: ", length(tiles.to.run))
  
  if (length(tiles.to.run) == 0) {
    message("All tiles already completed for model: ", nm)
    parallel::stopCluster(cl = c1)
    next
  }
  
  ptm = proc.time()
  foreach(j = tiles.to.run, .packages = c("Hmsc"), .errorhandling = "pass") %dopar% {#3:length(xmin)
    
    run.name <- sprintf("%06d", cell.sel.v[j])
    out.file <- paste0(dat.name, modelspec, run.name, ".Rdata")
    err.file <- paste0(dat.name, modelspec, run.name, "_ERROR.txt")
    
    if (file.exists(out.file)) {
      message("Skipping existing tile: ", run.name)
      return(NULL)
    }
    
    tryCatch({
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
      predY.loop.pa.int <- predict(pa, Gradient=Gradient.pa, expected=FALSE) ## TRUE gives probabilities instead of integer outcomes
      rm(Gradient.pa)

      #############################################################
      ## BIODIVERSITY METRICS FROM JOINT POSTERIOR
      ## Process posterior draw list directly to reduce memory use
      #############################################################
      n_draws <- length(predY.loop.pa.int)
      n_sites <- nrow(predY.loop.pa.int[[1]])
      
      richness_samples <- matrix(NA_real_, nrow = n_sites, ncol = n_draws)
      
      asc_mat <- matrix(FALSE, nrow = n_sites, ncol = n_draws)
      bry_mat <- matrix(FALSE, nrow = n_sites, ncol = n_draws)
      oct_mat <- matrix(FALSE, nrow = n_sites, ncol = n_draws)
      por_mat <- matrix(FALSE, nrow = n_sites, ncol = n_draws)
      
      for (s in seq_len(n_draws)) {
        y <- predY.loop.pa.int[[s]]
        richness_samples[, s] <- rowSums(y)
        asc_mat[, s] <- rowSums(y[, 1:7,   drop = FALSE]) > 0
        bry_mat[, s] <- rowSums(y[, 8:18,  drop = FALSE]) > 0
        oct_mat[, s] <- rowSums(y[, 20:28, drop = FALSE]) > 0
        por_mat[, s] <- rowSums(y[, 65:79, drop = FALSE]) > 0
      }
      
      richness_mean     <- rowMeans(richness_samples)
      richness_median   <- apply(richness_samples, 1, median)
      richness_se       <- apply(richness_samples, 1, se)
      richness_ci.lwr95 <- apply(richness_samples, 1, lwr95)
      richness_ci.upr95 <- apply(richness_samples, 1, upr95)
      
      summarise_group <- function(x) {
        list(
          mean     = rowMeans(x),
          median   = apply(x, 1, median),
          se       = apply(x, 1, se),
          ci.lwr95 = apply(x, 1, lwr95),
          ci.upr95 = apply(x, 1, upr95)
        )
      }
      
      asc <- summarise_group(asc_mat)
      bry <- summarise_group(bry_mat)
      oct <- summarise_group(oct_mat)
      por <- summarise_group(por_mat)

      # #############################################################
      # ## BIODIVERSITY METRICS FROM JOINT POSTERIOR (needs expected=FALSE), slower by ???
      # #############################################################
      # # Convert matrices
      # arr = simplify2array(predY.loop.pa.int)
      # # dimensions: [sites, species, samples]
      # ## ---- 1. SPECIES RICHNESS (ALL 83 SPECIES) ----
      # richness_samples <- apply(arr, c(1,3), sum)   # [sites, samples]
      # richness_mean   <- rowMeans(richness_samples)
      # richness_median <- apply(richness_samples, 1, median)
      # richness_se     <- apply(richness_samples, 1, se)
      # richness_ci.lwr95      <- apply(richness_samples, 1, lwr95)
      # richness_ci.upr95     <- apply(richness_samples, 1, upr95)
      # 
      # ## ---- 2. FUNCTION TO COMPUTE GROUP "ANY" ----
      # get_group_prob <- function(arr_subset) {
      #   any_mat <- apply(arr_subset, c(1,3), function(x) any(x == 1))  # [sites, samples]
      #   list(
      #     mean   = rowMeans(any_mat),
      #     median = apply(any_mat, 1, median),
      #     se     = apply(any_mat, 1, se),
      #     ci.lwr95     = apply(any_mat, 1, lwr95),
      #     ci.upr95    = apply(any_mat, 1, upr95)
      #   )
      # }
      # 
      # ## ---- 3. TAXA GROUP PROBABILITIES ----
      # # Ascidians (1–7)
      # asc <- get_group_prob(arr[, 1:7, , drop = FALSE])
      # # Bryozoans (8–18)
      # bry <- get_group_prob(arr[, 8:18, , drop = FALSE])
      # # Octocorals (20–28)
      # oct <- get_group_prob(arr[, 20:28, , drop = FALSE])
      # # Porifera (65–79)
      # por <- get_group_prob(arr[, 65:79, , drop = FALSE])

      save(richness_mean, richness_median, richness_se, richness_ci.lwr95, richness_ci.upr95,
           asc, bry, oct, por,
           global_row, cell_id.loop, sel.loop, XData.grid.loop, xy.grid.loop,
           file=out.file)
      rm(predY.loop.pa.int, richness_samples,
         asc_mat, bry_mat, oct_mat, por_mat, asc, bry, oct, por,
         richness_mean, richness_median, richness_se, richness_ci.lwr95, richness_ci.upr95)
      
      gc()
    }, error = function(e) {
      writeLines(c(paste0("model: ", nm), paste0("tile_index_j: ", j), paste0("tile_file_id: ", cell.sel.v[j]), paste0("error: ", conditionMessage(e))), con = err.file)
      return(NULL)
    })
  }
  computational.time = proc.time() - ptm
  message("Finished model: ", nm)
  parallel::stopCluster(cl = c1)
}


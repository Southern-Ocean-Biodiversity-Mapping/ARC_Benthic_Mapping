#############################################################
##### HMSC TILE → RASTER ASSEMBLY (UPDATED PIPELINE)
#############################################################

library(terra)

usr <- "VM"
source("0_SourceFile.R")

pred_dir  <- paste0(usr.main.dir, "4_model_prediction/pred_files")
model_dir <- paste0(usr.dropbox.dir, "data_products/modelling_files/circum_antarctic")
env_dir   <- paste0(usr.dropbox.dir, "data_environmental/derived")
map_dir <- paste0(usr.dropbox.dir, "data_products/predictive_maps/circum_antarctic")

res <- "2km"

thin    <- 10
samples <- 800
nChains <- 4

modelspec <- "_envonly"

model_ids <- c(
  "npp_cafe","npp_cbpm","npp_eppl","npp_vpmg",
  "fam_cafe","fam_cbpm","fam_eppl","fam_vpmg"
)

#############################################################
##### LOAD REFERENCE RASTER + GRID INDEXING
#############################################################

r.stack <- rast(file.path(env_dir,paste0("Circumpolar_EnvData_", res, "_shelf_mask_unscaled_variables.tif")))
r2 <- r.stack$depth
empty.ra <- rast(r2)
empty.ra[] <- NA

# Load lookup object
load(file.path(model_dir,paste0("4_model_prediction/hmsc_", res, "_model_cell_sel.Rdata")))
load(file.path(model_dir,paste0("4_model_prediction/hmsc_", res, "_model_cell_grid.Rdata")))

#############################################################
##### MAIN LOOP OVER MODELS
#############################################################

for (nm in model_ids) {
  message("==========================================")
  message("Assembling biodiversity maps for: ", nm)
  message("==========================================")
  
  #############################################################
  ##### IDENTIFY TILE FILES
  #############################################################
  pred.files <- list.files(pred_dir, pattern = paste0("model_cells_", nm, ".*\\.Rdata$"), full.names = TRUE)
  if (length(pred.files) == 0) {
    warning("No prediction files found for: ", nm)
    next
  }
  message("Number of tiles found: ", length(pred.files))
  
  #############################################################
  ##### INITIALISE STORAGE
  #############################################################
  all.sel.ra <- NULL
  all.pred.pa.mean   <- NULL
  all.pred.pa.median <- NULL
  all.pred.pa.se     <- NULL
  all.pred.pa.5      <- NULL
  all.pred.pa.95     <- NULL
  
  all.pred.ab.mean   <- NULL
  all.pred.ab.median <- NULL
  all.pred.ab.se     <- NULL
  all.pred.ab.5      <- NULL
  all.pred.ab.95     <- NULL
  
  #############################################################
  ##### LOOP THROUGH TILE FILES
  #############################################################
  for (i in seq_along(pred.files)) {
    if (i %% 50 == 0) message("Processing tile ", i)
    
    load(pred.files[i])  # loads tile outputs
    
    ## map tile rows → raster cell indices
    ## (which cells do we need to fill with data)
    # row.numbers <- as.numeric(rownames(xy.grid[sel.loop, ]))
    # this.sel.ra <- sel.not.na[sel[row.numbers]]
    this.sel.ra <- cell_id.loop
    
    #############################################################
    ## append indices
    #############################################################
    all.sel.ra <- c(all.sel.ra, this.sel.ra)
    
    #############################################################
    ## append PA
    #############################################################
    if (is.null(all.pred.pa.mean)) {
      all.pred.pa.mean   <- predY.pa.mean
      all.pred.pa.median <- predY.pa.median
      all.pred.pa.se     <- predY.pa.se
      all.pred.pa.5      <- predY.pa.5
      all.pred.pa.95     <- predY.pa.95
      all.pred.ab.mean   <- predY.ab.mean
      all.pred.ab.median <- predY.ab.median
      all.pred.ab.se     <- predY.ab.se
      all.pred.ab.5      <- predY.ab.5
      all.pred.ab.95     <- predY.ab.95
      
    } else {
      all.pred.pa.mean   <- rbind(all.pred.pa.mean, predY.pa.mean)
      all.pred.pa.median <- rbind(all.pred.pa.median, predY.pa.median)
      all.pred.pa.se     <- rbind(all.pred.pa.se, predY.pa.se)
      all.pred.pa.5      <- rbind(all.pred.pa.5, predY.pa.5)
      all.pred.pa.95     <- rbind(all.pred.pa.95, predY.pa.95)
      all.pred.ab.mean   <- rbind(all.pred.ab.mean,   data.frame(sp1=predY.ab.mean))
      all.pred.ab.median <- rbind(all.pred.ab.median, data.frame(sp1=predY.ab.median))
      all.pred.ab.se     <- rbind(all.pred.ab.se,     data.frame(sp1=predY.ab.se))
      all.pred.ab.5      <- rbind(all.pred.ab.5,      data.frame(sp1=predY.ab.5))
      all.pred.ab.95     <- rbind(all.pred.ab.95,     data.frame(sp1=predY.ab.95))
    }
    
    rm(predY.pa.mean, predY.pa.median, predY.pa.se, predY.pa.5, predY.pa.95,
       predY.ab.mean, predY.ab.median, predY.ab.se, predY.ab.5, predY.ab.95,
       sel.loop)
  }
  
  #############################################################
  ##### BUILD RASTERS — PRESENCE/ABSENCE
  #############################################################
  out_dir <- file.path(map_dir, paste0("hmsc_with_", nm))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # order everything by raster cell index
  ord <- order(all.sel.ra)
  all.sel.ra        <- all.sel.ra[ord]
  all.pred.pa.mean   <- all.pred.pa.mean[ord, , drop = FALSE]
  all.pred.pa.median <- all.pred.pa.median[ord, , drop = FALSE]
  all.pred.pa.se     <- all.pred.pa.se[ord, , drop = FALSE]
  all.pred.pa.5      <- all.pred.pa.5[ord, , drop = FALSE]
  all.pred.pa.95     <- all.pred.pa.95[ord, , drop = FALSE]
  
  # same for abundance
  all.pred.ab.mean   <- all.pred.ab.mean[ord, , drop = FALSE]
  all.pred.ab.median <- all.pred.ab.median[ord, , drop = FALSE]
  all.pred.ab.se     <- all.pred.ab.se[ord, , drop = FALSE]
  all.pred.ab.5      <- all.pred.ab.5[ord, , drop = FALSE]
  all.pred.ab.95     <- all.pred.ab.95[ord, , drop = FALSE]  

  message("Writing PA rasters...")
  
  for (j in seq_len(ncol(all.pred.pa.mean))) {
    if (j %% 10 == 0) message("Species ", j)
    pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
    
    valid <- !is.na(all.sel.ra)
    
    pred.ra[[1]][all.sel.ra[valid]] <- all.pred.pa.mean[valid, j]
    pred.ra[[2]][all.sel.ra[valid]] <- all.pred.pa.median[valid, j]
    pred.ra[[3]][all.sel.ra[valid]] <- all.pred.pa.se[valid, j]
    pred.ra[[4]][all.sel.ra[valid]] <- all.pred.pa.5[valid, j]
    pred.ra[[5]][all.sel.ra[valid]] <- all.pred.pa.95[valid, j]

    sp.name <- colnames(all.pred.pa.mean)[j]
    names(pred.ra) <- paste0(sp.name, c("_mean","_median","_se","_5","_95"))
    
    writeRaster(pred.ra,
      filename = file.path(out_dir, paste0("PA_", sp.name, ".tif")),
      overwrite = TRUE
    )
  }
  
  #############################################################
  ##### BUILD RASTERS — TOTAL ABUNDANCE
  #############################################################
  message("Writing abundance raster...")
  pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
  
  pred.ra[[1]][all.sel.ra] <- rowSums(all.pred.ab.mean)
  pred.ra[[2]][all.sel.ra] <- rowSums(all.pred.ab.median)
  pred.ra[[3]][all.sel.ra] <- rowSums(all.pred.ab.se)
  pred.ra[[4]][all.sel.ra] <- rowSums(all.pred.ab.5)
  pred.ra[[5]][all.sel.ra] <- rowSums(all.pred.ab.95)
  
  names(pred.ra) <- paste0("total_abundance", c("_mean","_median","_se","_5","_95"))
  
  writeRaster(pred.ra,
    filename = file.path(out_dir, "total_abundance.tif"),
    overwrite = TRUE
  )
  message("Finished model: ", nm)
  
  #############################################################
  ##### CLEAN MEMORY (IMPORTANT)
  #############################################################
  rm(all.sel.ra,
     all.pred.pa.mean, all.pred.pa.median, all.pred.pa.se, all.pred.pa.5, all.pred.pa.95,
     all.pred.ab.mean, all.pred.ab.median, all.pred.ab.se, all.pred.ab.5, all.pred.ab.95)
  
  gc()
}

message("ALL MODELS COMPLETE")


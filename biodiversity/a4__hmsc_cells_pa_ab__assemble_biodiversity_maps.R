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

for (nm in model_ids[2:8]) {
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
  out_dir <- file.path(map_dir, nm)
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











































































## fitting an hmsc using  Otsos book, course scripts and : https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13345&file=mee313345-sup-0002-AppendixS2.pdf
library(Hmsc)
library(terra)
library(bayesplot)
'%!in%' <- function(x,y)!('%in%'(x,y))

###############################
#res <- "500m"
res <- "2km"

## select which model to run
model.sel <- 1 ## full model
#model.sel <- 2 ## environment only model

if(model.sel==2){
  modelspec <- "_envonly"
} else modelspec <- ""

###############################

## specify model to load
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 800 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
nChains = 4

## presence absence model
modeltype = 1
model = 6
load(paste0("~/",res,"_model_cells_",model,"_pa_chains_4_thin_10_samples_800.Rdata"))
pa <- models[[model.sel]]
load(paste0("~/",res,"_model_cells_",model,"_pa_chains_4_thin_10_samples_800_MF.Rdata"))
pa.MF <- MF[[model.sel]]

## model fit
## cross validation
#load("~/2km_model_1_abundance_chains_4_thin_10_samples_800_2foldcv.Rdata")

## reorder species to alphabetical
sp.v <- order(pa$spNames)

#### environmental stuff
env.derived <- "/pvol/data_environmental/"
r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth
empty.ra <- rast(r2)
empty.ra[] <- NA

#### load model stuff
# m <- models[[1]]
load(paste0("/pvol/biodiversity_prediction/hmsc_",res,"_model_cell_sel.Rdata"))
load(paste0("/pvol/biodiversity_prediction/hmsc_",res,"_model_cell_grid.Rdata"))

# load(file=paste0("/pvol/biodiversity_prediction/",res,"_model_50km_cells_with_data.Rdata"))
load(file=paste0("/pvol/biodiversity_prediction/",res,"_model_100km_cells_with_data.Rdata"))

## identify prediction files that we need to load
pred.dir <- paste0("/pvol/biodiversity_prediction/pred_files",modelspec,"/",res,"/abundance/")
pred.files.raw <- list.files(pred.dir)
pred.files <- pred.files.raw[grep("Rdata", pred.files.raw)]

pred.files <- pred.files[grep("800_pred",pred.files)]
pred.files <- pred.files[grep("abundance",pred.files)]
# pred.files.pa <- pred.files[grep("pa_only",pred.files)]
pred.files.pa <- pred.files
pred.files.ab <- pred.files
# pred.files.ab <- pred.files[-grep("pa_only",pred.files)]

#########################################

###### PRESENCE/ABSENCE DATA:

## first define the character string to save the files:
base.str <- paste0("/pvol/biodiversity_prediction/",res,
                   "_model_cells_", as.character(model), "_",
                   c("pa","abundance")[modeltype], 
                   "_chains_",as.character(nChains),
                   "_thin_", ... = as.character(thin),
                   "_samples_", as.character(samples),modelspec,"_")

#### load all prediction files into single files, starting with the first file and then loop
#### THIS CAN BE DONE MUCH FASTER BY SPLITTING THE JOBS RATHER THAN ONE LONG RUN (3h)
load(paste0(pred.dir,pred.files.pa[1]))
## which cells do we need to fill with data
row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
all.sel.ra <- sel.not.na[sel[row.numbers]]
## fill objects for each value
all.predY.mean <- predY.pa.mean
all.predY.median <- predY.pa.median
all.predY.se <- predY.pa.se
all.predY.5 <- predY.pa.5
all.predY.95 <- predY.pa.95
## loop through all other species
for(i in 2:length(pred.files.pa)){
  print(i)
  load(paste0(pred.dir,pred.files.pa[i]))
  row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
  all.sel.ra <- c(all.sel.ra, sel.not.na[sel[row.numbers]])
  all.predY.mean <- rbind(all.predY.mean, predY.pa.mean)
  all.predY.median <- rbind(all.predY.median, predY.pa.median)
  all.predY.se <- rbind(all.predY.se, predY.pa.se)
  all.predY.5 <- rbind(all.predY.5, predY.pa.5)
  all.predY.95 <- rbind(all.predY.95, predY.pa.95)
}

## now fill raster cells with values from the prediction file, save mean and se output for each species
pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,1]
pred.ra[[2]][all.sel.ra]  <- all.predY.median[,1]
pred.ra[[3]][all.sel.ra]  <- all.predY.se[,1]
pred.ra[[4]][all.sel.ra]  <- all.predY.5[,1]
pred.ra[[5]][all.sel.ra]  <- all.predY.95[,1]
sp.nam <- colnames(all.predY.mean)[1]
names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
writeRaster(pred.ra, file=paste0(base.str,sp.nam,".tif"))
for(j in 2:ncol(all.predY.mean)){
  print(j)
  pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
  pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,j]
  pred.ra[[2]][all.sel.ra]  <- all.predY.median[,j]
  pred.ra[[3]][all.sel.ra]  <- all.predY.se[,j]
  pred.ra[[4]][all.sel.ra]  <- all.predY.5[,j]
  pred.ra[[5]][all.sel.ra]  <- all.predY.95[,j]
  sp.nam <- colnames(all.predY.mean)[j]
  names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
  writeRaster(pred.ra, file=paste0(base.str,sp.nam,".tif"))
}


###### ABUNDANCE DATA:

## abundance model
modeltype = 2
model = 8
load(paste0("~/",res,"_model_cells_",model,"_abundance_chains_4_thin_10_samples_800.Rdata"))
ab <- models[[model.sel]]
load(paste0("~/",res,"_model_cells_",model,"_abundance_chains_4_thin_10_samples_800_MF.Rdata"))
ab.MF <- MF[[model.sel]]

rm(models, MF)

## first define the character string to save the files:
base.str <- paste0("/pvol/biodiversity_prediction/",res,
                   "_model_cells_", as.character(model), "_",
                   c("pa","abundance")[modeltype], 
                   "_chains_",as.character(nChains),
                   "_thin_", ... = as.character(thin),
                   "_samples_", as.character(samples),modelspec,"_")

#### load all prediction files into single files, starting with the first file and then loop
#### THIS CAN BE DONE MUCH FASTER BY SPLITTING THE JOBS RATHER THAN ONE LONG RUN (3h)
load(paste0(pred.dir,pred.files.ab[1]))
## which cells do we need to fill with data
row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
all.sel.ra <- sel.not.na[sel[row.numbers]]
## fill objects for each value
all.predY.mean <- predY.mean
all.predY.median <- predY.median
all.predY.se <- predY.se
all.predY.5 <- predY.5
all.predY.95 <- predY.95
## loop through all other species
for(i in 2:length(pred.files.ab)){
  print(i)
  load(paste0(pred.dir,pred.files.ab[i]))
  row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
  all.sel.ra <- c(all.sel.ra, sel.not.na[sel[row.numbers]])
  all.predY.mean <- rbind(all.predY.mean, predY.mean)
  all.predY.median <- rbind(all.predY.median, predY.median)
  all.predY.se <- rbind(all.predY.se, predY.se)
  all.predY.5 <- rbind(all.predY.5, predY.5)
  all.predY.95 <- rbind(all.predY.95, predY.95)
}

## now fill raster cells with values from the prediction file, save mean and sd output for each species
pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,1]
pred.ra[[2]][all.sel.ra]  <- all.predY.median[,1]
pred.ra[[3]][all.sel.ra]  <- all.predY.se[,1]
pred.ra[[4]][all.sel.ra]  <- all.predY.5[,1]
pred.ra[[5]][all.sel.ra]  <- all.predY.95[,1]
sp.nam <- colnames(all.predY.mean)[1]
names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
writeRaster(pred.ra, file=paste0(base.str,sp.nam,".tif"))
for(j in 2:ncol(all.predY.mean)){
  print(j)
  pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
  pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,j]
  pred.ra[[2]][all.sel.ra]  <- all.predY.median[,j]
  pred.ra[[3]][all.sel.ra]  <- all.predY.se[,j]
  pred.ra[[4]][all.sel.ra]  <- all.predY.5[,j]
  pred.ra[[5]][all.sel.ra]  <- all.predY.95[,j]
  sp.nam <- colnames(all.predY.mean)[j]
  names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
  writeRaster(pred.ra, file=paste0(base.str,sp.nam,".tif"))
}

##### Setting up----
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
##############################################################################################################
##############################################################################################################
## load data (generated in "biodiversity_models_a_prep.R" in folder ARC_Benthic_Mapping/biodiversity)
biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")
load(file=paste0(biodiv.dir,"biodiversity_bio_dat.Rdata"))
load(file=paste0(biodiv.dir,"biodiversity_env_dat.Rdata"))

dat.hmsc <- cbind(dat_pa_clean, cell_metadata_env_clean[,c(21,22,43,69:85)])

##############################################################################################################
covars <- "depth + depth2 +
           waom4k_seafloorcurrents_mean+
           waom4k_seafloorcurrents_residual+
           waom4k_seafloortemperature+
           waom4k_seafloorsalinity+
           waom4k_test_susp08+
           waom4k_test_settle08+
           distance2canyons + distance2canyons2"
###############################################################################
## fitting an hmsc using Otsos course scripts
library(Hmsc)

## set up the data
Y <- 

## set up the model

## run MCMC and save the model





## script S3 - fitting a model
#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are

load(file = "models/unfitted_models") #models, modelnames

samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nChains = 4
for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  nm = length(models)
  for (model in 1:nm) {
    print(paste0("model = ",modelnames[model]))
    m = models[[model]]
    m = sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains) 
    models[[model]] = m
  }
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  save(models,modelnames,file=filename)
}










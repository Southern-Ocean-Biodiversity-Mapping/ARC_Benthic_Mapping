# ##### Setting up----
# library(raster)
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
## load data (generated in "biodiversity_models_a_prep.R" in folder ARC_Benthic_Mapping/biodiversity)
biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")
load(file=paste0(biodiv.dir,"biodiversity_bio_dat.Rdata"))
load(file=paste0(biodiv.dir,"biodiversity_env_dat.Rdata"))

# dat.hmsc <- cbind(dat_pa_clean, cell_metadata_env_clean[,c(21,22,43,69:ncol(cell_metadata_env_clean))])

##############################################################################################################

## fitting an hmsc using Otsos course scripts and : https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13345&file=mee313345-sup-0002-AppendixS2.pdf
library(Hmsc)

##### set up the data
Y <- dat_cov_pa ## species data
XData <- cell_metadata_env_clean_scaled[,c(21,22,43,69:76,78,81:87)] ## environmental data

##### set up the model
## study design
# studyDesign <- cell_metadata_env_clean_scaled[,c(10,11)] #18
# names(studyDesign) <- c("surveyID","transectID_full")
studyDesign <- data.frame(transectID_full=cell_metadata_env_clean_scaled[,11]) #18

## match up transect coordinates with transects in cells data
sel <- match(levels(cell_metadata_env_clean_scaled$cover_cells_transect1),transect.xy$transectID_full)
transect.xy.clean <- transect.xy[sel,]
transect.xy.clean$transectID_full <- factor(transect.xy.clean$transectID_full)

## random effect structure
#surveys <- levels(studyDesign$cover_cells_survey)

#each transect needs its own x-y value
# ntransects = length(transects)
# xy = matrix(0, nrow = ntransects, ncol = 2)
# for (i in 1:ntransects){
#   rows=studyDesign[,2]==transects[[i]]
#   xy[i,1] = mean(cell_metadata_env_clean_scaled[rows,]$proj_coord_x)
#   xy[i,2] = mean(cell_metadata_env_clean_scaled[rows,]$proj_coord_y)
# }

## spatial random effect
xy <- transect.xy.clean[,2:3]
colnames(xy) = c("x","y")
sRL = xy
rownames(sRL) = levels(transect.xy.clean[,1])
rL = HmscRandomLevel(sData=sRL)
rL$nfMin = 5
rL$nfMax = 10

## random effect
rL = HmscRandomLevel(units=studyDesign$transectID_full)

XFormula = ~ depth + depth2# + slope + distance2canyons

m = Hmsc(Y = Y, XData = XData, XFormula = XFormula,
         distr = "probit",
         studyDesign = studyDesign, ranLevels = list(transectID_full=rL))

## run MCMC and save the model
thin = 1
samples = 955
nChains = 4
set.seed(1)
ptm = proc.time()
m = sampleMcmc(m, samples = samples, thin = thin,
               adaptNf = rep(ceiling(0.4*samples*thin),1),
               transient = ceiling(0.5*samples*thin),
               nChains = nChains, #nParallel = nChains,
               initPar = "fixed effects")


computational.time = proc.time() - ptm
filename = file.path(ModelDir, paste("model_", as.character(model), "_",
                                     c("pa","abundance")[modeltype], "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples), ".Rdata", sep = ""))
save(m, file=filename, computational.time)


ModelDir <- "C:/Users/jjansen/Desktop/science/SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/"



























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










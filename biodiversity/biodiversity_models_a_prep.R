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

##### load biological and environmental data
load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Data.Rdata"))
## cell_metadata, count_cells, cover_cells
## image_metadata, count_images, cover_images

load(paste0(ARC_Data.dir,"annotation/Circumpolar_Annotation_Env_Data.RData"))
## cell_metadata_env
## image_metadata_env

##############################################################################################################
## CONSIDERING ONLY CELL DATA:
##############################################################################################################

##### remove NAs from data
waom.na.sel <- which(!is.na(cell_metadata_env$waom4k_seafloorcurrents_absolute))
cell_metadata_env_clean <- cell_metadata_env[waom.na.sel,]


##### setup environmental data
cell_metadata_env$depth2 <- poly(cell_metadata_env$depth,2)[,2]
cell_metadata_env$distance2canyons2 <- poly(cell_metadata_env$distance2canyons,2)[,2]
cell_metadata_env$waom4k_seafloorcurrents_absolute2 <- NA
cell_metadata_env$waom4k_seafloorcurrents_residual2 <- NA
cell_metadata_env$waom4k_seafloorcurrents_absolute2[waom.na.sel] <- poly(cell_metadata_env$waom4k_seafloorcurrents_absolute[waom.na.sel] ,2)[,2]
cell_metadata_env$waom4k_seafloorcurrents_residual2[waom.na.sel] <- poly(cell_metadata_env$waom4k_seafloorcurrents_residual[waom.na.sel] ,2)[,2]
cell_metadata_env$slope2 <- poly(cell_metadata_env$slope ,2)[,2]
cell_metadata_env$waom4k_test_settle082 <- poly(cell_metadata_env$waom4k_test_settle08 ,2)[,2]
cell_metadata_env$waom4k_test_susp082 <- poly(cell_metadata_env$waom4k_test_susp08 ,2)[,2]


##### scale environmental data
cell_metadata_env_clean_scaled <- cell_metadata_env_clean[,-c(1:20,77,78)]
scale.means <- NA
scale.vars <- NA
for(i in 1:ncol(cell_metadata_env_clean_scaled)){
  scale.means[i] <- mean(cell_metadata_env_clean_scaled[,i], na.rm=TRUE)
  scale.vars[i] <- var(cell_metadata_env_clean_scaled[,i], na.rm=TRUE)
  cell_metadata_env_clean_scaled[,i] <- (cell_metadata_env_clean_scaled[,i]-scale.means[i])/scale.vars[i]
}

##### setup biological data
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
sel_noid.cov <- grep("NoID",dataset.names)
sel_unsc.cov <- grep("Unscorable",dataset.names)

cover_cells_pa <- cover_cells
cover_cells_pa[cover_cells_pa>0] <- 1

cover_SF.prop <- rowSums(cover_cells[,sel_SF])/n_not_na
cover_SF <- rowSums(cover_cells[,sel_SF])
cover_SF_pa <- cover_SF
cover_SF_pa[cover_SF>0] <- 1
richness <- rowSums(cover_cells_pa[,-sel_sed])#/n_total
richness.l <- rowSums(cover_cells_pa[,-sel_sed])/log(n_total)
cover_all.prop <- rowSums(cover_cells[,-sel_sed])/n_not_na
cover_all <- rowSums(cover_cells[,-sel_sed])

## names of faunal groups for count_cells:
count.names <- names(count_cells)
## selector for each faunal class
sel_noid <- grep("NoID",count.names)
sel_echino <- grep("Echinoderms",count.names)
sel_crust <- grep("Crustacea",count.names)

count_mobile <- rowSums(count_cells[,-sel_noid]) ## remove NoIDs
count_echino <- rowSums(count_cells[,-sel_noid])
count_crust <- rowSums(count_cells[,-sel_crust])

######################################################################################################

##### individual species
dat_cov_species <- cover_cells[,-c(sel_sed,sel_noid.cov,sel_unsc.cov)]


##### large species groupings
dat_sum <- data.frame(cbind(cover_all, cover_SF, cover_SF_pa, richness, n_not_na))#, cell_metadata_env[,c(6:9,18:76)])
## count data
dat_sum$count_mobile <- dat_sum$count_echino <- dat_sum$count_crust <- NA
sel <- which(is.na(dat_sum$counts_N))
dat_sum$count_mobile[-sel] <- count_mobile
dat_sum$count_echino[-sel] <- count_echino
dat_sum$count_crust[-sel]  <- count_crust

dat_sum_clean <- dat_sum[waom.na.sel,]

##### presence-absence data
dat_pa <- cover_cells_pa[,-c(sel_sed,sel_noid.cov,sel_unsc.cov)]
dat_pa_clean <- dat_pa[waom.na.sel,]


######################################################################################################
biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")

save(dat, dat_clean, dat_pa, dat_pa_clean, file=paste0(biodiv.dir,"biodiversity_bio_dat.Rdata"))
save(cell_metadata_env_clean, cell_metadata_env_clean_scaled, scale.means, scale.vars, file=paste0(biodiv.dir,"biodiversity_env_dat.Rdata"))


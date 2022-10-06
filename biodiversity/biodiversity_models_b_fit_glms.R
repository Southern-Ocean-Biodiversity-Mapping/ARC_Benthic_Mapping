##### Setting up----
library(MASS)

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

dat.glm <- cbind(dat_cov_sum, cell_metadata_env_clean[,c(21,22,43,69:ncol(cell_metadata_env_clean))])

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

## Richness
summary(fit.r <- glm.nb(formula(paste0("richness~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.r <- update(fit.r, .~.-distance2canyons2))
summary(fit.r <- update(fit.r, .~.-waom4k_seafloorsalinity))

## total cover
summary(fit.c <- glm.nb(formula(paste0("cover_all~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))

## SF-abundance
summary(fit.sf <- glm.nb(formula(paste0("cover_SF~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.sf <- update(fit.sf, .~.-waom4k_seafloortemperature))

# ## mobile counts
# summary(fit.m <- glm.nb(formula(paste0("count_mobile~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
# summary(fit.m <- update(fit.m, .~.-waom4k_seafloortemperature))
# summary(fit.m <- update(fit.m, .~.-distance2canyons))
# summary(fit.m <- update(fit.m, .~.-waom4k_seafloorcurrents_residual))
# summary(fit.m <- update(fit.m, .~.-waom4k_seafloorsalinity))
# 
# ## echinoderms
# summary(fit.e <- glm.nb(formula(paste0("abund_echino~",covars,"+offset(log(n_not_na))")), data=dat.glm))
# summary(fit.e <- update(fit.e, .~.-waom4k_seafloortemperature))
# summary(fit.e <- update(fit.e, .~.-distance2canyons))
# summary(fit.e <- update(fit.e, .~.-waom4k_seafloorcurrents_residual))
# summary(fit.e <- update(fit.e, .~.-waom4k_seafloorsalinity))
# 
# ## crustacea
# summary(fit.cr <- glm.nb(formula(paste0("abund_crust~",covars,"+offset(log(n_not_na))")), data=dat.glm))
# summary(fit.cr <- update(fit.cr, .~.-distance2canyons))
# summary(fit.cr <- update(fit.cr, .~.-waom4k_seafloorsalinity))
# summary(fit.cr <- update(fit.cr, .~.-waom4k_seafloortemperature))
# summary(fit.cr <- update(fit.cr, .~.-waom4k_seafloorcurrents_residual))

#######################################################
cor(cell_metadata_env[,c(21,22,43,70:76)])
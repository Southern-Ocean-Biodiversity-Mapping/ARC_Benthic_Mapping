# ##### Setting up----
library(terra)
library(PerformanceAnalytics) ## plotting correlations
library(MASS)    ## glm.nb
library(glmmTMB) ## model fitting for glmm
library(DHARMa)  ## model diagnostics https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(mgcv)    ## gams

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

biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity")
distr.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/species_distributions")

##############################################################################################################

#res <- "500m"
res <- "2km"

##############################################################################################################
## load scaled environmental rasters:
env_stack_scaled <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled.tif"))

## load data
load(file=paste0(ARC_Data.dir,"Cell_level_env_",res,".Rdata"))
load(file=paste0(ARC_Data.dir,"Cell_level_bio_2pc_",res,".Rdata"))

r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth

#### check for NAs in the data:
## find NAs in waom data
waom.na.sel <- which(!is.na(cell_metadata_env_scaled$seafloorcurrents_absolute))
metadat <- cell_metadata_env_scaled[waom.na.sel,]
metadat$cellID <- factor(metadat$cellID)
cover_cells <- cover_mod[waom.na.sel,]
cover_groups <- cover_groupings[waom.na.sel,]
rm(cover_mod, cover_groupings, cell_metadata_env_scaled)

## presence absence data:
cov_pa.raw <- cover_cells[,-1]
cov_pa.raw[cov_pa.raw>0] <- 1

## remove rare species
cov_pa.raw2 <- cov_pa.raw[,-which(colSums(cov_pa.raw)<=18)]

## remove substrates, noid and unscorable
if(res=="2km") cov_pa <- cov_pa.raw2[,-c(grep("Sub",names(cov_pa.raw2)),grep("Unsco",names(cov_pa.raw2)),grep("NoID",names(cov_pa.raw2))[1:2])]

## remove NoIDs???


##############################################################################################################

dat.glm <- cbind(cover_groups, metadat)
dat.glm$image_quality_score <- dat.glm$image_quality_score-min(dat.glm$image_quality_score)+1
dat.glm$annotated_area <- dat.glm$cover_area * (dat.glm$cover_points_scorable/dat.glm$cover_points_total)


covars <- "depth + depth2 + tpi11 + logslope + tpi + seafloorcurrents_mean + seafloorcurrents_residual + seafloorcurrents_absolute + seafloortemperature + seafloorsalinity + test_settle08 + test_susp08 + distance2canyons + distance2canyons2"

## GLMs
summary(fit.glm.r <- glm.nb(formula(paste0("richness~",covars,"+log(annotated_area)")), data=dat.glm))
summary(fit.glm.r <- update(fit.glm.r, .~.-test_susp08))
summary(fit.glm.r <- update(fit.glm.r, .~.-seafloorcurrents_absolute))
summary(fit.glm.r <- update(fit.glm.r, .~.-seafloortemperature))
summary(fit.glm.r <- update(fit.glm.r, .~.-tpi11))
summary(fit.glm.r <- update(fit.glm.r, .~.-test_settle08))
summary(fit.glm.r <- update(fit.glm.r, .~.-tpi))

fit.sim <- simulateResiduals(fit.glm.r)
plot(fit.sim)

## GLMMs
summary(fit.glmm.r <- glmmTMB(richness ~ depth + depth2 + #tpi11 +
                                logslope + tpi +
                                seafloorcurrents_mean +
                                seafloorcurrents_residual +
                                #seafloorcurrents_absolute +
                                seafloortemperature +
                                #seafloorsalinity +
                                #test_settle08 +
                                #test_susp08 +
                                distance2canyons +
                                #distance2canyons2 +
                                log(annotated_area) +
                                (1|gear) + (1|cover_cells_transect1) + (1|cover_cells_survey),
                              data=dat.glm, family=nbinom1))
# summary(fit <- update(fit, .~.-seafloorcurrents_absolute))
# summary(fit <- update(fit, .~.-seafloorsalinity))
# summary(fit <- update(fit, .~.-tpi11))
# summary(fit <- update(fit, .~.-test_settle08))
# summary(fit <- update(fit, .~.-test_susp08))

fit.sim <- simulateResiduals(fit.glmm.r)
plot(fit.sim)

par(mfrow=c(3,3))
plotResiduals(fit.sim, dat.glm$depth[-c(70,71)])
plotResiduals(fit.sim, dat.glm$logslope[-c(70,71)])
plotResiduals(fit.sim, dat.glm$tpi11[-c(70,71)])
plotResiduals(fit.sim, dat.glm$waom4k_seafloorcurrents_mean[-c(70,71)])
plotResiduals(fit.sim, dat.glm$waom4k_seafloortemperature[-c(70,71)])
plotResiduals(fit.sim, dat.glm$waom4k_seafloorsalinity[-c(70,71)])
plotResiduals(fit.sim, dat.glm$waom4k_test_settle08[-c(70,71)])
plotResiduals(fit.sim, log(dat.glm$annotated_area)[-c(70,71)])

save(fit.glm.r, fit.glmm.r, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_richness.Rdata"))






###############

##


############
##### modelling richness
# plot(richness~annotated_area, data=dat.glm)
# m <- gam(richness~s(annotated_area), data=dat.glm)
# p_obj <- plot(m, residuals = TRUE)
# p_obj <- p_obj[[1]] # just one smooth so select the first component
# sm_df <- as.data.frame(p_obj[c("x", "se", "fit")])
# data_df <- as.data.frame(p_obj[c("raw", "p.resid")])
# ## plot
# ggplot(sm_df, aes(x = x, y = fit)) +
#   geom_rug(data = data_df, mapping = aes(x = raw, y = NULL),
#            sides = "b") +
#   geom_point(data = data_df, mapping = aes(x = raw, y = p.resid)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
#               alpha = 0.3) +
#   geom_line() +
#   labs(x = p_obj$xlab, y = p_obj$ylab)
# 
# plot(richness~cover_points_scorable, data=dat.glm)
# plot(richness~image_quality_score, data=dat.glm)
# plot(richness~cover_cells_survey, data=dat.glm)




##### modelling total cover

## total cover
summary(fit.glm.c <- stepAIC(glm.nb(formula(paste0("cover_all~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm)))

summary(fit.glmm.c <- glmmTMB(cover_all ~ #depth + #depth2 + 
                                #tpi11 +
                                logslope + #tpi +
                                #seafloorcurrents_mean +
                                #seafloorcurrents_residual +
                                #seafloorcurrents_absolute +
                                #seafloortemperature +
                                #seafloorsalinity +
                                test_settle08 +
                                #test_susp08 +
                                #distance2canyons +
                                #distance2canyons2 +
                                offset(log(cover_points_scorable)) +
                                (1|gear) + (1|cover_cells_transect1) + (1|cover_cells_survey),
                              data=dat.glm, family=nbinom1))
fit.sim <- simulateResiduals(fit.glmm.c)
plot(fit.sim)

save(fit.glm.c, fit.glmm.c, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_totalcover.Rdata"))


##### modelling SF cover

## SF-abundance
summary(fit.glm.sf <- glm.nb(formula(paste0("cover_SF~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-seafloorsalinity))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-tpi))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-tpi11))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-test_settle08))

summary(fit.glmm.sf <- glmmTMB(cover_SF ~ depth + depth2 + #tpi11 +
                                 logslope + #tpi +
                                 seafloorcurrents_mean +
                                 seafloorcurrents_residual +
                                 #seafloorcurrents_absolute +
                                 seafloortemperature +
                                 #seafloorsalinity +
                                 #test_settle08 +
                                 #test_susp08 +
                                 distance2canyons +
                                 #distance2canyons2 +
                                 offset(log(cover_points_scorable)) +
                                 (1|cover_cells_transect1) + (1|cover_cells_survey),
                               data=dat.glm, family=nbinom1))
fit.sim <- simulateResiduals(fit.glmm.sf)
plot(fit.sim)

save(fit.glm.sf, fit.glmm.sf, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_sfcover.Rdata"))


















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

##### modelling bryozoan cover

## B-abundance
summary(fit.glm.b <- glm.nb(formula(paste0("cover_B~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.glm.b <- update(fit.glm.b, .~.-waom4k_seafloorcurrents_absolute))
summary(fit.glm.b <- update(fit.glm.b, .~.-waom4k_seafloorsalinity))
summary(fit.glm.b <- update(fit.glm.b, .~.-tpi))

summary(fit.glmm.b <- glmmTMB(cover_B ~ depth + depth2 + logslope +
                                waom4k_seafloorcurrents_mean +
                                waom4k_test_susp08 +
                                waom4k_test_settle08 +
                                offset(log(cover_points_scorable)) +
                                (1|cover_cells_transect1) + (1|cover_cells_survey),
                              data=dat.glm, family=nbinom1)) #dat.glm[-c(70,71),]
fit.sim <- simulateResiduals(fit.glmm.b)
plot(fit.sim)

save(fit.glm.b, fit.glmm.b, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_bcover.Rdata"))

##### modelling sponge cover

## S-abundance
summary(fit.glm.s <- glm.nb(formula(paste0("cover_S~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.glm.s <- update(fit.glm.s, .~.-distance2canyons2))
summary(fit.glm.s <- update(fit.glm.s, .~.-waom4k_test_settle08))
summary(fit.glm.s <- update(fit.glm.s, .~.-waom4k_seafloorcurrents_absolute))

summary(fit.glmm.s <- glmmTMB(cover_S ~ depth + depth2 + logslope +
                                waom4k_seafloorcurrents_mean +
                                waom4k_seafloorcurrents_residual +
                                waom4k_seafloorsalinity +
                                waom4k_test_susp08 +
                                waom4k_test_settle08 +
                                distance2canyons +
                                offset(log(cover_points_scorable)) +
                                (1|cover_cells_transect1) + (1|cover_cells_survey),
                              data=dat.glm, family=nbinom1)) #dat.glm[-c(70,71),]
fit.sim <- simulateResiduals(fit.glmm.s)
plot(fit.sim)

save(fit.glm.s, fit.glmm.s, file=paste0(biodiv.dir,"/biodiversity_fit_",res,"_scover.Rdata"))

#######################################################
cor(cell_metadata_env[,c(21,22,43,70:76)])











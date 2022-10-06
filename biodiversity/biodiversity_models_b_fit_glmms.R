##### Setting up----
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
##############################################################################################################
##############################################################################################################
## load data (generated in "biodiversity_models_a_prep.R" in folder ARC_Benthic_Mapping/biodiversity)
biodiv.dir <- paste0(sci.dir,"SouthernOceanBiodiversityMapping/ARC_Benthic_Mapping/biodiversity/")
load(file=paste0(biodiv.dir,"biodiversity_bio_dat.Rdata"))
load(file=paste0(biodiv.dir,"biodiversity_env_dat.Rdata"))

## we can consider all non-correlated variables plus geomorph
## but that's a lot...:
length(sel.not.correlated)
names(cell_metadata_env_clean[,sel.not.correlated])
## but the cars data is patchy


#dat.glm <- cbind(dat_cov_sum, cell_metadata_env_clean[,c(8,10:11,18,20,78,79,sel.not.correlated,68,69)])
dat.glm <- cbind(dat_cov_sum, cell_metadata_env_clean_scaled[,c(8,10:11,18,20,78,79,sel.not.correlated,68,69)])
#dat.glm <- cbind(dat_cov_sum, cell_metadata_env_clean_scaled[,c(10:11,18,21:25,43,51,53,69:ncol(cell_metadata_env_clean_scaled))])
dat.glm$image_quality_score <- dat.glm$image_quality_score-min(dat.glm$image_quality_score)+1
dat.glm$annotated_area <- dat.glm$cover_area * (dat.glm$cover_points_scorable/dat.glm$cover_points_total)

##############################################################################################################
# covars <- names(dat.glm[,c(13:16,28:42)])
covars <- "depth + depth2 +
           tpi + tpi11 +
           logslope +
           waom4k_seafloorcurrents_mean +
           waom4k_seafloorcurrents_absolute +
           waom4k_seafloorcurrents_residual +
           waom4k_seafloortemperature +
           waom4k_seafloorsalinity +
           waom4k_test_susp08 +
           waom4k_test_settle08 +
           distance2canyons + distance2canyons2"

###############################################################################

##### modelling richness
# plot(richness~annotated_area, data=dat.glm)
m <- gam(richness~s(annotated_area), data=dat.glm)
p_obj <- plot(m, residuals = TRUE)
p_obj <- p_obj[[1]] # just one smooth so select the first component
sm_df <- as.data.frame(p_obj[c("x", "se", "fit")])
data_df <- as.data.frame(p_obj[c("raw", "p.resid")])
## plot
ggplot(sm_df, aes(x = x, y = fit)) +
  geom_rug(data = data_df, mapping = aes(x = raw, y = NULL),
           sides = "b") +
  geom_point(data = data_df, mapping = aes(x = raw, y = p.resid)) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 0.3) +
  geom_line() +
  labs(x = p_obj$xlab, y = p_obj$ylab)
# 
# plot(richness~cover_points_scorable, data=dat.glm)
# plot(richness~image_quality_score, data=dat.glm)
# plot(richness~cover_cells_survey, data=dat.glm)

## glm
summary(fit.glm.r <- glm.nb(formula(paste0("richness~",covars,"+log(annotated_area)")), data=dat.glm))
summary(fit.glm.r <- update(fit.glm.r, .~.-waom4k_seafloorcurrents_absolute))
summary(fit.glm.r <- update(fit.glm.r, .~.-waom4k_test_settle08))
summary(fit.glm.r <- update(fit.glm.r, .~.-tpi))
summary(fit.glm.r <- update(fit.glm.r, .~.-waom4k_seafloortemperature))
summary(fit.glm.r <- update(fit.glm.r, .~.-waom4k_test_susp08))

fit.sim <- simulateResiduals(fit.glm.r)
plot(fit.sim)


## fitting a glmm
## need to account for sampling effort
## can either do log(area) as a covariate, or use log(log(area+1)) as an offset
## Scott's suggestion was to use it as a covariate
# summary(fit <- glmmTMB(formula(paste0("richness~",paste(covars, collapse="+"),
#                                       "+offset(log(log(annotated_area+1)))",
#                                       # "+offset(log(image_quality_score))",
#                                       "+(1|gear)",
#                                       "+(1|cover_cells_transect1)",
#                                       "+(1|cover_cells_survey)"
# )), data=dat.glm[-c(70,71),], family=nbinom1))
# fit.sim <- simulateResiduals(fit)
# plot(fit.sim)
# summary(fit <- update(fit, .~.-waom4k_seafloorcurrents_residual))
# summary(fit <- update(fit, .~.-tpi))
# summary(fit <- update(fit, .~.-waom4k_seafloorcurrents_absolute))
# summary(fit <- update(fit, .~.-waom4k_test_susp08))
summary(fit.glmm.r <- glmmTMB(richness ~ depth + depth2 + tpi11 + logslope +
                         waom4k_seafloorcurrents_mean +
                         waom4k_seafloortemperature +
                         waom4k_seafloorsalinity +
                         waom4k_test_settle08 +
                         distance2canyons +
                         log(annotated_area) +
                         (1|gear) + (1|cover_cells_transect1) + (1|cover_cells_survey),
                       data=dat.glm[-c(70,71),], family=nbinom1))
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

save(fit.glm.r, fit.glmm.r, file=paste0(biodiv.dir,"biodiversity_fit_richness.Rdata"))



##### modelling total cover

## total cover
summary(fit.glm.c <- glm.nb(formula(paste0("cover_all~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.glm.c <- update(fit.glm.c, .~.-tpi))
summary(fit.glm.c <- update(fit.glm.c, .~.-tpi11))
summary(fit.glm.c <- update(fit.glm.c, .~.-waom4k_seafloorcurrents_absolute))

summary(fit.glmm.c <- glmmTMB(cover_all ~ depth + depth2 + tpi11 + logslope +
                         waom4k_seafloorcurrents_mean +
                         waom4k_seafloorcurrents_residual +
                         waom4k_seafloortemperature +
                         waom4k_seafloorsalinity +
                         waom4k_test_susp08 +
                         waom4k_test_settle08 +
                         distance2canyons + distance2canyons2 +
                         offset(log(cover_points_scorable)) +
                         (1|gear) + (1|cover_cells_transect1) + (1|cover_cells_survey),
                       data=dat.glm[-c(70,71),], family=nbinom1))
summary(fit.glmm.c <- update(fit.glmm.c, .~.-distance2canyons2))
summary(fit.glmm.c <- update(fit.glmm.c, .~.-distance2canyons))
summary(fit.glmm.c <- update(fit.glmm.c, .~.-waom4k_seafloorcurrents_residual))
summary(fit.glmm.c <- update(fit.glmm.c, .~.-waom4k_seafloortemperature))
summary(fit.glmm.c <- update(fit.glmm.c, .~.-tpi11))
summary(fit.glmm.c <- update(fit.glmm.c, .~.-waom4k_seafloorsalinity))

# summary(fit.c2 <- glmmTMB(formula(paste0("cover_all~",covars,
#                                          "+offset(log(cover_points_scorable))",
#                                          "+(1|gear) +(1|cover_cells_transect1) +(1|cover_cells_survey)")),
#                        data=dat.glm[-c(70,71),], family=nbinom1))
fit.sim <- simulateResiduals(fit.glmm.c)
plot(fit.sim)

save(fit.glm.c, fit.glmm.c, file=paste0(biodiv.dir,"biodiversity_fit_totalcover.Rdata"))


##### modelling SF cover

## SF-abundance
summary(fit.glm.sf <- glm.nb(formula(paste0("cover_SF~",covars,"+offset(log(cover_points_scorable))")), data=dat.glm))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-waom4k_seafloorsalinity))
summary(fit.glm.sf <- update(fit.glm.sf, .~.-waom4k_seafloorcurrents_absolute))

summary(fit.glmm.sf <- glmmTMB(cover_SF ~ depth + depth2 + logslope +
                                 waom4k_seafloorcurrents_mean +
                                 #waom4k_seafloorcurrents_absolute +
                                 waom4k_test_settle08 +
                                 offset(log(cover_points_scorable)) +
                                 (1|cover_cells_transect1) + (1|cover_cells_survey),
                               data=dat.glm[-c(70,71),], family=nbinom1))
fit.sim <- simulateResiduals(fit.glmm.sf)
plot(fit.sim)

save(fit.glm.sf, fit.glmm.sf, file=paste0(biodiv.dir,"biodiversity_fit_sfcover.Rdata"))

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

save(fit.glm.b, fit.glmm.b, file=paste0(biodiv.dir,"biodiversity_fit_bcover.Rdata"))

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

save(fit.glm.s, fit.glmm.s, file=paste0(biodiv.dir,"biodiversity_fit_scover.Rdata"))

#######################################################
cor(cell_metadata_env[,c(21,22,43,70:76)])


































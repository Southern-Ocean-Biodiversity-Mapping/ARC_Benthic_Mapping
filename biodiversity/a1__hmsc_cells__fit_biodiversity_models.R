
##############################################################################################################
library(Hmsc)
library(dplyr)
library(terra)

'%!in%' <- function(x,y)!('%in%'(x,y))

##############################################################################################################
###########################
##### SETTINGS ###########
###########################

usr <- "VM"
source("0_SourceFile.R")

bio.dir      <- paste0(usr.dropbox.dir, "data_biological/")
output_dir   <- paste0(usr.dropbox.dir, "data_products/modelling_files/circum_antarctic")


res <- "2km"

##############################################################################################################
###########################
##### LOAD NEW DATA ######
###########################

dat <- readRDS(file.path(output_dir, paste0("cover_modelling_inputs_", res, ".rds")))

# Extract key components
df <- dat$modelling_dataframe
resp_pa <- dat$responses$presence_absence |> select(-cell_id)
resp_counts <- dat$responses$counts |> select(-cell_id)

##############################################################################################################
###########################
##### CLEAN DATA #########
###########################

# Remove rows with NA predictors
complete_rows <- complete.cases(df)
df <- df[complete_rows, ]
resp_pa <- resp_pa[complete_rows, ]
resp_counts <- resp_counts[complete_rows, ]

##############################################################################################################
###########################
##### RESPONSE DATA ######
###########################

### Presence–absence (morphospecies)
cov_pa <- as.matrix(resp_pa)

# Remove rare species using prevalence threshold
prev_threshold <- 0.02
keep_sp <- colMeans(cov_pa) > prev_threshold
cov_pa <- cov_pa[, keep_sp]

### Total abundance of all biota
cov_ab <- df$total_abundance

### Richness
cov_rich <- df$richness_raw

##
#save(df, cov_pa, cov_ab, cov_rich, file=file.path(output_dir, paste0("2km_model_cells_data.Rdata")))


##############################################################################################################
###########################
##### PREDICTORS #########
###########################
# Manual predictor selection
model_vars_fixed <- c("depth","depth2","logslope","tpi",
                      "distance2canyons","distance2canyons2",
                      "seafloortemperature","seafloorcurrents_mean",
                      "seafloorcurrents_residual","seafloorsalinity",
                      "cover_points_scorable")
model_vars_swap_npp_mean <- c("cafe_mean", "cbpm_mean", "eppl_mean", "vpmg_mean")
model_vars_swap_npp_sd   <- c("cafe_sd", "cbpm_sd", "eppl_sd", "vpmg_sd")
model_vars_swap_fam_flx <- c("log.flux.mean.cafe",
                             "log.flux.mean.cbpm",
                             "log.flux.mean.eppl",
                             "log.flux.mean.vpmg")
model_vars_swap_fam_sed <- c("sed.mean.cafe",
                             "sed.mean.cbpm",
                             "sed.mean.eppl",
                             "sed.mean.vpmg")

model_vars <- c(model_vars_fixed, 
                model_vars_swap_npp_mean, model_vars_swap_npp_sd,
                model_vars_swap_fam_flx, model_vars_swap_fam_sed)
XData <- df[, model_vars]

## build pairs for formulas later
npp_pairs <- data.frame(
  mean = model_vars_swap_npp_mean,
  sd   = model_vars_swap_npp_sd,
  stringsAsFactors = FALSE
)
fam_pairs <- data.frame(
  flux = model_vars_swap_fam_flx,
  sed  = model_vars_swap_fam_sed,
  stringsAsFactors = FALSE
)

##############################################################################################################
###########################
##### STUDY DESIGN ########
###########################

# Keep structure — assumes these exist (update upstream if needed)
studyDesign <- data.frame(
  cellID      = factor(df$cell_id),
  surveyID    = factor(df$surveyID),
  transectID  = factor(df$transectID),
  gear        = factor(df$gear),
  year        = as.factor(df$year)
)

##############################################################################################################
##### RANDOM EFFECTS
##############################################################################################################
## survey
rL.s <- HmscRandomLevel(units = levels(studyDesign$surveyID))
## transect
rL.t <- HmscRandomLevel(units = levels(studyDesign$transectID))
## gear
rL.g <- HmscRandomLevel(units = levels(studyDesign$gear))
## year
rL.y <- HmscRandomLevel(units = levels(studyDesign$year))

##############################################################################################################
##### SPATIAL RANDOM EFFECT
##############################################################################################################
## spatial random effect
xy <- df[,c("proj_coord_x","proj_coord_y")]
colnames(xy) = c("x","y")
sRL = xy
rownames(sRL) = df$cell_id

### using the standard algorithm:
## 5min per iteration to fit the spatial model, which is way too long (~1 month for 10k iterations)
#rL = HmscRandomLevel(sData=sRL)

### trying NNGP:
## doesn't work, the error message is: "Failed updaters and their counts in chain 1  ( 15  attempts)"
# rL = HmscRandomLevel(sData=sRL, sMethod="NNPG")
# rL = setPriors(rL,nfMin=1,nfMax=1)

### trying GPP:
## ~ 40s for 10 iterations; 60s for 50 iterations, 84s for 100 iterations; using knots at 500km distance -> ~1h20min for 10k iterations
## ~ 3.5min for 10 iterations; 5min for 50 iterations; 6.6min for 100 iterations; using knots at 250km distance -> ~2h40min for 10k iterations
## ~ 7min for 2 iterations; 15min for 100 iterations using knots at 200km distance -> ~13h20min for 10k iterations
## BUT, on a 250km grid, with the full dataset, the predictions will take 41 days!!!
## first specifying knots on a grid
# xy.knots <- rbind(xy,c(2900000,0)) ## add a point to the right to allow mapping of East Antarctica
# Knots = constructKnots(xy.knots, knotDist = 250000, minKnotDist = 2500000)
# #Knots = constructKnots(xy, knotDist = 50000, minKnotDist = 2000000)
# plot(xy.knots[,1],xy.knots[,2],pch=18, asp=1)
# points(Knots[,1],Knots[,2],col='red',pch=18)

use_spatial <- TRUE
if (use_spatial) {
  ## knots at 200km distance, 250km min distance:
  ## add points between AP and Ross Sea
  xy.knots <- rbind(xy,
                    c(-2143647,498436),
                    c(-1916289,355285),
                    c(-2000000,100000),
                    c(-1983655,-368890),
                    c(-1739456,-638350),
                    c(-1621567,-975176),
                    c(-1360527,-1160431),
                    c(-900000,-1300000),
                    c(-627930,-1345685))
  ## add points to East Antarctica
  xy.knots <- rbind(xy.knots,
                    c(710952,-2154067),
                    c(1115144,-2288798),
                    c(2100359,-1724614),
                    c(2529812,-1017280),
                    c(2757170,-629930),
                    c(2824535,-318366),
                    c(2672963,-57326),
                    c(2656122,254237),
                    c(2487709,498436),
                    c(2386661,742635),
                    c(2268772,1256295),
                    c(2125621,1685748),
                    c(1645644,1812057),
                    c(820421,2056256),
                    c(1200000,2000000))
  ## add a point to Weddell Sea
  xy.knots <- rbind(xy.knots,
                    c(-1503678,1054199),
                    c(-1436312,1300000),
                    c(-1958393,1298398))
  Knots = constructKnots(xy.knots, knotDist = 200000, minKnotDist = 250000)
  rL = HmscRandomLevel(sData=sRL, sMethod='GPP', sKnot=Knots)
  rL$nfMax=10
  rL.s$nfMax=10
  rL.t$nfMax=10
  rL.g$nfMax=10
  rL.y$nfMax=10
}

##############################################################################################################
##### MODEL FORMULA
##############################################################################################################
build_formula <- function(vars_fixed, vars_swap) {
  as.formula(paste("~", paste(c(vars_fixed, vars_swap), collapse = " + ")))
}
## generate 8 formulas
formulas <- list()
## NPP models (mean + sd)
for (i in seq_len(nrow(npp_pairs))) {
  key <- paste0("npp_", gsub("_mean", "", npp_pairs$mean[i]))
  formulas[[key]] <- build_formula(
    model_vars_fixed,
    c(npp_pairs$mean[i], npp_pairs$sd[i])
  )
}
## Flux / sediment models
for (i in seq_len(nrow(fam_pairs))) {
  key <- paste0("fam_", gsub("log.flux.mean.", "", fam_pairs$flux[i]))
  
  formulas[[key]] <- build_formula(
    model_vars_fixed,
    c(fam_pairs$flux[i], fam_pairs$sed[i])
  )
}

##################################################
##### RUN ANALYSIS AND SAVE OUTPUTS (~5min per non-spatial model)
##################################################
all_models <- list()
all_times  <- list()

for (nm in names(formulas)[5:8]) {
  
  message("====================================")
  message("Building models for formula: ", nm)
  message("====================================")
  
  XFormula <- formulas[[nm]]
  
  models <- list()
  
  ##################################################
  ##### BUILD MODELS
  ##################################################
  
  models$mFULL <- Hmsc(
    Y = cov_pa,
    XData = XData,
    XFormula = XFormula,
    distr = "probit",
    studyDesign = studyDesign,
    ranLevels = list(cellID = rL, surveyID = rL.s)
  )
  
  models$mENV <- Hmsc(
    Y = cov_pa,
    XData = XData,
    XFormula = XFormula,
    distr = "probit",
    studyDesign = studyDesign,
    ranLevels = list(surveyID = rL.s)
  )
  
  models$mSPACE <- Hmsc(
    Y = cov_pa,
    XData = XData,
    XFormula = ~1,
    distr = "probit",
    studyDesign = studyDesign,
    ranLevels = list(cellID = rL, surveyID = rL.s)
  )
  
  models$mAB <- Hmsc(
    Y = cov_ab,
    XData = XData,
    XFormula = XFormula,
    distr = "lognormal poisson",
    studyDesign = studyDesign,
    ranLevels = list(surveyID = rL.s)
  )
  
  models$mRICH <- Hmsc(
    Y = cov_rich,
    XData = XData,
    XFormula = XFormula,
    distr = "lognormal poisson",
    studyDesign = studyDesign,
    ranLevels = list(surveyID = rL.s)
  )
  
  models$formula_id <- nm
  models$formula_text <- deparse(XFormula)
  
  ##################################################
  ##### RUN MCMC
  ##################################################
  
  thin      <- 10
  samples   <- 800
  transient <- ceiling(0.5 * samples * thin)
  nChains   <- 4
  
  set.seed(2)
  ptm <- proc.time()
  
  for (m in c("mENV","mAB")) {
    message("Sampling model: ", m)
    
    models[[m]] <- sampleMcmc(
      models[[m]],
      samples   = samples,
      thin      = thin,
      transient = transient,
      nChains   = nChains,
      nParallel = nChains,
      updater   = list(GammaEta = FALSE)
    )
  }
  
  print(runtime <- proc.time() - ptm)
  
  ##################################################
  ##### SAVE OUTPUT
  ##################################################
  
  filename.string <- paste0(
    res, "_model_cells_", nm,
    "_chains_", nChains,
    "_thin_", thin,
    "_samples_", samples
  )
  
  out_file <- file.path(
    "/pvol/2_fitting_and_running_models/",
    paste0(filename.string, ".Rdata")
  )
  
  save(
    models,
    XFormula,
    runtime,
    file = out_file
  )
  
  all_models[[nm]] <- models
  all_times[[nm]]  <- runtime
  
  message("Saved: ", out_file)
# }

##############################################################################################################
##### MODEL FIT & CV
##############################################################################################################

# model_ids <- names(formulas)
# 
# for (nm in model_ids) {
  
  message("====================================")
  message("Evaluating model: ", nm)
  message("====================================")
  
  model_file <- file.path(
    "/pvol/2_fitting_and_running_models/",
    paste0(
      res, "_model_cells_", nm,
      "_chains_", nChains,
      "_thin_", thin,
      "_samples_", samples,
      ".Rdata"
    )
  )
  
  load(model_file)  # loads: models, XFormula, runtime
  
  ##################################################
  ##### MODEL FIT
  ##################################################
  
  MF    <- list()
  preds <- list()
  
  for (m in c("mENV", "mAB")) {
    preds[[m]] <- computePredictedValues(models[[m]])
    MF[[m]]    <- evaluateModelFit(models[[m]], preds[[m]])
  }
  print("MF done")
  ##################################################
  ##### 5-FOLD CROSS-VALIDATION
  ##################################################
  ## This takes a long time, 5 days for one full model & 1h for an environment only model
  
  set.seed(2)
  partition <- createPartition(
    models$mENV,
    nfolds = 5,
    column = "transectID"
  )
  
  library(doParallel)
  UseCores <- 32
  cl <- makeCluster(UseCores, type = "FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
  registerDoParallel(cl)
  
  ptm <- proc.time()
  
  MF.cv    <- list()
  preds.cv <- list()
  
  for (m in c("mENV", "mAB")) {
    preds.cv[[m]] <- pcomputePredictedValues(
      models[[m]],
      partition = partition,
      nParallel = UseCores
    )
    
    MF.cv[[m]] <- evaluateModelFit(
      models[[m]],
      preds.cv[[m]]
    )
  }
  
  cv_time <- proc.time() - ptm
  stopCluster(cl)
  
  print("5-fold CV done")
  
  ##################################################
  ##### SAVE RESULTS
  ##################################################
  
  out_file <- file.path("/pvol/3_model_analysis/", paste0(res, "_model_cells_", nm, "_ENV_pa_ab_MF_CV.Rdata"))
  
  save(nm, XFormula, MF, MF.cv, preds, preds.cv, partition, runtime, cv_time, file = out_file)
  
  message("Saved evaluation results: ", out_file)
}







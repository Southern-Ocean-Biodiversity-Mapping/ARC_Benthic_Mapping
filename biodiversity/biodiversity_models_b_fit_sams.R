##### Setting up----
# library(raster)
# library(readxl)
# library(readr)
# library(dplyr)
# library(data.table)
# library(proj4)
# library(stringr)
# library(RColorBrewer)
# library(SOmap)
library(ecomix)

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

##### Select species with greater than ten occurrences across all sites.
spdata <- dat_cov_pa[,-which(colSums(dat_cov_pa)<10)][,1:20]

# samdat_10p <- cbind(spdata,cell_metadata_env[,c(21,22,43,70:76)])
#samdat_10p <- cbind(spdata, cell_metadata_env_clean_scaled[,c(21,22,43,69:ncol(cell_metadata_env_clean))])
samdat_10p <- cbind(spdata, cell_metadata_env_clean_scaled[,c(8,10:11,18,20,78,79,sel.not.correlated,68,69)])

##############################################################################################################
## Archetype formula
# archetype_formula <- as.formula(paste0(paste0('cbind(', 
#                                               paste(colnames(spdata),collapse=", "),
#                                               ") ~ ", 
#                                               "poly(depth, degree=2, raw=TRUE)+",
#                                               "poly(slope, degree=2, raw=TRUE)+",
#                                               "poly(waom4k_seafloorcurrents_mean, degree=2, raw=TRUE)+",
#                                               # "poly(waom4k_seafloorcurrents_residual, degree=2, raw=TRUE)+",
#                                               # "poly(waom4k_seafloortemperature, degree=2, raw=TRUE)+",
#                                               # "poly(waom4k_seafloorsalinity, degree=2, raw=TRUE)+",
#                                               "poly(waom4k_test_susp08, degree=2, raw=TRUE)+",
#                                               "poly(waom4k_test_settle08, degree=2, raw=TRUE)+",
#                                               "poly(distance2canyons, degree=2, raw=TRUE)"
#                                               # "poly(waom4k_seafloortemperature,degree=2,raw=TRUE)+
#                                               # poly(Oxygen,degree=2,raw=TRUE)+
#                                               # poly(Depth,degree=2,raw=TRUE)+
#                                               # poly(Productivity,degree=2,raw=TRUE)+
#                                               # Time"
# )))
archetype_formula <- as.formula(paste0(paste0('cbind(', 
                                              paste(colnames(spdata),collapse=", "),
                                              ") ~ ", 
                                              "depth + depth2 +",
                                              "slope + slope2 +",
                                              "waom4k_seafloorcurrents_absolute + waom4k_seafloorcurrents_absolute2 +",
                                              "waom4k_seafloorcurrents_residual + waom4k_seafloorcurrents_residual2 +",
                                              "waom4k_test_susp08 + waom4k_test_susp082 +",
                                              "waom4k_test_settle08 + waom4k_test_settle082 +",
                                              "distance2canyons + distance2canyons2"
)))

archetype_formula.null <- as.formula(paste0(paste0('cbind(', 
                                                   paste(colnames(spdata),collapse=", "),
                                                   ") ~ 1")))
archetype_formula.depth <- as.formula(paste0(paste0('cbind(', 
                                                    paste(colnames(spdata),collapse=", "),
                                                    ") ~ depth + depth2")))


## Species formula
species_formula <- ~ 1

## Fit a single model
sam_fit <- species_mix(archetype_formula = archetype_formula, # Archetype formula
                       species_formula = species_formula,    # Species formula
                       data = samdat_10p,            # Data
                       nArchetypes = 3,              # Number of groups (mixtures) to fit
                       family = 'bernoulli',         # Which family to use
                       control = list(quiet = TRUE))
plot(sam_fit,fitted.scale = 'logit')
plot(sam_fit,fitted.scale = 'logit', species="Bryozoan_Hard_BrAnt")

## multiple models
nArchetypes <- 6:9
sam_multifit <- species_mix.multifit(archetype_formula = archetype_formula, # Archetype formula
                                     # species_formula = species_formula,     # Species formula
                                     data = samdat_10p,                     # Data
                                     nArchetypes = nArchetypes,             # Number of groups (mixtures) to fit
                                     nstart = 10,                            # The number of fits per archetype.
                                     family = 'bernoulli',                  # Which family to use
                                     control = list(quiet = FALSE))
plot(sam_multifit,type="BIC")


## fit a single model multiple times:
sam_fit_list <- list()
for(i in 2:10){
  print(paste0(i," archetypes"))
  sam_fit_list[[i]] <- species_mix(archetype_formula = archetype_formula, # Archetype formula
                       species_formula = species_formula,    # Species formula
                       data = samdat_10p,            # Data
                       nArchetypes = i,              # Number of groups (mixtures) to fit
                       family = 'bernoulli',         # Which family to use
                       control = list(quiet = TRUE))
  message(paste0("BIC = ",round(sam_fit_list[[i]]$BIC,2)))
}



######################################## FROM RCP, NEED TO CHANGE TO SAM
##### stepwise variable selection
## specify null model
null_mod<-species_mix(archetype_formula=archetype_formula.null,
                      species_formula = species_formula,    # Species formula
                      nArchetypes=1,
                      data=samdat_10p,
                      family='bernoulli',
                      control=list(quiet=T))
null_mod$BIC

depth_mod<-species_mix(archetype_formula=archetype_formula.depth,
                       species_formula = species_formula,    # Species formula
                       nArchetypes=2,
                       data=samdat_10p,
                       family='bernoulli',
                       control=list(quiet=T))
depth_mod$BIC

env_vars <- c("distance2canyons","slope","waom4k_test_settle08")
weight <- rep(1,20)

## stepwise analysis
samdat_10p_r <- samdat_10p[,c(1:20,28,29,35,48)]
lin_step1<-SAM_fwd_step_linear(start_vars="depth", start_BIC=depth_mod$BIC,
                               add_vars=env_vars, species=colnames(spdata),
                               dist="bernoulli", data=samdat_10p_r,
                               nstarts=10, ## change to ~ 100 starts- this may take some time to run!
                               min.nSAM=2, max.nSAM=4,#weight=weight,
                               mc.cores=detectCores())   # if runing on VM van use more cores


##STEP BY STEP MANUALLY:
library(parallel)
start_vars="depth"
start_BIC=depth_mod$BIC
add_vars=env_vars
species=colnames(spdata)
dist="bernoulli"
data=samdat_10p
nstarts=50 ## change to ~ 100 starts- this may take some time to run!
min.nSAM=1
max.nSAM=3 #weight=weight,
mc.cores=detectCores()
init.sd=0.1

form.spp=NULL

#### Forward Selection-linear terms----
SAM_fwd_step_linear<-function(start_vars,             # names of variables to include in base model (usually from previous step)
                              start_BIC,              # value of BIC from best model in previous step
                              add_vars,               # names of variables to add (one at a time) in this step
                              species,                # names of species to model
                              dist= "Bernoulli",      # distribution for data model (see regimix())
                              nstarts=50,             # number of starts for regimix.multifit
                              form.spp=NULL,          # formula for species artifacts (see regimix())
                              data,                   # dataframe that contains all the data to fit regimix model
                              min.nSAM=1,             # minimum number of SAMs to consider
                              max.nSAM,               # maximum number of SAMs to consider
                              mc.cores=detectCores()#, # number of cores if parallel processing
                              # weight=NULL
)
{
  
  #set up results
  BICs<-as.data.frame(setNames(replicate((max.nSAM-min.nSAM+3),numeric(0), simplify = F), c("Var",paste0("SAM", rep(min.nSAM:max.nSAM)),"Start_BIC")))
  
  #loop through adding variables  
  for(j in 1:length(add_vars)){
    add<-add_vars[j]
    message(paste0("adding variable: ",add))
    temp_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~",           # 1+
                                paste0(paste(start_vars, collapse="+"), "+", paste(add, collapse = "+"))))
    #temp_form<-as.formula(paste("cbind(",paste(species, collapse=", "),")~",
    #                            paste0(paste(start_vars, collapse="+"), paste(add, collapse = "+"))))
    
    #run nSAMs  
    nSAMs_start<- list()
    for( ii in min.nSAM:max.nSAM){
      print(ii)
      nSAMs_start[[ii-diff(c(1,min.nSAM))]] <- species_mix.multifit(archetype_formula=temp_form,
                                                                    #species_formula=form.spp,
                                                                    data=data,
                                                                    nArchetypes=ii,
                                                                    #weights=weight,
                                                                    #inits="random2",
                                                                    nstart=nstarts,
                                                                    family=dist,
                                                                    control=list(quiet=T, init.sd=0.1, maxit=2000),
                                                                    mc.cores=mc.cores
      )
      
      print("fit done")
    }
    #get BICs
    SAM1_BICs <- sapply( nSAMs_start, function(x) sapply( x, function(y) y$BIC))
    # #Are any SAMs consisting of a small number of sites?  (A posteriori) If so remove.
    # RCP1_minPosteriorSites <- cbind( nrow(data), sapply( nRCPs_start[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))
    # RCP1_ObviouslyBad <- RCP1_minPosteriorSites < 2
    # RCP1_BICs[RCP1_ObviouslyBad] <- NA
    
    
    SAM1_minBICs <- apply( SAM1_BICs, 2, min, na.rm=TRUE)
    BICs[j,1:(ncol(BICs)-1)]<-c(paste0(add), round(SAM1_minBICs,1))
    
  }
  BICs$Start_BIC<-start_BIC
  return(BICs)
}






































matplot(t(lin_step1[,2:4]), type="l",col=1:11,lty=1, ylab="BIC", xlab="nRCP (-1)"); legend("topleft", legend=env_vars, lty=1, cex=0.7, col=1:11,bty="n") #cause have run 2-4 RCPs in code above, 1 on x-axis actually= RCP2
##
lin_step1
## choose best variable and run this, then repeat
lin_step2<-fwd_step_linear(start_vars="lflux", start_BIC=min(unlist(lin_step1)),
                           add_vars=env_vars[!env_vars %in% "lflux"],
                           species=species,  weight=weight,
                           dist="NegBin", data=sc_indat, nstarts=100, ## change to ~ 100 starts- this may take some time to run!
                           min.nRCP=2, max.nRCP =4, mc.cores=detectCores()) # if running on VM can use more cores

matplot(t(lin_step2[,2:4]), type="l",lty=1, col=1:11, ylab="BIC", xlab="nRCP (-1)"); legend("topleft", legend=env_vars[!env_vars %in% "lflux"], lty=1, cex=0.7, col=1:11,bty="n"); abline(h=min(unlist(lin_step1)), col="red", lty=2)





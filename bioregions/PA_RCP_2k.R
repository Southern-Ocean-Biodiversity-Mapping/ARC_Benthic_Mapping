####################################################
## Preliminary RCP analysis  @ 2km resolution     ##
## N.Hill Oct 2022                                ##
####################################################

#### 1) set up----
library(dplyr)
library(tidyr)
#library(raster)
library(terra)
#devtools::install_github('skiptoniam/ecomix')
library(ecomix)
library(ggplot2)
library(forcats)
library(ggpubr)

#comp= "nicole"
comp= "vm"

if(comp=="nicole"){
path<- "C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Benthic Diversity ARC\\Analysis\\ARC_Data\\"

# biological data at 2km cells and 2pc prevalence
load(paste0(path, "Cell_level_bio_2pc_2km.RData"))
load(paste0(path, "Cell_level_env_2km.Rdata"))
# add path for environmental data
}

if(comp=="vm"){
  path="/perm_storage/shared_space/BioMAS/"
  
  # biological data at 2km cells and 2pc prevalence
  load(paste0(path, "ARC_Data/Cell_level_bio_2pc_2km.RData"))
  load(paste0(path, "ARC_Data/Cell_level_env_2km.Rdata"))
  #load(paste0(path, "environmental_data/Circumpolar_EnvData_2km_shelf_mask_scaled_dataframe.Rdata"))
  load(paste0(path,"environmental_data/Circumpolar_Coastline.Rdata"))
}


#### 2) Format biological and environmental data ----
# try using presence-absence data as a starting point
# Need to work out how to do point count data (i.e. if binomial is possible)

# remove substrate points, noID and < 2pc prevalence
prev_2pc<-cover_mod [,-1]
prev_2pc[prev_2pc>0] <- 1
prev_2pc<-colSums(prev_2pc)
prev_2pc<-names(prev_2pc[prev_2pc>10])

cover_2pc<-cover_mod %>%
  dplyr::select(., all_of(prev_2pc)) %>%
  dplyr::select(., ! contains ("Sub_")) %>%
  dplyr::select( ., - NoID) %>%
  dplyr::select( ., - Unscorable)

# convert bio to presence-absence data
cover_2pc[cover_2pc>0] <- 1

# add CellID back into df
cover_2pc<-cbind(cover_mod[,"cellID"], cover_2pc)

# environmental factors to run
covars <- c("cellID","cover_N", "cover_cells_survey",
            "depth","depth2", 
            "logslope",
            "sst_mean",
            "sst_sd",
            "tpi11",
            "arag_mean",
            "no3_mean",
            "seafloorcurrents_mean",
            "test_settle08",
            "test_susp08")

#merge into one dataframe
bioenv<-left_join(cover_2pc, cell_metadata_env_scaled[,covars])

#check nas
colSums(is.na(bioenv))
# lots (31) for sst_mean and sd??? check with Jan- remove nas for now
bioenv<-na.omit(bioenv)

##### 3) Test run multifit RCPs and check time ----
# use survey as the sampling factor- will need to look at gamma penalty for this
# use Cover_N = number of images as the weighting factor

rcp_form<-as.formula(paste0("cbind(",paste(colnames(cover_2pc[-1]),
                                           collapse = ','),
                            ")~1+", paste0(covars[4:14], collapse= "+")))



start_time <- Sys.time() 

nstarts <-100
max.nRCP <- 8
nRCPs_samp <- list()
for( ii in 2:max.nRCP)
  nRCPs_samp[[ii]] <- regional_mix.multifit(rcp_formula = rcp_form,
                                            species_formula = spp_form,
                                            data = bioenv,
                                            nRCP=ii,
                                            family = "bernoulli",
                                            inits="random2",
                                            weights=cover_N,
                                            nstart=nstarts,
                                            mc.cores = 5)

end_time <- Sys.time() #2.9 hours

save(nRCPs_samp,file = paste0(path, "RCP/PA_Multi_100starts.RData"))

# get model BICs and plot
grps <- 2:max.nRCP
RCPsamp_BICs <- sapply( nRCPs_samp[-1], function(x) sapply( x, function(y) y$BIC))
RCPsamp_minPosteriorSites <- cbind(dim(bioenv)[1], sapply( nRCPs_samp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))

RCPsamp_minBICs <- apply( RCPsamp_BICs, 2, min)

df2a <- data.frame(grps=2:max.nRCP,bic=RCPsamp_minBICs)
df2b <- data.frame(grps=rep(2:max.nRCP, each=nrow( RCPsamp_BICs)),bic=as.numeric(RCPsamp_BICs))
gg1 <- ggplot(df2a,aes(x=grps,y=bic))+
  geom_point()+
  geom_line()+
  geom_point(data=df2b,aes(x=grps,y=bic))+
  scale_x_continuous("Number of Groups", labels = as.character(grps), breaks = grps)+
  ylab("BIC")
gg1

# 8 groups is not enough add additional RCPs to above

max.nRCP <- 14
for( ii in 9:max.nRCP)
  nRCPs_samp[[ii]] <- regional_mix.multifit(rcp_formula = rcp_form,
                                            species_formula = spp_form,
                                            data = bioenv,
                                            nRCP=ii,
                                            family = "bernoulli",
                                            inits="random2",
                                            weights=cover_N,
                                            nstart=nstarts,
                                            mc.cores = 6)

RCPsamp_BICs <- sapply( nRCPs_samp[-1], function(x) sapply( x, function(y) y$BIC))
RCPsamp_minPosteriorSites <- cbind(dim(bioenv)[1], sapply( nRCPs_samp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))

RCPsamp_minBICs <- apply( RCPsamp_BICs, 2, min)

df2a <- data.frame(grps=2:max.nRCP,bic=RCPsamp_minBICs)
df2b <- data.frame(grps=rep(2:max.nRCP, each=nrow( RCPsamp_BICs)),bic=as.numeric(RCPsamp_BICs))
gg1 <- ggplot(df2a,aes(x=grps,y=bic))+
  geom_point()+
  geom_line()+
  geom_point(data=df2b,aes(x=grps,y=bic))+
  scale_x_continuous("Number of Groups", labels = as.character(grps), breaks = grps)+
  ylab("BIC")
gg1


#hmmm BIC is still going down after 14 RCPs. Add more later- get rest of code working first
max.nRCP <- 21
for( ii in 15:max.nRCP)
  nRCPs_samp[[ii]] <- regional_mix.multifit(rcp_formula = rcp_form,
                                            species_formula = spp_form,
                                            data = bioenv,
                                            nRCP=ii,
                                            family = "bernoulli",
                                            inits="random2",
                                            weights=cover_N,
                                            nstart=nstarts,
                                            mc.cores = 6)
grps <- 2:max.nRCP
RCPsamp_BICs <- sapply( nRCPs_samp[-1], function(x) sapply( x, function(y) y$BIC))
RCPsamp_minPosteriorSites <- cbind(dim(bioenv)[1], sapply( nRCPs_samp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))

RCPsamp_minBICs <- apply( RCPsamp_BICs, 2, min)

df2a <- data.frame(grps=2:max.nRCP,bic=RCPsamp_minBICs)
df2b <- data.frame(grps=rep(2:max.nRCP, each=nrow( RCPsamp_BICs)),bic=as.numeric(RCPsamp_BICs))
gg1 <- ggplot(df2a,aes(x=grps,y=bic))+
  geom_point()+
  geom_line()+
  geom_point(data=df2b,aes(x=grps,y=bic))+
  scale_x_continuous("Number of Groups", labels = as.character(grps), breaks = grps)+
  ylab("BIC")
gg1



#### Get best model, run diagnostics, generate bootstraps ----
RCPsamp_goodun <- which.min( RCPsamp_BICs[,13])
control <- list( optimise=TRUE, quiet=TRUE)
RCPsamp_fin <- regional_mix(rcp_formula = rcp_form,
                            species_formula = spp_form,
                            data = bioenv,
                            nRCP=14,
                            family = "bernoulli",
                            weights=cover_N,
                            inits = unlist(nRCPs_samp[[14]][[RCPsamp_goodun]]$coef),
                            control = control,
                            titbits = TRUE)

plot(RCPsamp_fin, type="RQR", fitted.scale="log") 
#looks like some departure at one tail and increase in variance....
plot(RCPsamp_fin, type="RQR") 
#looks better not on log scale

#increase bootstraps at later date
rcpsamp_boots <- regional_mix.bootstrap(RCPsamp_fin, type="BayesBoot", nboot=50, mc.cores=5)
save(rcpsamp_boots,file = paste0(path, "RCP/PA_boots.RData"))

#### Secies Profiles----

sp_prof <- regional_mix.species_profile(RCPsamp_fin,rcpsamp_boots,type='response')

#format output for dotplot
sp_prof_dat<-as.data.frame(sp_prof$overall$mean) %>%
  tibble::rownames_to_column(., var="RCP") %>%
  pivot_longer(., cols= UBS_Sub:Crustaceans_Isopods_Serolidae, names_to= "species", values_to= "mean")

sp_prof_dat_sd<-as.data.frame(sp_prof$overall$sd) %>%
  tibble::rownames_to_column(., var="RCP") %>%
  pivot_longer(., cols=UBS_Sub:Crustaceans_Isopods_Serolidae , names_to= "species", values_to= "sd")

sp_prof_dat<-left_join(sp_prof_dat, sp_prof_dat_sd)
sp_prof_dat$low<-sp_prof_dat$mean- sp_prof_dat$sd
sp_prof_dat$upp<-sp_prof_dat$mean + sp_prof_dat$sd
#manually make upper and lower bounds 1,0
sp_prof_dat$low[sp_prof_dat$low <0] <- 0
sp_prof_dat$upp[sp_prof_dat$upp >1] <- 1

#set up to plot each RCP with 10 most prevalent species individually
RCPs<-paste0("RCP", 1:14)
sp_plot<-list()
for (i in 1: length(RCPs)) {
temp<-sp_prof_dat %>%
  filter(., RCP == RCPs[i]) %>%
  mutate(species = fct_reorder(species, mean))
  #arrange (., desc(mean))

sp_plot[[i]]<-ggplot(data= temp[1:10,])+ 
  geom_point(aes(x=mean, y=species)) +
  geom_linerange(aes(xmin= low,xmax=upp, y=   species))+
  xlim(0,1)+
  ggtitle(RCPs[i])+
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab("Occurrence")+
  ylab("Species")
} 

jpeg(filename=paste0(path,"ARC_Benthic_Mapping/bioregions/sp_prof_sel.jpg"), width=20, height=12, units="cm", res=300)
sp_prof_fig<-ggarrange(sp_plot[[1]],sp_plot[[4]], sp_plot[[14]],
                       ncol=2, nrow=2, align="h")
sp_prof_fig
dev.off()

#### spatial prediction ----

# get raster stack and turn into data frame for prediction. 
# Will need to modify this workflow as have already done hte prediction, but did not have the x,y values
env_stack<-rast(paste0(path, "environmental_data/Circumpolar_EnvData_2km_shelf_mask_scaled.tif"))
temp_pred<-as.data.frame(env_stack, xy=TRUE)
covars2<-c("x", "y", covars[4:14])
pred_stack.dat<-temp_pred[,covars2]
pred_stack.dat<-na.omit(pred_stack.dat)

#replace cover_cells_survey with most common survey - PS81
# think more about this in terms of which survey have more biota than others....
pred_stack.dat$cover_cells_survey<- factor("PS81", levels=levels(bioenv$cover_cells_survey))

# add line to remove na.s
#pred_stack.dat<-na.omit(pred_stack.dat[, covars[4:14]])

start_time <- Sys.time() 
RCPsamp_SpPreds <- predict(object=RCPsamp_fin, object2=rcpsamp_boots, newdata=pred_stack.dat, mc.cores=5)
end_time <- Sys.time() #18 minutes

save(RCPsamp_SpPreds,file = paste0(path, "RCP/RCPsamp_SpPreds.RData"))

#Hard class probabilities
hc_preds_df<-cbind(pred_stack.dat[,c("x", "y")], 
                   HC_RCP=apply(RCPsamp_SpPreds$ptPreds,1, which.max))
hc_preds_vec<-vect(hc_preds_df, geom= c("x", "y"), crs=crs(env_stack))
hc_pred_rast<-rasterize(hc_preds_vec, env_stack, field="HC_RCP")
plot(hc_pred_rast, type="class", col=grDevices::topo.colors(15))

#certainty for each hard class
hc_probs_df<-cbind(pred_stack.dat[,c("x", "y")], 
                   RCP_prob=apply(RCPsamp_SpPreds$ptPreds,1, max))
hc_probs_vec<-vect(hc_probs_df, geom= c("x", "y"), crs=crs(env_stack))
hc_probs_rast<-rasterize(hc_probs_vec, env_stack, field="RCP_prob")
plot(hc_probs_rast, col=brewer.pal(11, "Blues"))

#plot average RCP probabilities
av_preds_df<-cbind(pred_stack.dat[,c("x", "y")], RCPsamp_SpPreds$ptPreds)
av_preds_vec<-vect(av_preds_df, geom= c("x", "y"), crs=crs(env_stack))
# can directly input  RCPsamp_SpPreds$pt_preds here


# Work out how to do this more efficiently
RCPs<-paste0("RCP_", 1:14)
RCP_list<-list()
for (i in 1:length(RCPs)){
#for ( i in 1:2){
  RCP_list[[i]]<-rasterize(av_preds_vec, env_stack, field=RCPs[i])
}
par(mfrow=c(4,4))
for (i in 1:14){
  plot(RCP_list[[i]], main=RCPs[i])
}

RCP14<-rast(RCP_list)
names(test)<- paste0("RCP ", 1:14)
save(RCP_list, file = paste0(path, "RCP/RCP_rasters.RData"))

#### Plot partial response to environment
plot_vars<-covars[4:14]
part_effect_df<-effectPlotData(focal.predictors = plot_vars, mod=RCPsamp_fin)



par(mfrow=c(4,3), mar=c(5,4,2,2))
for(i in 1:length(part_effect_df)){
  temp<-as.data.frame(part_effect_df[[i]])
  
  pred<-predict(object=RCPsamp_fin, newdata=temp)
  matplot(x=as.vector(part_effect_df[[i]][,plot_vars[i]]), 
          y= pred, type='l', 
          xlab=plot_vars[i], ylab="Probability of Occurrence", ylim=c(0,1), lwd=3)
  legend(x="topright", legend=dimnames(pred)[[2]], col=1:length(dimnames(pred)[[2]]), lty=1, cex=0.9, bty="n")
}

  

####### Demonstration plots for discussion with Nicole Webster
library(RColorBrewer)
brewer.pal(brewer.pal.info["RdYlBu", "maxcolors"], "RdYlBu")
col1<-c("#A50026", "#D73027", "#F46D43", "#BF812D", "#FDAE61", "#FEE090", "#FFFFBF" ,"#E0F3F8",
        "#ABD9E9", "#74ADD1", "#35978F", "#4575B4" ,"#313695","#4D4D4D")
coast.proj<-vect(coast.proj)

plot(hc_pred_rast, type="class", col=col1,
     plg=list( title="Bioregion", cex=1.3), main= "Bioregions- Hard Class", cex.main=1.2)
lines(coast.proj)


#prob of hard class
plot(hc_probs_rast, col=brewer.pal(9, "Blues"),
     plg=list( title="Probability"), main="Most likely Bioregion", cex.main=1.2)
lines(coast.proj)

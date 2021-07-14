####################################################################
### Test GAM workflow with MI LL data -----                       ##
####################################################################

## 1) set-up
library(raster)
library(MASS)
library(mgcv)
library(zoo)

path<-"C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Toothfish FRDC\\Analysis\\"

# MI LL data
bioenv<-readRDS(paste0(path, "GAMS\\MI\\MI_ll_bioenv.rds"))

#one -ve cpue values??? exclude
bioenv<-bioenv[bioenv$cpue >0, ]
# keep data to Dec 2019- extent of BRAN data
bioenv<-bioenv[bioenv$Year <2020,]

#checking base standardisation/inclusion variables
table(bioenv$Vessel, bioenv$Year) #no real overlap between vessel and operating years (completely confounded??)
range(bioenv$cpue)
hist(bioenv$cpue)
summary(bioenv$CatchT)
hist(bioenv$CatchT)
table(bioenv$Year) # fewer hauls 2007-2010
table(bioenv$Year, bioenv$Month) #Aug- Apr. Initially Jul-Aug then spread, 2015 + all months
#factor or smooth?



##### setting up cross validation sets and use same for all comparisons/ model selection
#random selection of hauls, each year separately, all years except the last two...
#random- only keep average value across all folds
# yearly- want to look at values for each year.
CV_folds<- data.frame(Rand=sample(x=c(1:3),size=dim(bioenv)[1], replace=TRUE)
                      , Rand2=sample(x=c(1:3),size=dim(bioenv)[1], replace=TRUE))
#,Yearly=as.integer(as.factor(bioenv$Year)))
#,last2=ifelse(bioenv$Year<2018, 1, 2))

##### Set up some helper function for model selection and cross validation
## Measures of goodness of fit
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

MAE=function(m, o){
  mean(abs(m-o))
}


### Stepwise and cross val code
catch_stepwise<-function(dep_var,          #name of dependent variable
                         base_formula,     #character string containing elements of GAM formula to use as base model
                         start_vars=NULL,  #name of environmental predictor variables  in formula
                         add_vars,         #name of new environemntal variables to consider for new step
                         offset_var=NULL,  #names of offset variable
                         #family=tw,       #name of distribution for modelling
                         data,             #name of dataframe containing variables
                         cv_df             #name of dataframe containing cv folds
)
{
 
  #base_formula<-"s(FishingDepth) + s(Year,k=13, bs='cr') + as.factor(Month) + te(MidLon,MidLat, DaysSinceStart, d=c(2,1),bs=c('tp','cr'))"
  
  base_formula<-base_formula
  
  #loop through each add_var, create new formula:
  
  for(j in 1:length(add_vars)){
    res<-list()
    
    #if start_vars included  
    if(!is.null(start_vars)){
      temp_form<-as.formula(paste0(dep_var, "~" , base_formula, "+ offset(log(", offset_var, "))", 
                                   paste("+", paste("s(", start_vars,")", collapse = "+")),
                                   "+ s(", add_vars[j], ")"))
    }
    else{
      temp_form<-as.formula(paste0(dep_var, "~", base_formula, "+ s(", add_vars[j], ")")) 
    }
    
    # i) run GAM on entire dataset, extract goodness of fit metrics                     
    full_gam<-gam(temp_form, family=tw, data=data) 
    full_res<-data.frame(edf=round(sum(full_gam$edf),1), 
                         DevExpl= round(summary(full_gam)$dev.expl *100, 1),
                         AIC=round(AIC(full_gam),0),
                         r2=round(summary(full_gam)$r.sq, 3))
    
    # ii) run through cross validation options in cv dataframe, extract prediction metrics 
    cv_res<-list()
    for(k in 1:ncol(cv_df)){
      #set up dataframe for results
      cv_res[[k]]<-data.frame(RMSE= rep(NA, length(unique(cv_df[,k]))+1),
                              MAE= rep(NA, length(unique(cv_df[,k]))+1),
                              Spear= rep(NA, length(unique(cv_df[,k]))+1),
                              R2= rep(NA, length(unique(cv_df[,k]))+1))
      
      #run the actual cross validation for folds in column k                   
      for(l in 1:length(unique(cv_df[,k]))){
        #define datasets
        ind<-which(cv_df[,k]== l)
        train<-data[- ind, ]
        test<-data[ind,]
        
        #re-define  formula with number knots= equal to number of years in train data
        base_formula<-paste0("s(FishingDepth) + s(Year,k=", length(unique(train$Year)), ", bs='cr') + as.factor(Month) + te(MidLon,MidLat, DaysSinceStart, d=c(2,1),bs=c('tp','cr'))")
        
        #if start_vars included  
        if(!is.null(start_vars)){
          temp_form<-as.formula(paste0(dep_var, "~" , base_formula, "+ offset(log(", offset_var, "))", 
                                       paste("+", paste("s(", start_vars,")", collapse = "+")),
                                       "+ s(", add_vars[j], ")"))
        }
        else{
          temp_form<-as.formula(paste0(dep_var, "~", base_formula, "+ s(", add_vars[j], ")")) 
        }
        # run CV and extract results
        cv_gam<-gam(temp_form, family=tw, data=train)
        cv_pred<- predict.gam(cv_gam, newdata = test, type="response")
        cv_res[[k]]$RMSE[l]<-round(RMSE(cv_pred, test$CatchT),2)
        cv_res[[k]]$MAE[l]<-round(MAE(cv_pred, test$CatchT),2)
        cv_res[[k]]$Spear[l]<-round(cor(test$CatchT, cv_pred, method="spearman"),2)
        cv_res[[k]]$R2[l]<-round(cor(test$CatchT, cv_pred, method="pearson"),2)
      }
      cv_res[[k]][l+1,]<-round(colMeans(cv_res[[k]][1:l,], na.rm=TRUE) ,2) 
    }
    res[[j]]<-c(list(full_res), cv_res)
    #names(res[[j]])<-c("Full_mod", names(cv_df))
  }
  names(res)<-add_vars
  return(res)
}


## Base variables included in CPUE standardisation/ base model----
#choices: biomass or CPUE directly?
# month (factor or smooth)
# year (smooth TPS? to allow predcition outside current range. Check tails)
# if year is a term then seperate f(lat,lon)??
# depth * year interaction to account for local depletion??


##biomass----
# Null
catch_null<-gam(CatchT ~ offset(log(bioenv$Effort)), family=tw, data=bioenv)
AIC(catch_null) #7787

#year (cubic spline with 13 knots so are evenly spaced. Continuous so can predict to other years), month (factor because only 5 months Apr-Aug and not present all years),
#depth (TPS), space 
# year, month, depth, space * daily or monthly time
#daily time for short length correlation between fishing days
#monthly time for longer correlation cycles
bioenv$DaysSinceStart<-as.numeric(bioenv$Date - min(bioenv$Date))
bioenv$MonthSinceStart<-as.numeric((bioenv$MonthYear - min(bioenv$MonthYear)) *12)


# a) start with year as fixed factor and lat/long
catch_base_1<-gam(CatchT ~ offset(log(Effort))+ s(SoakTime) + s(Year, bs="cr") + as.factor(Month) + s(FishingDepth) +
                   te(MidLon,MidLat),family=tw, data=bioenv)
summary(catch_base_1)               
plot(catch_base_1, page=1)                


# b) look at days since start and low K for long-term trend



, DaysSinceStart, d=c(2,1),bs=c("tp","cr")), 
                family=tw, data=bioenv) 
summary(catch_base2) #DE= 36.2%, AIC= 6307
anova(catch_base2) #all factors significant
gam.check(catch_base2) #may need to increase k for year and space




#works running line by line- but first list of cv_res is NULL when running as function with multiple CV fold options



# add sst and then sst_anomaly
#also add select term that can drop terms out of model...
catch_sst<-gam(CatchT ~ offset(log(bioenv$Effort))+ s(Year,k=13, bs="cr") + as.factor(Month) + s(FishingDepth) +
                 te(MidLon,MidLat, DaysSinceStart, d=c(2,1),bs=c("tp","cr")) 
               + s(MI_sst_yr_mth2) +s(MI_sstanom_yr_mth), 
               select=TRUE, family=tw, data=bioenv) 
summary(catch_sst) #DE 37.4% (not much extra!)






#see if prediction works on data outside training- works
build_dat<-bioenv[bioenv$Year!= "2019",]
test_dat<-bioenv[bioenv$Year == "2019",]
test_gam<-gam(CatchT ~ offset(log(Effort))+ s(Year,k=12, bs="cr") + as.factor(Month) + s(FishingDepth) +
                           te(MidLon,MidLat, DaysSinceStart, d=c(2,1),bs=c("tp","cr")) 
                         + s(MI_sst_yr_mth2) +s(MI_sstanom_yr_mth), 
                         family=tw, data=build_dat) 

test_pred<-predict.gam(test_gam, newdata=test_dat, type="response")

#generate average yearly prediction from model and check against raw data....
#create standardised effort column to use for prediction
##cpue is in kg/hook (catchT/hooks*1000); model is in T so to get equivalent is per 1000 hooks?? 
pred_dat<-
sst_pred<-predict.gam()

# define different validation sets

###CPUE
# Null
cpue_null<-gam(cpue ~ 1, family=tw, data=bioenv)
AIC(cpue_null) #-1053

#year (cubic spline), month, depth, space
cpue_base<-gam(cpue ~ s(Year, k=13, bs="cr") + as.factor(Month) + s(FishingDepth) + s(MidLon,MidLat), 
                family=tw, data=bioenv) 
summary(cpue_base) #DE= 32.3%, AIC= -1942
anova(cpue_base) #depth not significant
gam.check(cpue_base) #may need to increase k for year and space

# year, month, depth, space * daily or monthly time
cpue_base2<-gam(cpue ~ s(Year,k=13, bs="cr") + as.factor(Month) + s(FishingDepth) +
                   te(MidLon,MidLat, DaysSinceStart, d=c(2,1),bs=c("tp","cr")), 
                 family=tw,  data=bioenv) 
summary(cpue_base2) #DE= 34.9%, AIC= -1983
anova(cpue_base2) #all factors significant
gam.check(cpue_base2) #may need to increase k for year and space

#### outputs/conclusions very similar


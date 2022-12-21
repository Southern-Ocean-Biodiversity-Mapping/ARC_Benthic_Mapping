library(tidyverse)
library(abind)
library(Hmsc)
library(sjSDM)
library(mistnet)
library(Metrics)

source('utils.R')

# Load data
load("../../data/modellingdata/modelling_data.RData")
Y = as.data.frame(df_bio_clean[, 2:ncol(df_bio_clean)])
Y_pa = as.data.frame((Y > 0) * 1)
X = as.data.frame(df_env_clean[, 4:ncol(df_env_clean)])

path_trained_models = "../trained_model"

list_algo = c("hmsc", "sjsdm", "mistnet")

df <- data.frame(algo=character(),
                 model_name=character(),
                 pred_r2_mean=double(), pred_r2_sd=double(),
                 pred_mae_mean=double(), pred_mae_sd=double(),
                 pred_rmse_mean=double(), pred_rmse_sd=double())

for (algo_ in list_algo) {
  print(paste0("Evaluating ", algo_, " models ..."))
  path_algo = paste0(path_trained_models, "/", algo_, "/abd/")
  
  lst_model_names = list.dirs(path = path_algo,
                              full.names = FALSE,
                              recursive = FALSE)
  
  if (algo_ == "hmsc") {
    path_best_pa = paste0(path_trained_models, "/", algo_, "/pa/20220927_model___hmsc")
    
    path_model_full = paste0(path_best_pa, "/full")
    load(path_model_full)
    
    pred_ = computePredictedValues(model_fit)
    
    pred = matrix(NA, nrow=model_fit$ny, ncol=model_fit$ns)
    sel = model_fit$distr[,1]==3
    if (sum(sel)>0){
      pred[,sel] = as.matrix(apply(abind(pred_[,sel,,drop=FALSE], along=3),
                                   c(1,2), median2))
    }
    sel = !model_fit$distr[,1]==3
    if (sum(sel)>0){
      pred[,sel] = as.matrix(apply(abind(pred_[,sel,,drop=FALSE], along=3),
                                   c(1,2), mean2))
    }
    
    lst_thr <- list()
    colnames(pred) <- colnames(Y_pa)
    for (morpho_taxon in colnames(Y_pa)) {
      roc_morpho_taxon <- roc(Y_pa[, morpho_taxon], pred[, morpho_taxon])
      lst_thr[morpho_taxon] = coords(roc_morpho_taxon, "best", ret = "threshold")
    }
    rm(model_fit)
    rm(pred)
  }
  
  for (model_name_ in lst_model_names) {
    print(model_name_)
    
    # Load model
    #path_model_full = paste0(path_algo, model_name_, "/full")
    #load(path_model_full)
    
    # EVALUATE PREDICTION POWER
    lst_pred_rmse <- lst_pred_r2 <- lst_pred_mae <- c()
    for (fold_ID in 1:5) {
      path_Kmodel = paste0(path_algo, model_name_, "/fold_", fold_ID)
      
      load(path_Kmodel)
      
      X_test = X[which(df_metadata_$fold == fold_ID), ]
      Y_test = Y[which(df_metadata_$fold == fold_ID), ]
      
      if (algo_ == "hmsc") {
        gradient = Hmsc::prepareGradient(hM = model_fit, XDataNew = as.data.frame(X_test))
        Y_pred_ = predict(model_fit,
                          Gradient=gradient,
                          expected = TRUE,
                          nParallel = length(model_fit$postList))
        
        Y_pred_abd_cod = apply(simplify2array(Y_pred_), 1:2, mean)
        #Y_pred = as.data.frame(Y_pred_)
        rm(model_fit)
        rm(gradient)
        
        path_Kmodel_pa = paste0(path_best_pa, "/fold_", fold_ID)
        load(path_Kmodel_pa)
        gradient = Hmsc::prepareGradient(hM = model_fit, XDataNew = as.data.frame(X_test))
        Y_pred_pa = predict(model_fit,
                          Gradient=gradient,
                          expected = TRUE,
                          nParallel = length(model_fit$postList))
        
        Y_pred_pa = apply(simplify2array(Y_pred_pa), 1:2, mean)
        
        Y_pred_pa_bin = Y_pred_pa * 0.
        for (morpho_taxon in names(lst_thr)) {
          Y_pred_pa_bin[, morpho_taxon] = (Y_pred_pa[, morpho_taxon] > lst_thr[morpho_taxon]) * 1
        }
        
        Y_pred_abd = Y_pred_abd_cod * Y_pred_pa_bin * 100.
        Y_pred_abd = as.data.frame(Y_pred_abd)
      } else if (algo_ == "sjsdm") {
        Y_pred = predict(model_fit, newdata = X_test)
        Y_pred_abd = Y_pred * 100
      } else if (algo_ == "mistnet") {
        Y_pred_ = predict(model_fit,
                          newdata=X_test,
                          n.importance.samples=model_fit$n.importance.samples)
        Y_pred = apply(Y_pred_, 1:2, mean)
        Y_pred_abd = Y_pred * 100
      } else if (algo_ == "brms") {
        Y_pred = data.frame()
        for (sp in names(m_species)) {
          if (nrow(Y_pred) == 0) {
            lst_pred = list()
            lst_pred[[sp]] = predict(m_species[[sp]], newdata = X_test)[, 1]
            Y_pred = rbind(Y_pred, lst_pred)
          } else {
            Y_pred[[sp]] = predict(m_species[[sp]], newdata = X_test)[, 1]
          }
        }
      }
      
      lst_pred_r2 = append(lst_pred_r2, mean(computeR2(Y_test, Y_pred_abd), na.rm=TRUE))
      lst_pred_rmse = append(lst_pred_rmse, mean(computeRMSE(Y_test, Y_pred_abd)))
      lst_pred_mae = append(lst_pred_mae, mae(actual=Y_test, predicted = as.matrix(Y_pred_abd)))
      
      rm(model_fit)
    }
    
    df = df %>% add_row(algo = algo_,
                        model_name = model_name_,
                        pred_r2_mean = round(mean(lst_pred_r2), 4),
                        pred_r2_sd = round(sd(lst_pred_r2), 4),
                        pred_mae_mean = round(mean(lst_pred_mae), 4),
                        pred_mae_sd = round(sd(lst_pred_mae), 4),
                        pred_rmse_mean = round(mean(lst_pred_rmse), 4),
                        pred_rmse_sd = round(sd(lst_pred_rmse), 4)
    )
  }
  print(df)
}

print(df)
#write.csv(df, "20221028_abd_hmsc.csv", row.names = FALSE)

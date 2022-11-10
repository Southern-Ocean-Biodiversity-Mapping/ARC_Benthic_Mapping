library(tidyverse)
library(abind)
library(caret)
#library(boral)
#library(Hmsc)
#library(sjSDM)
#library(mistnet)
library(brms)

source('utils.R')

# Load data
load("../04_data.RData")
Y = as.matrix(df_pa[, 2:ncol(df_pa)])
X = as.matrix(df_env_scaled[, 2:ncol(df_env_scaled)])

path_trained_models = "../trained_model"

list_algo = c("brms") #c("hmsc", "boral", "sjsdm", "mistnet", "brms")

df <- data.frame(algo=character(),
                 model_name=character(),
                 fit_rmse_mean=double(),
                 fit_auc_mean=double(),
                 pred_auc_mean=double(),
                 pred_sensitivity_mean=double(),
                 pred_specificity_mean=double(),
                 pred_precision_mean=double(),
                 pred_pcc_mean=double(),
                 fit_rmse_sd=double(),
                 fit_auc_sd=double(),
                 pred_auc_sd=double(),
                 pred_sensitivity_sd=double(),
                 pred_specificity_sd=double(),
                 pred_precision_sd=double(),
                 pred_pcc_sd=double())

for (algo_ in list_algo) {
  print(paste0("Evaluating ", algo_, " models ..."))
  path_algo = paste0(path_trained_models, "/", algo_, "/pa/")
  
  lst_model_names = list.dirs(path = path_algo,
                              full.names = FALSE,
                              recursive = FALSE)
  
  for (model_name_ in lst_model_names) {
    print(model_name_)
    
    # Load model
    path_model_full = paste0(path_algo, model_name_, "/full")
    load(path_model_full)
    
    # EVALUATE FIT
    rmse <- auc <- NULL
    if (algo_ == "boral") {
      pred = fitted(model_fit, est = "median", include.ranef = TRUE, linear.predictor = FALSE)$out

    } else if (algo_ == "sjsdm") {
      pred = predict(model_fit)
    } else if (algo_ == "mistnet") {
      pred = predict(model_fit,
                     newdata=X,
                     n.importance.samples=model_fit$n.importance.samples)
      pred = apply(pred, 1:2, mean)
    } else if (algo_ == "hmsc") {
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
    } else if (algo_ == "brms") {
      pred = data.frame()
      for (sp in names(m_species)) {
        if (nrow(pred) == 0) {
          lst_pred = list()
          lst_pred[[sp]] = predict(m_species[[sp]])[, 1]
          pred = rbind(pred, lst_pred)
        } else {
          pred[[sp]] = predict(m_species[[sp]])[, 1]
        }
      }
    }
    
    rmse_ = computeRMSE(Y, pred)
    auc_ = computeAUC(Y, pred)
    
    lst_thr <- list()
    colnames(pred) <- colnames(Y)
    for (morpho_taxon in colnames(Y)) {
      roc_morpho_taxon <- roc(Y[, morpho_taxon], pred[, morpho_taxon])
      lst_thr[morpho_taxon] = coords(roc_morpho_taxon, "best", ret = "threshold")
    }
    
    rm(model_fit)
    
    # EVALUATE PREDICTION POWER
    lst_pred_auc <- lst_pred_sensitivity <- lst_pred_specificity <- lst_pred_precision <- lst_pred_pcc <- c()
    for (fold_ID in 1:5) {
      path_Kmodel = paste0(path_algo, model_name_, "/fold_", fold_ID)

      load(path_Kmodel)
      
      X_test = X[which(df_metadata_$fold == fold_ID), ]
      Y_test = Y[which(df_metadata_$fold == fold_ID), ]
      
      if (algo_ == "boral") {
        Y_pred = predict(model_fit,
                         newX = X_test,
                         predict.type = "conditional",
                         scale = "response",
                         est = "median",
                         prob = 0.95,
                         return.alllinpred = FALSE)$linpred
      } else if (algo_ == "hmsc") {
        gradient = Hmsc::prepareGradient(hM = model_fit, XDataNew = as.data.frame(X_test))
        Y_pred_ = predict(model_fit,
                         Gradient=gradient,
                         expected = TRUE,
                         nParallel = length(model_fit$postList))
        
        Y_pred_ = apply(simplify2array(Y_pred_), 1:2, mean)
        Y_pred = as.data.frame(Y_pred_)
      } else if (algo_ == "sjsdm") {
        Y_pred = predict(model_fit, newdata = X_test)
      } else if (algo_ == "mistnet") {
        Y_pred_ = predict(model_fit,
                          newdata=X_test,
                          n.importance.samples=model_fit$n.importance.samples)
        Y_pred = apply(Y_pred_, 1:2, mean)
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
      
      pred_k_auc = mean(computeAUC(Y_test, Y_pred), na.rm = TRUE)
      lst_pred_auc = append(lst_pred_auc, pred_k_auc)
      
      colnames(Y_pred) <- colnames(Y_test)
      Y_pred_bin = Y_pred * 0.
      for (morpho_taxon in names(lst_thr)) {
        Y_pred_bin[, morpho_taxon] = (Y_pred[, morpho_taxon] > lst_thr[morpho_taxon]) * 1
      }
      
      conf_matrix = confusionMatrix(table(as.matrix(Y_pred_bin), as.matrix(Y_test)))
      lst_pred_sensitivity = append(lst_pred_sensitivity, conf_matrix$byClass["Sensitivity"][[1]])
      lst_pred_specificity = append(lst_pred_specificity, conf_matrix$byClass["Specificity"][[1]])
      lst_pred_precision = append(lst_pred_precision, conf_matrix$byClass["Precision"][[1]])
      lst_pred_pcc = append(lst_pred_pcc, conf_matrix$overall["Accuracy"][[1]])
      
      rm(model_fit)
    }

    df = df %>% add_row(algo = algo_,
                        model_name = model_name_,
                        fit_rmse_mean = round(mean(rmse_), 4),
                        fit_auc_mean = round(mean(auc_), 4),
                        pred_auc_mean = round(mean(lst_pred_auc, na.rm = TRUE), 4) * 100.,
                        pred_sensitivity_mean = round(mean(lst_pred_sensitivity, na.rm = TRUE), 4) * 100.,
                        pred_specificity_mean = round(mean(lst_pred_specificity, na.rm = TRUE), 4) * 100.,
                        pred_precision_mean = round(mean(lst_pred_precision, na.rm = TRUE), 4) * 100.,
                        pred_pcc_mean = round(mean(lst_pred_pcc, na.rm = TRUE), 4) * 100.,
                        fit_rmse_sd = round(sd(rmse_), 4),
                        fit_auc_sd = round(sd(auc_), 4),
                        pred_auc_sd = round(sd(lst_pred_auc, na.rm = TRUE), 4) * 100.,
                        pred_sensitivity_sd = round(sd(lst_pred_sensitivity, na.rm = TRUE), 4) * 100.,
                        pred_specificity_sd = round(sd(lst_pred_specificity, na.rm = TRUE), 4) * 100.,
                        pred_precision_sd = round(sd(lst_pred_precision, na.rm = TRUE), 4) * 100.,
                        pred_pcc_sd = round(sd(lst_pred_pcc, na.rm = TRUE), 4) * 100.
    )
  }
  print(df)
}

print(df)
#write.csv(df, "20221028_fit_withThr_mistnet_brms.csv", row.names = FALSE)

library(caret)
library(Hmsc)
library(abind)

source('../evaluate/utils.R')

# Load data
load("../../data/modellingdata/modelling_data.RData")
Y = as.data.frame(df_bio_clean[, 2:ncol(df_bio_clean)])
Y_pa = as.data.frame((Y > 0) * 1)

# Load model
load("../../data/modellingdata/model_PA_final")

pred_ = computePredictedValues(model_fit)
Y_pred_ = apply(simplify2array(pred_), 1:2, mean)
Y_pred = as.data.frame(Y_pred_)

rmse_ = computeRMSE(Y_pa, Y_pred)
auc_ = computeAUC(Y_pa, Y_pred)

mean_rmse = mean(rmse_)
sd_rmse = sd(rmse_)
print(paste0('RMSE: ', round(mean_rmse, 2), ' +/- ', round(sd_rmse, 2)))
mean_auc = mean(auc_)
sd_auc = sd(auc_)
print(paste0('AUC: ', round(mean_auc, 2), ' +/- ', round(sd_auc, 2)))

colnames(Y_pred) <- colnames(Y_pa)
Y_pred_bin = Y_pred * 0.
for (morpho_taxon in colnames(Y_pa)) {
  roc_morpho_taxon <- roc(Y_pa[, morpho_taxon], Y_pred[, morpho_taxon])
  thr = coords(roc_morpho_taxon, "best", ret = "threshold")
  Y_pred_bin[, morpho_taxon] = (Y_pred[, morpho_taxon] > thr[[1]]) * 1
}

conf_matrix = confusionMatrix(table(as.matrix(Y_pred_bin), as.matrix(Y_pa)))
print('Confusion Matrix')
conf_matrix
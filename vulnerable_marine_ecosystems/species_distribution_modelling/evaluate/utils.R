library(pROC)


computeRMSE = function(Y, predY){
  ns = dim(Y)[2]
  RMSE = rep(NA,ns)
  for (i in 1:ns){
    RMSE[i] = sqrt(mean((Y[,i]-predY[,i])^2, na.rm=TRUE))
  }
  return(RMSE)
}

computeR2 = function(Y, predY, method="pearson"){
  ns = dim(Y)[2]
  R2 = rep(NA,ns)
  for (i in 1:ns){
    co = cor(Y[,i], predY[,i], method=method, use='pairwise')
    R2[i] = sign(co)*co^2
  }
  return(R2)
}

computeAUC = function(Y, predY){
  ns = dim(Y)[2]
  AUC = rep(NA,ns)
  ## take care that Y has only levels {0,1} as specified in auc() below
  Y <- ifelse(Y > 0, 1, 0)
  for (i in 1:ns){
    sel = !is.na(Y[,i])
    if(length(unique(Y[sel,i]))==2)
      AUC[i] = pROC::auc(Y[sel,i],predY[sel,i], levels=c(0,1),direction="<")
  }
  return(AUC)
}

median2 = function(x){return (median(x,na.rm=TRUE))}
mean2 = function(x){return (mean(x,na.rm=TRUE))}
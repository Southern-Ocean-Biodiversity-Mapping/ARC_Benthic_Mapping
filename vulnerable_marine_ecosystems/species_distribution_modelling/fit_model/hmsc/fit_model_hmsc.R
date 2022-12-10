fit_hmsc = function(folder_name, Y, X, currency_data, distribution, n_parallel, n_chain, n_sample, thin) {
  # Save params
  save(n_chain, n_parallel, n_sample, thin, currency_data, distribution,
       file=paste0(folder_name, "/params.Rdata"))
  
  # Full model
  fname = paste0(folder_name, "/full")
  if (!file.exists(fname)) {
    print("Full model ...")
    if (distribution == "normal") {
      model_raw <- Hmsc(Y=Y, XData = X, distr=distribution, YScale = TRUE)
    } else {
      model_raw <- Hmsc(Y=Y, XData = X, distr=distribution)
    }
    
    model_fit = sampleMcmc(model_raw,
                           samples = n_sample,
                           thin = thin,
                           transient = ceiling(0.5 * n_sample * thin),
                           nChains = n_chain,
                           nParallel = n_parallel)
    
    save(model_fit, file=fname)
  }

  # kFold models
  for (fold_ID in unique(df_metadata_$fold)) {
    print(paste0('Fitting model #', fold_ID))
    
    Y_fold = Y[which(df_metadata_$fold != fold_ID), ]
    X_fold = X[which(df_metadata_$fold != fold_ID), ]
    
    fname <- file.path(folder_name, paste0("fold_", fold_ID))
    if (!file.exists(fname)) {
      if (distribution == "normal") {
        model_raw <- Hmsc(Y=Y_fold, XData = X_fold, distr=distribution, YScale = TRUE)
      } else {
        model_raw <- Hmsc(Y=Y_fold, XData = X_fold, distr=distribution)
      }
      
      model_fit = sampleMcmc(model_raw,
                             samples = n_sample,
                             thin = thin,
                             transient = ceiling(0.5 * n_sample * thin),
                             nChains = n_chain,
                             nParallel = n_parallel)
      
      save(model_fit, file=fname)
    }
  }
}







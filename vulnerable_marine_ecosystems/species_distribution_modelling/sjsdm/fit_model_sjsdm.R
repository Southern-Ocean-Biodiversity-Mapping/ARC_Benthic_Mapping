fit_sjsdm = function(folder_name, X, Y, hyperparams) {
  dir.create(folder_name)
  
  # Save params
  save(hyperparams,
       file=paste0(folder_name, "/params.Rdata"))
  
  # Full model
  fname = paste0(folder_name, "/full")
  if (!file.exists(fname)) {
    print("Full model ...")
    if (hyperparams$method_env_ == "DNN") {
      env_ = DNN(X, hidden = c(10L, 10L, 10L))
    } else {
      env_ = linear(data = X)
    }
    model_fit <- sjSDM(Y = Y,
                       env = env_,
                       se = FALSE,
                       family=binomial("probit"),
                       iter=hyperparams$n_iter_,
                       sampling = hyperparams$n_sampling_)
    fname = paste0(folder_name, "/", "full")
    save(model_fit, file=fname)
  }
  
  # kFold models
  for (fold_ID in unique(df_metadata_$fold)) {
    print(paste0('Fitting model #', fold_ID))
   
    Y_fold = Y[which(df_metadata_$fold != fold_ID), ]
    X_fold = X[which(df_metadata_$fold != fold_ID), ]
    
    fname <- file.path(folder_name, paste0("fold_", fold_ID))
    if (!file.exists(fname)) {
      Y_fold = Y[which(df_metadata_$fold != fold_ID), ]
      X_fold = X[which(df_metadata_$fold != fold_ID), ]
      
      if (hyperparams$method_env_ == "DNN") {
        env_ = DNN(X_fold, hidden = c(10L, 10L, 10L))
      }
      else {
        env_ = linear(data = X_fold)
      }
      model_fit <- sjSDM(Y = Y_fold,
                         env = env_,
                         se = FALSE,
                         family=binomial("probit"),
                         iter=hyperparams$n_iter_,
                         sampling = hyperparams$n_sampling_)
      fname = paste0(folder_name, "/fold_", fold_ID)
      save(model_fit, file=fname)
    }
  }

}

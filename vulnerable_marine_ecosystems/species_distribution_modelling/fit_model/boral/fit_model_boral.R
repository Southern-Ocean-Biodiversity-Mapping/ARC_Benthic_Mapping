fit_boral <- function(folder_name,
                      Y,
                      X,
                      currency_data,
                      distribution,
                      n_burning,
                      n_iter,
                      n_thin,
                      seed) {
  n <- nrow(Y)
  p <- ncol(Y)
  
  # Save params
  save(n_burning, n_iter, n_thin, currency_data, distribution,
       file=paste0(folder_name, "/params.Rdata"))
  
  mcmc_control <- list(n.burnin = n_burning,
                       n.iteration = n_iter,
                       n.thin = n_thin,
                       seed = seed)
  
  fname_model <- file.path(folder_name, "full.txt")
  if (!file.exists(fname_model)) {
    print('Fitting model # Full')
    model_fit <- boral(y = Y,
                       X = X,
                       family = distribution,
                       row.eff = "none",
                       mcmc.control = mcmc_control,
                       save.model = TRUE,
                       model.name = fname_model)
  
    # Save model
    save(model_fit,
         file=paste0(folder_name, "/full"))
  }
  
  # kFold models
  #for (fold_ID in unique(df_metadata_$fold)) {
  #  print(paste0('Fitting model #', fold_ID))
  #  
  #  Y_fold = Y[which(df_metadata_$fold != fold_ID), ]
  #  X_fold = X[which(df_metadata_$fold != fold_ID), ]
  # 
  #  fname_model <- file.path(folder_name, paste0("fold_", fold_ID, ".txt"))
  #  if (!file.exists(fname_model)) {
  #    model_fit <- boral(y = Y_fold,
  #                       X = X_fold,
  #                       family = distribution,
  #                       row.eff = "none",
  #                       mcmc.control = mcmc_control,
  #                       save.model = TRUE,
  #                       model.name = fname_model)
  #    
  #    fname = paste0(folder_name, "/fold_", fold_ID)
  #    save(model_fit, file=fname)
  #  }
  #}
}
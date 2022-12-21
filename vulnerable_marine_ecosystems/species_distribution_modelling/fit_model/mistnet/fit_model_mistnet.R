fit_model = function(X, Y, hyperparams){

  MNet_mod = mistnet(
    x = X,
    y = Y,
    layer.definitions = list(
      defineLayer(
        nonlinearity = rectify.nonlinearity(),
        size = hyperparams$n_layer1_,
        prior = gaussian.prior(mean = 0, sd = .1)
      ),
      defineLayer(
        nonlinearity = sigmoid.nonlinearity(),
        size = ncol(Y),
        prior = gaussian.prior(mean = 0, sd = .1)
      )
    ),
    loss=bernoulliLoss(),
    updater = adagrad.updater(learning.rate = hyperparams$learning_rate_),
    sampler = gaussian.sampler(ncol = 10L, sd=1),
    n.importance.samples = hyperparams$n_importance_samples_,
    n.minibatch = hyperparams$n_minibatch_,
    training.iterations = 0,
    initialize.biases = TRUE,
    initialize.weights = TRUE
  )
  MNet_mod$layers[[1]]$biases[] = 1 # First layer biases equal 1
  
  start.time = Sys.time()
  while(
    difftime(Sys.time(), start.time, units = "secs") < hyperparams$fit_seconds_
  ){
    MNet_mod$fit(100)
    cat(".")
    # Update prior variance
    for(layer in MNet_mod$layers){
      layer$prior$update(
        layer$weights, 
        update.mean = FALSE, 
        update.sd = TRUE,
        min.sd = .01
      )
    }
    # Update mean for final layer
    MNet_mod$layers[[2]]$prior$update(
      layer$weights, 
      update.mean = TRUE, 
      update.sd = FALSE,
      min.sd = .01
    )
  } # End while
  
  MNet_mod
}

fit_mistnet = function(folder_name, X, Y, hyperparams){
  
  dir.create(folder_name, showWarnings = FALSE)
  
  # Save params
  save(hyperparams,
       file=paste0(folder_name, "/params.Rdata"))
  
  fname = paste0(folder_name, "/full")
  if (!file.exists(fname)) {
    print('Fitting model # Full')
    model_fit <- fit_model(X = X, Y=Y, hyperparams = hyperparams)
    # Save model
    save(model_fit, file=fname)
  }
  
  # kFold models
  for (fold_ID in unique(df_metadata_$fold)) {
    fname = paste0(folder_name, "/fold_", fold_ID)
    if (!file.exists(fname)) {
      print(paste0('Fitting model #', fold_ID))
      
      Y_fold = Y[which(df_metadata_$fold != fold_ID), ]
      X_fold = X[which(df_metadata_$fold != fold_ID), ]
      
      model_fit <- fit_model(X = X_fold, Y = Y_fold, hyperparams = hyperparams)

      save(model_fit, file=fname)
    }
  }
}
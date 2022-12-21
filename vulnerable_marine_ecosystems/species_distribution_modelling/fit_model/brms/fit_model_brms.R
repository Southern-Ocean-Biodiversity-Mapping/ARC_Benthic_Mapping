fit_brms = function(folder_name, X, Y, distribution, currency, trials){
  
  dir.create(folder_name, showWarnings = FALSE)
  
  if (currency == "abd") {
    n_chains = 2
  } else {
    n_chains = 4
  }
  
  # Save params
  save(distribution, currency,
       file=paste0(folder_name, "/params.Rdata"))
  
  lst_species = names(Y)
  
  fname = paste0(folder_name, "/full")
  if (!file.exists(fname)) {
    print('Fitting model # Full')
    m_species = list()
    for (species in lst_species) {
      print(species)
      bio = data.frame(bio_data=Y[, species])
      dat_ = cbind(bio, X)
      
      if (trials) {
        m <- brm(bio_data | trials(1) ~ depth+depth_sq+slope_log+ice_sd+ssh_mean+ssh_sd+sst_sd+tpi11+arag_mean+o2_mean+waom4k_seafloorcurrents_mean+waom4k_seafloorcurrents_residual+waom4k_seafloortemperature+waom4k_seafloorsalinity+waom4k_test_settle08+waom4k_test_susp08,
                   data = dat_, family = distribution, chains=n_chains)
      } else {
        if (currency == "abd") {
          m <- brm(bio_data ~ depth+slope_log+ssh_mean+tpi11+arag_mean+o2_mean+waom4k_seafloorcurrents_mean+waom4k_seafloortemperature+waom4k_seafloorsalinity,
                   data = dat_, family = distribution, chains=n_chains)
        } else {
          m <- brm(bio_data ~ depth+depth_sq+slope_log+ice_sd+ssh_mean+ssh_sd+sst_sd+tpi11+arag_mean+o2_mean+waom4k_seafloorcurrents_mean+waom4k_seafloorcurrents_residual+waom4k_seafloortemperature+waom4k_seafloorsalinity+waom4k_test_settle08+waom4k_test_susp08,
                   data = dat_, family = distribution, chains=n_chains)
        }
      }

      m_species[[species]] = m
    }
    # Save model
    save(m_species, file=fname)
  }
  
  # kFold models
  for (fold_ID in unique(df_metadata_$fold)) {
    fname = paste0(folder_name, "/fold_", fold_ID)
    if (!file.exists(fname)) {
      print(paste0('Fitting model #', fold_ID))
      
      Y_fold = Y[which(df_metadata_$fold != fold_ID), ]
      X_fold = X[which(df_metadata_$fold != fold_ID), ]
      
      m_species = list()
      for (species in lst_species) {
        bio = data.frame(bio_data=Y_fold[, species])
        dat_ = cbind(bio, X_fold)
        
        if (trials) {
          m <- brm(bio_data | trials(1) ~ depth+depth_sq+slope_log+ice_sd+ssh_mean+ssh_sd+sst_sd+tpi11+arag_mean+o2_mean+waom4k_seafloorcurrents_mean+waom4k_seafloorcurrents_residual+waom4k_seafloortemperature+waom4k_seafloorsalinity+waom4k_test_settle08+waom4k_test_susp08,
                   data = dat_, family = distribution, chains=n_chains)
        } else {
          if (currency == "abd") {
            m <- brm(bio_data ~ depth+slope_log+ssh_mean+tpi11+arag_mean+o2_mean+waom4k_seafloorcurrents_mean+waom4k_seafloortemperature+waom4k_seafloorsalinity,
                     data = dat_, family = distribution, chains=n_chains)
          } else {
            m <- brm(bio_data ~ depth+depth_sq+slope_log+ice_sd+ssh_mean+ssh_sd+sst_sd+tpi11+arag_mean+o2_mean+waom4k_seafloorcurrents_mean+waom4k_seafloorcurrents_residual+waom4k_seafloortemperature+waom4k_seafloorsalinity+waom4k_test_settle08+waom4k_test_susp08,
                     data = dat_, family = distribution, chains=n_chains)
          }
        }
        
        m_species[[species]] = m
      }
      # Save model
      save(m_species, file=fname)
    }
  }
}
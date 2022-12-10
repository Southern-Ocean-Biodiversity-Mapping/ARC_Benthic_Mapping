library(brms)

seed_ = 7109
set.seed(seed_)

source('fit_model_brms.R')

load("../../data/modellingdata/modelling_data.RData")
X = as.data.frame(df_env_clean[, 4:ncol(df_env_clean)])
Y = as.data.frame(df_bio_clean[, 2:ncol(df_bio_clean)])

# Params
currency_data = "abd" #"abd"

if (currency_data == "pa") {
  folder_root = "trained_models/brms/pa"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame((Y > 0) * 1)
  
  distr = c("bernoulli", "binomial")
} else if (currency_data == "abd") {
  folder_root = "trained_models/brms/abd"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame(sapply(Y, as.integer))
  
  distr = c("zero_inflated_negbinomial"
            #"zero_inflated_poisson"
            #"hurdle_poisson"
            #"hurdle_negbinomial"
            )
}

n_param = length(distr)

for (i in 1:n_param) {
  rdn_number = sample(1:10000, 1)
  folder_name_cur = paste0(folder_root, '/', str_pad(rdn_number, 5, pad = "0"))
  print(folder_name_cur)
  
  trials = FALSE
  if (distr[[i]] == 'binomial' | distr[[i]] == 'zero_inflated_binomial') {
    trials = TRUE
  }
  
  model_fit = fit_brms(folder_name = folder_name_cur,
                       X = X,
                       Y = Y,
                       distribution = distr[[i]],
                       currency=currency_data,
                       trials)
}
library(mistnet)

seed_ = 7109
set.seed(seed_)

source('fit_model_mistnet.R')

load("../../data/modellingdata/modelling_data.RData")
X = as.data.frame(df_env_clean[, 4:ncol(df_env_clean)])
Y = as.data.frame(df_bio_clean[, 2:ncol(df_bio_clean)])

# Params
currency_data = "abd" #"pa"
n_minibatch = rep(25, 6)
sampler_size = rep(c(5, 10), 3)
n_importance_samples = rep(25, 6)
n_layer1 = c(8, 8, 16, 16, 24, 24)
learning_rate = rep(c(0.1, 0.01), 3)
fit_seconds = rep(90, 6)

if (currency_data == "pa") {
  folder_root = "trained_models/mistnet/pa"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame((Y > 0) * 1)
} else if (currency_data == "abd") {
  folder_root = "trained_models/mistnet/abd"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame(Y / 100)
}

n_param = length(n_minibatch)

for (i in 1:n_param) {
  rdn_number = sample(1:10000, 1)
  folder_name_cur = paste0(folder_root, '/', str_pad(rdn_number, 5, pad = "0"))
  print(folder_name_cur)
  
  hyperparam_cur = list(currency_data_=currency_data,
                        n_minibatch_=n_minibatch[[i]],
                        sampler_size_=sampler_size[[i]],
                        n_importance_samples_=n_importance_samples[[i]],
                        n_layer1_=n_layer1[[i]],
                        learning_rate_=learning_rate[[i]],
                        fit_seconds_=fit_seconds[[i]])
  
  model_fit = fit_mistnet(folder_name = folder_name_cur,
                          X = X,
                          Y = Y,
                          hyperparams = hyperparam_cur)
}
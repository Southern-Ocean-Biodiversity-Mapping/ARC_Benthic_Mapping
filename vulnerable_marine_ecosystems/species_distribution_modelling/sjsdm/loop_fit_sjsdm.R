library(sjSDM)
library(stringr)

seed_ = 7109
set.seed(seed_)

source('fit_model_sjSDM.R')

load("../04_data.RData")

X = as.matrix(df_env_scaled[, 2:ncol(df_env_scaled)])

# Params
currency_data = "abd" #"pa"
method_env = c(rep("linear", 4), rep("DNN", 4))
n_sampling = rep(c(100L, 500L), 4)
n_iter = rep(c(500L, 500L, 1000L, 1000L), 2)

if (currency_data == "pa") {
  folder_root = "../trained_model/sjsdm/pa"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.matrix(df_pa[, 2:ncol(df_pa)])
} else if (currency_data == "abd") {
  folder_root = "../trained_model/sjsdm/abd"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.matrix(df_abd[, 2:ncol(df_abd)]) / 100
}

n_param = length(n_iter)

for (i in 1:n_param) {
  rdn_number = sample(1:10000, 1)
  folder_name_cur = paste0(folder_root, '/', str_pad(rdn_number, 5, pad = "0"))
  print(folder_name_cur)
  
  hyperparam_cur = list(currency_data_=currency_data,
                        method_env_=method_env[[i]],
                        n_sampling_=n_sampling[[i]],
                        n_iter_=n_iter[[i]])
  
  fit_sjsdm(folder_name = folder_name_cur,
            X = X,
            Y = Y,
            hyperparams = hyperparam_cur)
}
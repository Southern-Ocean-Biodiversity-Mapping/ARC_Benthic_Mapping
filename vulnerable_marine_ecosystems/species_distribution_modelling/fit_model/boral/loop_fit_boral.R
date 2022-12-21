packages <- c("boral")
install.packages(setdiff(packages, rownames(installed.packages())))

library(boral)
library(stringr)

seed_ = 7109
set.seed(seed_)

source('fit_model_boral.R')

load("../../data/modellingdata/modelling_data.RData")
X = as.data.frame(df_env_clean[, 4:ncol(df_env_clean)])
Y = as.data.frame(df_bio_clean[, 2:ncol(df_bio_clean)])

# Params
currency_data_ = "abd" #"pa"
lst_n_burning = c(10, 10, 10) #c(10000, 20000, 10000)
lst_n_iter = c(100, 100, 100) #c(80000, 40000, 40000)
lst_n_thin = c(30, 30, 60)
distr = "ztnegative.binomial" #"binomial"

if (currency_data_ == "pa") {
  folder_root = "trained_models/boral/pa"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame((Y > 0) * 1)
} else if (currency_data_ == "abd") {
  folder_root = "trained_models/boral/abd"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame(Y / 100)
  Y = sapply(Y, as.integer)
  
  Y[Y == 0] = NA
}

n_param = length(lst_n_burning)
n_param = 1

for (i in 1:n_param) {
  rdn_number = sample(1:10000, 1)
  folder_name_cur = paste0(folder_root, "/", str_pad(rdn_number, 5, pad = "0"), "_boral")
  print(folder_name_cur)
  dir.create(folder_name_cur, showWarnings = FALSE)
  
  fit_boral(folder_name = folder_name_cur,
            Y = Y,
            X = X,
            currency_data = currency_data_,
            distribution = distr,
            n_burning = lst_n_burning[[i]],
            n_iter = lst_n_iter[[i]],
            n_thin = lst_n_thin[[i]],
            seed = seed_
            )
}
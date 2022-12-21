library(Hmsc)
library(stringr)

seed_ = 7109
set.seed(seed_)

source('fit_model_hmsc.R')

load("../../data/modellingdata/modelling_data.RData")
X = as.data.frame(df_env_clean[, 4:ncol(df_env_clean)])
Y = as.data.frame(df_bio_clean[, 2:ncol(df_bio_clean)])

# Params
currency_data_ = "abd" # pa"
distr = "probit"
nChains = c(2, 2, 2, 2, 4, 4, 4, 4)
nParallel = rep(1, 8)
nSamples = rep(c(1000, 500), 4)
thinning = rep(c(100, 100, 10, 10), 2)

if (currency_data_ == "pa") {
  folder_root = "trained_models/hmsc/pa"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame((Y > 0) * 1)
} else if (currency_data_ == "abd") {
  folder_root = "trained_models/hmsc/abd"
  dir.create(folder_root, showWarnings = FALSE)
  
  Y = as.data.frame(Y / 100)
}

if (is.numeric(as.matrix(Y)) || is.logical(as.matrix(Y)) && is.finite(sum(Y, na.rm=TRUE))) {
  print("Y looks OK")
} else {
  print("Y should be numeric and have finite values")	}

if (any(is.na(X))) {
  print("X has NA values - not allowed for")
} else {
  print("X looks ok")	}

n_param = length(nParallel)

for (i in 1:n_param) {
  rdn_number = sample(1:10000, 1)
  folder_name_cur = paste0(folder_root, "/", str_pad(rdn_number, 5, pad = "0"), "_hmsc")
  print(folder_name_cur)
  dir.create(folder_name_cur, showWarnings = FALSE)
  
  fit_hmsc(folder_name = folder_name_cur,
           Y = Y,
           X = X,
           currency_data = currency_data_,
           distribution = distr,
           n_parallel=nParallel[[i]],
           n_chain=nChains[[i]],
           n_sample=nSamples[[i]],
           thin=thinning[[i]]
  )
}

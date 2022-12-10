library(Hmsc)
library(stringr)

seed_ = 7109
set.seed(seed_)

load("../preprocessing/modelling_data.RData")

X = as.data.frame(df_env_clean[, 4:ncol(df_env_clean)])
Y = as.data.frame(df_bio_clean[, 2:ncol(df_bio_clean)])
Y_pa = as.data.frame((Y > 0) * 1)
  
# Params
distribution = "probit"
nChains = 2
nParallel = 1
nSamples = 500
thinning = 100

if (is.numeric(as.matrix(Y)) || is.logical(as.matrix(Y)) && is.finite(sum(Y, na.rm=TRUE))) {
  print("Y looks OK")
} else {
  print("Y should be numeric and have finite values")	}

if (any(is.na(X))) {
  print("X has NA values - not allowed for")
} else {
  print("X looks ok")	}

model_raw <- Hmsc(Y=Y, XData = X, distr=distribution)

model_fit = sampleMcmc(model_raw,
                       samples = nSamples,
                       thin = thinning,
                       transient = ceiling(0.5 * nSamples * thinning),
                       nChains = nChains,
                       nParallel = nParallel)

fname = "model_PA_final"
save(model_fit, file=fname)

# Evaluate convergence
mpost = convertToCodaObject(model_fit, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
# Gelman diagnostics, i.e. the Potential scale reduction factors
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
ma = psrf.beta[,1]
quantile(ma)

library(rstan)

full_data <- readRDS("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/full_data.rds")

y0 <- as.numeric(full_data[,"outcome"])
x0 <- as.matrix(full_data[,-"outcome"])

fileName <- "./log_reg.stan"
stanCode <- readChar(fileName, file.info(fileName)$size)

dat <- list(
  n      = nrow(x0),
  p      = ncol(x0),
  ntrail = 1,
  x      = x0,
  y      = y0
)

rtime <- proc.time()
resStan <- stan(
  model_code = stanCode,
  data = dat,
  chains = 1,
  iter = 10000,
  warmup = 5000,
  thin = 5
)

res <- resStan@sim$samples[[1]][1:ncol(x0)]
rtime <- proc.time() - rtime

fname <- paste0("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/full_res.rds")
saveRDS(list(res = res, time = rtime[3]), fname)
  

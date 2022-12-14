##
library(rstan)
library(glue)

## Get index number
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  i <- 1
} else {
  i <- args[[1]]
}

start_time <- Sys.time()

## Load the subset of the data
load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/partition{i}.rds"))

## Checks
# sum(y0)
# table(y0, x0[,4])

## Call individual submission
fileName <- "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/stoc_log_reg.stan"
stanCode <- readChar(fileName, file.info(fileName)$size)

dat <- list(n      = nrow(x0),
            p      = ncol(x0),
            ntrail = 1,
            nrep   = nrep0,
            x      = x0,
            y      = y0)

resStan <- stan(model_code = stanCode, data = dat,
                chains = 1, iter = 10000, warmup = 5000, thin = 5)

## Check
## library(bayestestR)
## describe_posterior(as.data.frame(resStan))


res <- do.call(cbind, resStan@sim$samples[[1]][1:ncol(x0)])
## describe_posterior(as.data.frame(res))

fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/stan/res{i}.rds")
save(res, file = fname)

end_time <- Sys.time()
save(end_time - start_time, file = glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/stan/time{i}.rds"))


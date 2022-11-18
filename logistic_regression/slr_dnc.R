## Get index number
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  i <- 1
} else {
  i <- args[[1]]
}

## Load the subset of the data
load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/partition{i}.rds"))

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

res <- do.call(cbind, resStan@sim$samples[[1]][1:ncol(x0)])

fname <- paste0("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/res{i}.rds")
save(res, file = fname)
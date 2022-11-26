## Run on login node
## If needed to remove previous partitions:
## rm /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/*

## Submit the partitioning script:
## qsub -pe smp -2 -cwd -e /dev/null -o /dev/null /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partition.job


## Submit the job scripts:
## stan (currently thinking we get rid of this)
## qsub -pe smp -2 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out -t 1-25 /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/slr.job

## mh
## SAVE TIME HERE, COMPARE TO THE FINAL COMPLETION TIME FOR TOTAL RUN TIME
## date > /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/start_time.txt
## qsub -q BIOSTAT -pe smp -2 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out -t 1-25 /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/mh_dnc.job

################################################################################
####### LOG INTO COMPUTE NODE BEFORE PROCEEDING ################################
################################################################################

################################################################################
#### Combine the results from MH Runs  #########################################
################################################################################

########################
#### Libraries #########
########################
library(glue)
library(magrittr)
library(bayestestR)
library(coda)
library(dplyr)

###############################
#### Check acceptance rate ####
###############################

all_counts <- numeric(25)
for (j in 1:25) {
  load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/acc_counts/mh_large/acc_counts{j}.rds"))
  all_counts[j] <- acc_counts
}

min(all_counts) / 10000
max(all_counts) / 10000
mean(all_counts) / 10000

###############################
#### Check average run time ###
###############################

times <- numeric(25)
for (j in 1:25) {
  load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/time{j}.rds"))
  times[j] <- end_time
}

max(times) ## time in seconds
as.Date((max(times) / (3600*24)), origin = "1970-01-01")

library(lubridate)
library(stringr)
with_tz(as_datetime(max(times)), "America/Chicago")
  
date_string <- read_lines("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/start_time.txt")
split_date <- str_split(date_string, " ")[[1]]
month_num <- which(month.abb == split_date[2])
as_datetime(paste0(split_date[6], "-", month_num, "-", split_date[3], " ", split_date[4]), tz = "America/Chicago")

with_tz(as_datetime(max(times)), "America/Chicago") - as_datetime(paste0(split_date[6], "-", month_num, "-", split_date[3], " ", split_date[4]), tz = "America/Chicago")

## "%a %b %d %h:%s:%m %Z %Y"
as.POSIXct(read_lines("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/start_time.txt"), format = "%a %b %d %h:%s:%m")
######################################
#### Put all draws in single list ####
######################################

start_time <- Sys.time()
results <- list()
for (j in 1:25) {
  load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/res{j}.rds"))
  results[[j]] <- res
}

########################
#### Remove burn-in ####
########################
remove_burnin <- function(x, burnin) {
  
  x[-(1:burnin),]
  
}

results <- lapply(results, remove_burnin, 1000)
## TIME END HERE
end_time <- Sys.time()
tdiff <- end_time - start_time
fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/read_time.rds")
save(tdiff, file = fname)

################################################################################
#### Recenter Draws using simple method ########################################
################################################################################
## TIME START HERE
start_time <- Sys.time()
K <- 25
subset_mean <- t(sapply(1:K, function(i) colMeans(results[[i]]))) # rows indicate subset
global_mean <- colMeans(bind_rows(lapply(results, data.frame)))
recenter <- t(sapply(1:K, function(i) subset_mean[i, ] - global_mean)) # rows indicate subset
results_recentered <- lapply(1:K,function(i) results[[i]] - matrix(recenter[i, ],
                                                                   nrow = nrow(results[[i]]),
                                                                   ncol = ncol(results[[i]]),
                                                                   byrow = T))
full_data_draws <- do.call(rbind, results_recentered)
## TIME END HERE
end_time <- Sys.time()
tdiff <- end_time - start_time
fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/simple_center_time.rds")
save(tdiff, file = fname)

#####################
#### Diagnostics ####
#####################
full_data_draws <- full_data_draws %>% as.mcmc()
heidel.diag(full_data_draws)
raftery.diag(full_data_draws)

#################
#### Summary ####
#################
tmp <- describe_posterior(as.data.frame(full_data_draws))
tmp$Median
names(tmp)

head(full_data_draws)

tmp_carrier <- exp(full_data_draws[,5])
res_carrier <- describe_posterior(as.data.frame(tmp_carrier))
res_carrier$Median
res_carrier$CI_low
res_carrier$CI_high
################################################################################
###### Combine using WASP approximation ########################################
################################################################################
## TIME START HERE
start_time <- Sys.time()
computeBarycenter <- function (meanList, covList) {
  library(matrixStats)
  library(expm)
  library(MASS)
  
  ncomp <- length(covList)
  ndim <- nrow(covList[[1]])
  wts <- rep(1, ncomp) / ncomp
  
  baryMean <- rowMeans(do.call(cbind, meanList))
  baryCov <- diag(1.0, ndim)
  barySd <- sqrtm(baryCov)
  err <- baryCov
  cnt <- 1
  while ((norm(err, type = "F") > 1e-6) & cnt < 500) {
    if (cnt %% 10 == 0)  cat("iter: ", cnt, "\n")
    
    ssj <- matrix(0.0, nrow = ndim, ncol = ndim)
    for (ii in 1:ncomp) {
      ssj <- ssj + sqrtm(baryCov %*% covList[[ii]])
    }
    
    tmp <- solve(sqrtm(baryCov))
    baryCovNew <- tmp %*% tcrossprod(ssj / ncomp) %*% tmp
    err <- baryCov - baryCovNew
    baryCov <- baryCovNew
    barySd <- sqrtm(baryCov)
    cnt <- cnt + 1
  }
  
  list(mean = baryMean, cov = baryCov, sqrt = barySd, iter = cnt)
}

sampleBetas <- function (betaList) {
  
  library(matrixStats)
  library(expm)
  
  npart <- length(betaList)
  meanBetas <- list()
  covBetas <- list()
  
  meanBetas <- lapply(betaList, function(x) colMeans(x))
  covBetas <- lapply(betaList, function(x) cov(x))
  
  resBary <- computeBarycenter(meanBetas, covBetas)
  muBetas <- resBary$mean
  sigBetas <- resBary$cov
  sqrtSigBetas <- resBary$sqrt
  
  baryList <- list()
  for (ii in 1:npart) {
    tmp <- as(chol2inv(chol(covBetas[[ii]])), "symmetricMatrix")
    tmp1 <- matrix(meanBetas[[ii]], nrow = nrow(betaList[[ii]]), ncol = ncol(betaList[[ii]]), byrow = TRUE)
    centScaledSamps <- sqrtm(tmp) %*% (t(betaList[[ii]] - tmp1))
    baryList[[ii]] <- t(muBetas + sqrtSigBetas %*% centScaledSamps)
  }
  
  list(
    wasp = do.call(rbind, baryList)
  )
}


bary <- sampleBetas(results)
bary_res <- bind_rows(bary)
## TIME END HERE
end_time <- Sys.time()
tdiff <- end_time - start_time
fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/bary_center_time.rds")
save(tdiff, file = fname)

################################################################################
## Look at the results from the two methods ####################################
################################################################################
describe_posterior(as.data.frame(full_data_draws))
describe_posterior(as.data.frame(as.mcmc(bary_res)))

tmp_carrier <- exp(as.data.frame(as.mcmc(bary_res))[,5])
res_carrier <- describe_posterior(as.data.frame(tmp_carrier))
res_carrier$Median
res_carrier$CI_low
res_carrier$CI_high

## Times
load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/bary_center_time.rds")


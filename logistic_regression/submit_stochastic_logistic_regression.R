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
  qsub -q BIOSTAT -pe smp -2 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out -t 1-25 /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/mh_dnc.job

################################################################################
####### LOG INTO COMPUTE NODE BEFORE PROCEEDING ################################
################################################################################

################################################################################
#### Combine the results from MH Runs  #########################################
################################################################################
all_counts <- numeric(25)

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
  times[j] <- tdiff
}

mean(times) ## Minutes

######################################
#### Put all draws in single list ####
######################################

## TIME START HERE
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


################################################################################
#### Recenter Draws using simple method ########################################
################################################################################
## TIME START HERE
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


#####################
#### Diagnostics ####
#####################
full_data_draws <- full_data_draws %>% as.mcmc()
heidel.diag(full_data_draws)
raftery.diag(full_data_draws)

#################
#### Summary ####
#################
describe_posterior(as.data.frame(full_data_draws))

################################################################################
###### Combine using WASP approximation ########################################
################################################################################
## TIME START HERE
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

################################################################################
## Look at the results from the two methods ####################################
################################################################################
describe_posterior(as.data.frame(full_data_draws))
describe_posterior(as.data.frame(as.mcmc(bary_res)))

library(rstan)
library(glue)
library(dplyr)


## clean up folder to store partitions in
system("rm /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/*", intern = FALSE)
system("rm /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/*", intern = FALSE)

### Will want to partition and submit here
### Also include combining data in this file

## Number of partitions
npart <- 25

## Submit the partitioning script
sub_command <- glue(
  "#!/bin/bash 
  qsub -pe smp -2 -cwd -e /dev/null -o /dev/null /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partition.job"
)

system(sub_command, intern = FALSE)

# continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/")) < npart
# i <- 0
# 
# while (continue & i < 100) {
#   
#   Sys.sleep(10)
#   i <- 1 + 1
#   continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/")) < npary
#   
# }

## submit command
# sub_command <- glue(
#   "#!/bin/bash 
#   qsub -pe smp -2 -cwd -e /dev/null -o /dev/null -t 1-{npart} /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/slr.job"
# )
# 
# system(sub_command, intern = FALSE)


## RUN THIS INSTEAD ON LOGIN NODE
## qsub -pe smp -2 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out -t 1-25 /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/slr.job

## RUN THIS INSTEAD ON LOGIN NODE
## qsub -q BIOSTAT -pe smp -2 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out -t 1-25 /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/mh_dnc.job

## Check every minute if all the files are there before combining
# continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/")) < npart
# i <- 0
# 
# while (continue & i < 100) {
#   
#   Sys.sleep(10)
#   i <- 1 + 1
#   continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/")) < npary
#   
# }


####### LOG INTO COMPUTE NODE BEFORE PROCEEDING ################################




###### Combine everything
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

npart <- 25

results <- list()
for (i in 1:npart) {
  fname <- paste0("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/res{i}.rds")
  load(fname)
  results[[sid]] <- res[-(1:1000), ]
}
bary <- sampleBetas(results)

rname <- paste0("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/final_results/results_250.rds")
save(bary, file = rname)


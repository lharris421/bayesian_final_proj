library(parallel)
library(coda)
library(bayestestR)
library(adaptMCMC)
library(mvtnorm)
library(tidyverse)
library(magrittr)
library(glue)
library(matrixStats)
library(expm)
library(MASS)

########################
#### WASP functions ####
########################

computeBarycenter <- function (meanList, covList) {
  
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

#######################
#### Log Posterior ####
#######################

log_post_fun <- function(param, X, y, ratio, sigma, mu) {
  (ratio)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

#######################
#### Main function ####
#######################

inner_draws <- function(j, NN = 1e4) {
  
  # load data
  load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/large/partition{j}.rds"))
  
  #priors
  mu <- c(-1.61, rep(0, ncol(x0)))
  sd_temp <- sqrt(diag(var(x0)))
  scale <- c(1, sd_temp[2], 1, 1, sd_temp[5])  # need to make sure we only scale continuous vars
  sigma <- diag(c(40^2, 3^2 * scale)) 
  
  x0 <- cbind(rep(1, nrow(x0)), x0)
  colnames(x0)[1] <- "intercept"
  
  temp <- glm(y0 ~x0-1, family = "binomial")
  
  set.seed(2022)
  full_data_draws <-  MCMC(p = log_post_fun,
                           n = NN,
                           init = coef(temp),
                           acc.rate = 0.234,
                           X = x0, 
                           y = y0, 
                           ratio = nrep0,
                           mu = mu,
                           sigma = sigma)
  
  full_data_draws$samples
}

#########################
#### Run in parallel ####
#########################

start.time <- Sys.time()
K <- 25
cl <- makeCluster(min(detectCores(), 25))
clusterExport(
  cl,
  c("inner_draws", "log_post_fun")
)
clusterEvalQ(cl, {
  library(MASS)
  library(glue)
  library(adaptMCMC)
})

results <-  parLapply(
  cl,
  1:K,
  function(x) inner_draws(j = x, NN = 1e4)
)
stopCluster(cl)

########################
#### Remove burn-in ####
########################

remove_burnin <- function(x, burnin) {
  x[-(1:burnin),]
}

results <- lapply(results, remove_burnin, 1000)

##########################################
#### Combine using WASP approximation ####
##########################################

bary <- sampleBetas(results)
bary_res <- bind_rows(bary)
full_data_draws <- bary_res$wasp
colnames(full_data_draws) <- colnames(results[[1]])

## Summary
summary <- describe_posterior(as.data.frame(full_data_draws))
hpd <- hdi(as.data.frame(full_data_draws))

end.time <- Sys.time()
elapsed <- end.time-start.time

results <- tibble(var = summary[,1],
                  coef = summary[,2],
                  hpd_lower = hpd[,3],
                  hdp_upper = hpd[,4])

#####################
#### Diagnostics ####
#####################

full_data_draws <- full_data_draws %>% as.mcmc()
heidel.diag(full_data_draws)
raftery.diag(full_data_draws)

##############
#### Save ####
##############

out <- list(result_table = results,
            diagnostics = list(heidel = heidel.diag(full_data_draws),
                               raftery = raftery.diag(full_data_draws)),
            comp_time = elapsed)
save(out,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/normal_approx_dnc_WASP.Rdata")

temp <- out$result_table
temp %>% mutate(across(coef:hdp_upper, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  mutate(nice = paste0(coef, " (", hpd_lower, ", ", hdp_upper, ")"))



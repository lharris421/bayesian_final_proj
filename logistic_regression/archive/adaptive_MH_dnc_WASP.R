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

args <- commandArgs(trailingOnly = TRUE)
seed_in <- args[1]

set.seed(seed_in) 
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
elapsed <- difftime(end.time, start.time, units = "mins")

results <- tibble(var = summary[,1],
                  coef = summary[,2],
                  central_lower = summary[, 4],
                  central_upper= summary[, 5],
                  hpd_lower = hpd[,3],
                  hdp_upper = hpd[,4])

#####################
#### Diagnostics ####
#####################

full_data_draws <- full_data_draws %>% as.mcmc()
heidel.diag(full_data_draws)
raftery.diag(full_data_draws)

######################
#### Bayes Factor ####
######################

load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/full_data_2.rds")

score <- mod_dat %>% select_at(vars(CHF:Depression)) %>% rowSums()
full_data <- mod_dat %>%
  dplyr::select(proc, smoker, age, sex, carrier) %>%
  mutate(score = scale(score))

y <- full_data$proc * 1
X <- as.matrix(full_data)[, -1]

mu <- c(-1.61, rep(0, ncol(X)))
sd_temp <- sqrt(diag(var(X)))
scale <- c(1, sd_temp[2], 1, 1, sd_temp[5])  # need to make sure we only scale continuous vars
sigma <- diag(c(40^2, 3^2 * scale)) 

X <- cbind(rep(1, nrow(X)), X)
colnames(X)[1] <- "intercept"

set.seed(666)
BF <- bayesfactor_parameters(posterior = as.data.frame(full_data_draws[,5]), # Need posterior draws
                             prior = data.frame(rnorm(nrow(full_data_draws), mu[5], sqrt(sigma[5,5]))), # also need prior draws
                             null = 0) # did not converge
##############
#### Save ####
##############

out <- list(result_table = results,
            BF = NA, #as.numeric(BF),
            diagnostics = list(heidel = heidel.diag(full_data_draws),
                               raftery = raftery.diag(full_data_draws)),
            comp_time = elapsed)

save(out, full_data_draws,
     file = paste0("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/adaptive_MH_dnc_WASP_", 
                   seed_in, ".Rdata"))

temp <- out$result_table
temp %>% mutate(across(coef:central_upper, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  dplyr::select(var:central_upper) %>% 
  mutate(nice = paste0(coef, " (", central_lower, ", ", central_upper, ")"))



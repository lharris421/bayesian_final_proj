library(parallel)
library(coda)
library(bayestestR)
library(adaptMCMC)
library(mvtnorm)
library(tidyverse)
library(magrittr)
library(glue)


#######################
#### Log Posterior ####
#######################

log_post_fun <- function(param, X, y, ratio, sigma, mu) {
  (ratio)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

########################
#### Optim function ####
########################

optim_fun <- function(init, X, y, ratio, sigma, mu) {
  optim(par = init,
        fn = log_post_fun,
        method = "BFGS",
        control = list(fnscale= -1,
                       maxit = 1e6),
        hessian = T,
        sigma = sigma, 
        mu = mu,
        X = X,
        y = y,
        ratio = ratio)
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
  
  Opt <- optim_fun(init = rep(0, ncol(x0)),
                   X = x0, 
                   y = y0,
                   ratio = nrep0,
                   sigma = sigma,
                   mu = mu)
  
  params <- Opt$par
  names(params) <- colnames(x0)
  SigNew <- chol2inv(chol(-Opt$hessian))
  
  ## Draw from multivariate normal
  set.seed(666)
  MASS::mvrnorm(NN, mu =  params, Sigma = SigNew)
}

#########################
#### Run in parallel ####
#########################

start.time <- Sys.time()
K <- 25
cl <- makeCluster(min(detectCores(), 25))
clusterExport(
  cl,
  c("inner_draws", "log_post_fun",
    "optim_fun")
)
clusterEvalQ(cl, {
  library(MASS)
  library(glue)
})

results <-  parLapply(
  cl,
  1:K,
  function(x) inner_draws(j = x, NN = 1e4)
)
stopCluster(cl)

########################
#### Recenter Draws ####
########################

subset_mean <- t(sapply(1:K, function(i) colMeans(results[[i]]))) # rows indicate subset
## This step was slightly off... it assumed perfect balance
## global_mean <- colMeans(subset_mean)
global_mean <- colMeans(bind_rows(lapply(results, data.frame)))
recenter <- t(sapply(1:K, function(i) subset_mean[i, ] - global_mean)) # rows indicate subset
## Should this be minus?
results_recentered <- lapply(1:K,function(i) results[[i]] - matrix(recenter[i, ],
                                                                   nrow = nrow(results[[i]]),
                                                                   ncol = ncol(results[[i]]),
                                                                   byrow = T))
full_data_draws <- do.call(rbind, results_recentered)

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
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/adaptive_MH_dnc_WASP.Rdata")

temp <- out$result_table
temp %>% mutate(across(coef:hdp_upper, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  mutate(nice = paste0(coef, " (", hpd_lower, ", ", hdp_upper, ")"))



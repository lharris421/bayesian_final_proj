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

inner_draws <- function(j) {
  
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
  
  return(list(params = params,
              Sig = SigNew))
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
  function(x) inner_draws(j = x)
)
stopCluster(cl)

##########################################
#### Combine using WASP approximation ####
##########################################

mean_list <- lapply(1:K, function(i) results[[i]]$params)
Sigs_list <- lapply(1:K, function(i) results[[i]]$Sig)

bary <- computeBarycenter(mean_list, Sigs_list)
global_mean <- bary$mean
SigNew <- bary$cov

end.time <- Sys.time()
elapsed <- end.time-start.time

## Summary
res <- tibble()
for (i in 1:length(global_mean)){
  tmp <- global_mean[i]+c(-1,1)*qnorm(0.975)*sqrt(SigNew[i,i])
  res <-  bind_rows(res,
                    tibble(param = global_mean[i],
                           lower_CI = tmp[1], 
                           upper_CI = tmp[2]))
}

results <- tibble(parameter = names(global_mean)) %>% 
  bind_cols(res)

##############
#### Save ####
##############

out <- list(result_table = results,
            comp_time = elapsed)

save(out,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/normal_approx_dnc_WASP_V2.Rdata")

temp <- out$result_table
temp %>% mutate(across(param:upper_CI, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  dplyr::select(parameter:upper_CI) %>% 
  mutate(nice = paste0(param, " (", lower_CI, ", ", upper_CI, ")"))



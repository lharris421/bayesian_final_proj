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
  MASS::mvrnorm(NN, mu =  params, Sigma = SigNew)
}

############################
#### Parallel function #####
############################

par_fun <- function(){
  K <- 25
  cl <- makeCluster(min(detectCores(), 25))
  clusterExport(
    cl,
    c("inner_draws", "log_post_fun",
      "optim_fun"),
    envir = environment()
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
  gc()
  
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
  
  
  return(full_data_draws)
}

set.seed(666)
full_data_draws <- par_fun()

time <- microbenchmark::microbenchmark(par_fun(),
                                       unit = "s",
                                       times = 10)

elapsed <- tibble(summary(time)) %>% 
  dplyr::select(-expr, -neval) %>% 
  mutate(across(everything(), ~./60)) #convert to minutes

#################
#### Summary ####
#################

summary <- describe_posterior(as.data.frame(full_data_draws))
hpd <- hdi(as.data.frame(full_data_draws))

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
                             null = 0)
##############
#### Save ####
##############

out <- list(result_table = results,
            BF = as.numeric(BF),
            diagnostics = list(heidel = heidel.diag(full_data_draws),
                               raftery = raftery.diag(full_data_draws)),
            comp_time = elapsed)

save(out, full_data_draws,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/normal_approx_dnc.Rdata")

temp <- out$result_table
temp %>% mutate(across(coef:central_upper, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  dplyr::select(var:central_upper) %>% 
  mutate(nice = paste0(coef, " (", central_lower, ", ", central_upper, ")"))


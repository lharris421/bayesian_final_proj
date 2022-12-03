library(parallel)
library(coda)
library(bayestestR)
library(adaptMCMC)
library(mvtnorm)
library(tidyverse)
library(magrittr)

###################
#### Load Data ####
###################

load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/full_data_2.rds")

####################################
#### Augment data for algorithm ####
####################################

score <- mod_dat %>% select_at(vars(CHF:Depression)) %>% rowSums()
full_data <- mod_dat %>%
  dplyr::select(proc, smoker, age, sex, carrier) %>%
  mutate(score = scale(score))


y <- full_data$proc * 1
X <- as.matrix(full_data)[, -1]
N <- 1e4 ## Number of draws

########################################################
#### Set prior values / set up results df ##############
########################################################

mu <- c(-1.61, rep(0, ncol(X)))
sd_temp <- sqrt(diag(var(X)))
scale <- c(1, sd_temp[2], 1, 1, sd_temp[5])  # need to make sure we only scale continuous vars
sigma <- diag(c(40^2, 3^2 * scale)) 

##############################
#### Add intercept column ####
##############################

X <- cbind(rep(1, nrow(X)), X)
colnames(X)[1] <- "intercept"

#######################
#### Log Posterior ####
#######################

log_post_fun <- function(param) {
  
  sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

###########################
#### Run Normal Approx ####
###########################

run_norm <- function(){
  Opt <- optim(par = rep(0,ncol(X)),
               fn = log_post_fun, 
               method = "BFGS",
               control = list(fnscale= -1,
                              maxit = 1e6),
               hessian = T)
  
  params <- Opt$par
  names(params) <- colnames(X)
  SigNew <-  chol2inv(chol(-Opt$hessian))
  
  ## Draw from multivariate normal
  beta_draws <- MASS::mvrnorm(N, mu = params, Sigma = SigNew)
  return(beta_draws)
  
}

set.seed(666) 
beta_draws <- run_norm()

time <- microbenchmark::microbenchmark(run_norm(),
                                       unit = "s",
                                       times = 10)

elapsed <- tibble(summary(time)) %>% 
  dplyr::select(-expr, -neval) %>% 
  mutate(across(everything(), ~./60)) #convert to minutes

#################
#### Summary ####
#################

summary <- describe_posterior(as.data.frame(beta_draws))
hpd <- hdi(as.data.frame(beta_draws))

results <- tibble(var = summary[,1],
                  coef = summary[,2],
                  central_lower = summary[, 4],
                  central_upper= summary[, 5],
                  hpd_lower = hpd[,3],
                  hdp_upper = hpd[,4])

#####################
#### Diagnostics ####
#####################

beta_draws <- beta_draws %>% as.mcmc()
heidel.diag(beta_draws)
raftery.diag(beta_draws)

######################
#### Bayes Factor ####
######################

set.seed(666)
BF <- bayesfactor_parameters(posterior = as.data.frame(beta_draws[,5]), # Need posterior draws
                       prior = data.frame(rnorm(nrow(beta_draws), mu[5], sqrt(sigma[5,5]))), # also need prior draws
                       null = 0)
##############
#### Save ####
##############

out <- list(result_table = results,
            BF = as.numeric(BF),
            diagnostics = list(heidel = heidel.diag(beta_draws),
                               raftery = raftery.diag(beta_draws)),
            comp_time = elapsed)
save(out, beta_draws,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/normal_approx_full.Rdata")

temp <- out$result_table
temp %>% mutate(across(coef:central_upper, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  dplyr::select(var:central_upper) %>% 
  mutate(nice = paste0(coef, " (", central_lower, ", ", central_upper, ")"))





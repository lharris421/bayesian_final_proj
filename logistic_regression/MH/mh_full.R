## source("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/MH/mh_full.R")

## Start the timer!!!
start_time <- Sys.time()

######################################
#### Get info from qsub ##############
######################################
args = commandArgs(trailingOnly=TRUE)
str(args)
seed <- as.numeric(args)

#######################
#### Libraries ########
#######################
## rm(list = ls())
library(tidyverse)
library(magrittr)
library(glue)
library(bayestestR)
library(coda)

#######################
#### Load Data ########
#######################

load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/full_data_2.rds")

##################################################
#### Consideration for prior on intercept ########
##################################################

score <- mod_dat %>% select_at(vars(CHF:Depression)) %>% rowSums()
full_data <- mod_dat %>%
  dplyr::select(proc, smoker, age, sex, carrier) %>%
  mutate(score = scale(score))

####################################
#### Augment data for algorithm ####
####################################

y <- full_data$proc
X <- as.matrix(full_data)[, -1]
N <- 1e4 ## Number of draws

########################################################
#### Set prior values / set up results df ##############
########################################################
mu <- c(-1.61, rep(0, ncol(X)))
sd_temp <- sqrt(diag(var(X)))
scale <- c(1, sd_temp[2], 1, 1, sd_temp[5])  # need to make sure we only scale continuous vars
sigma <- diag(c(40^2, 3^2 * scale))

beta_draws_mh <- matrix(0.0, N, ncol(X) + 1)

glm_data <- bind_cols(y = y, X)
glm_fit <- glm(y ~ ., glm_data, family = "binomial")
beta_draws_mh[1, ] <- coef(glm_fit)

acc_count <- 0
sd_prop <- .25
proposal_cov <- diag(summary(glm_fit)$coefficients[,2]^2)*sd_prop

##############################
#### Add intercept column ####
##############################

X <- cbind(rep(1, nrow(X)), X)

#######################
#### Log Posterior ####
#######################

log_post_fun <- function(param) {
  
  sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

#######################
#### Run MH ###########
#######################

set.seed(seed)
for(i in 2:N) {
  
  beta_draws_mh[i,] <- beta_draws_mh[i - 1,]
  
  ## Draw beta
  beta_proposal <- MASS::mvrnorm(1, mu = beta_draws_mh[i,], Sigma = proposal_cov)
  
  acc_prob <- exp(log_post_fun(beta_proposal) - log_post_fun(beta_draws_mh[i,]))
  
  if (runif(1) < acc_prob) {
    
    beta_draws_mh[i,] <- beta_proposal
    acc_count <- acc_count + 1
    
  }
  
}

fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/res_full_{seed}.rds")
save(beta_draws_mh, file = fname)

## Save the time it took to run
end_time <- Sys.time()
tdiff <- difftime(end_time, start_time, units = "mins")
fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/mh_full_time_{seed}.rds")
save(tdiff, file = fname)


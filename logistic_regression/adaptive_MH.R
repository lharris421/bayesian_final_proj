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
NN <- 1e4 ## Number of draws

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

log_post_fun <- function(param, X, y) {
  
  sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

#########################
#### Run adaptive MH ####
#########################

start.time <- Sys.time()

#load glm full data
load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/glm.Rdata")

set.seed(2022)
beta_draws <-  MCMC(p = log_post_fun,
                    n = NN,
                    init = out$result_table$coef,
                    acc.rate = 0.234,
                    X = X, 
                    y = y)
rm(out)

beta_draws <- beta_draws$samples

########################
#### Remove burn-in ####
########################

remove_burnin <- function(x, burnin) {
  x[-(1:burnin),]
}

beta_draws <- remove_burnin(beta_draws, 1000)

#################
#### Summary ####
#################

summary <- describe_posterior(as.data.frame(beta_draws))
hpd <- hdi(as.data.frame(beta_draws))

end.time <- Sys.time()

elapsed <- end.time-start.time

results <- tibble(var = summary[,1],
                  coef = summary[,2],
                  hpd_lower = hpd[,3],
                  hdp_upper = hpd[,4])


#####################
#### Diagnostics ####
#####################

heidel.diag(beta_draws)
raftery.diag(beta_draws)

##############
#### Save ####
##############

out <- list(result_table = results,
            diagnostics = list(heidel = heidel.diag(beta_draws),
                               raftery = raftery.diag(beta_draws)),
            comp_time = elapsed)
save(out,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/adaptive_MH_full.Rdata")

temp <- out$result_table
temp %>% mutate(across(coef:hdp_upper, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  mutate(nice = paste0(coef, " (", hpd_lower, ", ", hdp_upper, ")"))

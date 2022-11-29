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

start.time <- Sys.time()

Opt <- optim(par = rep(0,ncol(X)),
             fn = log_post_fun, 
             method = "BFGS",
             control = list(fnscale= -1,
                            maxit = 1e6),
             hessian = T)

params <- Opt$par
names(params) <- colnames(X)
SigNew <-  chol2inv(chol(-Opt$hessian))

end.time <- Sys.time()

elapsed <- end.time-start.time

## Summary
res <- tibble()
for (i in 1:length(params)){
  tmp <- params[i]+c(-1,1)*qnorm(0.975)*sqrt(SigNew[i,i])
  res <-  bind_rows(res,
                    tibble(Mean = params[i],
                           lower_CI = tmp[1], 
                           upper_CI = tmp[2]))
  
}

results <- tibble(parameter = colnames(X)) %>% 
  bind_cols(res)

##############
#### Save ####
##############

out <- list(result_table = results,
            comp_time = elapsed)
save(out,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/normal_approx_full_V2.Rdata")

temp <- out$result_table
temp %>% mutate(across(Mean:upper_CI, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  dplyr::select(parameter:upper_CI) %>% 
  mutate(nice = paste0(Mean, " (", lower_CI, ", ", upper_CI, ")"))





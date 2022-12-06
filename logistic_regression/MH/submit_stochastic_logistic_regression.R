## Run on login node
## If needed to remove previous partitions:
## rm /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/*

## Submit the partitioning script:
## qsub -pe smp -2 -cwd -e /dev/null -o /dev/null /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partition.job


## Submit the job scripts:
## stan (currently thinking we get rid of this)
## qsub -pe smp -2 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out -t 1-25 /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/slr.job

## mh
## SAVE TIME HERE, COMPARE TO THE FINAL COMPLETION TIME FOR TOTAL RUN TIME
date > /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/start_time_10101.txt
qsub -q BIOSTAT -pe smp -2 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out -t 1-25 /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/MH/mh_dnc.job 10101

################################################################################
####### LOG INTO COMPUTE NODE BEFORE PROCEEDING ################################
################################################################################

################################################################################
#### Combine the results from MH Runs  #########################################
################################################################################

########################
#### Libraries #########
########################
library(glue)
library(magrittr)
library(bayestestR)
library(coda)
library(dplyr)
library(lubridate)
library(stringr)
library(readr)

###############################
#### Check acceptance rate ####
###############################
all_counts <- numeric(25)
seed <- 10101
for (j in 1:25) {
  load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/acc_counts/mh_large/acc_counts{j}_{seed}.rds"))
  all_counts[j] <- acc_counts
}

min(all_counts) / 10000
max(all_counts) / 10000
mean(all_counts) / 10000

###############################
#### Check average run time ###
###############################
source("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/MH/mh_dnc_timing.R")

summary(time_simple)
summary(time_wasp)


#####################
#### Load in a run ##
#####################
load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/simple_{seed}.rds"))
load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/wasp_{seed}.rds"))

#####################
#### Diagnostics ####
#####################
full_data_draws <- full_data_draws %>% as.mcmc()
heidel.diag(full_data_draws)
raftery.diag(full_data_draws)

#################
#### Summary ####
#################
tmp <- describe_posterior(as.data.frame(full_data_draws))

tmp_carrier <- exp(full_data_draws[,5])
res_carrier <- describe_posterior(as.data.frame(tmp_carrier))
res_carrier$Median
res_carrier$CI_low
res_carrier$CI_high

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

################################################################################
###### WASP approximation ######################################################
################################################################################
bary_res %<>% as.mcmc()
# BF <- bayesfactor_parameters(posterior = as.data.frame(bary_res[,5]), # Need posterior draws
#                              prior = data.frame(rnorm(nrow(bary_res), mu[5], sqrt(sigma[5,5]))), # also need prior draws
#                              null = 0)

################################################################################
## Look at the results from the two methods ####################################
################################################################################
describe_posterior(as.data.frame(full_data_draws))
describe_posterior(as.data.frame(as.mcmc(bary_res)))

tmp_carrier <- exp(as.data.frame(as.mcmc(bary_res))[,5])
res_carrier <- describe_posterior(as.data.frame(tmp_carrier))
res_carrier$Median
res_carrier$CI_low
res_carrier$CI_high



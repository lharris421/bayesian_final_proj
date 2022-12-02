qsub -q BIOSTAT -pe smp -8 -e /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/err -o /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/out /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/mh_full.job 1004

#######################
#### Timing ###########
#######################
## Load in all the timing files and look at the distribution



#######################
#### Results ##########
#######################
load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/res_full.rds")
acc_count / N ## Acceptance rate

tmp_carrier <- exp(as.data.frame(beta_draws_mh)[,5])
res_carrier <- describe_posterior(as.data.frame(tmp_carrier))
res_carrier$Median
res_carrier$CI_low
res_carrier$CI_high

#####################
#### Diagnostics ####
#####################
# full_data_draws <- beta_draws_mh %>% as.mcmc()
# heidel.diag(full_data_draws)
# raftery.diag(full_data_draws)

######################
#### Bayes Factor ####
######################
# load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/full_data_2.rds")
# 
# score <- mod_dat %>% select_at(vars(CHF:Depression)) %>% rowSums()
# full_data <- mod_dat %>%
#   dplyr::select(proc, smoker, age, sex, carrier) %>%
#   mutate(score = scale(score))
# 
# y <- full_data$proc * 1
# X <- as.matrix(full_data)[, -1]
# 
# mu <- c(-1.61, rep(0, ncol(X)))
# sd_temp <- sqrt(diag(var(X)))
# scale <- c(1, sd_temp[2], 1, 1, sd_temp[5])  # need to make sure we only scale continuous vars
# sigma <- diag(c(40^2, 3^2 * scale)) 
# 
# X <- cbind(rep(1, nrow(X)), X)
# colnames(X)[1] <- "intercept"
# 
# set.seed(666)
# BF <- bayesfactor_parameters(posterior = as.data.frame(full_data_draws[,5]), # Need posterior draws
#                              prior = data.frame(rnorm(nrow(full_data_draws), mu[5], sqrt(sigma[5,5]))), # also need prior draws
#                              null = 0)

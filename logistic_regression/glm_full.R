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

glm_fit <- glm(y ~ X, family = "binomial")

time <- microbenchmark::microbenchmark(glm(y ~ X, family = "binomial"),
                                       unit = "s",
                                       times = 10)

elapsed <- tibble(summary(time)) %>% 
  dplyr::select(-expr, -neval) %>% 
  mutate(across(everything(), ~./60)) #convert to minutes

results <- tibble(var = c("intercept", colnames(X)),
                          coef = coef(glm_fit)) %>%
  bind_cols(broom::confint_tidy(glm_fit, func = confint.default))  

######################
#### Save Results ####
######################

out <- list(result_table = results,
            comp_time = elapsed)
save(out,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/glm.Rdata")

temp <- out$result_table
temp %>% mutate(across(coef:conf.high, ~format(round(exp(.), 3), nsmall = 3))) %>% 
  mutate(nice = paste0(coef, " (", conf.low, ", ", conf.high, ")"))

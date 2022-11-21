library(parallel)
library(coda)
library(bayestestR)
library(adaptMCMC)
library(mvtnorm)
library(tidyverse)
library(magrittr)

  
#######################
#### Simulate Data ####
#######################

set.seed(666)
N <- 1e4
x1 <- rnorm(N, sd = 3)           # some continuous variables 
x2 <- rnorm(N , sd = 2)
x3 <-rnorm(N)
beta <- c(-3, 3.8, 1.1, 2.3) # True Betas
X <- cbind(1, x1, x2, x3)
eta <- X %*% beta      
pr <- 1 / (1 + exp(-eta))         # pass through an inv-logit function
y <- rbinom(N, 1, pr)  

#######################
#### Log Posterior ####
#######################

log_post_fun <- function(param, X, y, N, sigma, mu) {
  sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

#########################
#### Run adaptive MH ####
#########################
#priors
sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(X[,-1])))))
mu <- rep(0,4)

temp <- glm(y ~ x1 + x2 + x3, family = "binomial")
NN <- 1e4 # posterior draws

set.seed(2022)
full_data_draws <-  MCMC(p = log_post_fun,
                         n = NN,
                         init = coef(temp),
                         acc.rate = 0.234,
                         X = X, 
                         y = y, 
                         N = N,
                         mu = mu,
                         sigma = sigma)

#####################
#### Diagnostics ####
#####################

# Looks bad
full_data_draws <- as.mcmc(full_data_draws$samples)
heidel.diag(full_data_draws)
raftery.diag(full_data_draws)
plot(full_data_draws)

#################
#### Summary ####
#################
describe_posterior(as.data.frame(full_data_draws))

## Truth
beta




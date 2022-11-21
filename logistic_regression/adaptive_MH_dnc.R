library(parallel)
library(coda)
library(bayestestR)
library(adaptMCMC)
library(mvtnorm)
library(tidyverse)
library(magrittr)


## Divide and conqure with multivariate normal prior
  
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

log_post_fun <- function(param, X, y, N, m, sigma, mu) {
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

#######################
#### Main function ####
#######################

inner_draws <- function(i, X, y, N, NN = 1e4) {
  fold_data_y <- y[fold_idx == i]
  fold_data_X <- X[fold_idx == i,]
  
  #priors
  temp <- glm(fold_data_y ~ fold_data_X-1, family = "binomial")
  sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(fold_data_X[,-1]))))) # need to make sure we only scale continuous vars
  mu <- rep(0,4)
  
  full_data_draws <-  MCMC(p = log_post_fun,
                           n = NN,
                           init = coef(temp),
                           acc.rate = 0.234,
                           X = X, 
                           y = y, 
                           N = N,
                           m = sum(fold_idx == i),
                           mu = mu,
                           sigma = sigma)
  
  full_data_draws$samples
}

########################
#### Partition Data ####
########################
# Create folds
K <- 4
## Why create folds this way? Why not force a balance? YOU ARE CORRECT!
## i.e. (sample(1:N, replace = FALSE)) %% 4 + 1
# fold_idx <- sample(1:K, size = N, replace = TRUE)
fold_idx <- (sample(1:N, replace = FALSE)) %% 4 + 1

cl <- makeCluster(min(detectCores(), K))
clusterExport(
  cl,
  c("X", "y", "fold_idx", "N",
    "inner_draws", "log_post_fun")
)
clusterEvalQ(cl, {
  library(MASS)
  library(adaptMCMC)
})

results <-  parLapply(
  cl,
  1:K,
  function(x) inner_draws(i = x, X = X, y = y, N = N, NN = 1e4)
)
stopCluster(cl)

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

#####################
#### Diagnostics ####
#####################

# Looks bad
full_data_draws <- full_data_draws %>% as.mcmc()
plot(full_data_draws)
heidel.diag(full_data_draws)
raftery.diag(full_data_draws)

# Apply individually 
check_heidel <- function(x) {
  
  x %>% as.mcmc() %>% heidel.diag()
  
}

check_raftery <- function(x) {
  
  x %>% as.mcmc() %>% raftery.diag()
  
}

lapply(results_recentered, check_heidel)
lapply(results_recentered, check_raftery)

#################
#### Summary ####
#################
describe_posterior(as.data.frame(full_data_draws))

## Truth
beta




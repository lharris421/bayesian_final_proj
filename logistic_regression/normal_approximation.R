library(parallel)
library(coda)
library(bayestestR)
library(adaptMCMC)
library(mvtnorm)
library(tidyverse)
library(magrittr)

  

## Simulate Data
set.seed(666)
N <- 1e4
x1 <- rnorm(N)           # some continuous variables 
x2 <- rnorm(N)
x3 <-rnorm(N)
beta <- c(-2, 0.11,1.34,2.3)
X <-  cbind(1, x1,x2,x3)
eta <- X %*% beta      
pr <- 1/(1+exp(-eta))         # pass through an inv-logit function
y <- rbinom(N,1,pr)  


## Flat prior
  
log_post_fun <- function(param){
  sum(y*X%*%param - log(1 + exp(X%*%param)))
}  

Opt <- optim(par = c(10,0,2,1),
             fn = log_post_fun, 
             # gr = gradient_fun,
             method = "BFGS",
             control = list(fnscale= -1,
                            maxit = 1e6),
             hessian = T)

(Opt$par)


## Multivariate normal prior
mu <- rep(0,4)
sigma <- diag(c(40^2, 3^2 *sqrt(diag(var(X[,-1])))))

log_post_fun <- function(param) {
  sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

Opt <- optim(par = rep(0,4),
             fn = log_post_fun, 
             method = "BFGS",
             control = list(fnscale= -1,
                            maxit = 1e6),
             hessian = T)

params <- Opt$par
SigNew <-  chol2inv(chol(-Opt$hessian))

## Draw from multivariate normal
beta_draws <- MASS::mvrnorm(1e4, mu = params, Sigma = SigNew)

## 
beta
describe_posterior(as.data.frame(beta_draws))


## Divide and conqure with multivariate normal prior
  
#######################
#### Simulate Data ####
#######################

set.seed(666)
N <- 1e4
x1 <- rnorm(N)           # some continuous variables 
x2 <- rnorm(N)
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

########################
#### Optim function ####
########################

optim_fun <- function(init, X, y, N, m, sigma, mu) {
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
        N = N,
        m = m)
}

#######################
#### Main function ####
#######################

inner_draws <- function(i, X, y, N, NN = 1e4) {
  fold_data_y <- y[fold_idx == i]
  fold_data_X <- X[fold_idx == i,]
  
  #priors
  sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(fold_data_X[,-1])))))
  mu <- rep(0,4)
  
  Opt <- optim_fun(init = rep(0, ncol(fold_data_X)),
                   X = fold_data_X, 
                   y = fold_data_y,
                   m = sum(fold_idx == i),
                   N = N,
                   sigma = sigma,
                   mu = mu)
  
  params <- Opt$par
  SigNew <- chol2inv(chol(-Opt$hessian))
  
  ## Draw from multivariate normal
  ## beta_draws <- MASS::mvrnorm(NN, mu =  params, Sigma = SigNew)
  ## return(beta_draws)
  set.seed(666) # Added seed
  MASS::mvrnorm(NN, mu =  params, Sigma = SigNew)
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
    "inner_draws", "log_post_fun",
    "optim_fun")
)
clusterEvalQ(cl, {
  library(MASS)
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




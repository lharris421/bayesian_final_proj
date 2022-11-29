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
  sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(fold_data_X[,-1]))))) # need to make sure we only scale continuous vars
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
  
  return(list(params = params,
              Sig = SigNew))
}

########################
#### Partition Data ####
########################
# Create folds
K <- 4
## Why create folds this way? Why not force a balance? YOU ARE CORRECT!
## i.e. (sample(1:N, replace = FALSE)) %% 4 + 1
# fold_idx <- sample(1:K, size = N, replace = TRUE)
fold_idx <- (sample(1:N, replace = FALSE)) %% K + 1

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

global_mean <- colMeans(t(sapply(1:K, function(i) rbind(results[[i]]$params))))
Sigs <- lapply(1:K, function(i) results[[i]]$Sig)
SigNew <- Reduce('+', Sigs)
SigNew <- SigNew/K

#################
#### Summary ####
#################

res <- tibble()
for (i in 1:length(global_mean)){
  tmp <- global_mean[i]+c(-1,1)*qnorm(0.975)*sqrt(SigNew[i,i])
  res <-  bind_rows(res,
                    tibble(param = global_mean[i],
                           lower_CI = tmp[1], 
                           upper_CI = tmp[2]))
         
}

## Truth
beta


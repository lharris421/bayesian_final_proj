#######################
#### Libraries ########
#######################
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

## Prior values
mu <- rep(0,4)
sigma <- diag(c(40^2, 3^2 *sqrt(diag(var(X[,-1])))))  # need to make sure we only scale continuous vars

log_post_fun <- function(param) {
  sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

beta_draws_mh <- matrix(0.0, N, ncol(X),
                       dimnames = list(NULL,colnames(X)))

glm_fit <- glm(y ~ x1 + x2 + x3, family = "binomial")
beta_draws_mh[1, ] <- coef(glm_fit)

acc_count <- 0
sd_prop <- .85
proposal_cov <- diag(summary(glm_fit)$coefficients[,2]^2)*sd_prop
set.seed(666)
pb <- txtProgressBar(0, N, style=3)

for(i in 2:N){
  
  beta_draws_mh[i,] <- beta_draws_mh[i - 1,]
    
  ## Draw beta
  beta_proposal <- MASS::mvrnorm(1, mu = beta_draws_mh[i,], Sigma = proposal_cov)
      
  acc_prob = exp(log_post_fun(beta_proposal) - log_post_fun(beta_draws_mh[i,]))
      
  if(runif(1) < acc_prob) {
    
    beta_draws_mh[i,] <- beta_proposal
    acc_count <- acc_count + 1
  
  }
  
  setTxtProgressBar(pb, i)
  
}

acc_count / N
beta
describe_posterior(as.data.frame(beta_draws_mh))

#######################
#### D & C ############
#######################
sd_prop <- .25

log_post_fun <- function(param, X, y, N, m, sigma, mu) {
  
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

inner_draws <- function(i, X, y, N, NN = 1e4) {
  
  fold_data_y <- y[fold_idx == i]
  fold_data_X <- X[fold_idx == i,]
  
  #priors
  sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(fold_data_X[,-1])))))
  mu <- rep(0,4)
  
  beta_draws_mh <- matrix(0.0, NN, ncol(X), dimnames = list(NULL,colnames(X)))
  
  glm_data <- bind_cols(y = fold_data_y, fold_data_X[,-1])
  glm_fit <- glm(y ~ ., glm_data, family = "binomial")
  beta_draws_mh[1, ] <- coef(glm_fit)
  
  acc_count <- 0
  
  proposal_cov <- diag(summary(glm_fit)$coefficients[,2]^2)*sd_prop
  
  for(i in 2:NN){
    
    beta_draws_mh[i,] <- beta_draws_mh[i - 1,]
    
    ## Draw beta
    beta_proposal <- MASS::mvrnorm(1, mu = beta_draws_mh[i,], Sigma = proposal_cov)
    
    acc_prob = exp(
      log_post_fun(beta_proposal, fold_data_X, fold_data_y, N, m = length(fold_data_y), sigma, mu) - 
      log_post_fun(beta_draws_mh[i,], fold_data_X, fold_data_y, N, m = length(fold_data_y), sigma, mu)
    )
    
    if(runif(1) < acc_prob) {
      
      beta_draws_mh[i,] <- beta_proposal
      acc_count <- acc_count + 1
      
    }
    
  }
  
  return(list(beta_draws_mh, acc_count))
  
}


# Create folds
set.seed(666)
K <- 4
fold_idx <- sample(1:K, size = N, replace = TRUE)

cl <- makeCluster(min(detectCores(), K))
clusterExport(
  cl,
  c("X", "y", "fold_idx", "N",
    "inner_draws", "log_post_fun", "sd_prop")
)
clusterEvalQ(cl, {
  library(MASS)
  library(dplyr)
})

results <-  parLapply(
  cl,
  1:K,
  function(x) inner_draws(i = x, X = X, y = y, N = N, NN = 1e4)
)
stopCluster(cl)

str(results)

sapply(results, function(x) x[[2]]) / 1e4

results <- lapply(results, function(x) x[[1]])

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


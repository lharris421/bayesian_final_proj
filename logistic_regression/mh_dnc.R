##
library(glue)
library(dplyr)

## Get index number
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  j <- 1
} else {
  j <- args[[1]]
}

start_time <- Sys.time()

## Load the subset of the data
load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/partition{j}.rds"))

## Call individual submission
## Need to make this into a function that submits to argon
sd_prop <- .012
N <- round(nrep0 * nrow(x0))

log_post_fun <- function(param, X, y, N, m, sigma, mu) {
  
  X <- cbind(rep(1, nrow(X)), X)
  
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

inner_draws <- function(X, y, N, NN = 1e4) {
  
  #priors
  sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(X)))))
  mu <- rep(0, ncol(X) + 1)
  
  beta_draws_mh <- matrix(0.0, NN, (ncol(X) + 1))
  
  glm_data <- bind_cols(y = y, X)
  glm_fit <- glm(y ~ ., glm_data, family = "binomial")
  beta_draws_mh[1, ] <- coef(glm_fit)
  
  acc_count <- 0
  
  proposal_cov <- diag(summary(glm_fit)$coefficients[,2]^2)*sd_prop
  
  for(i in 2:NN){
    
    beta_draws_mh[i,] <- beta_draws_mh[i - 1,]
    
    ## Draw beta
    beta_proposal <- MASS::mvrnorm(1, mu = beta_draws_mh[i,], Sigma = proposal_cov)
    
    acc_prob = exp(
      log_post_fun(beta_proposal, X, y, N, m = length(y), sigma, mu) - 
        log_post_fun(beta_draws_mh[i,], X, y, N, m = length(y), sigma, mu)
    )
    
    if(runif(1) < acc_prob) {
      
      beta_draws_mh[i,] <- beta_proposal
      acc_count <- acc_count + 1
      
    }
    
    if (i %% 100 == 0) {print(i)}
    
  }
  
  return(list(beta_draws_mh, acc_count))
  
}


res_full <- inner_draws(x0, y0, N)

## res_full[[2]] / 10000

## Check
## library(bayestestR)
## describe_posterior(as.data.frame(resStan))
res <- res_full[[1]]
acc_counts <- res_full[[2]]

## describe_posterior(as.data.frame(res))

fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh/res{j}.rds")
save(res, file = fname)

fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/acc_counts/acc_counts{j}.rds")
save(acc_counts, file = fname)

end_time <- Sys.time()
tdiff <- end_time - start_time
fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh/time{j}.rds")
save(tdiff, file = fname)

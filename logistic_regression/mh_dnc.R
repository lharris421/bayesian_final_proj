##
library(glue)

## Get index number
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  i <- 1
} else {
  i <- args[[1]]
}

## Load the subset of the data
load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/partition{i}.rds"))

## Call individual submission
## Need to make this into a function that submits to argon
sd_prop <- .25

log_post_fun <- function(param, X, y, N, m, sigma, mu) {
  
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

inner_draws <- function(X, y, N, NN = 1e5) {
  
  #priors
  sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(X[,-1])))))
  mu <- rep(0,4)
  
  beta_draws_mh <- matrix(0.0, NN, ncol(X), dimnames = list(NULL,colnames(X)))
  
  glm_data <- bind_cols(y = y, X[,-1])
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
    
  }
  
  return(list(beta_draws_mh, acc_count))
  
}




## Check
## library(bayestestR)
## describe_posterior(as.data.frame(resStan))


res <- do.call(cbind, resStan@sim$samples[[1]][1:ncol(x0)])
## describe_posterior(as.data.frame(res))

fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/res{i}.rds")
save(res, file = fname)

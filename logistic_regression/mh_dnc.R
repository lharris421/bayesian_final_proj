######################
#### Start timing ####
######################
## start_time <- Sys.time()

###################
#### Libraries ####
###################
library(glue)
library(dplyr)
library(stringr)

######################################
#### Get index from qsub #############
######################################
args = commandArgs(trailingOnly=TRUE)
str(args)
args <- as.numeric(str_split(str_remove(args, "^--args\s"), " "))
print(args)
print(j)
print(seed)
str(args)
if (length(args) == 0){
  j <- 1
} else {
  j <- args[1]
  seed <- args[2]
}

######################################
#### Load the subset of the data #####
######################################
load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/large/partition{j}.rds"))


sd_prop <- .012 ## Controls step size of MH algorithm
N <- round(nrep0 * nrow(x0)) ## Full dataset size, needed for log posterior

################################
#### Log posterior function ####
################################

log_post_fun <- function(param, X, y, N, m, sigma, mu) {
  
  X <- cbind(rep(1, nrow(X)), X)
  
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

#####################
#### Code for MH ####
#####################

inner_draws <- function(X, y, N, NN = 1e4) {
  
  # priors
  mu <- c(-1.61, rep(0, ncol(X)))
  sd_temp <- sqrt(diag(var(X)))
  scale <- c(1, sd_temp[2], 1, 1, sd_temp[5])  # need to make sure we only scale continuous vars
  sigma <- diag(c(40^2, 3^2 * scale)) 
  
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
    
    ## if (i %% 100 == 0) {print(i)}
    
  }
  
  return(list(beta_draws_mh, acc_count))
  
}


######################
#### Run MH ##########
######################
set.seed(seed)
res_full <- inner_draws(x0, y0, N)

######################
#### Save results ####
######################
res <- res_full[[1]]
acc_counts <- res_full[[2]]

fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/res{j}_{seed}.rds")
save(res, file = fname)

fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/acc_counts/mh_large/acc_counts{j}_{seed}.rds")
save(acc_counts, file = fname)

end_time <- Sys.time()
## tdiff <- end_time - start_time

## Only save the end time
fname <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/time{j}_{seed}.rds")
save(end_time, file = fname)

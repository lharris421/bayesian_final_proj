library(parallel)
library(coda)
library(bayestestR)
library(mvtnorm)
library(tidyverse)
library(magrittr)
library(microbenchmark)

#######################
#### Simulate Data ####
#######################

set.seed(666)
N <- 1e5
x1 <- rnorm(N)           # some continuous variables 
x2 <- rnorm(N)
x3 <-rnorm(N)
beta <- c(-3, 3.8, 1.1, 2.3) # True Betas
X <- cbind(1, x1, x2, x3)
eta <- X %*% beta      
pr <- 1 / (1 + exp(-eta))         # pass through an inv-logit function
y <- rbinom(N, 1, pr)  

sd_prop <- .012 ## Controls step size of MH algorithm

################################
#### Log posterior function ####
################################

log_post_fun <- function(param, X, y, N, m, sigma, mu) {
  
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

#####################
#### Code for MH ####
#####################

inner_draws <- function(fold_idx, i, x, y, N, NN = 1e4) {
  
  y <- y[fold_idx == i]
  X <- x[fold_idx == i,]
  
  # priors
  mu <- rep(0, ncol(X))
  sd_temp <- sqrt(diag(var(X[,-1])))
  sigma <- diag(c(40^2, 3^2 * sd_temp)) 
  
  beta_draws_mh <- matrix(0.0, NN, ncol(X))
  
  glm_data <- bind_cols(y = y, X)
  glm_fit <- glm(y ~ . -1, glm_data, family = "binomial")
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
    
    if (runif(1) < acc_prob) {
      
      beta_draws_mh[i,] <- beta_proposal
      acc_count <- acc_count + 1
      
    }
    
    if (i %% 100 == 0) {print(i)}
    
  }
  
  ## return(list(beta_draws_mh, acc_count))
  beta_draws_mh
  
}

#######################
#### Rep function #####
#######################

rep_fun <- function(X, y, fold_idx, N){
  
  cl <- makeCluster(min(detectCores(), K))
  
  clusterExport(
    cl,
    c("X", "y", "fold_idx", "N",
      "inner_draws", "log_post_fun", "sd_prop"),
    envir = environment()
  )
  
  clusterEvalQ(cl, {
    library(MASS)
    library(dplyr)
  })
  
  results <-  parLapply(
    cl,
    1:K,
    function(x) inner_draws(fold_idx = fold_idx,
                            i = x, x = X, y = y, N = N, NN = 1e4)
  )
  
  stopCluster(cl)
  
  ########################
  #### Recenter Draws ####
  ########################
  
  subset_mean <- t(sapply(1:K, function(i) colMeans(results[[i]]))) # rows indicate subset
  global_mean <- colMeans(bind_rows(lapply(results, data.frame)))
  recenter <- t(sapply(1:K, function(i) subset_mean[i, ] - global_mean)) # rows indicate subset
  results_recentered <- lapply(1:K,function(i) results[[i]] - matrix(recenter[i, ],
                                                                     nrow = nrow(results[[i]]),
                                                                     ncol = ncol(results[[i]]),
                                                                     byrow = T))
  full_data_draws <- do.call(rbind, results_recentered)
  
  return(full_data_draws)
}


########################
#### Partition Data ####
########################

out <- tibble()
K <- 1

repeat{
  # Create folds
  fold_idx <- (sample(1:N, replace = FALSE)) %% K + 1
  m_j <- max(table(fold_idx))
  
  if (m_j<500){
    break
  }

  full_data_draws <- rep_fun(X = X, y = y, 
                             fold_idx = fold_idx, N = N)
  
  gc()
  
  time <- microbenchmark::microbenchmark(rep_fun(X = X, y = y, 
                                                 fold_idx = fold_idx, N = N),
                                         unit = "s",
                                         times = 10)  
  
  gc()
  
  #################
  #### Summary ####
  #################
  summary <- describe_posterior(as.data.frame(full_data_draws))
  
  tmp <- tibble(K = K, 
                m_j = m_j, 
                coef = summary[2,2],
                lower_ci = summary[2,4],
                upper_ci = summary[2,5])
  
  tmp <- bind_cols(tmp, 
                   tibble(summary(time)) %>% 
                     dplyr::select(-expr))
  
  out <- bind_rows(out, tmp)
  
  K <- K+1
  
  if(K %% 5 == 0){
    print(K)
  }
}

save(out,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/mh_sim_plot_data.Rdata")

##############
#### Plot ####
##############

load("/Volumes/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/mh_sim_plot_data.Rdata")

N <- 1e5

# Plot of estimate and CI by K
as_tibble(out) %>% 
  select(-m_j) %>% 
  left_join(as_tibble(out) %>% select(K, M_J=m_j) %>% 
              slice(c(1, seq(0, nrow(out), by = 5))) %>% 
              mutate(m_j = paste0(M_J, " (", trimws(format(round(M_J/N *100, 2), nsmall =2)), "%)")),
            by = "K") %>% 
  mutate(m_j = ifelse(is.na(m_j), "", m_j)) %>% 
  mutate(K = as.factor(K)) %>% 
  ggplot(aes(K, coef, label = m_j))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), colour="black", width=.1)+
  geom_point(size=3, shape=21, fill="blue")+
  geom_hline(yintercept = 3.8, col = "red") +
  geom_text(aes(y = upper_ci), angle = 90, nudge_y = +0.025,
            fontface = "bold", size = 3.5)+
  scale_x_discrete(labels = c(1, seq(0, nrow(out), by = 5)),
                   breaks = c(1, seq(0, nrow(out), by = 5)))+
  annotate("text", x = 197, y=3.805, label = "Truth", col = "red",
           fontface = "bold", size = 3.5)+ 
  expand_limits(x = -2, y = 3.7)+
  expand_limits(x = 202, y = 3.7)+
  ylim(c(3.7,4.05))+
  ylab("Estimate")+
  ggtitle("Estimate and 95% Cred Int by K - Metropolis Hastings")

# Plot of computation time by K
as_tibble(out) %>% 
  select(-m_j) %>% 
  mutate(K = as.factor(K)) %>% 
  ggplot(aes(K, mean))+
  geom_point(size=3, shape=21, fill="blue")+
  geom_errorbar(aes(ymin=min, ymax=max), colour="black", width=.1)+
  scale_x_discrete(labels = c(1, seq(0, nrow(out), by = 5)),
                   breaks = c(1, seq(0, nrow(out), by = 5)))+
  expand_limits(x = -2, y = 2.5)+
  expand_limits(x = 202, y = 2.5)+
  ylab("Computation time (seconds)")+
  ggtitle("Computation time by K - Metropolis Hastings")






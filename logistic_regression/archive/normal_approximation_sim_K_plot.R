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
N <- 1e5
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

log_post_fun <- function(param, X, y, N, m
                         #, sigma, mu
) {
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param)))
  # -
  #   drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu)) # flat prior
}

########################
#### Optim function ####
########################

optim_fun <- function(init, X, y, N, m
                      #, sigma, mu
) {
  optim(par = init,
        fn = log_post_fun,
        method = "BFGS",
        control = list(fnscale= -1,
                       maxit = 1e6),
        hessian = T,
        # sigma = sigma, 
        # mu = mu,
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
  
  #priors - make it flat prir
  # sigma <- diag(c(40^2, 3^2 * sqrt(diag(var(fold_data_X[,-1]))))) # need to make sure we only scale continuous vars
  # mu <- rep(0,4)
  
  Opt <- optim_fun(init = rep(0, ncol(fold_data_X)),
                   X = fold_data_X, 
                   y = fold_data_y,
                   m = sum(fold_idx == i),
                   N = N
                   # ,
                   # sigma = sigma,
                   # mu = mu
  )
  
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

out <- matrix(0, nrow = 1, ncol = 6,
              dimnames = list(rep("", 1),
                              c("K", "m_j",
                                "coef",
                                "lower_ci",
                                "upper_ci",
                                "comp_time")))
K <- 1
repeat{
  # Create folds
  fold_idx <- (sample(1:N, replace = FALSE)) %% K + 1
  m_j <- max(table(fold_idx))
  if (m_j<500){
    break
  }
  start.time <- Sys.time()
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
  global_mean <- colMeans(bind_rows(lapply(results, data.frame)))
  recenter <- t(sapply(1:K, function(i) subset_mean[i, ] - global_mean)) # rows indicate subset
  results_recentered <- lapply(1:K,function(i) results[[i]] - matrix(recenter[i, ],
                                                                     nrow = nrow(results[[i]]),
                                                                     ncol = ncol(results[[i]]),
                                                                     byrow = T))
  full_data_draws <- do.call(rbind, results_recentered)
  end.time <- Sys.time()
  elapsed <- end.time-start.time
  gc()
  #################
  #### Summary ####
  #################
  summary <- describe_posterior(as.data.frame(full_data_draws))
  
  if(K==1){
    out[K, ] <- c(K, m_j, summary[2,2], summary[2,4], summary[2,5], elapsed)
  } else{
    out <- rbind(out, c(K, m_j, summary[2,2], summary[2,4], summary[2,5], elapsed))
  }
  K <- K+1
}

save(out,
     file = "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/normal_approx_sim_plot_data.Rdata")

##############
#### Plot ####
##############

load("/Volumes/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/normal_approx_sim_plot_data.Rdata")

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
  expand_limits(x = -5, y = 3.7)+
  expand_limits(x = 205, y = 3.7)+
  ylim(c(3.7,4.05))+
  ylab("Estimate")

# Plot of computation time by K
as_tibble(out) %>% 
  select(-m_j) %>% 
  left_join(as_tibble(out) %>% select(K, M_J=m_j) %>% 
              slice(c(1, seq(0, nrow(out), by = 5))) %>% 
              mutate(m_j = paste0(M_J, " (", trimws(format(round(M_J/N *100, 2), nsmall =2)), "%)")),
            by = "K") %>% 
  mutate(m_j = ifelse(is.na(m_j), "", m_j)) %>% 
  mutate(K = as.factor(K)) %>% 
  ggplot(aes(K, comp_time, group =1))+
  geom_line(col = "blue") +
  geom_point(size=3, shape=21, fill="blue")+
  scale_x_discrete(labels = c(1, seq(0, nrow(out), by = 5)),
                   breaks = c(1, seq(0, nrow(out), by = 5)))+
  expand_limits(x = -2, y = 2.5)+
  expand_limits(x = 202, y = 2.5)+
  ylab("Computation time (seconds)")



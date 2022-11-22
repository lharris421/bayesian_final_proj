Normal Approximation
================
2022-11-18

**Likelihood (non-D&C)**

$$
L_i(\boldsymbol{\beta}) = \Big(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{y_i}\Big(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{1-y_i}
$$

**Multivariate normal prior on **β****

$$
\pi({\boldsymbol{\beta}}) \propto \exp\Big\\{ -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\Big\\}
$$

**Posterior (non-D&C)**

$$
\pi({\boldsymbol{\theta}}\|y_i) \propto \Big\[\Big(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{ y_i}\Big(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{1- y_i} \Big\] \cdot \Big\[\exp\Big\\{ -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\Big\\}  \Big\]
$$

*Log-posterior (non-D&C)*

$$
\begin{aligned}
\log(\pi({\boldsymbol{\theta}}\|y_i) ) &\propto y_i \log(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}})+(1-y_i) \log(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}) -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
&\propto  y_i\log(\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})- y_i\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})+ y_i \log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
&\propto y_i\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}} -\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
\end{aligned}
$$

Therefore,

$$
\begin{aligned}
\log(\pi({\boldsymbol{\theta}}\|\mathrm{\bf{y}}) ) 
= \sum\_{i=1}^n (y_i\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}} - \log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})) - \frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu}) + C\\\\
\end{aligned}
$$

for subset j

$$
\begin{aligned}
\log(\pi({\boldsymbol{\theta}}\|\mathrm{y}\_j) ) 
= \Big(\frac{n}{m_j}\Big) \sum\_{i=1}^{m_j} (y_i\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}} - \log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})) - \frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu}) + C\\\\
\end{aligned}
$$

**Solve for MLE**

``` r
## Simulate Data
set.seed(666)
N <- 1e4
x1 <- rnorm(N, sd = 3)           # some continuous variables 
x2 <- rnorm(N, sd = 10)
x3 <-rnorm(N)
beta <- c(-2, 0.11,1.34,2.3)
X <-  cbind(1, x1,x2,x3)
eta <- X %*% beta      
pr <- 1/(1+exp(-eta))         # pass through an inv-logit function
y <- rbinom(N,1,pr)  
```

**Posterior with flat prior check to see if we get beta from optim**

``` r
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
```

    ## [1] -2.2037901  0.1233374  1.4126835  2.4478958

**Posterior with multivariate normal prior (non-D&C)**

``` r
mu <- rep(0,4)

# need to only scale continuous vars in design matrix
# one way to identify dummy var columns in design matrix is to count unique in column and if > 2
# then it is continuous but could be computational expensive
# especially if n is large,

scale_vars <- sapply(1:ncol(X), function(i) length(unique(X[,i])) > 2)
var_scale_vars <- diag(var(X[,scale_vars]))
scaled <- rep(1, ncol(X))
scaled[which(scale_vars)] <- var_scale_vars
sigma <- diag(c(40^2, rep(3^2, ncol(X)-1)) *sqrt(scaled))

# Orto make it faster we can just hard code priors which is what we shoud do
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

## Summary
beta
```

    ## [1] -2.00  0.11  1.34  2.30

``` r
describe_posterior(as.data.frame(beta_draws))
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |         95% CI |   pd |          ROPE | % in ROPE
    ## ----------------------------------------------------------------------
    ## V1        |  -2.20 | [-2.39, -2.02] | 100% | [-0.10, 0.10] |        0%
    ## V2        |   0.12 | [ 0.08,  0.16] | 100% | [-0.10, 0.10] |    11.18%
    ## V3        |   1.41 | [ 1.32,  1.50] | 100% | [-0.10, 0.10] |        0%
    ## V4        |   2.45 | [ 2.25,  2.64] | 100% | [-0.10, 0.10] |        0%

**Posterior with multivariate normal prior (D&C)**

``` r
#######################
#### Simulate Data ####
#######################

set.seed(666)
N <- 1e4
x1 <- rnorm(N, sd = 3)           # some continuous variables 
x2 <- rnorm(N, sd = 10)
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
  scaling <- var(fold_data_X[,-1])
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
## Why create folds this way? Why not force a balance?
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
```

    ## [[1]]
    ## [1] "MASS"      "stats"     "graphics"  "grDevices" "utils"     "datasets" 
    ## [7] "methods"   "base"     
    ## 
    ## [[2]]
    ## [1] "MASS"      "stats"     "graphics"  "grDevices" "utils"     "datasets" 
    ## [7] "methods"   "base"     
    ## 
    ## [[3]]
    ## [1] "MASS"      "stats"     "graphics"  "grDevices" "utils"     "datasets" 
    ## [7] "methods"   "base"     
    ## 
    ## [[4]]
    ## [1] "MASS"      "stats"     "graphics"  "grDevices" "utils"     "datasets" 
    ## [7] "methods"   "base"

``` r
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

full_data_draws <- full_data_draws %>% as.mcmc()
heidel.diag(full_data_draws)
```

    ##                                    
    ##      Stationarity start     p-value
    ##      test         iteration        
    ## var1 passed       1         0.973  
    ## var2 passed       1         0.988  
    ## var3 passed       1         0.987  
    ## var4 passed       1         1.000  
    ##                               
    ##      Halfwidth Mean  Halfwidth
    ##      test                     
    ## var1 passed    -3.00 0.001171 
    ## var2 passed     3.89 0.001357 
    ## var3 passed     1.14 0.000403 
    ## var4 passed     2.22 0.000933

``` r
raftery.diag(full_data_draws)
```

    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                        
    ##  Burn-in  Total Lower bound  Dependence
    ##  (M)      (N)   (Nmin)       factor (I)
    ##  2        3710  3746         0.99      
    ##  1        3755  3746         1.00      
    ##  2        3771  3746         1.01      
    ##  2        3810  3746         1.02

``` r
# Apply individually 
check_heidel <- function(x) {
  
  x %>% as.mcmc() %>% heidel.diag()
  
}

check_raftery <- function(x) {
  
  x %>% as.mcmc() %>% raftery.diag()
  
}

lapply(results_recentered, check_heidel)
```

    ## [[1]]
    ##                                    
    ##      Stationarity start     p-value
    ##      test         iteration        
    ## var1 passed       1         0.516  
    ## var2 passed       1         0.588  
    ## var3 passed       1         0.590  
    ## var4 passed       1         0.827  
    ##                               
    ##      Halfwidth Mean  Halfwidth
    ##      test                     
    ## var1 passed    -3.00 0.002248 
    ## var2 passed     3.89 0.002457 
    ## var3 passed     1.14 0.000715 
    ## var4 passed     2.22 0.001861 
    ## 
    ## [[2]]
    ##                                    
    ##      Stationarity start     p-value
    ##      test         iteration        
    ## var1 passed       1         0.521  
    ## var2 passed       1         0.586  
    ## var3 passed       1         0.595  
    ## var4 passed       1         0.798  
    ##                               
    ##      Halfwidth Mean  Halfwidth
    ##      test                     
    ## var1 passed    -3.00 0.002615 
    ## var2 passed     3.89 0.003088 
    ## var3 passed     1.14 0.000934 
    ## var4 passed     2.22 0.002203 
    ## 
    ## [[3]]
    ##                                    
    ##      Stationarity start     p-value
    ##      test         iteration        
    ## var1 passed       1         0.520  
    ## var2 passed       1         0.586  
    ## var3 passed       1         0.592  
    ## var4 passed       1         0.826  
    ##                               
    ##      Halfwidth Mean  Halfwidth
    ##      test                     
    ## var1 passed    -3.00 0.002538 
    ## var2 passed     3.89 0.002970 
    ## var3 passed     1.14 0.000865 
    ## var4 passed     2.22 0.002080 
    ## 
    ## [[4]]
    ##                                    
    ##      Stationarity start     p-value
    ##      test         iteration        
    ## var1 passed       1         0.517  
    ## var2 passed       1         0.585  
    ## var3 passed       1         0.589  
    ## var4 passed       1         0.805  
    ##                               
    ##      Halfwidth Mean  Halfwidth
    ##      test                     
    ## var1 passed    -3.00 0.002441 
    ## var2 passed     3.89 0.002750 
    ## var3 passed     1.14 0.000802 
    ## var4 passed     2.22 0.002062

``` r
lapply(results_recentered, check_raftery)
```

    ## [[1]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                        
    ##  Burn-in  Total Lower bound  Dependence
    ##  (M)      (N)   (Nmin)       factor (I)
    ##  2        3710  3746         0.99      
    ##  2        3771  3746         1.01      
    ##  2        3710  3746         0.99      
    ##  2        3771  3746         1.01      
    ## 
    ## 
    ## [[2]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                        
    ##  Burn-in  Total Lower bound  Dependence
    ##  (M)      (N)   (Nmin)       factor (I)
    ##  2        3741  3746         0.999     
    ##  2        3710  3746         0.990     
    ##  2        3680  3746         0.982     
    ##  2        3771  3746         1.010     
    ## 
    ## 
    ## [[3]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                        
    ##  Burn-in  Total Lower bound  Dependence
    ##  (M)      (N)   (Nmin)       factor (I)
    ##  2        3680  3746         0.982     
    ##  2        3710  3746         0.990     
    ##  2        3710  3746         0.990     
    ##  2        3802  3746         1.010     
    ## 
    ## 
    ## [[4]]
    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                        
    ##  Burn-in  Total Lower bound  Dependence
    ##  (M)      (N)   (Nmin)       factor (I)
    ##  2        3741  3746         0.999     
    ##  2        3710  3746         0.990     
    ##  2        3710  3746         0.990     
    ##  2        3741  3746         0.999

``` r
#################
#### Summary ####
#################
describe_posterior(as.data.frame(full_data_draws))
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |         95% CI |   pd |          ROPE | % in ROPE
    ## ----------------------------------------------------------------------
    ## V1        |  -3.00 | [-3.25, -2.76] | 100% | [-0.10, 0.10] |        0%
    ## V2        |   3.89 | [ 3.61,  4.18] | 100% | [-0.10, 0.10] |        0%
    ## V3        |   1.14 | [ 1.06,  1.23] | 100% | [-0.10, 0.10] |        0%
    ## V4        |   2.22 | [ 2.01,  2.43] | 100% | [-0.10, 0.10] |        0%

``` r
## Truth
beta
```

    ## [1] -3.0  3.8  1.1  2.3

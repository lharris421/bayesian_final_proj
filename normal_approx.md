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
\pi({\boldsymbol{\theta}}\|\mathrm{\bf{y}})\_i \propto \Big\[\Big(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{ y_i}\Big(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{1- y_i} \Big\] \cdot \Big\[\exp\Big\\{ -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\Big\\}  \Big\]
$$

*Log-posterior (non-D&C)*

$$
\begin{aligned}
\ell_i(\pi({\boldsymbol{\theta}}\|\mathrm{\bf{y}}) ) &\propto y_i \log(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}})+(1-y_i) \log(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}) -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
&\propto  y_i\log(\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})- y_i\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})+ y_i \log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
&\propto y_i\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}} -\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
\end{aligned}
$$

**Solve for MLE**

``` r
## Simulate Data
set.seed(666)
N <- 1e4
x1 <- rnorm(N)           # some continuous variables 
x2 <- rnorm(N)
x3 <-rnorm(N)
beta <- c(1, 0.11,1.34,2.3)
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

    ## [1] 1.00181368 0.08753796 1.38046846 2.26045203

**Posterior with multivariate normal prior (non-D&C)**

``` r
mu <- rep(0,4)
sigma <- diag(c(40^2, 3^2 *sqrt(diag(var(X[,-1])))))

log_post_fun <- function(param){
  sum(y*X%*%param - log(1 + exp(X%*%param)))-
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

Opt <- optim(par = rep(0,4),
             fn = log_post_fun, 
             method = "BFGS",
             control = list(fnscale= -1,
                            maxit = 1e6),
             hessian = T)
(params <- Opt$par)
```

    ## [1] 1.00192414 0.08760457 1.37917357 2.25989564

``` r
SigNew <-  chol2inv(chol(-Opt$hessian))

## Draw from multivariate normal
beta_draws <- MASS::mvrnorm(1e4, mu = params, Sigma = SigNew)
```

**Posterior with multivariate normal prior (D&C)**

``` r
#######################
#### Simulate Data ####
#######################

set.seed(666)
N <- 1e4
x1 <- rnorm(N)           # some continuous variables 
x2 <- rnorm(N)
x3 <-rnorm(N)
beta <- c(1,3.8,1.1,2.3)
X <-  cbind(1, x1,x2,x3)
eta <- X %*% beta      
pr <- 1/(1+exp(-eta))         # pass through an inv-logit function
y <- rbinom(N,1,pr)  

#######################
#### Log Posterior ####
#######################

log_post_fun <- function(param, X, y, N, m, sigma, mu){
  (N/m)*sum(y*X%*%param - log(1 + exp(X%*%param)))-
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

########################
#### Optim function ####
########################

optim_fun <- function(init, X, y, N, m, sigma, mu){
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

inner_draws <- function(i, X, y, N, NN = 1e4){
  fold_data_y <- y[fold_idx == i]
  fold_data_X <- X[fold_idx == i, ]
  
  #priors
  sigma <- diag(c(40^2, 3^2 *sqrt(diag(var(fold_data_X[,-1])))))
  mu <- rep(0,4)
  
  Opt <- optim_fun(init = rep(0, ncol(fold_data_X)),
                   X = fold_data_X, 
                   y = fold_data_y,
                   m = sum(fold_idx == i),
                   N = N,
                   sigma = sigma,
                   mu = mu)
  params <- Opt$par
  SigNew <-  chol2inv(chol(-Opt$hessian))
  
  ## Draw from multivariate normal
 beta_draws <- MASS::mvrnorm(NN, mu =  params, Sigma = SigNew)
 
 return(beta_draws)
}

########################
#### Partition Data ####
########################
# Create folds
K <- 4
fold_idx <- sample(1:K, size = N, replace = TRUE)

cl <- makeCluster(min(detectCores(), K))
clusterExport(cl,
              c("X", "y", "fold_idx", "N",
                "inner_draws", "log_post_fun",
                "optim_fun"))
clusterEvalQ(cl,{
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
results <-  parLapply(cl,
                      1:K,
                      function(x){inner_draws(i = x,
                                              X = X, 
                                              y = y,
                                              N = N,
                                              NN = 1e4)})
stopCluster(cl)

########################
#### Recenter Draws ####
########################

subset_mean <- t(sapply(1:K,function(i) colMeans(results[[i]]))) # rows indicate subset
global_mean <- colMeans(subset_mean)
recenter <- t(sapply(1:K,function(i) subset_mean[i, ] - global_mean)) # rows indicate subset
results_recentered <- lapply(1:K,function(i) results[[i]] + matrix(recenter[i, ],
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
```

    ##                                    
    ##      Stationarity start     p-value
    ##      test         iteration        
    ## var1 passed       16001     0.05849
    ## var2 passed           1     0.08317
    ## var3 failed          NA     0.00458
    ## var4 passed           1     0.10207
    ##                               
    ##      Halfwidth Mean  Halfwidth
    ##      test                     
    ## var1 passed    0.911 0.0664   
    ## var2 failed    3.774 0.6269   
    ## var3 <NA>         NA     NA   
    ## var4 passed    2.238 0.1807

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
    ##  24       32672 3746         8.72      
    ##  27       33462 3746         8.93      
    ##  12       17284 3746         4.61      
    ##  12       18604 3746         4.97

``` r
#################
#### Summary ####
#################
describe_posterior(as.data.frame(full_data_draws))
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |       95% CI |   pd |          ROPE | % in ROPE
    ## --------------------------------------------------------------------
    ## V1        |   0.94 | [0.71, 1.28] | 100% | [-0.10, 0.10] |        0%
    ## V2        |   3.72 | [3.04, 4.64] | 100% | [-0.10, 0.10] |        0%
    ## V3        |   1.07 | [0.87, 1.53] | 100% | [-0.10, 0.10] |        0%
    ## V4        |   2.27 | [1.92, 2.51] | 100% | [-0.10, 0.10] |        0%
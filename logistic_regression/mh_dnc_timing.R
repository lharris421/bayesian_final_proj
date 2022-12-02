## WASP Functions
computeBarycenter <- function (meanList, covList) {
  library(matrixStats)
  library(expm)
  library(MASS)
  
  ncomp <- length(covList)
  ndim <- nrow(covList[[1]])
  wts <- rep(1, ncomp) / ncomp
  
  baryMean <- rowMeans(do.call(cbind, meanList))
  baryCov <- diag(1.0, ndim)
  barySd <- sqrtm(baryCov)
  err <- baryCov
  cnt <- 1
  while ((norm(err, type = "F") > 1e-6) & cnt < 500) {
    if (cnt %% 10 == 0)  cat("iter: ", cnt, "\n")
    
    ssj <- matrix(0.0, nrow = ndim, ncol = ndim)
    for (ii in 1:ncomp) {
      ssj <- ssj + sqrtm(baryCov %*% covList[[ii]])
    }
    
    tmp <- solve(sqrtm(baryCov))
    baryCovNew <- tmp %*% tcrossprod(ssj / ncomp) %*% tmp
    err <- baryCov - baryCovNew
    baryCov <- baryCovNew
    barySd <- sqrtm(baryCov)
    cnt <- cnt + 1
  }
  
  list(mean = baryMean, cov = baryCov, sqrt = barySd, iter = cnt)
}

sampleBetas <- function (betaList) {
  
  library(matrixStats)
  library(expm)
  
  npart <- length(betaList)
  meanBetas <- list()
  covBetas <- list()
  
  meanBetas <- lapply(betaList, function(x) colMeans(x))
  covBetas <- lapply(betaList, function(x) cov(x))
  
  resBary <- computeBarycenter(meanBetas, covBetas)
  muBetas <- resBary$mean
  sigBetas <- resBary$cov
  sqrtSigBetas <- resBary$sqrt
  
  baryList <- list()
  for (ii in 1:npart) {
    tmp <- as(chol2inv(chol(covBetas[[ii]])), "symmetricMatrix")
    tmp1 <- matrix(meanBetas[[ii]], nrow = nrow(betaList[[ii]]), ncol = ncol(betaList[[ii]]), byrow = TRUE)
    centScaledSamps <- sqrtm(tmp) %*% (t(betaList[[ii]] - tmp1))
    baryList[[ii]] <- t(muBetas + sqrtSigBetas %*% centScaledSamps)
  }
  
  list(
    wasp = do.call(rbind, baryList)
  )
}

seeds <- 1001:1010
time_seeds <- numeric(length(seeds))
wasp_diff <- numeric(length(seeds))
for (seed in seeds) {
  
  times <- numeric(25)
  for (j in 1:25) {
    load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/time{j}_{seed}.rds"))
    times[j] <- end_time
  }
  
  date_string <- read_lines(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/mh_large/start_time_{seed}.txt"))
  split_date <- str_split(date_string, " ")[[1]]
  month_num <- which(month.abb == split_date[2])
  ## as_datetime(paste0(split_date[6], "-", month_num, "-", split_date[3], " ", split_date[4]), tz = "America/Chicago")
  
  ## Need to be careful with the units here
  time_seeds[seed - 1000] <- difftime(with_tz(as_datetime(max(times)), "America/Chicago"), as_datetime(paste0(split_date[6], "-", month_num, "-", split_date[3], " ", split_date[4]), tz = "America/Chicago"), "mins")
  
  start_time <- Sys.time()
  results <- list()
  for (j in 1:25) {
    load(glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/res{j}_{seed}.rds"))
    results[[j]] <- res
  }
  
  remove_burnin <- function(x, burnin) {
    
    x[-(1:burnin),]
    
  }
  
  results <- lapply(results, remove_burnin, 1000)
  ## TIME END HERE
  end_time <- Sys.time()
  tdiff <- difftime(end_time, start_time, units = "mins")
  time_seeds[seed - 1000] <- time_seeds[seed - 1000] + tdiff
  
  start_time <- Sys.time()
  K <- 25
  subset_mean <- t(sapply(1:K, function(i) colMeans(results[[i]]))) # rows indicate subset
  global_mean <- colMeans(bind_rows(lapply(results, data.frame)))
  recenter <- t(sapply(1:K, function(i) subset_mean[i, ] - global_mean)) # rows indicate subset
  results_recentered <- lapply(1:K,function(i) results[[i]] - matrix(recenter[i, ],
                                                                     nrow = nrow(results[[i]]),
                                                                     ncol = ncol(results[[i]]),
                                                                     byrow = T))
  full_data_draws <- do.call(rbind, results_recentered)
  
  end_time <- Sys.time()
  tdiff <- difftime(end_time, start_time, units = "mins")
  time_seeds[seed - 1000] <- time_seeds[seed - 1000] + tdiff
  
  start_time <- Sys.time()
  bary <- sampleBetas(results)
  bary_res <- bind_rows(bary)
  end_time <- Sys.time()
  wasp_diff[seed - 1000] <- (difftime(end_time, start_time, units = "mins") - tdiff)
  
  ## Save the draws
  save(full_data_drawms, file = glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/simple_{seed}.rds"))
  save(bary_res, file = glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/results/wasp_{seed}.rds"))

}

time_simple <- time_seed
time_wasp <- time_seed + wasp_diff

save(time_simply, file = glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/mh_large/simple_{seed}.rds"))
save(time_wasp, file = glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/times/results/wasp_{seed}.rds"))



library(rstan)
library(glue)
library(dplyr)


## clean up folder to store partitions in
system("rm /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/*", intern = FALSE)
system("rm /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/*", intern = FALSE)

### Will want to partition and submit here
### Also include combining data in this file

## Number of partitions
npart <- 250

## Submit the partitioning script
sub_command <- glue(
  "#!/bin/bash 
  qsub -pe smp -2 -cwd -e /dev/null -o /dev/null /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partition.job"
)

system(sub_command, intern = FALSE)

continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/")) < npart
i <- 0

while (continue & i < 100) {
  
  Sys.sleep(10)
  i <- 1 + 1
  continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/")) < npary
  
}

## submit command
sub_command <- glue(
  "#!/bin/bash 
  qsub -pe smp -2 -cwd -e /dev/null -o /dev/null -t 1-{npart} /Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/slr.job"
)

system(sub_command, intern = FALSE)

## Check every minute if all the files are there before combining
continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/")) < npart
i <- 0

while (continue & i < 100) {
  
  Sys.sleep(10)
  i <- 1 + 1
  continue <- length(list.files("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/results/")) < npary
  
}


#######




###### Combine everything

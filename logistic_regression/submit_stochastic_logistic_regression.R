library(rstan)
library(glue)

### Will want to partition and submit here
### Also include combining data in this file

## May need to likewise submit this

### Partition and save the data temporarily
npart <- 20

load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/full_data_2.rds")
full_data <- mod_dat
rm("mod_dat")

partitions <- (sample(1:nrow(full_data), replace = FALSE) %% npart) + 1

for (i in 1:npart) {
  
  y0 <- as.numeric(unlist(full_data[partitions == i,"proc"]))
  x0 <- as.matrix(full_data)[partitions == i, -"proc"]
  nrep0 <- length(y)
  name0 <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/partition{i}.rds")
  save(
    y0,
    x0,
    nrep0,
    file = name
  )
  
}

## submit command
sub_command <- glue(
  "#!/bin/bash 
  qsub -pe smp -2 -cwd -e /dev/null -o /dev/null -t 1-{npart} mh.job"
)

system(sub_command, intern = FALSE)

## Check every minute if all the files are there before combining


#######




###### Combine everything

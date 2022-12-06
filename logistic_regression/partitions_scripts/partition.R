### 
library(glue)
library(dplyr)

### Partition and save the data temporarily
load("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/full_data_2.rds")

## head(mod_dat)
## table(mod_dat$carrier, mod_dat$proc)

score <- mod_dat %>% select_at(vars(CHF:Depression)) %>% rowSums()
full_data <- mod_dat %>%
  dplyr::select(proc, smoker, age, sex, carrier) %>%
  mutate(score = scale(score))

nrow(full_data)
table(full_data$carrier)
table(full_data$proc)
table(full_data$sex)
table(full_data$smoker)
mean(full_data$age)

# carriers <- full_data %>%
#   filter(carrier == TRUE)
# 
# controls <- full_data %>%
#   filter(carrier == FALSE)
# 
# set.seed(1234)
# controls <- controls[sample(1:10, nrow(controls), replace = TRUE) == 1,]
# 
# full_data <- bind_rows(carriers, controls)

## check summary(glm(proc ~ ., full_data, family = "binomial"))

set.seed(1234)
npart <- 25
partitions <- (sample(1:nrow(full_data), replace = FALSE) %% npart) + 1

for (i in 1:npart) {
  
  y0 <- as.numeric(unlist(full_data[partitions == i,"proc"]))
  x0 <- as.matrix(full_data)[partitions == i, -1]
  nrep0 <- nrow(full_data) / nrow(x0)
  name0 <- glue("/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/logistic_regression/partitions/large/partition{i}.rds")
  save(
    y0,
    x0,
    nrep0,
    file = name0
  )
  
}
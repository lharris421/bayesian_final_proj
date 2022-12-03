library(tidyverse)
library(magrittr)
library(glue)

###################
#### Load Data ####
###################

main_path <- "/Shared/Statepi_Marketscan/aa_lh_bayes/bayesian_final_proj/data/"
results <- tibble()

## Normal approx ---------------------------------------------------------------

## non d&c
load(paste0(main_path, "normal_approx_full.Rdata"))
results <- bind_rows(results, 
                     out$comp_time %>% 
                       mutate(type =  "Non-D&C",
                              method = "Normal Approximation"))
rm(out, beta_draws)

## D&C recentering
load(paste0(main_path, "normal_approx_dnc.Rdata"))
results <- bind_rows(results, 
                     out$comp_time %>% 
                       mutate(type =  "D&C Recentering",
                              method = "Normal Approximation"))
rm(out, full_data_draws)

## D&C WASP
load(paste0(main_path, "normal_approx_dnc_WASP.Rdata"))
results <- bind_rows(results, 
                     out$comp_time %>% 
                       mutate(type =  "D&C WASP",
                              method = "Normal Approximation"))
rm(out, full_data_draws)

## Adaptive MH -----------------------------------------------------------------

## non d&c
tmp <- list.files(main_path)[str_detect(list.files(main_path), "adaptive_MH_full_")]
time_tmp <- c()

for (i in tmp){
  load(paste0(main_path, i))
  time_tmp <- c(time_tmp, out$comp_time)
  rm(out)
}

results <- bind_rows(results, 
                     tibble(min = min(time_tmp),
                            lq = quantile(time_tmp, probs = 0.25),
                            mean = mean(time_tmp),
                            median = median(time_tmp),
                            uq = quantile(time_tmp, probs = 0.75),
                            max = max(time_tmp),
                            type =  "Non-D&C",
                            method = "Adaptive - MH"))

## D&C recentering
load(paste0(main_path, "adaptive_MH_dnc.Rdata"))
results <- bind_rows(results, 
                     out$comp_time %>% 
                       mutate(type =  "D&C Recentering",
                              method = "Adaptive - MH"))
rm(out, full_data_draws)

## D&C WASP
load(paste0(main_path, "adaptive_MH_dnc_WASP.Rdata"))
results <- bind_rows(results, 
                     out$comp_time %>% 
                       mutate(type =  "D&C WASP",
                              method = "Adaptive - MH"))
rm(out, full_data_draws)


## MH --------------------------------------------------------------------------
time_data <- data.frame(
  method = c("MH", "MH", "MH", "Adaptive - MH", "Adaptive - MH", "Adaptive - MH", "Normal", "Normal", "Normal", "GLM"),
  group = c("Non-D&C", "D&C Recentering", "D&C WASP", "Non-D&C", " D&C Recentering", "D&C WASP", "Non-D&C", "D&C Recentering", "D&C WASP", "Non-D&C"),
  median = c(391.20, 14.67, 14.67, 203.13, 12.06, 13.90, 8.12, 0.71, 0.70, 0.72), 
  q1 = c(387.90, 13.72, 13.73, 199.74, 11.03, 11.98, 7.06, 0.65, 0.64, 0.67), 
  q3 = c(398.0, 15.52, 15.52, 210.73, 14.10, 16.10, 10.00, 0.98, 0.97, 1.00)
) 

## non d&c
tmp <- as_tibble(time_data) %>% 
  filter(group == "Non-D&C" & method == "MH") %>% 
  mutate(mean = NA,
         min = NA,
         max = NA) %>% 
  rename(lq=q1,
         uq=q3,
         type = group)

results <- bind_rows(results, tmp)

## D&C recentering
tmp <- as_tibble(time_data) %>% 
  filter(group == "D&C Recentering" & method == "MH") %>% 
  mutate(mean = NA,
         min = NA,
         max = NA) %>% 
  rename(lq=q1,
         uq=q3,
         type = group)

results <- bind_rows(results, tmp)

## D&C WASP
tmp <- as_tibble(time_data) %>% 
  filter(group == "D&C WASP" & method == "MH") %>% 
  mutate(mean = NA,
         min = NA,
         max = NA) %>% 
  rename(lq=q1,
         uq=q3,
         type = group)

results <- bind_rows(results, tmp)

## GLM --------------------------------------------------------------------------
load(paste0(main_path, "glm.Rdata"))
results <- bind_rows(results, 
                     out$comp_time %>% 
                       mutate(type =  "Non-D&C",
                              method = "GLM"))
rm(out)

##############
#### Save ####
##############

save(results, 
     file = paste0(main_path,  "time_plot_data.Rdata"))

##############
#### Plot ####
##############

load(paste0(main_path, "time_plot_data.Rdata"))

results <- results %>% 
  mutate(type = factor(type, levels = c("Non-D&C", "D&C Recentering", "D&C WASP")),
         method = factor(method, levels = c("MH", "Adaptive - MH", "Normal Approximation", "GLM"))) %>% 
  rename(Method = method)

ggplot(results, aes(x=type, y=median, fill=Method)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c("red",
                             "darkorange",
                             "forestgreen",
                             "darkblue")) +
  geom_errorbar(aes(ymin=lq, ymax=uq), width=0.2, size=1,
                position=position_dodge(.9)) +
  geom_text(data=results, aes(x=type, y=25, label=format(round(median, 2), nsmall = 2)),
            position=position_dodge(0.9),
            fontface = "bold", size = 6)+
  labs(title="Results - Comparison of Computation Time", x="", y = "Computation Time (Minutes)") +
  theme_bw(base_size = 22) +
  theme(legend.position="bottom",
        legend.title = element_text(face = "bold"))


ggplot(results, aes(x=type, y=median, fill=Method)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c("red",
                             "darkorange",
                             "forestgreen",
                             "darkblue")) +
  geom_errorbar(aes(ymin=lq, ymax=uq), width=0.2, size=1,
                position=position_dodge(0.9)) +
  geom_text(data=results, aes(x=type, y=1.5, label=format(round(median, 2), nsmall = 2)),
            position=position_dodge(0.9),
            fontface = "bold", size = 6)+
  labs(title="Results - Comparison of Computation Time", x="", y = "Computation Time (Minutes)") +
  theme_bw(base_size = 22) +
  coord_cartesian(ylim=c(0, 20)) +
  theme(legend.position="bottom",
        legend.title = element_text(face = "bold"))

ggplot(results, aes(x=type, y=median, fill=Method)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values=c("red",
                             "darkorange",
                             "forestgreen",
                             "darkblue")) +
  geom_errorbar(aes(ymin=lq, ymax=uq), width=0.2, size=1,
                position=position_dodge(0.9)) +
  geom_text(data=results, aes(x=type, y=0.8, label=format(round(median, 2), nsmall = 2)),
            position=position_dodge(0.9),
            fontface = "bold", size = 6)+
  labs(title="Results - Comparison of Computation Time", x="", y = "Computation Time (Minutes)") +
  theme_bw(base_size = 22) +
  coord_cartesian(ylim=c(0, 2)) +
  theme(legend.position="bottom",
        legend.title = element_text(face = "bold"))



library(ggplot2)
time_data <- data.frame(
  method = c("MH", "MH", "MH", "AMH", "AMH", "AMH", "Normal", "Normal", "Normal", "GLM"),
  group = c("Non-D&C", "Simple", "Wasp", "Non-D&C", "Simple", "Wasp", "Non-D&C", "Simple", "Wasp", "Non-D&C"),
  median = c(391.20, 14.67, 14.67, 203.13, 12.06, 13.90, 8.12, 0.71, 0.70, 0.72), 
  q1 = c(387.90, 13.72, 13.73, 199.74, 11.03, 11.98, 7.06, 0.65, 0.64, 0.67), 
  q3 = c(398.0, 15.52, 15.52, 210.73, 14.10, 16.10, 10.00, 0.98, 0.97, 1.00)
) %>%
  mutate(method = factor(method, levels = c("MH", "AMH", "Normal", "GLM"), ordered = TRUE))

ggplot(time_data, aes(x=group, y=median, fill=method)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=q1, ymax=q3), width=.2,
                position=position_dodge(.9)) +
  labs(title="Results - Comparison of Computation Time", x="", y = "Computation Time (Minutes)") +
  theme_bw() 

ggplot(time_data, aes(x=group, y=median, fill=method)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=q1, ymax=q3), width=.2,
                position=position_dodge(.9)) +
  labs(title="Results - Comparison of Computation Time", x="", y = "Computation Time (Minutes)") +
  theme_bw() +
  coord_cartesian(ylim=c(0, 20))

result_table <- read.csv("SuppFigure4/bmi1-bmi8-bmia.csv")

library(ggplot2)

desired_order <- c(
  "1-year-old BMI -> Adult BMI (total)",
  "1-year-old BMI -> Adult BMI (indirect)",
  "1-year-old BMI -> Adult BMI (direct)",
  "1-year-old BMI -> 8-year-old BMI (direct)",
  "8-year-old BMI -> Adult BMI (direct)"
)

df$quantity <- factor(df$quantity, levels = desired_order)

ggplot(df, aes(x = beta_hat, y = quantity, 
               xmin = lower, xmax = upper)) + 
  geom_point() +  # Mean estimates as points
  geom_errorbarh(height = 0.1) +  # Horizontal error bars
  theme_bw() +
  labs(x = expression(hat(beta)), y = NULL, 
       title = "Estimates from FLOW-MR") +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 0, linetype = "dashed")

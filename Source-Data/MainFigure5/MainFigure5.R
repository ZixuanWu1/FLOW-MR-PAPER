library(ggforce)
library(ggplot2)

df = read.csv("MainFigure5/stroke_result.csv")




df4 = df[df$outcome == "Stroke", ]
df3 = df[df$outcome == "SBP (A)", ]
df2 = df[df$outcome == "BMI (A)", ]
df1 = df[df$outcome == "LDL-C (A)", ]

# Plot
p1 = ggplot(df1, aes(y = exposure, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0, linetype = "dashed") + xlab("Estimate") + ylab("Exposure -> Outcome")

p2 = ggplot(df2, aes(y = exposure, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0, linetype = "dashed") + xlab("Estimate") + ylab("Exposure -> Outcome")


p3 = ggplot(df3, aes(y = exposure, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0,linetype = "dashed") + xlab("Estimate") + ylab("Exposure -> Outcome")

p4 = ggplot(df4, aes(y = exposure, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0,linetype = "dashed") + xlab("Estimate") + ylab("Exposure -> Outcome") 


wrap_plots(p1+ labs(title = "Direct Effects on LDL-C (A)"), p2+ labs(title = "Direct Effects on BMI (A)"), p3+ labs(title = "Direct Effects on SBP (A)"),  p4+ labs(title = "Direct Effects on Stroke"), nrow = 2) +
  plot_layout(guides = "collect")  & 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),    axis.text.x = element_text(size = 14),   # Increase X-axis tick text size
        axis.text.y = element_text(size = 14)  , plot.title = element_text(size = 18,  hjust = 0.5))



library(ggforce)
library(ggplot2)

df = read.csv("MainFigure5/stroke_result.csv")




df1 = df[df$outcome == "Stroke", ]
df2 = df[df$outcome == "SBP (A)", ]
df3 = df[df$outcome == "BMI (A)", ]

# Plot
p1 = ggplot(df1, aes(y = exposure_outcome, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0, linetype = "dashed") + xlab("Estimate") + ylab("Exposure -> Outcome")

p2 = ggplot(df2, aes(y = exposure_outcome, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0, linetype = "dashed") + xlab("Estimate") + ylab("Exposure -> Outcome")


p3 = ggplot(df3, aes(y = exposure_outcome, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0,linetype = "dashed") + xlab("Estimate") + ylab("Exposure -> Outcome")

wrap_plots(p1+ labs(title = "Direct Effects on Stroke"), p2+ labs(title = "Direct Effects on SBP (A)"), p3+ labs(title = "Direct Effects on BMI (A)")) +
  plot_layout(guides = "collect")  & 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())



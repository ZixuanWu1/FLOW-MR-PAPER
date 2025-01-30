library(ggforce)
library(ggplot2)

df1 = read.csv("MainFigure4/breast_k=4_direct.csv")
df2 = read.csv("MainFigure4/breast_k=4_indirect.csv")


df1$exposure_outcome = factor(df1$exposure_outcome, levels =  c("1-year-old BMI -> 8-year-old BMI", 
                                                              "1-year-old BMI -> Adult BMI", "8-year-old BMI -> Adult BMI",  "1-year-old BMI -> Breast Cancer",  
                                                              "8-year-old BMI -> Breast Cancer", "Adult BMI -> Breast Cancer" )   )

p1 = ggplot(df1, aes(y = exposure_outcome, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0) + xlab("Estimate") + ylab("Exposure -> Outcome")

p2 = ggplot(df2, aes(y = exposure_outcome, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(aes(x = mean), size = 2) +  # Add points for mean
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) + geom_vline(xintercept = 0) + xlab("Estimate") + ylab("Exposure -> Outcome")


wrap_plots(p1 + labs(title = "Direct Effects"), p2 + labs(title = "Indirect Effects")) +plot_annotation(
  theme = theme(plot.title = element_text(hjust = 0.5))
) +
  plot_layout(guides = "collect") & 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())

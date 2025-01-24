library(patchwork)
library(dplyr)
library(ggpubr)
library(ggplot2)
# Summarize information from data to data frame



df1 = read.csv("SuppFigure2.1/k=3_time.csv")
df2 = read.csv("SuppFigure2.1/multivariate_time.csv")

df1$pval = factor((df1$pval), levels = c(1e-8, 1e-6, 1e-4, 1e-2))
df2$pval = factor((df2$pval), levels = c(1e-8, 1e-6, 1e-4, 1e-2))


c1 = ggplot(df1 , aes(x = nsnp, y = coverage, color = Method)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(16)) +
  geom_point(position = position_dodge(16)) +
  labs(x = "Number of SNPs", y = "Runtime (log10 seconds)") +
  theme_bw() + ggtitle("K = 3") + theme(legend.position = "bottom") 

c2 = ggplot(df2 , aes(x = nsnp, y = coverage, color = Method)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(16)) +
  geom_point(position = position_dodge(16)) +
  labs(x = "Number of SNPs", y = "Runtime (log10 seconds)") +
  theme_bw() + ggtitle("Multivariate") +    theme(legend.position = "none") + 
  guides(color = "none", fill = "none", linetype = "none") 


wrap_plots(c1, c2, guides = "collect") +
  plot_layout(guides = "collect",axis_titles = "collect") &    theme(legend.position = "bottom")# Ensures only one legend is shown

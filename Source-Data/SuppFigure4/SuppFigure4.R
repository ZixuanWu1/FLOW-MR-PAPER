library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)

df1 = read.csv("k=3_bodysize_real.csv")
df2 = read.csv("k=3_bmi_real.csv")
 
df1$pval = factor((df1$pval), levels = c(1e-8, 1e-6, 1e-4, 1e-2))
df1$Trait = factor((df1$Trait), levels = c(1, 2))

df2$pval = factor((df2$pval), levels = c(1e-8, 1e-6, 1e-4, 1e-2))
df2$Trait = factor((df2$Trait), levels = c(1, 2))




df1 = df1[df1$Method != "FLOW-MR-indirect", ]
ggplot(df1, aes(x = pval, y = estimate, color = Trait)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `lower`, ymax = `upper`), width = 0.2, position = position_dodge(width = 0.3)) +
  labs(x = "p-value Threshold", y = expression(hat(beta))
  ) +
  theme_bw() +
  facet_wrap(~ Method, nrow = 2, scales = "free_x") + 
  scale_color_manual(values = c("1" = "#F8766D", "2" = "#00BFC4"),
                     labels = c("1" = "early-life Body Size (UK Biobank)", "2" = "Adult BMI")) + 
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +  theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



df2 = df2[df2$Method !="FLOW-MR-indirect", ]
ggplot(df2, aes(x = pval, y = estimate, color = Trait)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `lower`, ymax = `upper`), width = 0.2, position = position_dodge(width = 0.3)) +
  labs(x = "p-value Threshold", y = expression(hat(beta))
  ) +
  theme_bw() +
  facet_wrap(~ Method, nrow = 2, scales = "free_x") + 
  ylim(c(-1, 1)) +   scale_color_manual(values = c("1" = "#F8766D", "2" = "#00BFC4"),
                                        labels = c("1" = "8-year-old BMI", "2" = "Adult BMI"))+ 
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +  theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 





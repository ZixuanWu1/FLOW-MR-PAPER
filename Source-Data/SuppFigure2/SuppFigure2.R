library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)

df = read.csv("SuppFigure2/k=4_sim.csv")


df$pval = factor(unique(df$pval), levels = c(1e-8, 1e-6, 1e-4, 1e-2))
K = 3
n = 100
# Compute the coverage proportion
true_val = c(0.2,0.2, 1)
covers = list()
for (trait in 1:K){
  df_sub = df[df$Trait == trait, ]
  df_sub = df_sub[df_sub$Method != "FLOW-MR-indirect", ]
  df_sub$result <- with(df_sub, ifelse(`lower` <= true_val[trait] & `upper` >= true_val[trait], TRUE, FALSE))
  
  coverage_df_sub <- df_sub %>%
    group_by(Method, Trait, pval) %>%
    summarize(coverage = mean(result)) %>%
    ungroup()
  coverage_df_sub$lower = 0
  coverage_df_sub$upper = 0
  
  for(i in 1:nrow(coverage_df_sub)){
    coverage_df_sub$lower[i] =  as.vector( 
      prop.test(x = coverage_df_sub$coverage[i] * n, n = n)$conf.int)[1]
    coverage_df_sub$upper[i] =  as.vector( 
      prop.test(x = coverage_df_sub$coverage[i] * n, n = n)$conf.int)[2]
  }
  
  covers[[trait]] = ggplot(coverage_df_sub, aes(x = pval, y = coverage, color = Method)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.5)) +
    geom_point(position = position_dodge(0.5)) +
    labs(x = "P-value threshold", y = "Coverage", title = paste("Coverage for Exp", toString(trait))) +
    theme_bw() + geom_abline(slope = 0, intercept = .95, linetype = "dashed")+   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}





# Compute the Length
lengths = list()
for (trait in 1:K){
  df_sub = df[df$Trait == trait, ]
  df_sub = df_sub[df_sub$Method != "Bayes-indirect", ]
  df_sub = df_sub[df_sub$length < 3, ]
  
  
  y_min <- quantile(df_sub$length, 0.01)  # Lower 5% quantile
  y_max <- quantile(df_sub$length, 0.99)  # Upper 95% quantile
  
  
  lengths[[trait]] =  ggplot(df_sub, aes(x = pval, y = length, fill = Method)) +
    geom_boxplot(outlier.shape = NA,size = .1, position = position_dodge()) +  # Bar plot with dodged bars
    labs(x = "P-value threshold", y = "Length",  title = paste("Length for Exp", toString(trait), "indirect") ) +
    theme_bw() +    theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
    coord_cartesian(ylim = c(y_min, y_max))  # Automatical
  
}

# indirect effect

covers_ind = list()
true_ind = c(0.83, .7)
for (trait in 1:(K - 1)){
  df_sub = df[df$Method == "FLOW-MR-indirect", ]
  df_sub = df_sub[df_sub$Trait == trait, ]
  df_sub$result <- with(df_sub, ifelse(`lower` <= true_ind[trait] & `upper` >= true_ind[trait], TRUE, FALSE))
  
  coverage_df_sub <- df_sub %>%
    group_by(Method, Trait, pval) %>%
    summarize(coverage = mean(result)) %>%
    ungroup()
  coverage_df_sub$lower = 0
  coverage_df_sub$upper = 0
  
  for(i in 1:nrow(coverage_df_sub)){
    coverage_df_sub$lower[i] =  as.vector( 
      prop.test(x = coverage_df_sub$coverage[i] * n, n = n)$conf.int)[1]
    coverage_df_sub$upper[i] =  as.vector( 
      prop.test(x = coverage_df_sub$coverage[i] * n, n = n)$conf.int)[2]
  }
  
  covers_ind[[trait]] = ggplot(coverage_df_sub, aes(x = pval, y = coverage, color = Method)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.5)) +
    geom_point(position = position_dodge(0.5)) +
    labs(x = "P-value threshold", y = "Coverage", title = paste("Coverage for Exp", toString(trait), "indirect")) +
    theme_bw()+ theme(legend.position = "none") + geom_abline(slope = 0, intercept = .95, linetype = "dashed")+   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}


lengths_ind = list()
for (trait in 1:(K - 1)){
  df_sub = df[df$Method == "FLOW-MR-indirect", ]
  
  df_sub = df_sub[df_sub$Trait == trait, ]
  df_sub = df_sub[df_sub$length < 2, ]
  
  y_min <- quantile(df_sub$length, 0.01)  # Lower 5% quantile
  y_max <- quantile(df_sub$length, 0.99)  # Upper 95% quantile
  
  
  
  lengths_ind[[trait]] =  ggplot(df_sub, aes(x = pval, y = length, fill = Method)) +
    geom_boxplot(outlier.shape = NA, size = .1, position = position_dodge()) +  # Bar plot with dodged bars
    labs(x = "P-value threshold", y = "Length",  title = paste("Length for Exp", toString(trait), "indirect") ) +
    theme_bw()+ theme(legend.position = "none")+   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_cartesian(ylim = c(y_min, y_max))  # Automatical
}







wrap_plots(list(covers[[1]]  + ggtitle(expression(X[1] ~ "→" ~ Y ~ "direct")),
                covers[[2]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[2] ~ "→" ~ Y ~ "direct")),
                covers[[3]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[3] ~ "→" ~ Y ~ "direct")),
                covers_ind[[1]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[1] ~ "→" ~ Y ~ "indirect")) + 
                  guides(color = "none", fill = "none", linetype = "none"),
                covers_ind[[2]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[2] ~ "→" ~ Y ~ "indirect")) + 
                  guides(color = "none", fill = "none", linetype = "none")
                
), nrow = 1) +
  plot_layout(guides = 'collect', widths = c(1.5, 1.5, 1.5,1,1),axis_titles = "collect") &
  theme(legend.position = "bottom")


wrap_plots(list(lengths[[1]]  + ggtitle(expression(X[1] ~ "→" ~ Y ~ "direct")),
                lengths[[2]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[2] ~ "→" ~ Y ~ "direct")),
                lengths[[3]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[3] ~ "→" ~ Y ~ "direct")),
                lengths_ind[[1]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[1] ~ "→" ~ Y ~ "indirect")) + 
                  guides(color = "none", fill = "none", linetype = "none"),
                lengths_ind[[2]] + theme(axis.title.y = element_blank()) + ggtitle(expression(X[2] ~ "→" ~ Y ~ "indirect")) + 
                  guides(color = "none", fill = "none", linetype = "none")
                
), nrow = 1) +
  plot_layout(guides = 'collect', widths = c(1.5, 1.5, 1.5,1,1),axis_titles = "collect") &
  theme(legend.position = "bottom")

df = read.csv("SuppFigure6/k=3_sim_cor_first.csv")
method_colors <- c(
  "IVW" = "#8494FF",  
  "FLOW-MR" = "#F8766D",  
  "GRAPPLE" = "#B39D2F",  
  "MR-Horse" = "#00A9FF",  
  "MR-Egger" = "#00BFC4",
  "MR-Robust" = "#FF61CC"
)

df$pval = factor((df$pval), levels = c(1e-8, 1e-6, 1e-4, 1e-2))
df$Trait = factor((df$Trait), levels = c(1, 2))
K = 1
n = 100
true_val = c(0.7)
covers = list()
for (trait in 1:K){
  df_sub = df[df$Trait == trait, ]
  df_sub = df[df$Method != "MVMR-cML", ]
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
    labs(x = "P-value threshold", y = "Coverage", title = paste("X1 -> Y")) +
    theme_bw() + ylim(NA,1)+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_abline(slope = 0, intercept = .95, linetype = "dashed")
  # Disable the legend in p1

  
}



# Compute the Length
lengths = list()
for (trait in 1:K){
  df_sub = df[df$Trait == trait, ]
  df_sub = df_sub[df_sub$length < 2, ]
  df_sub = df
  
  lengths[[trait]] =  ggplot(df_sub, aes(x = pval, y = length, fill = Method)) +
    geom_boxplot(outlier.shape = NA,size = .1, position = position_dodge()) +  # Bar plot with dodged bars
    labs(x = "P-value threshold", y = "Interval Width",  title = paste("X1 -> Y t")) +
    theme_bw() +   theme(legend.position = "none") +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +     scale_fill_manual(values = method_colors) 

  # Disable the legend in p1

  
}




## Means


means = list()
for (trait in 1:K){
  df_sub = df[df$Trait == trait, ]
  df_sub = df_sub[df_sub$Method != "Bayes-indirect", ]
  df_sub = df_sub[,c(1,2, 3, 4)]
  
  
  
  means[[trait]] = ggplot(df_sub, aes(x = pval, y = estimate, color = Method)) +
    stat_summary(fun = mean, geom = "point", size = 2,               position = position_dodge(width = 0.5))  +
    labs(x = "P-value threshold", y = "Mean Estimate",title = paste("X1 -> Y ")) +
    theme_bw() + geom_abline(slope = 0, intercept = .7)  +     scale_color_manual(values = method_colors) 
  
}

wrap_plots(list( covers[[1]]  + ggtitle(expression(X[1] ~ "→" ~ X[2] )) + 
                   guides(color = "none", fill = "none", linetype = "none"),
               lengths[[1]]  + ggtitle(expression(X[1] ~ "→" ~ X[2] )) , means[[1]]+ ggtitle(expression(X[1] ~ "→" ~ X[2] )) + 
                 guides(color = "none", fill = "none", linetype = "none") ),
            , nrow = 1) +
  plot_layout(guides = 'collect', widths = c(1,1,1),axis_titles = "collect")  &
  theme(legend.position = "bottom")



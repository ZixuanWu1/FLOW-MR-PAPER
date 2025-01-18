df = read.csv("k=3_sim_cor.csv")

df$pval = factor((df$pval), levels = c(1e-8, 1e-6, 1e-4, 1e-2))
df$Trait = factor((df$Trait), levels = c(1, 2))


# Compute the coverage proportion
true_val = c(0.3, 1)
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
    labs(x = "pval", y = "Coverage", title = paste("Coverage for Exp", toString(trait))) +
    theme_bw()  + geom_abline(slope = 0, intercept = .95, linetype = "dashed")+   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}



# Compute the Length
lengths = list()
for (trait in 1:K){
  df_sub = df[df$Trait == trait, ]
  df_sub = df_sub[df_sub$Method != "FLOW-MR-indirect", ]
  df_sub = df_sub[df_sub$length < 2, ]
  
  
  lengths[[trait]] =  ggplot(df_sub, aes(x = pval, y = length, fill = Method)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(), size = .1) +  # Bar plot with dodged bars
    labs(x = "pval", y = "Interval Width",  title = paste("Length for Exp", toString(trait)) ) +
    theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

# indirect effect
df_sub = df[df$Method == "FLOW-MR-indirect", ]


covers_ind = list()
true_ind = c(0.7)
for (trait in 1:(K - 1)){
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
    labs(x = "pval", y = "Coverage", title = paste("Coverage for Exp", toString(trait), "indirect")) +
    theme_bw() + theme(legend.position = "none") +   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

lengths_ind = list()
for (trait in 1:(K - 1)){
  df_sub = df_sub[df_sub$Trait == trait, ]
  df_sub = df_sub[df_sub$length < 2, ]
  
  
  lengths_ind[[trait]] =  ggplot(df_sub, aes(x = pval, y = length, fill = Method)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(), size = .1) +  # Bar plot with dodged bars
    labs(x = "pval", y = "Interval Width",  title = paste("Length for Exp", toString(trait), "indirect") ) +
    theme_bw() + theme(legend.position = "none")+   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}



wrap_plots(list(
  covers[[2]] + ggtitle(expression(X[2] ~ "→" ~ Y ~ "direct")),
  covers[[1]] + ggtitle(expression(X[1] ~ "→" ~ Y ~ "direct")) + theme(axis.title.y = element_blank()),
  covers_ind[[1]] + 
    theme(axis.title.y = element_blank()) + 
    ggtitle(expression(X[1] ~ "→" ~ Y ~ "indirect")) + 
    guides(color = "none", fill = "none", linetype = "none") # Suppress legend elements
), nrow = 1) +
  plot_layout(guides = 'collect', widths = c(1.5, 1.5, 1), axis_titles = "collect") &
  theme(legend.position = "bottom")


wrap_plots(list(                lengths[[2]]  + ggtitle(expression(X[2] ~ "→" ~ Y ~ "direct")) , 
                                lengths[[1]]  + ggtitle(expression(X[1] ~ "→" ~ Y ~ "direct"))+ theme(axis.title.y = element_blank()),
                                lengths_ind[[1]] + theme(axis.title.y = element_blank()) + 
                                  ggtitle(expression(X[1] ~ "→" ~ Y ~ "indirect")) +  guides(color = "none", fill = "none", linetype = "none")
), nrow = 1) +
  plot_layout(guides = 'collect', widths = c(1.5, 1.5, 1),axis_titles = "collect")  &
  theme(legend.position = "bottom")
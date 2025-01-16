library(gridExtra)
library(ggplot2)
library(ggforce)

bodysize3 = read.csv("k=3_real_bodysize.csv")
bodysize3$threshold <- factor(bodysize3$threshold, levels = c(1e-8, 1e-6, 1e-4, 1e-2))
bodysize3$exposure[bodysize3$exposure =="Childhood Body Size"] = "early-life Body Size"
bodysize3$exposure = factor(bodysize3$exposure, levels = c("early-life Body Size", "Adult BMI"))

bmi3 = read.csv("k=3_real_bmi.csv")
bmi3$threshold <- factor(bmi3$threshold, levels = c(1e-8, 1e-6, 1e-4, 1e-2))
bmi3$exposure[bmi3$exposure =="Childhood BMI"] = "8-year-old BMI"

bmi4 = read.csv("k=4_real.csv")
bmi4$threshold <- factor(bmi4$threshold, levels = c(1e-8, 1e-6, 1e-4, 1e-2))



plot_mv <-function(table, ylimit  = NA, methods = c("MVMR-IVW",
                                                   "GRAPPLE",
                                                   "FLOW-MR"),
                  facet = F, facety = c(-3, 3), facetl= F, facetly = c(0, 2), K = 3){
  table = table[which(table$method %in% methods),]
  mmm =  c("MVMR-IVW",
           "GRAPPLE",
           "FLOW-MR")
  if(K == 2){
    mmm[1] = "IVW"
  }
  tables = list()
  tables[[1]] = table[which(table$method == "MVMR-IVW"),]
  tables[[2]] = table[which(table$method == "GRAPPLE"),]
  tables[[3]] = table[which(table$method == "FLOW-MR"),]
  p = list()
  for(i in 1:3){
    
    df = tables[[i]]
    p[[i]]<- ggplot(df, aes(x=threshold, y=beta_hat, fill=exposure)) + 
      geom_point(aes( color = factor(exposure)),
                 position=position_dodge(.9), show.legend = FALSE) +
      geom_errorbar(aes(ymin=lower, ymax=upper, color = factor(exposure)), width=.2,
                    position=position_dodge(.9))  + theme_bw() +
      labs( x =  "p-value Threshold", y = expression(hat(beta)),
            title = mmm[i] ) + geom_hline(yintercept=0,linetype=2)+scale_x_discrete(labels= c(1e-8,                                                                                           1e-6,
                                                                                              1e-4,
                                                                                              1e-2)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))  + theme(legend.title = element_blank()) + scale_fill_manual(values = c("red", "blue", "green"))
    if(i == 2 & facet == T){
      p[[i]] = p[[i]] + facet_zoom(ylim = facety)
    }
    if(!is.na(ylimit[1])){
      p[[i]] = p[[i]] + ylim(ylimit)
    }
  }
  
  for(i in 1:3){
    df = tables[[i]]
    p[[i + 3]]<- ggplot(df, aes(x=threshold, y= (upper - lower), fill=exposure)) + 
      geom_point(aes( color = factor(exposure)),
                 position=position_dodge(.9), show.legend = FALSE) +
      geom_errorbar(aes(ymin=0, ymax=(upper - lower), color = factor(exposure)), 
                    width=.2, position=position_dodge(.9)) +
      theme_bw() + labs( x= "p-value Threshold", y = "Length", 
                         title = mmm[i]) + scale_x_discrete(labels= c(1e-8,
                                                                      1e-6,
                                                                      1e-4,
                                                                      1e-2)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + theme(legend.title = element_blank())
    if(i == 2 & facetl == T){
      p[[i + 3]] = p[[i + 3]] + facet_zoom(ylim = facetly)
    }
    
    
  }
  
  wrap_plots(p[[1]]  , p[[2]],p[[3]],  nrow = 1)   + plot_layout(guides = "collect", axis_titles = "collect") &  theme(legend.position = "bottom")

}


plot_mv(bodysize3)
plot_mv(bmi3, facet = T, facety = c(-2,2))
plot_mv(bmi4)


flow_4 = readRDS( "k=4_flow.RData")
flow_4

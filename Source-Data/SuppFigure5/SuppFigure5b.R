library(ggpubr)
library(ggplot2)
library(gridExtra)
df = read.csv("SuppFigure5/sbpc_bmia_sbpa.csv")
df$threshold = factor(df$threshold)
plot_mv_K <- function(table, ylimit  = NA, methods = c("MVMR-IVW",
                                          "GRAPPLE",
                                          "FLOW-MR"), pvals = c(1e-8, 1e-6, 1e-4, 1e-2)){
  table = table[which(table$method %in% methods),]
  mmm = c("MVMR-IVW",
          "GRAPPLE",
          "FLOW-MR")
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
      labs( x =  "log threshold", y = expression(hat(beta)),
            title = mmm[i] ) + geom_hline(yintercept=0,linetype=2)+scale_x_discrete(labels= pvals)
    if(!is.na(ylimit[1])){
      p[[i]] = p[[i]] + ylim(ylimit)
    }
  }
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  grid.arrange(p[[1]]  + theme(legend.position="none") , p[[2]]  + theme(legend.position="none") , p[[3]] + theme(legend.position="none"),
               mylegend = g_legend(p[[1]]), nrow = 1)
  

}

plot_mv_K(df, methods = c("MVMR-IVW", "GRAPPLE", "FLOW-MR"))


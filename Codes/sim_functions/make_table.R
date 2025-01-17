require(gridExtra)
library(ggplot2)

maketable <- function(result, CI_level = 0.95,   pvals = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)){

  methods = c("IVW", "Egger", "Weighted_median", "NOME", "Grapple", "Bayes_first",
              "Bayes_second")
  n = length(pvals)
  table = matrix(nrow = n * 7, ncol = 7 )
  for(i in 1:n){
    for(j in 1:4){
      table[ ((i - 1) *7 + j ), 1] = methods[j]
      table[ ((i - 1) *7 + j ), 2] = pvals[i]
      res =  result[[methods[j]]][[i]]
      table[ ((i - 1) *7 + j ), 3] = res$b 
      table[ ((i - 1) *7 + j ), 4] = res$se
      table[ ((i - 1) *7 + j ), 5] = res$b - qnorm((CI_level + 1)/2) * res$se
      table[ ((i - 1) *7 + j ), 6] = res$b + qnorm((CI_level + 1)/2) * res$se
      table[ ((i - 1) *7 + j ), 7] = res$pval
    }
    for(j in 5:5){
      table[ ((i - 1) *7 + j ), 1] = methods[j]
      table[ ((i - 1) *7 + j ), 2] = pvals[i]
      grapple =  result[[methods[j]]][[i]]
      table[ ((i - 1) *7 + j ), 3] = grapple$beta.hat
      table[ ((i - 1) *7 + j ), 4] = sqrt(grapple$beta.var)
      table[ ((i - 1) *7 + j ), 5] = grapple$beta.hat - 1.96 * sqrt(grapple$beta.var)
      table[ ((i - 1) *7 + j ), 6] = grapple$beta.hat + 1.96 * sqrt(grapple$beta.var)
      table[ ((i - 1) *7 + j ), 7] = grapple$beta.p.value
    }
    
    for(j in 6:7){
      table[ ((i - 1) *7 + j ), 1] = methods[j]
      table[ ((i - 1) *7 + j ), 2] = pvals[i]
      bayes =  result[[methods[j]]][[i]]
      table[ ((i - 1) *7 + j ), 3] = bayes[1,1]
      table[ ((i - 1) *7 + j ), 4] = bayes[1,3]
      table[ ((i - 1) *7 + j ), 5] =  bayes[1,4]
      table[ ((i - 1) *7 + j ), 6] = bayes[1, 6]
      table[ ((i - 1) *7 + j ), 7] = NA
    }
  }
  
  
  
  
  colnames(table) = c("method", "threshold", "beta_hat", "se", "lower", "upper", "p-value")
  table <- as.data.frame(table)
  for(i in c(3, 4,5,6,7)){
    table[, i] = as.numeric(table[, i])
  }
  table[,2] = factor(table[,2], levels = pvals)
  return(table)
}



maketable_MV <- function(result, exps = c("BMI child", "BMI adult") ){
  mr_result = result$MVMR
  mv_grapple = result$Grapple
  mv_first= result$first
  mv_second = result$second
  pvals = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
  table = matrix(nrow = 56, ncol = 8 )
  for(i in 1:7){
    table[(8 * i - 7) : (8 * i), 1] = rep(exps, 4)
    
    table[(8 * i - 7): (8 * i), 2] = rep(pvals[i], 8)
    
    table[(8 * i - 7): (8 * i), 3] = c("MVMR", 
                                       "MVMR", "Grapple", 
                                       "Grapple", "Bayes_first", "Bayes_first",
                                       "Bayes_second", "Bayes_second")
    
    #Estimate
    table[c( (8 * i - 7), (8 * i - 6)     ), 4] = mr_result[[i]][, 1]
    table[c( (8 * i - 5), (8 * i -4)     ), 4] = mv_grapple[[i]]$beta.hat
    table[c( (8 * i - 3), (8 * i - 2)     ), 4] = mv_first[[i]][c(2, 1), 1]
    table[c( (8 * i - 1), (8 * i - 0)     ), 4] = mv_second[[i]][c(2, 1), 1]
    #SE
    table[c( (8 * i - 7), (8 * i - 6)     ), 5] = mr_result[[i]][, 2]
    table[c( (8 * i - 5), (8 * i - 4)     ), 5] = 
      sqrt(c( mv_grapple[[i]]$beta.var[1,1],mv_grapple[[i]]$beta.var[2,2] ))
    table[c( (8 * i - 3), (8 * i - 2)     ), 5] = mv_first[[i]][c(2, 1), 3]
    table[c( (8 * i - 1), (8 * i - 0)     ), 5] = mv_second[[i]][c(2, 1), 3]
    #LOWER
    table[c( (8 * i - 7), (8 * i -6)     ), 6] = mr_result[[i]][, 1] -
      1.96 *  mr_result[[i]][, 2]
    
    table[c( (8 * i - 5), (8 * i - 4)     ), 6] =  mv_grapple[[i]]$beta.hat - 1.96 * 
      sqrt(c( mv_grapple[[i]]$beta.var[1,1],mv_grapple[[i]]$beta.var[2,2] ))
    
    table[c( (8 * i - 3), (8 * i - 2)     ), 6] = mv_first[[i]][c(2, 1), 4]
    table[c( (8 * i - 1), (8 * i - 0)     ), 6] = mv_second[[i]][c(2, 1), 4]
    #UPPER
    
    table[c( (8 * i - 7), (8 * i - 6)     ), 7] = mr_result[[i]][, 1] +
      1.96 *  mr_result[[i]][, 2]
    
    table[c( (8 * i - 5), (8 * i - 4)     ), 7] =  mv_grapple[[i]]$beta.hat + 1.96 * 
      sqrt(c( mv_grapple[[i]]$beta.var[1,1],mv_grapple[[i]]$beta.var[2,2] ))
    
    table[c( (8 * i - 3), (8 * i - 2)     ), 7] = mv_first[[i]][c(2,1), 6]
    table[c( (8 * i - 1), (8 * i - 0)     ), 7] = mv_second[[i]][c(2,1), 6]
    
    #P-value
    table[c( (8 * i - 7), (8 * i - 6)     ), 8] = mr_result[[i]][, 4]
    table[c( (8 * i - 5), (8 * i - 4)     ), 8] =  mv_grapple[[i]]$beta.p.value
    table[(8 * i - 3):(8 * i - 0)    , 8] = rep(NA, 4)
  }
  
  colnames(table) = c("exposure", "threshold", "method", "beta_hat", 
                      "se", "lower", "upper", "p-value")
  table <- as.data.frame(table)
  for(i in c(4,5,6,7, 8)){
    table[, i] = as.numeric(table[, i])
  }
  table[,2] = factor(table[,2], levels = pvals)
  return(table)
}

plot_uni <-function(table, ylimit = NA, methods =  c("IVW", "Egger", "Weighted_median", "NOME", "Grapple", "Bayes_first",
                                                     "Bayes_second")
                    
                    ){
  
  table = table[which(table$method %in% methods),]
  p1<- ggplot(table, aes(x=threshold, y=beta_hat, fill=method)) + 
    geom_point(aes(color = factor(method)),
               position=position_dodge(.9), show.legend = FALSE) +
    geom_errorbar(aes(ymin=lower, ymax=upper, color = factor(method)), width=.2,
                  position=position_dodge(.9))  + theme_bw() +
    labs( x= "threshold", y = expression(hat(beta))) + geom_hline(yintercept=0,linetype=2)
  if(!is.na(ylimit)){
    p1 = p1 + ylim(ylimit)
  }
  
  
  p2<- ggplot(table, aes(x=threshold, y= (upper - lower), fill=method)) + 
    geom_point(aes(color = factor(method)),
               position=position_dodge(.9), show.legend = FALSE) +
    geom_errorbar(aes(ymin=0, ymax=(upper - lower), color = factor(method)), width=.2,
                  position=position_dodge(.9)) +
    theme_bw() +
    labs( x= "threshold", y = "Length") + geom_hline(yintercept=0,linetype=1)
  
  if(!is.na(ylimit)){
    p2 = p2 + ylim(ylimit)
  }
  
  grid.arrange(p1, p2,nrow = 2)

  
  
}

plot_mv <-function(table, ylimit  = NA, methods = c("MVMR",
                                                    "Grapple",
                                                    "Bayes_first",
                                                    "Bayes_second")){
  table = table[which(table$method %in% methods),]
  mmm = c("MVMR",
          "Grapple",
          "Bayes")
  tables = list()
  tables[[1]] = table[which(table$method == "MVMR"),]
  tables[[2]] = table[which(table$method == "Grapple"),]
  tables[[3]] = table[which(table$method == "Bayes_first"),]
  tables[[4]] = table[which(table$method == "Bayes_second"),]
  p = list()
  for(i in 1:3){
    df = tables[[i]]
    p[[i]]<- ggplot(df, aes(x=threshold, y=beta_hat, fill=exposure)) + 
      geom_point(aes( color = factor(exposure)),
                 position=position_dodge(.9), show.legend = FALSE) +
      geom_errorbar(aes(ymin=lower, ymax=upper, color = factor(exposure)), width=.2,
                    position=position_dodge(.9))  + theme_bw() +
      labs( x =  "log threshold", y = expression(hat(beta)),
      title = mmm[i] ) + geom_hline(yintercept=0,linetype=2)+scale_x_discrete(labels= -8:-1)
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
  
  for(i in 1:3){
    df = tables[[i]]
    p[[i]]<- ggplot(df, aes(x=threshold, y= (upper - lower), fill=exposure)) + 
      geom_point(aes( color = factor(exposure)),
                 position=position_dodge(.9), show.legend = FALSE) +
      geom_errorbar(aes(ymin=0, ymax=(upper - lower), color = factor(exposure)), 
                    width=.2, position=position_dodge(.9)) +
      theme_bw() + labs( x= "log threshold", y = "Length", 
            title = mmm[i]) + scale_x_discrete(labels= -8:-1)
    

  }
  
  grid.arrange(p[[1]]  + theme(legend.position="none") , p[[2]]  + theme(legend.position="none") , p[[3]] + theme(legend.position="none"),
               mylegend = g_legend(p[[1]]), nrow = 1)
}


plot_mv2 <-function(table, ylimit  = NA, methods = c("MVMR",
                                                    "Grapple",
                                                    "Bayes_first",
                                                    "Bayes_second"), th = -8:-1){
  table = table[which(table$method %in% methods),]
  table = table[ ((th[1] + 8) * 8 + 1):56,]
  tables = list()
  tables[[1]] = table[which(table$method == "MVMR"),]
  tables[[2]] = table[which(table$method == "Grapple"),]
  tables[[3]] = table[which(table$method == "Bayes_first"),]
  tables[[4]] = table[which(table$method == "Bayes_second"),]
  p = list()
  for(i in 1:4){
    df = tables[[i]]
    p[[i]]<- ggplot(df, aes(x=threshold, y=beta_hat, fill=exposure)) + 
      geom_point(aes( color = factor(exposure)),
                 position=position_dodge(.9), show.legend = FALSE) +
      geom_errorbar(aes(ymin=lower, ymax=upper, color = factor(exposure)), width=.2,
                    position=position_dodge(.9))  + theme_bw() +
      labs( x =  "log threshold", y = expression(hat(beta)),
            title = toString(tables[[i]]$method[1]) ) + geom_hline(yintercept=0,linetype=2)+scale_x_discrete(labels= th)
    
  }
  
  grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 2, nrow = 2)
  
  for(i in 1:4){
    df = tables[[i]]
    p[[i]]<- ggplot(df, aes(x=threshold, y= (upper - lower), fill=exposure)) + 
      geom_point(aes( color = factor(exposure)),
                 position=position_dodge(.9), show.legend = FALSE) +
      geom_errorbar(aes(ymin=0, ymax=(upper - lower), color = factor(exposure)), 
                    width=.2, position=position_dodge(.9)) +
      theme_bw() + labs( x= "log threshold", y = "Length", 
                         title = toString(tables[[i]]$method[1])) + scale_x_discrete(labels= th)
    
    if(!is.na(ylimit)){
      p[[i]] = p[[i]] + ylim(ylimit)
    }
  }
  
  grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 2, nrow = 2)
}





maketable_MV3 <- function(result, exps = c("1-year old BMI", " 8-year-old BMI", "Adult BMI") ){
  mr_result = result$MVMR
  mv_grapple = result$Grapple
  mv_first= result$first
  mv_second = result$second
  pvals = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
  table = matrix(nrow = 84, ncol = 8 )
  for(i in 1:7){
    table[(12 * i - (11)) : (12 * i), 1] = rep(exps, 4)
    
    table[(12 * i - 11): (12 * i), 2] = rep(pvals[i], 12)
    
    table[(12 * i - 11): (12 * i), 3] = c("MVMR", 
                                       "MVMR","MVMR",  "Grapple",  "Grapple", 
                                       "Grapple", "Bayes_first", "Bayes_first",
                                       "Bayes_first","Bayes_second", "Bayes_second", "Bayes_second")
    
    #Estimate
    table[c( (12 * i - 11): (12 * i - 9)     ), 4] = mr_result[[i]][, 1]
    table[c( (12 * i - 8): (12 * i -6)     ), 4] = mv_grapple[[i]]$beta.hat
    table[c( (12 * i - 5): (12 * i - 3)     ), 4] = mv_first[[i]][c(3, 2, 1), 1]
    table[c( (12 * i - 2): (12 * i - 0)     ), 4] = mv_second[[i]][c(3, 2, 1), 1]
    #SE
    table[c( (12 * i - 11): (12 * i - 9)     ), 5]= mr_result[[i]][, 2]
    table[c( (12 * i - 8): (12 * i -6)     ), 5] = 
      sqrt(c( mv_grapple[[i]]$beta.var[1,1],mv_grapple[[i]]$beta.var[2,2], mv_grapple[[i]]$beta.var[3,3] ))
    table[c( (12 * i - 5): (12 * i - 3)     ), 5] = mv_first[[i]][c(3, 2, 1), 3]
    table[c( (12 * i - 2): (12 * i - 0)     ), 5] = mv_second[[i]][c(3, 2, 1), 3]
    #LOWER
    table[c( (12 * i - 11): (12 * i - 9)     ), 6] = mr_result[[i]][, 1] -
      1.96 *  mr_result[[i]][, 2]
    
    table[c( (12 * i - 8): (12 * i -6)     ), 6] =  mv_grapple[[i]]$beta.hat - 1.96 * 
      sqrt(c( mv_grapple[[i]]$beta.var[1,1],mv_grapple[[i]]$beta.var[2,2], mv_grapple[[i]]$beta.var[3,3] ))
    
    table[c( (12 * i - 5): (12 * i - 3)     ), 6]= mv_first[[i]][c(3, 2, 1), 4]
    table[c( (12 * i - 2): (12 * i - 0)     ), 6] = mv_second[[i]][c(3, 2, 1), 4]
    #UPPER
    
    table[c( (12 * i - 11): (12 * i - 9)     ), 7] = mr_result[[i]][, 1] +
      1.96 *  mr_result[[i]][, 2]
    
    table[c( (12 * i - 8): (12 * i -6)     ), 7] =  mv_grapple[[i]]$beta.hat + 1.96 * 
      sqrt(c( mv_grapple[[i]]$beta.var[1,1],mv_grapple[[i]]$beta.var[2,2], mv_grapple[[i]]$beta.var[3,3] ))
    
    table[c( (12 * i - 5): (12 * i - 3)     ), 7]= mv_first[[i]][c(3, 2, 1), 6]
    table[c( (12 * i - 2): (12 * i - 0)     ), 7] = mv_second[[i]][c(3, 2, 1), 6]
    
    #P-value
    table[c( (12 * i - 11): (12 * i - 9)     ), 8] = mr_result[[i]][, 4]
    table[c( (12 * i - 8): (12 * i - 6)     ), 8] =  mv_grapple[[i]]$beta.p.value
    table[(8 * i - 5):(8 * i - 0)    , 8] = rep(NA, 6)
  }
  
  colnames(table) = c("exposure", "threshold", "method", "beta_hat", 
                      "se", "lower", "upper", "p-value")
  table <- as.data.frame(table)
  for(i in c(4,5,6,7, 8)){
    table[, i] = as.numeric(table[, i])
  }
  table[,2] = factor(table[,2], levels = pvals)
  return(table)
}




maketable_MV_K <- function(result, exps = c("1-year old BMI", " 8-year-old BMI", "Adult BMI"),
                           pvals= c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2) ){
  K = length(exps)
  mr_result = result$MVMR
  mv_grapple = result$Grapple
  mv_first= result$first
  n_pval = length(pvals)
  table = matrix(nrow = 3 * K * n_pval, ncol = 8 )
  M = 3 * K
  for(i in 1:n_pval){
    table[(M * i - (M - 1)) : (M * i), 1] = rep(exps, 3)
    
    table[(M * i - (M - 1)): (M * i), 2] = rep(pvals[i], M)
    
    table[(M * i - (M - 1)): (M * i), 3] = c(rep("MVMR", K), rep("Grapple", K), rep("Bayes", K))
    
    #Estimate
    table[c( (M * i - (M - 1)): (M * i - (M - K))     ), 4] = mr_result[[i]][, 1]
    table[c( (M * i - (M - K - 1)): (M * i - (M - 2 * K))     ), 4] = mv_grapple[[i]]$beta.hat
    table[c( (M * i - (M - 2 * K - 1)): (M * i - 0)     ), 4] = mv_first[[i]][seq(K, 1, by = -1), 1]
    #SE
    table[c( (M * i - (M - 1)): (M * i - (M - K))     ), 5]= mr_result[[i]][, 2]
    vars = c()
    if( K == 1){
      vars = c(vars, mv_grapple[[i]]$beta.var)
    }
    else{
      for (j in 1:K){
        vars = c(vars, mv_grapple[[i]]$beta.var[j, j])
      }
    }

    table[c( (M * i - (M - K - 1)): (M * i - (M - 2 * K))     ), 5] = sqrt(vars)
    table[c( (M * i - (M - 2 * K - 1)): (M * i - 0)     ), 5] = mv_first[[i]][seq(K, 1, by = -1), 3]
    #LOWER
    table[c( (M * i - (M - 1)): (M * i - (M - K))     ), 6] = mr_result[[i]][, 1] -
      1.96 *  mr_result[[i]][, 2]
    table[c( (M * i - (M - K - 1)): (M * i - (M - 2 * K))     ), 6] =  mv_grapple[[i]]$beta.hat - 1.96 * sqrt(vars)
    table[c( (M * i - (M - 2 * K - 1)): (M * i - 0)     ), 6]= mv_first[[i]][seq(K, 1, by = -1), 4]
    #UPPER
    
    table[c( (M * i - (M - 1)): (M * i - (M - K))     ), 7] = mr_result[[i]][, 1] +
      1.96 *  mr_result[[i]][, 2]
    
    table[c( (M * i - (M - K - 1)): (M * i - (M - 2 * K))     ), 7] =  mv_grapple[[i]]$beta.hat + 1.96 * sqrt(vars)
    
    table[c( (M * i - (M - 2 * K - 1)): (M * i - 0)     ), 7]= mv_first[[i]][seq(K, 1, by = -1), 6]

    #P-value
    table[c( (M * i - (M - 1)): (M * i - (M - K))     ), 8] = mr_result[[i]][, 4]
    table[c( (M * i - (M - K - 1)): (M * i - (M - 2 * K))     ), 8] =  mv_grapple[[i]]$beta.p.value
    table[c( (M * i - (M - 2 * K - 1)): (M * i - 0)     ) , 8] = rep(NA, K)
  }
  
  colnames(table) = c("exposure", "threshold", "method", "beta_hat", 
                      "se", "lower", "upper", "p-value")
  table <- as.data.frame(table)
  for(i in c(4,5,6,7, 8)){
    table[, i] = as.numeric(table[, i])
  }
  table[,2] = factor(table[,2], levels = pvals)
  return(table)
}


plot_mv_K <-function(table, ylimit  = NA, methods = c("MVMR",
                                                    "Grapple",
                                                    "Bayes"), pvals = c(1e-8, 1e-6, 1e-4)){
  table = table[which(table$method %in% methods),]
  mmm = c("MVMR",
          "Grapple",
          "Bayes")
  tables = list()
  tables[[1]] = table[which(table$method == "MVMR"),]
  tables[[2]] = table[which(table$method == "Grapple"),]
  tables[[3]] = table[which(table$method == "Bayes"),]
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
  
  for(i in 1:3){
    df = tables[[i]]
    p[[i]]<- ggplot(df, aes(x=threshold, y= (upper - lower), fill=exposure)) + 
      geom_point(aes( color = factor(exposure)),
                 position=position_dodge(.9), show.legend = FALSE) +
      geom_errorbar(aes(ymin=0, ymax=(upper - lower), color = factor(exposure)), 
                    width=.2, position=position_dodge(.9)) +
      theme_bw() + labs( x= "log threshold", y = "Length", 
                         title = mmm[i]) + scale_x_discrete(labels= -8:-1)
    
    
  }
  
  grid.arrange(p[[1]]  + theme(legend.position="none") , p[[2]]  + theme(legend.position="none") , p[[3]] + theme(legend.position="none"),
               mylegend = g_legend(p[[1]]), nrow = 1)
}


library(extraDistr)
library(MVMR)
library(GRAPPLE)
library(FLOWMR)
library(extraDistr)
library(patchwork)
library(doParallel)
source("t1d_table.R")
source("gibbs_vector.R")

runMVMR3 <- function(data, 
                     p_value =  c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
                     init = "Random",iter = 16000, warmup = 8000, bayes = TRUE, cor_mat = NULL, K  = 3, alpha0 = 3, alpha1 = 3){
  nn <- colnames(data)
  result_grapple = list()
  result_bayes = list()
  result = list()
  indirect = list()
  
  total = list()
  
  n_snp = rep(0, length(p_value))
  F_stat = matrix(0, length(p_value), 2)
  for(i in 1:length(p_value)){
    dat = data[which(data$selection_pvals < p_value[i]),]
    
    if("SNP" %in% nn){
      F.data <- format_mvmr(BXGs = data[data$selection_pvals < p_value[i], grep("gamma_exp", nn)],
                            BYG = data[data$selection_pvals < p_value[i],"gamma_out1"],
                            seBXGs = data[data$selection_pvals < p_value[i], grep("se_exp", nn)],
                            seBYG = data[data$selection_pvals < p_value[i], "se_out1"],
                            RSID = data[data$selection_pvals < p_value[i], "SNP"])
      
    }
    else{
      F.data <- format_mvmr(BXGs = data[data$selection_pvals < p_value[i], grep("gamma_exp", nn)],
                            BYG = data[data$selection_pvals < p_value[i],"gamma_out1"],
                            seBXGs = data[data$selection_pvals < p_value[i], grep("se_exp", nn)],
                            seBYG = data[data$selection_pvals < p_value[i], "se_out1"])
      
    }
    
    result[[i]] = ivw_mvmr(r_input = F.data)
    
    n_snp[i] = dim(dat)[1]
    str = strength_mvmr(r_input = F.data, gencov = 0)
    F_stat[i, ] = c(str$exposure1, str$exposure2)
    
    result_grapple[[i]] = 
      GRAPPLE::grappleRobustEst(dat, plot.it = FALSE, cor.mat = cor_mat)
    
    if(bayes){
      if(K == 3){
        new_order = 3:1
        Gamma_hat =rbind(dat$gamma_out1, 
                         dat$gamma_exp2,
                         dat$gamma_exp1)
        
        Sd_hat = rbind(dat$se_out1, 
                       dat$se_exp2,
                       dat$se_exp1)
        result1 = BayesMediation(Gamma_hat, Sd_hat, inv = TRUE,
                                 init = init, iter = iter, warmup  = warmup, , cor = cor_mat[new_order, new_order], Raw = T,
                                 indirect = T, alpha_0 = alpha0, alpha_1 = alpha1)
      }
      
      result_bayes[[i]] = result1$summary
      total[[i]] = FLOWMR::total_effect( result1$raw, 3)
      indirect[[i]] = FLOWMR::indirect_effect( result1$raw, 3, path = "all")
      
    }
    
    
  }
  
  if(bayes){
    
    
    
    return(list(MVMR = result, Grapple = result_grapple, 
                first = result_bayes, indirect = indirect, total = total, n_snp = n_snp, F_stat = F_stat))
  }
  else{
    return( list(MVMR = result, Grapple = result_grapple
    )  )
  }
}


# BMI Adult Analysis

sel.file <- c("BMIchild-egg15.csv" )
exp.file <- c("BMI1year-moba19", "BMI8year-moba19")
out.file <- "BMIadult-ukb.csv"
plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"
data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01,
                               plink_exe = "plink_mac_20210606/plink", cal.cor = T)


filtered.data <-data.list$data;
mv_breast_bmi = runMVMR3(filtered.data, p_value =c(1e-8, 1e-6, 1e-4, 1e-3), iter = 6000,
                         warmup = 3000, cor_mat = data.list$cor.mat)

# SBP Adult Analysis



sel.file <- c( "SBP-gera17.csv", "BMI-giant17eu.csv")
exp.file <- c("sbp_child.csv", "BMIadult-ukb.csv" )
out.file <-"SBP-ukb.csv"
plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"


data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01,
                               plink_exe = "plink_mac_20210606/plink", cal.cor = T)
filtered.data <-data.list$data;
mv_sbp_bmi = runMVMR3(filtered.data, p_value =c(1e-8, 1e-6, 1e-4, 1e-3), iter = 16000,
                         warmup = 8000, cor_mat = data.list$cor.mat, alpha0 = 5, alpha1= 5)
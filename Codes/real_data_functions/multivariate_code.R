library(doParallel)
library(FLOWMR)
source("gibbs_vector.R")
library(MVMR)
library(GRAPPLE)
library(MASS)
library(TwoSampleMR)

runMVMR7 <- function(data, 
                     p_value =  c(1e-8,  1e-6,  1e-4,  1e-2),
                     init = "Random",iter = 16000, warmup = 8000, bayes =
                       TRUE, cor_mat = NULL, n_traits = c(1,3,3),  alpha_0 = 1.1, alpha_1 = 1.1, indirect_est = FALSE){
  nn <- colnames(data)
  result_grapple = list()
  result_bayes = list()
  result = list()
  indirect = list()
  for(i in 1:length(p_value)){
    dat = data[which(  data$selection_pvals  < p_value[i]),]
    
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
    
    
    result_grapple[[i]] = 
      GRAPPLE::grappleRobustEst(dat, plot.it = FALSE, cor.mat = cor_mat)
    
    if(bayes){
      Gamma_hat =rbind(dat$gamma_out1, 
                       dat$gamma_exp6,
                       dat$gamma_exp5,
                       dat$gamma_exp4,
                       dat$gamma_exp3,
                       dat$gamma_exp2,
                       dat$gamma_exp1)
      
      Sd_hat = rbind(dat$se_out1, 
                     dat$se_exp6,
                     dat$se_exp5,
                     dat$se_exp4,
                     dat$se_exp3,
                     dat$se_exp2,
                     dat$se_exp1)
      
      new_order = c(7, 6,5,4,3,2,1)
      result1 = BayesMediation(Gamma_hat, Sd_hat, inv = TRUE,
                               init = init, iter = iter, warmup  = warmup, cor = cor_mat[new_order,new_order], Raw = F,
                               indirect = indirect_est, n_traits = n_traits, alpha_0 = alpha_0,
                               alpha_1 = alpha_1)
      result_bayes[[i]] = result1$summary
      indirect[[i]] = result1$indirect_effect
    }
    
    
  }
  
  if(bayes){
    
    return(list(MVMR = result, Grapple = result_grapple, 
                first = result_bayes, indirect = indirect))
  }
  else{
    return( list(MVMR = result, Grapple = result_grapple
    )  )
  }
}


sel.file <- c("LDL-gera18.csv", "BMI-giant17eu.csv", "SBP-gera17.csv", "bmi_child_agg.csv")
exp.file <- c("ldl_child.csv",  "bmi_child.csv", "sbp_child.csv",
              "LDL-glgc13.csv",  "BMI_adult.csv", "SBP-ukb.csv")
out.file <- "AS-Malik18EU.csv"
plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"


data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =1e-4,
                               plink_exe = "plink_mac_20210606/plink", cal.cor = T)

filtered.data <-data.list$data;
data = filtered.data
cor_mat = data.list$cor.mat

result = runMVMR7(data = data, cor_mat = cor_mat, p_value = c(1e-8, 1e-6, 1e-4), 
                  n_traits = c(1,1, 2, 1, 2), alpha_0 = 3, alpha_1 = 3, indirect_est= T, iter = 16000, warmup = 8000)



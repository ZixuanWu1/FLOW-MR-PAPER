source("sim_functions.R")

runMVMR7 <- function(data, 
                     p_value =  c(1e-8,  1e-6,  1e-4,  1e-2),
                     init = "Random",iter = 16000, warmup = 8000, cor_mat = NULL, n_traits = c(1,3,3),  alpha_0 =6, alpha_1 = 6){
  nn <- colnames(data)
  result_flow = list()
  n_snp = rep(0, length(p_value))
  F_stat = matrix(0, length(p_value), 6)
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
    
    
    n_snp[i] = dim(dat)[1]
    print(n_snp)
    str = strength_mvmr(r_input = F.data, gencov = 0)
    F_stat[i, ] = c(str$exposure1, str$exposure2, str$exposure3, str$exposure4, str$exposure5, str$exposure6)
    
    
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
                             init = init, iter = iter, warmup  = warmup, cor = cor_mat[new_order, new_order], Raw = T,
                             n_traits = n_traits, alpha_0 = alpha_0,
                             alpha_1 = alpha_1)
    result_flow[[i]] = result1
    
    
    
    
  }
  
  return(list(result_flow = result_flow, n_snp = n_snp, F_stat = F_stat))
  
}




sel.file <- c("LDL-gera18.csv", "BMI-giant17eu.csv", "SBP-gera17.csv", "BMIchild-egg15.csv")
exp.file <- c("ldl_child.csv",  "bmi_child.csv", "sbp_child.csv",
              "LDL-glgc13.csv",  "BMIadult-ukb.csv", "SBP-ukb.csv")
out.file <- "AS-Malik18EU.csv"
plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"


data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =1e-4,
                               plink_exe = "plink_mac_20210606/plink", cal.cor = T)

filtered.data <-data.list$data;
data = filtered.data
cor_mat = data.list$cor.mat

result = runMVMR7(data = data,  p_value = c(1e-8, 1e-6), cor_mat = data.list$cor.mat,
                  n_traits = c(1,1, 2, 1, 2), alpha_0 = 10, alpha_1 =10, iter = 8000, warmup =3000)

result = runMVMR7(data = data,  p_value = c(1e-4), cor_mat = data.list$cor.mat,
                  n_traits = c(1,1, 2, 1, 2), alpha_0 = 10, alpha_1 =10, iter = 8000, warmup =3000)

data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =1e-3,
                               plink_exe = "plink_mac_20210606/plink", cal.cor = T)

filtered.data <-data.list$data;
data = filtered.data
cor_mat = data.list$cor.mat

result = runMVMR7(data = data,  p_value = c(1e-3), cor_mat = data.list$cor.mat,
                  n_traits = c(1,1, 2, 1, 2), alpha_0 = 10, alpha_1 =10, iter = 8000, warmup =3000)


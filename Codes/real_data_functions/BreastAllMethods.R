load("BreastAllMethods.RData")
source("sim_functions.R")
source("t1d_table.R")


runMVMR3 <- function(data, 
                     p_value =  c(1e-8,  1e-6,  1e-4,  1e-2),
                     init = "Random",iter = 16000, warmup = 8000, bayes =
                       TRUE, cor_mat = NULL, n_traits = c(1,3,3),  alpha_0 = 1.1, alpha_1 = 1.1, indirect_est = FALSE){
  nn <- colnames(data)
  result_grapple = list()
  result_bayes = list()
  result_egger = list()
  result_robust = list()
  result_cml = list()
  result_horse = list()
  
  result = list()
  indirect = list()
  
  time_ivw = list()
  time_grapple = list()
  time_bayes = list()
  time_egger = list()
  time_robust = list()
  time_horse = list()
  time_cml = list()
  for(i in 1:length(p_value)){
    dat = data[which(  data$selection_pvals  < p_value[i]),]
    
    if("SNP" %in% nn){
      F.data <- format_mvmr(BXGs = data[data$selection_pvals < p_value[i], grep("gamma_exp", nn)],
                            BYG = data[data$selection_pvals < p_value[i],"gamma_out1"],
                            seBXGs = data[data$selection_pvals < p_value[i], grep("se_exp", nn)],
                            seBYG = data[data$selection_pvals < p_value[i], "se_out1"],
                            RSID = data[data$selection_pvals < p_value[i], "SNP"])
      
    } else{
      F.data <- format_mvmr(BXGs = data[data$selection_pvals < p_value[i], grep("gamma_exp", nn)],
                            BYG = data[data$selection_pvals < p_value[i],"gamma_out1"],
                            seBXGs = data[data$selection_pvals < p_value[i], grep("se_exp", nn)],
                            seBYG = data[data$selection_pvals < p_value[i], "se_out1"])
      
    }
    t1 = Sys.time()
    
    result[[i]] = ivw_mvmr(r_input = F.data)
    t2 = Sys.time()
    time_ivw[[i]] = t2 - t1
    
    t1 = Sys.time()
    result_grapple[[i]] = 
      GRAPPLE::grappleRobustEst(dat, plot.it = FALSE, cor.mat = cor_mat, diagnosis = F)
    t2 = Sys.time()
    time_grapple[[i]] = t2 - t1
    
    
    bx_matrix = cbind(
      dat$gamma_exp2,
      dat$gamma_exp1)
    
    sebx_matrix = cbind( 
      dat$se_exp2,
      dat$se_exp1)
    
    input = (mr_mvinput(bx = bx_matrix, bxse = sebx_matrix,
                        by = dat$gamma_out1,
                        byse = dat$se_out1))
    
    
    t1 = Sys.time()
    
    result_egger[[i]] = mr_mvegger(input)
    t2 = Sys.time()
    time_egger[[i]] = t2 - t1
    
    t1 = Sys.time()
    result_robust[[i]] = mvmr_robust(bx_matrix, dat$gamma_out1,
                                     dat$se_out1, k.max = 1000, maxit.scale = 1000)
    
    t2 = Sys.time()
    time_robust[[i]] = t2 - t1
    
    
    t1 = Sys.time()
    
    Sig_inv_l = invcov_mvmr(se_bx=sebx_matrix,se_by=dat$se_out1,rho_mat = cor_mat)
    
    
    MVcML_res = MVMRcML::MVmr_cML_DP(bx_matrix, matrix(dat$gamma_out1), sebx_matrix, Sig_inv_l, n = 4207)
    
    MVcML_BIC_SE = MVMRcML::MVcML_SdTheta(bx_matrix, matrix(dat$gamma_out1), Sig_inv_l,
                                          theta=MVcML_res$BIC_theta,
                                          zero_ind = setdiff(1:length(dat$gamma_out1),MVcML_res$BIC_invalid))
    
    MVcMLBIC_pval = pnorm(-abs(MVcML_res$BIC_theta/MVcML_BIC_SE))*2
    result_cml[[i]] = list( MVcML_res, MVcML_BIC_SE, MVcMLBIC_pval )
    
    t2 = Sys.time()
    time_cml[[i]] = t2 - t1
    
    input = data.frame(cbind(cbind( dat$gamma_out1,
                                    dat$gamma_exp2,
                                    dat$gamma_exp1), cbind( dat$se_out1, 
                                                            dat$se_exp2,
                                                            dat$se_exp1)))
    
    colnames(input) = c("betaY",  "betaX2", "betaX1",
                        "betaYse",  "betaX2se", "betaX1se")
    t1 = Sys.time()
    
    result_horse[[i]] = mvmr_horse(input)
    
    t2 = Sys.time()
    time_horse[[i]] = t2 - t1
    if(bayes){
      new_order <- c(3, 1, 2)
      Gamma_hat =rbind(dat$gamma_out1, 
                       dat$gamma_exp2,
                       dat$gamma_exp1)
      
      Sd_hat = rbind(dat$se_out1, 
                     dat$se_exp2,
                     dat$se_exp1)
      t1 = Sys.time()
      result1 = BayesMediation(Gamma_hat, Sd_hat, inv = TRUE,
                               init = init, iter = iter, warmup  = warmup, cor = cor_mat[new_order, new_order], Raw = F,
                               indirect = indirect_est, n_traits = n_traits, alpha_0 = alpha_0,
                               alpha_1 = alpha_1)
      result_bayes[[i]] = result1$summary
      indirect[[i]] = result1$indirect_effect
      
      t2 = Sys.time()
      time_bayes[[i]] = t2 - t1
      
    }
    
    t4 = Sys.time()
  }
  
  if(bayes){
    
    return(list(MVMR = result, Grapple = result_grapple, 
                first = result_bayes, indirect = indirect,
                Egger = result_egger, Robust = result_robust, CML = result_cml, Horse = result_horse,
                time_ivw = time_ivw,
                time_grapple = time_grapple,
                time_bayes = time_bayes,
                time_egger = time_egger ,
                time_robust = time_robust,
                time_horse = time_horse,
                time_cml = time_cml
    ))
  }
  else{
    return( list(MVMR = result, Grapple = result_grapple
    )  )
  }
} 


runMVMR4 <- function(data, 
                     p_value =  c(1e-8,  1e-6,  1e-4,  1e-2),
                     init = "Random",iter = 16000, warmup = 8000, bayes =
                       TRUE, cor_mat = NULL, n_traits = c(1,3,3),  alpha_0 = 1.1, alpha_1 = 1.1, indirect_est = FALSE){
  nn <- colnames(data)
  result_grapple = list()
  result_bayes = list()
  result_egger = list()
  result_robust = list()
  result_cml = list()
  result_horse = list()
  
  result = list()
  indirect = list()
  
  time_ivw = list()
  time_grapple = list()
  time_bayes = list()
  time_egger = list()
  time_robust = list()
  time_horse = list()
  time_cml = list()
  for(i in 1:length(p_value)){
    dat = data[which(  data$selection_pvals  < p_value[i]),]
    
    if("SNP" %in% nn){
      F.data <- format_mvmr(BXGs = data[data$selection_pvals < p_value[i], grep("gamma_exp", nn)],
                            BYG = data[data$selection_pvals < p_value[i],"gamma_out1"],
                            seBXGs = data[data$selection_pvals < p_value[i], grep("se_exp", nn)],
                            seBYG = data[data$selection_pvals < p_value[i], "se_out1"],
                            RSID = data[data$selection_pvals < p_value[i], "SNP"])
      
    } else{
      F.data <- format_mvmr(BXGs = data[data$selection_pvals < p_value[i], grep("gamma_exp", nn)],
                            BYG = data[data$selection_pvals < p_value[i],"gamma_out1"],
                            seBXGs = data[data$selection_pvals < p_value[i], grep("se_exp", nn)],
                            seBYG = data[data$selection_pvals < p_value[i], "se_out1"])
      
    }
    
    result[[i]] = ivw_mvmr(r_input = F.data)
    
    
    result_grapple[[i]] = 
      GRAPPLE::grappleRobustEst(dat, plot.it = FALSE, cor.mat = cor_mat, diagnosis = F)
    
    bx_matrix = cbind(
      dat$gamma_exp3,
      dat$gamma_exp2,
      dat$gamma_exp1)
    
    sebx_matrix = cbind( 
      dat$se_exp3,
      dat$se_exp2,
      dat$se_exp1)
    
    
    input = (mr_mvinput(bx = bx_matrix, bxse = sebx_matrix,
                        by = dat$gamma_out1,
                        byse = dat$se_out1))
    
    
    t1 = Sys.time()
    
    result_egger[[i]] = mr_mvegger(input)
    t2 = Sys.time()
    time_egger[[i]] = t2 - t1
    
    t1 = Sys.time()
    result_robust[[i]] = mvmr_robust(bx_matrix, dat$gamma_out1,
                                     dat$se_out1, k.max = 1000, maxit.scale = 1000)
    
    t2 = Sys.time()
    time_robust[[i]] = t2 - t1
    
    
    t1 = Sys.time()
    
    Sig_inv_l = invcov_mvmr(se_bx=sebx_matrix,se_by=dat$se_out1,rho_mat = cor_mat)
    
    
    MVcML_res = MVMRcML::MVmr_cML_DP(bx_matrix, matrix(dat$gamma_out1), sebx_matrix, Sig_inv_l, n = 4207)
    
    MVcML_BIC_SE = MVMRcML::MVcML_SdTheta(bx_matrix, matrix(dat$gamma_out1), Sig_inv_l,
                                          theta=MVcML_res$BIC_theta,
                                          zero_ind = setdiff(1:length(dat$gamma_out1),MVcML_res$BIC_invalid))
    
    MVcMLBIC_pval = pnorm(-abs(MVcML_res$BIC_theta/MVcML_BIC_SE))*2
    result_cml[[i]] = list( MVcML_res, MVcML_BIC_SE, MVcMLBIC_pval )
    
    t2 = Sys.time()
    time_cml[[i]] = t2 - t1
    
    
    t2 = Sys.time()
    
    input = data.frame(cbind(cbind( dat$gamma_out1, 
                                    dat$gamma_exp3,
                                    dat$gamma_exp2,
                                    dat$gamma_exp1), cbind( dat$se_out1, 
                                                            dat$se_exp3,
                                                            dat$se_exp2,
                                                            dat$se_exp1)))
    
    colnames(input) = c("betaY", "betaX3", "betaX2", "betaX1",
                        "betaYse",  "betaX3se", "betaX2se", "betaX1se")
    t1 = Sys.time()
    
    result_horse[[i]] = mvmr_horse(input)
    
    t2 = Sys.time()
    time_horse[[i]] = t2 - t1
    
    
    if(bayes){
      new_order = c(4,1,2,3)
      Gamma_hat =rbind(dat$gamma_out1, 
                       dat$gamma_exp3,
                       dat$gamma_exp2,
                       dat$gamma_exp1)
      
      Sd_hat = rbind(dat$se_out1, 
                     dat$se_exp3,
                     dat$se_exp2,
                     dat$se_exp1)
      t1 = Sys.time()
      result1 = BayesMediation(Gamma_hat, Sd_hat, inv = TRUE,
                               init = init, iter = iter, warmup  = warmup, cor = cor_mat[new_order, new_order], Raw = F,
                               indirect = indirect_est, n_traits = n_traits, alpha_0 = alpha_0,
                               alpha_1 = alpha_1)
      t2 = Sys.time()
      time_bayes[[i]] = t2 - t1
      result_bayes[[i]] = result1$summary
      indirect[[i]] = result1$indirect_effect
    }
    
    t4 = Sys.time()
  }
  
  if(bayes){
    
    return(list(MVMR = result, Grapple = result_grapple, 
                first = result_bayes, indirect = indirect,
                Egger = result_egger, Robust = result_robust, CML = result_cml, Horse = result_horse,
                time_ivw = time_ivw,
                time_grapple = time_grapple,
                time_bayes = time_bayes,
                time_egger = time_egger ,
                time_robust = time_robust,
                time_horse = time_horse,
                time_cml = time_cml
    ))
  }
  else{
    return( list(MVMR = result, Grapple = result_grapple
    )  )
  }
}


# 
# sel.file <- c("BMI-giant17eu.csv", "bmi_child_agg.csv" )
# exp.file <- c( "childhood_body_size.csv","BMI_adult.csv" )
# out.file <- "Breast-Micha17erp.csv"
# plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"
# 
# data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01,
#                                plink_exe = "plink_mac_20210606/plink", cal.cor = T)
# 
# filtered.data1 <-data.list$data;
# cor_mat1 = data.list$cor.mat
# 
# 
# sel.file <- c("bmi_child_agg.csv", "BMI-giant17eu.csv" )
# exp.file <- c("childhood_bmi_nat_comm_2019_8years", "BMI_adult.csv")
# out.file <- "Breast-Micha17erp.csv"
# plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"
# data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01,
#                                plink_exe = "plink_mac_20210606/plink", cal.cor = T)
# cor_mat2 = data.list$cor.mat
# 
# 
# filtered.data2 <-data.list$data
# 
# 
# sel.file <- c("bmi_child_agg.csv", "BMI-giant17eu.csv" )
# exp.file <- c("childhood_bmi_nat_comm_2019_1year", "childhood_bmi_nat_comm_2019_8years", "BMI_adult.csv")
# out.file <- "Breast-Micha17erp.csv"
# plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"
# data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01,
#                                plink_exe = "plink_mac_20210606/plink", cal.cor = T)
# 
# 
# filtered.data3 <-data.list$data;
# cor_mat3 = data.list$cor.mat

results1 = runMVMR3(data = filtered.data1, cor_mat = cor_mat1, p_value = c(1e-3, 1e-4, 1e-6, 1e-8),
                    n_traits = NULL , alpha_0 = 2, alpha_1 = 2, indirect_est= TRUE, iter = 6000,
                    warmup = 3000)

results2 = runMVMR3(data = filtered.data2, cor_mat = cor_mat2, p_value = c(1e-3, 1e-4, 1e-6, 1e-8),
                    n_traits = NULL , alpha_0 = 2, alpha_1 = 2, indirect_est= TRUE, iter = 6000,
                    warmup = 3000)

results3 = runMVMR4(data = filtered.data3, cor_mat = cor_mat3, p_value = c(1e-3, 1e-4, 1e-6, 1e-8),
                    n_traits = NULL , alpha_0 = 2, alpha_1 = 2, indirect_est= TRUE, iter = 6000,
                    warmup = 3000)
save.image("/home/zixuanwu/mr_mediation/real_result_new.RData")

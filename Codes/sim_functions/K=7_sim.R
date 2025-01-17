setwd("/home/zixuanwu/mr_mediation/mediation_simulation")
load("K=7_dat.RData")
source("sim_functions.R")
source("make_table.R")



# sel.file <- c("LDL-gera18.csv", "HDL-gera18.csv", "TG-gera18.csv")
# exp.file <- c("ldl_child.csv", "hdl_child.csv", "tg_child.csv",
#               "LDL-glgc13.csv", "HDL-glgc13.csv","TG-glgc13.csv")
# out.file <- "CAD-Nelson17.csv"
# plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"
# 
# 
# data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =1e-2,
#                                plink_exe = "plink_mac_20210606/plink", cal.cor = F)
# 
# filtered.data <-data.list$data;
# cor_mat = data.list$cor.mat
# 
# # Generate Pleiotropy
# set.seed(123)
# taus = list()
# 
# out = sample(threshold(filtered.data$gamma_out1, 0.01) ) 
# exps = list()
# for (i in 6:1){
#   name = paste("gamma_exp", toString(i), sep = "")
#   exps[[7 - i]] =  sample(threshold(filtered.data[[name]], 0.01)  )
# }
# 
# taus[[1]] = out
# for(i in 2:7){
#   taus[[i]] = exps[[i - 1]]
# }
# 
# # Get covariance matrix
# Sd_hat = rbind(filtered.data$se_out1,filtered.data$se_exp6,filtered.data$se_exp5, 
#                filtered.data$se_exp4,filtered.data$se_exp3,
#                filtered.data$se_exp2, filtered.data$se_exp1)
# 
# # Set B
# 
# B = matrix(0, 7, 7)
# B[4, 7] = 1
# B[3, 6] = 1
# B[2, 5] = 1
# B[1, 4] = 0.8
# 
# #Generate data
# results = list()


args <- commandArgs(trailingOnly=TRUE)
start <- as.numeric(args[1])


for(i in c(start)){
  set.seed(i)
  dat = generate_fixed_true(sd = Sd_hat, taus = taus, cor_mat = diag(7), B)
  
  Gamma_hat = dat$hat
  Gamma = dat$true
  
  zs1 = Gamma[7,] / Sd_hat[7,]
  zs2 = Gamma[6,] / Sd_hat[6,]
  zs3 = Gamma[5,] / Sd_hat[5,]
  zs4 = Gamma[4,] / Sd_hat[4,]
  zs5 = Gamma[3,] / Sd_hat[3,]
  zs6 = Gamma[2,] / Sd_hat[2,]
  
  ps1 = (1 - 2 * abs(pnorm(zs1) - 0.5))
  ps2 = (1 - 2 * abs(pnorm(zs2) - 0.5))
  ps3 = (1 - 2 * abs(pnorm(zs3) - 0.5))
  ps4 = (1 - 2 * abs(pnorm(zs4) - 0.5))
  ps5 = (1 - 2 * abs(pnorm(zs5) - 0.5))
  ps6 = (1 - 2 * abs(pnorm(zs6) - 0.5))
  
  

  
  df = data.frame(gamma_out1 = Gamma_hat[1,],
                  gamma_exp1 = Gamma_hat[7,],
                  gamma_exp2 = Gamma_hat[6,],
                  gamma_exp3 = Gamma_hat[5,],
                  gamma_exp4 = Gamma_hat[4,],
                  gamma_exp5 = Gamma_hat[3,],
                  gamma_exp6 = Gamma_hat[2,],
                  se_out1 = Sd_hat[1,], 
                  se_exp1 = Sd_hat[7,],
                  se_exp2 = Sd_hat[6,],
                  se_exp3 = Sd_hat[5,],
                  se_exp4 = Sd_hat[4,],
                  se_exp5 = Sd_hat[3,],
                  se_exp6 = Sd_hat[2,],
                  selection_pvals = pmin(ps1, ps2, ps3, ps4, ps5, ps6) * 6)
  
  results[[i]] = runMVMR7(data = df, cor_mat = cor_mat, p_value = c(1e-2, 1e-4, 1e-6,1e-8), 
                          n_traits = c(1,3, 3), alpha_0 = 2, alpha_1 = 2, indirect_est= TRUE, iter = 6000,
                          warmup = 3000)
  print(i)
}

saveRDS(results, paste( "/home/zixuanwu/mr_mediation/results/result", toString(start), ".RData", sep = ""))

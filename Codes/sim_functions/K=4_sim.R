load("K=4_dat.RData")
source("sim_functions.R")
source("make_table.R")


# sel.file <- c("bmi_child_agg.csv", "BMI-giant17eu-F.csv" )
# exp.file <- c("childhood_bmi_nat_comm_2019_1year", "childhood_bmi_nat_comm_2019_8years", "BMI_adult.csv")
# out.file <- "T2D-diagram12-M.csv"
# plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"
# 
# data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01,
#                                plink_exe = "plink_mac_20210606/plink", cal.cor = T)
# 
# filtered.data <-data.list$data;
# 
# N = length(filtered.data$SNP)
# Sd_hat = rbind(filtered.data$se_out1, filtered.data$se_exp3,
#                filtered.data$se_exp2, filtered.data$se_exp1)
# B = rbind(c(0, 1, 0.2, 0.2), c(0, 0, 0.7, 0.2 ), c(0, 0, 0, 0.7), c(0, 0, 0, 0))
# 
# set.seed(123)
# out = sample(threshold(filtered.data$gamma_out1, 0.01) )
# exp3 = sample(threshold(filtered.data$gamma_exp3, 0.01)  )
# exp2 = sample(threshold(filtered.data$gamma_exp2, 0.01 ) )
# exp1 =sample(threshold(filtered.data$gamma_exp1, 0.01)  )
# taus = list(out, exp3, exp2, exp1)
# results = list()


args <- commandArgs(trailingOnly=TRUE)
start <- as.numeric(args[1])


for(i in 1:5){
  set.seed(5 * (start-1) + i)
  dat = generate_fixed_true(sd = Sd_hat, taus = taus, cor_mat = diag(4), B)
  
  
  Gamma_hat = dat$hat
  Gamma = dat$true
  
  zs1 = Gamma[4,] / Sd_hat[4,]
  zs2 = Gamma[3,] / Sd_hat[3,]
  zs3 = Gamma[2,] / Sd_hat[2,]
  
  
  ps1 = (1 - 2 * abs(pnorm(zs1) - 0.5))
  ps2 = (1 - 2 * abs(pnorm(zs2) - 0.5))
  ps3 = (1 - 2 * abs(pnorm(zs3) - 0.5))
  
  df = data.frame(gamma_out1 = Gamma_hat[1,],
                  gamma_exp1 = Gamma_hat[4,],
                  gamma_exp2 = Gamma_hat[3,],
                  gamma_exp3 = Gamma_hat[2,],
                  se_out1 = Sd_hat[1,], 
                  se_exp1 = Sd_hat[4,],
                  se_exp2 = Sd_hat[3,],
                  se_exp3 = Sd_hat[2,],
                  selection_pvals = 3 * pmin(ps1, ps2, ps3))
  
  results[[i]] = runMVMR4(data = df, cor_mat = diag(4), p_value = c(1e-2, 1e-4, 1e-6, 1e-8),
                          n_traits = NULL , alpha_0 = 2, alpha_1 = 2, indirect_est= TRUE, iter = 6000,
                          warmup = 3000)
  print(i)
}

saveRDS(results, paste( "/home/zixuanwu/mr_mediation/results/result", toString(start), ".RData", sep = ""))

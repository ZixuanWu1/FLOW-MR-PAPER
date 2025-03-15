load("K=3_dat.RData")
source("sim_functions.R")
source("make_table.R")

# ## Simulation
# sel.file <- c("BMIchild-egg15.csv", "BMI-giant17eu-F.csv" )
# exp.file <- c("BMI8year_moba19", "BMIadult-ukb.csv")
# out.file <- "T2D-diagram12-M.csv"
# plink_refdat <- "data_maf0.01_rs_ref/data_maf0.01_rs_ref"
# 
# data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01,
#                                plink_exe = "plink_mac_20210606/plink", cal.cor = F)
# 
# filtered.data <-data.list$data;
# 
# # Fix pleiotropy. Only change the noise.
# N = length(filtered.data$SNP)
# Sd_hat = rbind(filtered.data$se_out1, filtered.data$se_exp2, filtered.data$se_exp1)
# B = matrix(c(1,0,0,1,1,0, .3, 0.7, 1), nrow = 3) - diag(3)
# 
# 
# 
# set.seed(1234)
# out = sample(threshold(filtered.data$gamma_out1, 0.01) )
# exp2 = sample(threshold(filtered.data$gamma_exp2, 0.01 ) )
# exp1 =sample(threshold(filtered.data$gamma_exp1, 0.01)  )
# taus = list(out, exp2, exp1)
# 
# results = list()

args <- commandArgs(trailingOnly=TRUE)
start <- as.numeric(args[1])


for(i in 1:10){
  set.seed(10 * (start-1) + i)
  dat = generate_fixed_true(sd = Sd_hat, taus = taus, cor_mat = diag(3), B)
  
  
  Gamma_hat = dat$hat
  Gamma = dat$true
  
  zs1 = Gamma[3,] / Sd_hat[3,]
  zs2 = Gamma[2,] / Sd_hat[2,]
  
  ps1 = (1 - 2 * abs(pnorm(zs1) - 0.5))
  ps2 = (1 - 2 * abs(pnorm(zs2) - 0.5))
  
  df = data.frame(gamma_out1 = Gamma_hat[1,],
                  gamma_exp1 = Gamma_hat[3,],
                  gamma_exp2 = Gamma_hat[2,],
                  se_out1 = Sd_hat[1,], 
                  se_exp1 = Sd_hat[3,],
                  se_exp2 = Sd_hat[2,],
                  selection_pvals = 2 * pmin(ps1, ps2))
  
  results[[i]] = runMVMR3(data = df, cor_mat = diag(3), p_value = c(1e-2, 1e-4, 1e-6, 1e-8),
                         n_traits = NULL , alpha_0 = 2, alpha_1 = 2, indirect_est= TRUE, iter = 6000,
                          warmup = 3000)
print(i)
}

saveRDS(results, paste( "/home/zixuanwu/mr_mediation/results/result", toString(start), ".RData", sep = ""))

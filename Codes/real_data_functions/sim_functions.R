library(doParallel)
library(MrMediation)
source("gibbs_vector.R")
library(MVMR)
library(GRAPPLE)
library(MASS)
library(MendelianRandomization)
library(MVMRcML)
source("Other_methods.R")

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
    
    Sig_inv_l = invcov_mvmr(se_bx=sebx_matrix,se_by=dat$se_out1,rho_mat = diag(3))
    

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
      Gamma_hat =rbind(dat$gamma_out1, 
                       dat$gamma_exp2,
                       dat$gamma_exp1)
      
      Sd_hat = rbind(dat$se_out1, 
                     dat$se_exp2,
                     dat$se_exp1)
      t1 = Sys.time()
      result1 = BayesMediation(Gamma_hat, Sd_hat, inv = TRUE,
                               init = init, iter = iter, warmup  = warmup, cor = cor_mat, Raw = F,
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
    
    Sig_inv_l = invcov_mvmr(se_bx=sebx_matrix,se_by=dat$se_out1,rho_mat = diag(4))
    
    
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
                               init = init, iter = iter, warmup  = warmup, cor = cor_mat, Raw = F,
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


runMVMR7 <- function(data, 
                     p_value =  c(1e-8,  1e-6,  1e-4,  1e-2),
                     init = "Random",iter = 16000, warmup = 8000, bayes =
                       TRUE, cor_mat = NULL, n_traits = c(1,3,3),  alpha_0 = 1.1, alpha_1 = 1.1, indirect_est = FALSE){
  nn <- colnames(data)
  result_grapple = list()
  result_bayes = list()
  result_egger = list()
  result_robust = list()
  result_horse = list()
  
  result = list()
  indirect = list()
  
  time_ivw = list()
  time_grapple = list()
  time_bayes = list()
  time_egger = list()
  time_robust = list()
  time_horse = list()

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
    
    bx_matrix = cbind( dat$gamma_exp6,
                       dat$gamma_exp5,
                       dat$gamma_exp4,
                       dat$gamma_exp3,
                       dat$gamma_exp2,
                       dat$gamma_exp1)
    
    sebx_matrix = cbind( dat$se_exp6,
                         dat$se_exp5,
                         dat$se_exp4,
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
    
    
    
    
    
    
    input = data.frame(cbind(cbind( dat$gamma_out1, dat$gamma_exp6,
                                    dat$gamma_exp5,
                                    dat$gamma_exp4,
                                    dat$gamma_exp3,
                                    dat$gamma_exp2,
                                    dat$gamma_exp1), cbind( dat$se_out1, dat$se_exp6,
                                                            dat$se_exp5,
                                                            dat$se_exp4,
                                                            dat$se_exp3,
                                                            dat$se_exp2,
                                                            dat$se_exp1)))
    
    colnames(input) = c("betaY", "betaX6", "betaX5", "betaX4", "betaX3", "betaX2", "betaX1",
                        "betaYse", "betaX6se", "betaX5se", "betaX4se", "betaX3se", "betaX2se", "betaX1se")
    t1 = Sys.time()
    
    result_horse[[i]] = mvmr_horse(input)
    
    t2 = Sys.time()
    time_horse[[i]] = t2 - t1
    
    
    
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
      t1 = Sys.time()
      result1 = BayesMediation(Gamma_hat, Sd_hat, inv = TRUE,
                               init = init, iter = iter, warmup  = warmup, cor = cor_mat, Raw = F,
                               indirect = indirect_est, n_traits = n_traits, alpha_0 = alpha_0,
                               alpha_1 = alpha_1)
      t2 = Sys.time()
      time_bayes[[i]] = t2 - t1
      
      result_bayes[[i]] = result1$summary
      indirect[[i]] = result1$indirect_effect
    }
    
  }
  
  if(bayes){
    
    return(list(MVMR = result, Grapple = result_grapple, 
                first = result_bayes, indirect = indirect,
                Egger = result_egger, Robust = result_robust, Horse = result_horse,
                time_ivw = time_ivw,
                time_grapple = time_grapple,
                time_bayes = time_bayes,
                time_egger = time_egger ,
                time_robust = time_robust,
                time_horse = time_horse
    ))
  }
  else{
    return( list(MVMR = result, Grapple = result_grapple
    )  )
  }
}
threshold <- function(X, epsilon){
  newX = abs(X) - epsilon
  newX[newX < 0] = 0
  return(sign(X) * newX )
}

generate_fixed_true <- function( sd,  taus, cor_mat, B, inv = FALSE){
  N = dim(sd)[2]
  K = dim(sd)[1]
  if(inv == FALSE){
    B = solve(diag(K) - B) - diag(K)
  }
  #Generate the he pleiotropy effects
  alpha = matrix(0, K, N)
  for(i in 1:K){
    alpha[i, ] = taus[[i]]
  }
  
  #Compute Gamma
  Gamma = (diag(K) + B)%*% alpha
  #Add noise to Gamma
  Gamma_hat = Gamma
  
  for(j in 1:N){
    Gamma_hat[, j] = Gamma_hat[, j] + mvrnorm(n = 1, mu = rep(0, K), Sigma = diag(sd[, j]) %*% cor_mat %*% diag(sd[, j]))
  }
  
  
  return(list(true = Gamma, hat = Gamma_hat))
}


dIVW <- function(Gamma_hat, Sd_hat, cor_mat, exp_names = NULL){
  K = dim(Gamma_hat)[1]
  P = dim(Gamma_hat)[2]
  
  W = diag(1/Sd_hat[1,]^2) / sum(1/Sd_hat[1,]^2)
  
  Gamma_hat = Gamma_hat[K:1, ]
  Sd_hat = Sd_hat[K:1, ]
  cor_mat = cor_mat[K:1, K:1]
  
  result = Gamma_hat %*% W %*% t(Gamma_hat)
  
  for(i in 1:P){
    result = result - diag(Sd_hat[, i]) %*% cor_mat %*% diag(Sd_hat[, i]) * diag(W)[i]
  }
  
  L = as.matrix(expand(lu(result))$L)
  
  for(k in 1:K){
    L[, k] = L[, k] / L[k, k]
  }
  L = diag(K) - solve(L)
  output = L[K:1, K:1]
  
  
  
  Ls = array(dim = c(K, K, 10000))
  
  for(j in 1:10000){
    
    success <- FALSE
    
    while (!success) {
      result <- tryCatch({
        # Code that might produce an error
        # Example:
        samples = sample(1:P, P, TRUE)
        Gamma_sub = Gamma_hat[, samples]
        W_sub = W[samples, samples]
        W_sub = W_sub / sum(diag(W_sub))
        
        
        result = Gamma_sub %*% W_sub %*% t(Gamma_sub)
        
        for(i in 1:P){
          result = result - diag(Sd_hat[, samples[i]]) %*% cor_mat %*% diag(Sd_hat[, samples[i]]) * diag(W_sub)[i]
        }
        
        L = as.matrix(expand(lu(result))$L)
        
        
        for(k in 1:K){
          L[, k] = L[, k] / L[k, k]
        }
        L = diag(K) - solve(L)
        Ls[,,j] = L[K:1, K:1]
        
        # If no error, set success to TRUE
        success <- TRUE
      }, error = function(e) {
        NULL
      })
    }
    if (j %% 1000 == 0){
      print(j)
      
    }
    
    
    
    
    
  }
  lower = matrix(0,K, K)
  upper = matrix(0,K, K)
  median = matrix(0,K, K)
  
  for(i in 1:K){
    for(j in 1:K){
      lower[i, j] = quantile(Ls[i, j,], 0.025)
      upper[i, j] = quantile(Ls[i, j,], 0.975)
      
      median[i, j] = quantile(Ls[i, j,], 0.5)
    }
  }
  
  lower = round(lower, 2)
  upper = round(upper, 2)
  median = round(median, 2)
  sig = which(upper < 0 | lower > 0, arr.ind = T)
  if (! is.null(exp_names)){
    for (i in 1:nrow(sig)){
      for(j in 1:2){
        sig[i, j] = exp_names[K + 1 - as.numeric(sig[i, j])]
      }
    }
  }
  sig = sig[, 2:1]
  return(list(output = output, lower= lower, upper = upper, median = median, sig_ind = sig))
}

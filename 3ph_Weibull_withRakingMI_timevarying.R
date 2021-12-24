## Weibull

library(survey)
library(mvtnorm)
library(sandwich)
library(lubridate)
library(tidyverse)
#source('func_3phaseC.R')
#source('all_funcs.R')

######################################
args   = commandArgs(TRUE)
TASKID = as.numeric(args[2])
njob   = TASKID
## Using data_emrb, depois de quebrar em mes, retornar p/ escala original

set.seed(2108 + njob)

N <- n1 <- 15000
n2 <- 1250
n3 <- 250
Nsim <- 100
beta <- c(.5,.5); beta0 <- -3
Nimpute <- 3
setting <- 'setting 1'

beta0_chart <- c(-5, 4.75)
beta0_emr <- c(-3.5, 3)
sx_emr   <- 1
sx_chart <- .01

coefRak <- coefRak_3p <- varsRak <- varsRak_3p <- varsRakb_3p <- coefIPW <- coefIPWb <- varsIPW <- varsIPWb <- matrix(NA, nrow=Nsim, ncol=length(beta)+1)
coefMI <- coefMI_2p <- varsMI <- varsMI_2p <- coefCox <- coefglm <- coefglm2 <- coefglm2a <- varsCox <- varsglm <- matrix(NA, nrow=Nsim, ncol=length(beta)+1)
coefRak_3pMI <- coefRak_2phMI <- varsRak_3pMI <- varsRak_2phMI <- matrix(NA, nrow=Nsim, ncol=length(beta)+1)
RubinSE <- RubinSE_2p <- Estimate <- Estimate_2p <- AnalysisN <- AnalysisN_2p <- NULL
RobinsWangSE <- RobinsWangSE_2p <- RobinsWangSEd <- matrix(NA, nrow=Nsim, ncol=length(beta)+1)

file_est <- paste0('Coefs_unbiased_beta_0_chart', beta0_chart[1], '_sx_chart_', sx_chart[1], '_', N, n1, n2,'.txt')
file_var <- paste0('Vars_unbiased_beta_0_chart', beta0_chart[1], '_sx_chart_', sx_chart[1], '_', N, n1, n2,'.txt')

gamma_x1 <- FALSE
if (gamma_x1) {
  file_est <- paste0('Coefs_unbiased_gammaX1_beta_0_chart', beta0_chart[1], '_sx_chart_', sx_chart[1], '_', N, n1, n2,'.txt')
  file_var <- paste0('Vars_unbiased_gammaX1_beta_0_chart', beta0_chart[1], '_sx_chart_', sx_chart[1], '_', N, n1, n2,'.txt')
}

nr <- 3
n_months <- 6
sum_ifs_2ph_MI_mat <- matrix(NA, ncol=length(beta)+2, nrow=Nsim)
run_3phase <- TRUE
for (i in 1:Nsim){
  
  ## Data generate
  data_true <- simulWeib(N=N, lambda=0.0015, rho=1, beta=beta, rateC=0.2, nr=6, maxT=n_months,
                         beta0_emr=beta0_emr, beta0_chart=beta0_chart, gamma_x1=gamma_x1)
  data_emr <- data_chart <- data_phone <- data_true
  fit_cox <- coxph(Surv(time, y) ~ X1 + X2, data=data_true)
  fit <- glm(y ~ X1 + X2, family=binomial, data=data_true)
  
  coefCox[i,] <- c(NA, coef(fit_cox))
  varsCox[i,] <- c(NA, summary(fit_cox)$coef[,2]^2)
  coefglm[i,] <- coef(fit)
  varsglm[i,] <- summary(fit)$coef[,2]^2

#  data_samp <- fun_samp(data_emr, data_chart, data_phone, n1, n2, biased=FALSE)
#  data_chart <- data_samp$data_chart
#  data_phone <- data_samp$data_phone
#  data_extend <- fun_extend(data_emr, data_chart, data_phone)

  data_emr <- data_true %>%
    mutate(uid = paste0(id, X2)) %>%
    group_by(uid) %>%
    mutate(y_emr_cum = cumsum(y_emr)) %>%
    filter(y_emr == y_emr_cum) %>%
    ungroup() %>%
    mutate(uid = paste0(id, X1, X2))
  
  samp1 <- aggregate(data_emr$y_emr, list(data_emr$id), max)
  R1.1 <- sample(samp1$Group.1[samp1$x == 1], n2)
  R1.0 <- sample(samp1$Group.1[samp1$x == 0], n2)
  R1 <- c(R1.1, R1.0)
  prob1.1 <- n2/length(samp1$Group.1[samp1$x == 1])
  prob1.0 <- n2/length(samp1$Group.1[samp1$x == 0])
  prob1 <- data.frame(id = samp1$Group.1, probs = NA)
  prob1[prob1$id %in% samp1$Group.1[samp1$x == 1],]$probs <- prob1.1
  prob1[prob1$id %in% samp1$Group.1[samp1$x == 0],]$probs <- prob1.0

  data_emr <- data_emr %>%
    mutate(R1 = ifelse(id %in% R1, 1, 0)) %>%
    left_join(prob1, by='id')
  
  data_chartb <- data_true %>%
    mutate(uid = paste0(id, X2)) %>%
    group_by(uid) %>%
    mutate(y_chart_cum = cumsum(y_chart)) %>%
    filter(y_chart == y_chart_cum) %>%
    ungroup() %>%
    mutate(uid = paste0(id, X1, X2)) %>%
    filter(id %in% R1) %>%
    dplyr::select(uid, id, y_chart, X1_chart)

  if (setting == 'setting 1'){
    samp2 <- samp1[samp1$Group.1 %in% R1,]
    R2.1 <- sample(samp2$Group.1[samp2$x == 1], n3)
    R2.0 <- sample(samp2$Group.1[samp2$x == 0], n3)
    R2 <- c(R2.1, R2.0)
    prob2.1 <- n3/length(samp2$Group.1[samp2$x == 1])
    prob2.0 <- n3/length(samp2$Group.1[samp2$x == 0])
    prob2 <- data.frame(id = samp2$Group.1, probs2a = NA)
    prob2[prob2$id %in% samp2$Group.1[samp2$x == 1],]$probs2a <- prob2.1
    prob2[prob2$id %in% samp2$Group.1[samp2$x == 0],]$probs2a <- prob2.0
  }
  if (setting == 'setting 2') {
    samp2b <- aggregate(data_chartb$y_chart, list(data_chartb$id), max)
    samp2b <- samp2b[samp2b$Group.1 %in% R1,]
    R2.1 <- sample(samp2$Group.1[samp2$x == 1], n3)
    R2.0 <- sample(samp2$Group.1[samp2$x == 0], n3)
    R2 <- c(R2.1, R2.0)
  }
  
  data_chartb <- data_chartb %>%
    mutate(R2 = ifelse(id %in% R2, 1, 0)) %>%
    dplyr::select(-id)
    
  data_phoneb <- data_true %>%
    filter(id %in% R2) %>%
    mutate(uid = paste0(id, X2)) %>%
    group_by(uid) %>%
    mutate(y_cum = cumsum(y)) %>%
    filter(y == y_cum) %>%
    ungroup() %>%
    mutate(uid = paste0(id, X1, X2)) %>%
    left_join(prob2, by='id') %>%
    dplyr::select(uid, y, X1, probs2a)


    
  data_emr2b <- data_emr %>%
    dplyr::select(-y_chart, -y, -X1_chart, -X1) %>%
    left_join(data_chartb, by='uid') %>%
    left_join(data_phoneb, by='uid') %>%
    mutate(R1 = ifelse(is.na(y_chart), 0, R1),
           #R2 = ifelse(is.na(y) & R2==1, NA, R2),
           probs2 = ifelse(is.na(y), NA, probs*probs2a)) #%>%
  data_emr2b[is.na(data_emr2b$y),]$R2 <- NA
    #filter(!(R1==1 & is.na(y_chart))) %>%
    #filter(!(R2==1 & is.na(y)) | is.na(R2))
  
  fit <- glm(y ~ X1 + X2, family=binomial, data=data_emr2b)
  coefglm2[i,] <- coef(fit)
  fit <- glm(y ~ X1 + X2, family=binomial, data=data_emr)
  coefglm2a[i,] <- coef(fit)
  
    
  dd_emr <- data_emr2b
  ifs3ph_step1 <- ifs3ph_step2 <- ifs2ph_step1 <- list()
  
  res_rak_ipw <- fun_raking(data_emr2b=data_emr2b)
  (coefRak[i,] <- coef(res_rak_ipw$resCAL_2ph))
  (coefIPW[i,] <- coef(res_rak_ipw$resIPW_2ph))
  
  #(varsRak[i,] <- summary(res_rak_ipw$resCAL_2ph)$coef[,2]^2)
  #(varsIPW[i,] <- summary(res_rak_ipw$resIPW_2ph)$coef[,2]^2)
  
  
  res_rak3ph <- fun_raking_3ph(data_emr2b=dd_emr,
                               totals1=res_rak_ipw$totals1,
                               totals2=res_rak_ipw$totals1,
                               datcsv1=res_rak_ipw$datcsv1,
                               datcsv2=res_rak_ipw$datcsv2)
  
  
  (coefRak_3p[i,] <- summary(res_rak3ph$resIPW_2ph)$coef[,1])
  
  
  newCoef <- newVar <- matrix(NA, nrow=Nimpute, ncol=length(beta)+1)
  newCoef_2p <- newVar_2p <- matrix(NA, nrow=Nimpute, ncol=length(beta)+1)

#  print(i)
#}

  dd_emr <- dd_emr %>%
    group_by(id) %>%
    mutate(X1_chart_M0 = dplyr::lag(X1_chart),
           X1_chart_M1 = dplyr::lead(X1_chart),
           X1_emr_M0 = dplyr::lag(X1_emr),
           X1_emr_M1 = dplyr::lead(X1_emr),
           X1_M0 = dplyr::lag(X1),
           X1_M1 = dplyr::lead(X1),
           y_M0 = dplyr::lag(y),
           y_M1 = dplyr::lead(y),
           y_chart_M0 = dplyr::lag(y_chart),
           y_chart_M1 = dplyr::lead(y_chart),
           y_emr_M0 = dplyr::lag(y_emr),
           y_emr_M1 = dplyr::lead(y_emr),
           y_M0 = ifelse(is.na(y_M0) & !is.na(y), y, y_M0),
           y_M1 = ifelse(is.na(y_M1) & !is.na(y), y, y_M1),
           y_emr_M0 = ifelse(is.na(y_emr_M0) & !is.na(y_emr), y_emr, y_emr_M0),
           y_emr_M1 = ifelse(is.na(y_emr_M1) & !is.na(y_emr), y_emr, y_emr_M1),
           y_chart_M0 = ifelse(is.na(y_chart_M0) & !is.na(y_chart), y_chart, y_chart_M0),
           y_chart_M1 = ifelse(is.na(y_chart_M1) & !is.na(y_chart), y_chart, y_chart_M1),
           X1_emr_M0 = ifelse(is.na(X1_emr_M0) & !is.na(X1_emr), X1_emr, X1_emr_M0),
           X1_emr_M1 = ifelse(is.na(X1_emr_M1) & !is.na(X1_emr), X1_emr, X1_emr_M1),
           X1_chart_M0 = ifelse(is.na(X1_chart_M0) & !is.na(X1_chart), X1_chart, X1_chart_M0),
           X1_chart_M1 = ifelse(is.na(X1_chart_M1) & !is.na(X1_chart), X1_chart, X1_chart_M1),
           X1_emr_M0 = ifelse(is.na(X1_emr_M0) & !is.na(X1_emr), X1_emr, X1_emr_M0),
           X1_emr_M1 = ifelse(is.na(X1_emr_M1) & !is.na(X1_emr), X1_emr, X1_emr_M1),
           X1_M0 = ifelse(is.na(X1_M0) & !is.na(X1), X1, X1_M0),
           X1_M1 = ifelse(is.na(X1_M1) & !is.na(X1), X1, X1_M1)) %>%
    ungroup()
  
  for (k in 1:Nimpute){
    if (k == 10) {
      cat('\n Starting Imputation, step = ', k)
      cat('\n 3-phase MI')
    }
    
    if (run_3phase){
      dd_emr_3ph <- func_imp3ph_YX(dd_emr)
#      dd_emr_3phb <- dd_emr_3ph %>% mutate(uidA = paste0(id, X1_emr, X2)) %>% group_by(uidA) %>%
#        mutate(ycum_3ph  = cumsum(y_step2),
#               ycum2_3ph = cumsum(ycum_3ph)) %>% filter(ycum2_3ph <= 1)
      (newMod_3ph <- glm(y_step2 ~ X1_step2 + X2, data=dd_emr_3ph, family=binomial))
      newCoef[k,] <- newMod_3ph$coeff
      newVar[k,]  <- diag(vcov(newMod_3ph))
      
      dd_emr_3ph1 <- dd_emr_3ph#func_imp3ph_YX_raking(dd_emr)
      
 #     dd_emr_3phb1 <- dd_emr_3ph1 %>% filter(R1==1) %>%
#        mutate(uidA = paste0(id, X1_emr, X2)) %>% group_by(uidA) %>%
#        mutate(ycum_3ph  = cumsum(y_step1_imputed),
#               ycum2_3ph = cumsum(ycum_3ph)) %>% filter(ycum2_3ph <= 1)

      #      X1_temp <- aggregate(dd_emr_3phb1$X1_step1_imputed, list(dd_emr_3phb1$id), mean)
      #      colnames(X1_temp) <- c('id', 'X1_new_imputed_avg')
      #      dd_emr_3phb1 <- dd_emr_3phb1 %>%
      #        group_by(id) %>%
      #        arrange(desc(y_step1_imputed)) %>%
      #        slice(1) %>%
      #        ungroup() %>%
      #        left_join(X1_temp, by='id') %>%
      #        as.data.frame()
      
      fit_3ph_step1 <- glm(y_step1_imputed ~ X1_step1_imputed + X2,
                           family=binomial, data=dd_emr_3ph1[!is.na(dd_emr_3ph1$y_step1_imputed),])
      fit_3ph_step1b <- influence_func(fit_3ph_step1)
      fit_3ph_step1b <- aggregate(fit_3ph_step1b, list(dd_emr_3ph1$id[!is.na(dd_emr_3ph1$y_step1_imputed)]), sum)
      colnames(fit_3ph_step1b) <- c('id', paste0('if', 1:(ncol(fit_3ph_step1b)-1)))
      ifs3ph_step1[[k]]  <- fit_3ph_step1b

#      dd_emr_3phb2 <- dd_emr_3ph1 %>%
#        mutate(uidA = paste0(id, X1_emr, X2)) %>%
#        group_by(uidA) %>%
#        mutate(ycum_3ph  = cumsum(y_step2_imputed),
#               ycum2_3ph = cumsum(ycum_3ph)) %>% filter(ycum2_3ph <= 1)
 
      #      X1_temp <- aggregate(dd_emr_3phb2$X1_step2_imputed, list(dd_emr_3phb2$id), mean)
      #      colnames(X1_temp) <- c('id', 'X1_new_imputed_avg')
      #      dd_emr_3phb2 <- dd_emr_3phb2 %>%
      #        group_by(id) %>%
      #        arrange(desc(y_step2_imputed)) %>%
      #        slice(1) %>%
      #        ungroup() %>%
      #        left_join(X1_temp, by='id') %>%
      #        as.data.frame()
      
      fit_3ph_step2 <- glm(y_step2_imputed ~ X1_step2_imputed + X2,
                           family=binomial, data=dd_emr_3phb2)
      fit_3ph_step2b <- influence_func(fit_3ph_step2)
      fit_3ph_step2b <- aggregate(fit_3ph_step2b, list(dd_emr_3phb2$id), sum)
      colnames(fit_3ph_step2b) <- c('id', paste0('if', 1:(ncol(fit_3ph_step2b)-1)))
      ifs3ph_step2[[k]]  <- fit_3ph_step2b
      
      dd_emr_3ph <- dd_emr_3phb
    }
    
    ## 2-phase     
    dd_emr_2ph <- func_imp2ph_YX(dd_emr)
    #dd_emr_2phb <- dd_emr_2ph %>% mutate(uidA = paste0(id, X1_emr, X2)) %>% group_by(uidA) %>%
    #  mutate(ycum_2ph  = cumsum(y_new),
    #         ycum2_2ph = cumsum(ycum_2ph)) %>% filter(ycum2_2ph <= 1)
    (newMod_2ph <- glm(y_new ~ X1_new + X2, data=dd_emr_2ph, family=binomial))
    newCoef_2p[k,] <- newMod_2ph$coeff
    newVar_2p[k,]  <- diag(vcov(newMod_2ph))
    
    dd_emr_2ph1 <- dd_emr_2ph#func_imp2ph_YX_raking(dd_emr)   
#    dd_emr_2ph1 <- dd_emr_2ph1 %>% mutate(uidA = paste0(id, X1_emr, X2)) %>% group_by(uidA) %>%
#      mutate(ycum_2ph  = cumsum(y_new_imputed),
#             ycum2_2ph = cumsum(ycum_2ph)) %>% filter(ycum2_2ph <= 1)
    
    #    X1_temp <- aggregate(dd_emr_2ph1$X1_new_imputed, list(dd_emr_2ph1$id), mean)
    #    colnames(X1_temp) <- c('id', 'X1_new_imputed_avg')
    #    dd_emr_2ph1 <- dd_emr_2ph1 %>%
    #      group_by(id) %>%
    #      arrange(desc(y_new_imputed)) %>%
    #      slice(1) %>%
    #      ungroup() %>%
    #      left_join(X1_temp, by='id') %>%
    #      as.data.frame()
    
    fit_2ph_step1  <- glm(y_new_imputed ~ X1_new_imputed + X2, family=binomial, data=dd_emr_2ph1)
    fit_2ph_step1b <- influence_func(fit_2ph_step1)
    fit_2ph_step1b <- aggregate(fit_2ph_step1b, list(dd_emr_2ph1$id), sum)
    colnames(fit_2ph_step1b) <- c('id', paste0('if', 1:(ncol(fit_2ph_step1b)-1)))
    ifs2ph_step1[[k]]  <- fit_2ph_step1b
  }
  
  coefMI_2p[i,] <- colMeans(newCoef_2p)
  varsMI_2p[i,] <- colMeans(newVar_2p)+(Nimpute+1)*apply(newCoef_2p, 2, var)/Nimpute
  
  res_rak_2phMI <- fun_rak2ph_MI(dd_emr, ifs2ph_step1)
  (coefRak_2phMI[i,] <- coef(res_rak_2phMI))
  (varsRak_2phMI[i,] <- summary(res_rak_2phMI)$coef[,2]^2)
  
  if (run_3phase){
    coefMI[i,] <- colMeans(newCoef)
    varsMI[i,] <- colMeans(newVar)+(Nimpute+1)*apply(newCoef, 2, var)/Nimpute
    
    res_rak_3phMI <- fun_rak3ph_MI(dd_emr, ifs3ph_step1, ifs3ph_step2)
    (coefRak_3pMI[i,] <- coef(res_rak_3phMI))
  }
  
  
  
  
  if (i > 1){
    cat('\n step = ',i ,'\n')
    cat('\n Coef Cox = ',   round(colMeans(coefCox[1:i,]), 3))
    cat('\n EmpVar Cox = ', round(apply(coefCox[1:i,], 2, var), 3))
    cat('\n EVar Cox = ',   round(colMeans(varsCox[1:i,]), 4))
    
    cat('\n step = ',i ,'\n')
    cat('\n Coef Cox = ',   round(colMeans(coefglm2[1:i,]), 3))
    cat('\n EmpVar Cox = ', round(apply(coefglm2[1:i,], 2, var), 3))
    
    cat('\n step = ',i ,'\n')
    cat('\n MI 3-phase = ', round(colMeans(coefMI[1:i,]), 3))
    cat('\n EmpVar MI = ', round(apply(coefMI[1:i,], 2, var), 4))
    cat('\n EVar MI = ',   round(colMeans(RobinsWangSE[1:i,]), 4))
    
    cat('\n step = ',i ,'\n')
    cat('\n MI 2-phase = ', round(colMeans(coefMI_2p[1:i,]), 3))
    cat('\n EmpVar MI = ', round(apply(coefMI_2p[1:i,], 2, var), 4))
    cat('\n EVar MI = ',   round(colMeans(RobinsWangSE_2p[1:i,]), 4))
    
    cat('\n step = ',i ,'\n')
    cat('\n Raking = ', round(colMeans(coefRak[1:i,]), 3))
    cat('\n EmpVar Rak = ', round(apply(coefRak[1:i,], 2, var), 3))
    cat('\n EVar Rak = ',   round(colMeans(varsRak[1:i,]), 4))
    
    cat('\n step = ',i ,'\n')
    cat('\n IPW = ', round(colMeans(coefIPW[1:i,]), 3))
    cat('\n EmpVar IPW = ', round(apply(coefIPW[1:i,], 2, var), 3))
    cat('\n EVar IPW = ',   round(colMeans(varsIPW[1:i,]), 4))
    
    cat('\n step = ',i ,'\n')
    cat('\n Raking 3-phase = ', round(colMeans(coefRak_3p[1:i,]), 3))
    cat('\n EmpVar Rak = ', round(apply(coefRak_3p[1:i,], 2, var), 3))
    cat('\n EVar Rak = ',   round(colMeans(varsRak_3p[1:i,]), 4))
    
    cat('\n step = ',i ,'\n')
    cat('\n Raking 2-phase-MI = ', round(colMeans(coefRak_2phMI[1:i,]), 3))
    cat('\n EmpVar Rak = ', round(apply(coefRak_2phMI[1:i,], 2, var), 3))
    cat('\n EVar Rak = ',   round(colMeans(varsRak_2phMI[1:i,]), 4))
    
    cat('\n step = ',i ,'\n')
    cat('\n Raking 3-phase-MI = ', round(colMeans(coefRak_3pMI[1:i,]), 3))
    cat('\n EmpVar Rak = ', round(apply(coefRak_3pMI[1:i,], 2, var), 3))
    cat('\n EVar Rak = ',   round(colMeans(varsRak_3pMI[1:i,]), 4))
    
  }
  
  ddest <- matrix(c(coefCox[i,], coefIPW[i,], coefRak[i,], coefRak_3p[i,],
                    coefRak_2phMI[i,], coefRak_3pMI[i,],
                    coefMI_2p[i,], coefMI_2p[i,], coefMI[i,], coefMI[i,]), nrow=1)
  ddvar <- matrix(c(varsCox[i,], varsIPW[i,], varsRak[i,], varsRak_3p[i,],
                    varsRak_2phMI[i,], varsRak_3pMI[i,],
                    RobinsWangSE_2p[i,], varsMI_2p[i,], RobinsWangSE[i,], varsMI[i,]), nrow=1)  
#  if (i == 1) {
#    save2 <- paste('N = ', N, 'n1 = ', n1, 'n2 = ', n2, 'beta = ', paste(beta, collapse=" "), ' njob = ', njob)
#    if (!file.exists(file_est)) {
#      file.create(file_est)
#      write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
#      file.create(file_var)
#      write.table(save2, file = file_var,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
#    }
#  }
#  
#  write.table(ddest, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)
#  write.table(ddvar, file = file_var, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE) 
  
}

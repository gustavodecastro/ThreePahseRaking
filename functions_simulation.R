fun_samp <- function(data_emr, data_chart, data_phone, n1, n2, biased=FALSE){
  ## Get sampling probabilities
  ids_y_step1 <- aggregate(data_emr$y_emr, list(data_emr$id), max)
  ids_y1 <- c(sample(ids_y_step1[ids_y_step1$x>0,]$Group.1, n1), sample(ids_y_step1[ids_y_step1$x==0,]$Group.1, n1))
  ids_y_step1 <- ids_y_step1 %>% mutate(R = ifelse(Group.1 %in% ids_y1, 1, 0))
  ids_y_step1$probs <- predict(glm(R ~ x, data=ids_y_step1, family=binomial), type='response')
  ids_y_step1$probs <- ifelse(ids_y_step1$R==1, ids_y_step1$probs, 1-ids_y_step1$probs)
  ids_y_step1a <- ids_y_step1 %>% mutate(id = Group.1) %>% dplyr::select(id, probs)
  id_R1 <- ids_y_step1$Group.1[ids_y_step1$R==1]
  data_chart  <- data_chart %>% filter(id %in% id_R1) %>% left_join(ids_y_step1a, by='id')
  #glm(failed ~ X1 + X2 + offset(log(py)), data=data_chart, family=binomial)
  
  ## Get sampling probabilities
  if (biased){
    ids_y_step2  <- aggregate(data_chart$y_chart, list(data_chart$id), max)
  } else {
    ids_y_step2  <- aggregate(data_chart$y_emr, list(data_chart$id), max)  
  }
  ids_y2 <- c(sample(ids_y_step2[ids_y_step2$x>0,]$Group.1, n2), sample(ids_y_step2[ids_y_step2$x==0,]$Group.1, n2))
  ids_y_step2  <- ids_y_step2 %>% mutate(R = ifelse(Group.1 %in% ids_y2, 1, 0))
  ids_y_step2$probs2a <- predict(glm(R ~ x, data=ids_y_step2, family=binomial), type='response')
  ids_y_step2$probs2a <- ifelse(ids_y_step2$R==1, ids_y_step2$probs2a, 1-ids_y_step2$probs2a)
  ids_y_step2b <- ids_y_step2 %>% mutate(id = Group.1) %>% dplyr::select(id, probs2a)
  id_R2 <- ids_y_step2$Group.1[ids_y_step2$R==1]
  data_phone <- data_phone %>% filter(id %in% id_R2) %>% left_join(ids_y_step1a, by='id') %>% left_join(ids_y_step2b, by='id') %>%
    mutate(probs2 = probs2a*probs)
  
  return(list(data_chart=data_chart, data_phone=data_phone))
}

fun_extend <- function(data_emr, data_chart, data_phone) {
  data_chart_temp <- data_chart %>% dplyr::select(id, uid)
  data_emr <- data_emr %>%
    mutate(R1 = ifelse(id %in% data_chart_temp$id, 1, 0)) %>%
    mutate(R2 = ifelse(!(id %in% data_chart$id), NA,
                       ifelse(id %in% data_phone$id, 1, 0)))  
  data_emr2 <- data_emr[rep(1:nrow(data_emr), data_emr$n_followup),]
  
  data_emr3 <- data_emr2 %>%
    group_by(id) %>%
    mutate(nobs = 1:n(), uidA = paste0(y, X1, nobs)) %>% ungroup() %>%
    group_by(uid) %>%
    mutate(y = ifelse(row_number()==n(), y[n()], 0)) %>% ungroup()
  
  data_emr2 <- data_emr2 %>%
    group_by(id) %>%
    mutate(nobs = 1:n(), uidA = paste0(y, X1, nobs)) %>% ungroup() %>%
    group_by(uid) %>%
    mutate(y_emr = ifelse(row_number()==n(), y_emr[n()], 0)) %>% ungroup()
  
  data_chart2 <- data_chart[rep(1:nrow(data_chart), data_chart$n_followup),]
  data_chart2 <- data_chart2 %>%
    group_by(id) %>%
    mutate(nobs = 1:n(), uidA = paste0(y, X1, nobs)) %>% ungroup() %>%
    group_by(uid) %>%
    mutate(y_chart = ifelse(row_number()==n(), y_chart[n()], 0)) %>% ungroup()
  
  data_phone2 <- data_phone[rep(1:nrow(data_phone), data_phone$n_followup),]
  data_phone2 <- data_phone2 %>%
    group_by(id) %>%
    mutate(nobs = 1:n(), uidA = paste0(y, X1, nobs)) %>% ungroup() %>%
    group_by(uid) %>%
    mutate(y = ifelse(row_number()==n(), y[n()], 0)) %>% ungroup()
  
  data_chart2_temp <- data_chart2 %>% dplyr::select(uidA, y_chart, X1_chart, probs)
  data_phone2_temp <- data_phone2 %>% dplyr::select(uidA, y, X1, probs2, probs2a)
  data_emr2b <- data_emr2 %>% dplyr::select(-y, -X1, -y_chart, -X1_chart) %>%
    left_join(data_chart2_temp, by='uidA') %>% left_join(data_phone2_temp, by='uidA')    
  
  return(list(d1=data_emr2b, d2=data_emr3))
}

func_imp3ph_YX <- function(dd_emr){
  imp_y_step1 <- glm(y ~ X1_chart + X1_chart_M0 + X1_chart_M1 +
                       X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_chart + y_chart_M0 + y_chart_M1 +
                       y_emr + y_emr_M0 + y_emr_M1,
                     family = binomial, data = dd_emr[which(dd_emr$R2==1),])
  desmatrix<-desmatrixfn(imp_y_step1,summary(imp_y_step1),Nobs=1)
  linpred<-as.vector(desmatrix%*%t(cbind(1, dd_emr[which(dd_emr$R1==1),
                                                   names(imp_y_step1$coefficients[-1])])))
  linpred<-ifelse(linpred > 20,20,linpred)
  linpred<-ifelse(linpred < -20,-20,linpred)
  y_step1<-rbinom(length(linpred),1,exp(linpred)/(1+exp(linpred)))
  dd_emr$y_step1 <- NA
  dd_emr$y_step1[which(dd_emr$R1==1)] <- c(y_step1)
  dd_emr$y_step1_imputed <- NA
  dd_emr$y_step1_imputed <- dd_emr$y_step1
  dd_emr$y_step1[which(dd_emr$R2==1)] <- dd_emr$y[which(dd_emr$R2==1)]
  
  imp_x1_step1 <- lm(X1 ~ X1_chart + X1_chart_M0 + X1_chart_M1 +
                       X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_step1 + y_emr + y_emr_M0 + y_emr_M1,
                     data = dd_emr[which(dd_emr$R2==1),])
  estSigma_step1 <- sqrt(summary(imp_x1_step1)$sigma^2*rchisq(Nimpute, summary(imp_x1_step1)$df[2], ncp=0)/summary(imp_x1_step1)$df[2])
  estPar_step1   <- rmvnorm(1, mean=imp_x1_step1$coeff, sigma=summary(imp_x1_step1)$sigma^2*summary(imp_x1_step1)$cov)
  dd_emr$X1_step1 <- NA
  
  dd_emr$X1_step1[which(dd_emr$R1==1)] <- c(t(estPar_step1[1,] %*%
                                                t(cbind(1, dd_emr[which(dd_emr$R1==1),
                                                                  names(imp_x1_step1$coefficients[-1])])) +
                                                rnorm(sum(dd_emr$R1==1), 0, estSigma_step1[1])))
    
  dd_emr$X1_step1_imputed <- NA
  dd_emr$X1_step1_imputed <- dd_emr$X1_step1
  dd_emr$X1_step1[which(dd_emr$R2==1)] <- dd_emr$X1[which(dd_emr$R2==1)]
  
  imp_y_step2 <- glm(y_step1 ~ X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_emr + y_emr_M0 + y_emr_M1,
                     family = binomial, data = dd_emr[which(dd_emr$R1==1),])
  desmatrix<-desmatrixfn(imp_y_step2,summary(imp_y_step2),Nobs=1)
  linpred<-as.vector(desmatrix%*%t(cbind(1, dd_emr[, names(imp_y_step2$coefficients[-1])])))
  linpred<-ifelse(linpred > 20,20,linpred)
  linpred<-ifelse(linpred < -20,-20,linpred)
  y_step2<-rbinom(length(linpred),1,exp(linpred)/(1+exp(linpred)))
  dd_emr$y_step2 <- c(y_step2)
  dd_emr$y_step2_imputed <- dd_emr$y_step2
  dd_emr$y_step2[which(dd_emr$R1==1)] <- dd_emr$y_step1[which(dd_emr$R1==1)]
  
  imp_x1_step2 <- lm(X1_step1 ~ X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_step2 + y_emr + y_emr_M0 + y_emr_M1,
                     data = dd_emr[which(dd_emr$R1==1),])
  estSigma_step2 <- sqrt(summary(imp_x1_step2)$sigma^2*rchisq(Nimpute, summary(imp_x1_step2)$df[2], ncp=0)/summary(imp_x1_step2)$df[2])
  estPar_step2   <- rmvnorm(1, mean=imp_x1_step2$coeff, sigma=summary(imp_x1_step2)$sigma^2*summary(imp_x1_step2)$cov)
  dd_emr$X1_step2 <- c(t(estPar_step2[1,] %*% t(cbind(1, dd_emr[, names(imp_x1_step2$coefficients[-1])])) + rnorm(nrow(dd_emr), 0, estSigma_step2[1])))
  #dd_emr$X1_step2 <- c(t(estPar_step2[1,] %*% t(cbind(1, dd_emr$X1_emr, dd_emr$X2, dd_emr$y_step2, dd_emr$y_emr)) + rnorm(nrow(dd_emr), 0, estSigma_step2[1])))
  dd_emr$X1_step2_imputed <- dd_emr$X1_step2
  dd_emr$X1_step2[which(dd_emr$R1==1)] <- dd_emr$X1_step1[which(dd_emr$R1==1)]
  
  dd_emr
}

func_imp3ph_YX_raking <- function(dd_emr){
  imp_y_step1 <- glm(y ~ X1_chart + X1_chart_M0 + X1_chart_M1 + 
                       X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_chart + y_chart_M0 + y_chart_M1 + 
                       y_emr + y_emr_M0 + y_emr_M1,
                     family = binomial, data = dd_emr[which(dd_emr$R2==1),])
  desmatrix<-desmatrixfn(imp_y_step1,summary(imp_y_step1),Nobs=1)
  linpred<-as.vector(desmatrix%*%t(cbind(1, dd_emr[which(dd_emr$R1==1),
                                                   names(imp_y_step1$coefficients[-1])])))
  linpred<-ifelse(linpred > 20,20,linpred)
  linpred<-ifelse(linpred < -20,-20,linpred)
  y_step1<-rbinom(length(linpred),1,exp(linpred)/(1+exp(linpred)))
  dd_emr$y_step1 <- NA
  dd_emr$y_step1[which(dd_emr$R1==1)] <- c(y_step1)
  dd_emr$y_step1_imputed <- NA
  dd_emr$y_step1_imputed <- dd_emr$y_step1
  
  imp_x1_step1 <- lm(X1 ~ X1_chart + X1_chart_M0 + X1_chart_M1 + 
                       X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_chart + y_chart_M0 + y_chart_M1 + 
                       y_emr + y_emr_M0 + y_emr_M1 + y_step1_imputed,
                     data = dd_emr[which(dd_emr$R2==1),])
  estSigma_step1 <- sqrt(summary(imp_x1_step1)$sigma^2*rchisq(Nimpute, summary(imp_x1_step1)$df[2], ncp=0)/summary(imp_x1_step1)$df[2])
  estPar_step1   <- rmvnorm(1, mean=imp_x1_step1$coeff, sigma=summary(imp_x1_step1)$sigma^2*summary(imp_x1_step1)$cov)
  dd_emr$X1_step1 <- NA
  dd_emr$X1_step1[which(dd_emr$R1==1)] <- c(t(estPar_step1[1,] %*% t(cbind(1, dd_emr[which(dd_emr$R1==1), names(imp_x1_step1$coefficients[-1])])) +
                                                rnorm(sum(dd_emr$R1==1), 0, estSigma_step1[1])))
  dd_emr$X1_step1_imputed <- NA
  dd_emr$X1_step1_imputed <- dd_emr$X1_step1
  
  imp_y_step2 <- glm(y_step1_imputed ~ X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_emr + y_emr_M0 + y_emr_M1,
                     family = binomial, data = dd_emr[which(dd_emr$R1==1),])
  desmatrix<-desmatrixfn(imp_y_step2,summary(imp_y_step2),Nobs=1)
  linpred<-as.vector(desmatrix%*%t(cbind(1, dd_emr[, names(imp_y_step2$coefficients[-1])])))
  linpred<-ifelse(linpred > 20,20,linpred)
  linpred<-ifelse(linpred < -20,-20,linpred)
  y_step2<-rbinom(length(linpred),1,exp(linpred)/(1+exp(linpred)))
  dd_emr$y_step2 <- c(y_step2)
  dd_emr$y_step2_imputed <- dd_emr$y_step2
  
  imp_x1_step2 <- lm(X1_step1_imputed ~ X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                       y_emr + y_emr_M0 + y_emr_M1 + y_step2_imputed,
                     data = dd_emr[which(dd_emr$R1==1),])
  estSigma_step2 <- sqrt(summary(imp_x1_step2)$sigma^2*rchisq(Nimpute, summary(imp_x1_step2)$df[2], ncp=0)/summary(imp_x1_step2)$df[2])
  estPar_step2   <- rmvnorm(1, mean=imp_x1_step2$coeff, sigma=summary(imp_x1_step2)$sigma^2*summary(imp_x1_step2)$cov)
  dd_emr$X1_step2 <- c(t(estPar_step2[1,] %*% t(cbind(1, dd_emr[, names(imp_x1_step2$coefficients[-1])])) + rnorm(nrow(dd_emr), 0, estSigma_step2[1])))
  dd_emr$X1_step2_imputed <- dd_emr$X1_step2
  
  dd_emr
}

func_imp2ph_YX <- function(dd_emr){
  imp_y1_2ph <- glm(y ~ X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                      y_emr + y_emr_M0 + y_emr_M1,
                    family = binomial, data = dd_emr[which(dd_emr$R2==1),])
  desmatrix<-desmatrixfn(imp_y1_2ph,summary(imp_y1_2ph),Nobs=1)
  linpred<-as.vector(desmatrix%*%t(cbind(1, dd_emr[, names(imp_y1_2ph$coefficients[-1])])))
  linpred<-ifelse(linpred > 20,20,linpred)
  linpred<-ifelse(linpred < -20,-20,linpred)
  dd_emr$y_new <- rbinom(length(linpred),1,exp(linpred)/(1+exp(linpred)))
  dd_emr$y_new_imputed <- dd_emr$y_new
  dd_emr$y_new[which(dd_emr$R2==1)] <- dd_emr$y[which(dd_emr$R2==1)]
  
  imp_x1_2ph <- lm(X1 ~ X1_emr + X1_emr_M0 + X1_emr_M1 + X2 +
                     y_new + y_emr + y_emr_M0 + y_emr_M1,
                   data = dd_emr[which(dd_emr$R2==1),])
  estSigma_2ph <- sqrt(summary(imp_x1_2ph)$sigma^2*rchisq(Nimpute, summary(imp_x1_2ph)$df[2], ncp=0)/summary(imp_x1_2ph)$df[2])
  estPar_2ph   <- rmvnorm(1, mean=imp_x1_2ph$coeff, sigma=summary(imp_x1_2ph)$sigma^2*summary(imp_x1_2ph)$cov)
  dd_emr$X1_new <- c(t(estPar_2ph[1,] %*%
                         t(cbind(1, dd_emr[, names(imp_x1_2ph$coefficients[-1])])) +
                         rnorm(nrow(dd_emr), 0, estSigma_2ph[1])))
  dd_emr$X1_new_imputed <- dd_emr$X1_new 
  dd_emr$X1_new[which(dd_emr$R2==1)] <- dd_emr$X1[which(dd_emr$R2==1)]
  
  dd_emr
}

func_imp2ph_YX_raking <- function(dd_emr){
  imp_y1_2ph <- glm(y ~ X1_emr + X2 + y_emr, family = binomial, data = dd_emr[which(dd_emr$R2==1),])
  desmatrix<-desmatrixfn(imp_y1_2ph,summary(imp_y1_2ph),Nobs=1)
  linpred<-as.vector(desmatrix%*%t(cbind(1, dd_emr[, names(imp_y1_2ph$coefficients[-1])])))
  linpred<-ifelse(linpred > 20,20,linpred)
  linpred<-ifelse(linpred < -20,-20,linpred)
  dd_emr$y_new <- rbinom(length(linpred),1,exp(linpred)/(1+exp(linpred)))
  dd_emr$y_new_imputed <- dd_emr$y_new
  
  imp_x1_2ph <- lm(X1 ~ X1_emr + X2 + y_new_imputed + y_emr, data = dd_emr[which(dd_emr$R2==1),])
  estSigma_2ph <- sqrt(summary(imp_x1_2ph)$sigma^2*rchisq(Nimpute, summary(imp_x1_2ph)$df[2], ncp=0)/summary(imp_x1_2ph)$df[2])
  estPar_2ph   <- rmvnorm(1, mean=imp_x1_2ph$coeff, sigma=summary(imp_x1_2ph)$sigma^2*summary(imp_x1_2ph)$cov)
  dd_emr$X1_new <- c(t(estPar_2ph[1,] %*% t(cbind(1, dd_emr$X1_emr, dd_emr$X2, dd_emr$y_new_imputed, dd_emr$y_emr)) + rnorm(nrow(dd_emr), 0, estSigma_2ph[1])))
  dd_emr$X1_new_imputed <- dd_emr$X1_new 
  
  dd_emr    
}

desmatrixfn<-function(mod,summ.mod,Nobs=Nimpute){
  desmatrix<-rmvnorm(Nobs, mean=mod$coeff[rownames(summ.mod$coef)], sigma=(summ.mod$dispersion*summ.mod$cov.unscaled[rownames(summ.mod$coef),rownames(summ.mod$coef)]))
}


fun_raking <- function(data_emr2b=dd_emr) {
  fit_ph1 <- glm(y_emr ~ X1_emr + X2, family=binomial, data=data_emr2b)
  ifs1a <- influence_func(fit_ph1)
  ifs   <- aggregate(ifs1a, list(data_emr2b$id), sum)
  colnames(ifs) <- c('id', paste0('if', 1:ncol(ifs1a)))
  totals1 <- c(`(Intercept)`=length(unique(data_emr2b$id)), colSums(ifs)[-1])
  
  fit_ph2 <- glm(y_chart ~ X1_chart + X2, family=binomial, data=data_emr2b[which(data_emr2b$R1==1),])
  ifs2a <- influence_func(fit_ph2)
  ifs2  <- aggregate(ifs2a, list(data_emr2b$id[which(data_emr2b$R1==1)]), sum)
  colnames(ifs2) <- c('id', paste0('if', 1:ncol(ifs2a)))
  datcsv2 <- data_emr2b[which(data_emr2b$R2==1),] %>% left_join(ifs2, by='id') %>% dplyr::select(id, if1, if2, if3) %>% group_by(id) %>% slice(1) %>% ungroup()
  totals2 <- c(`(Intercept)`=length(unique(datcsv2$id)), colSums(ifs2)[-1])
  
  data_phoneb <- data_emr2b[which(data_emr2b$R2==1),] %>% left_join(ifs, by='id')
  if_designR <- survey::svydesign(id = ~id, probs = ~ probs2, data = data_phoneb)
  if_cal2ph  <- survey::calibrate(design=if_designR, calfun='raking', formula=~if1+if2+if3, population = totals1)
  resCAL_2ph <- svyglm(y ~ X1 + X2, design=if_cal2ph, family=binomial)
  resIPW_2ph <- svyglm(y ~ X1 + X2, design=if_designR, family=binomial)
  
  return(list(resCAL_2ph=resCAL_2ph, resIPW_2ph=resIPW_2ph, datcsv1=ifs, totals1=totals1, datcsv2=datcsv2, totals2=totals2))
}




fun_raking_3ph <- function(data_emr2b, totals1, totals2, datcsv1, datcsv2) {
  data_chartb <- data_emr2b %>% filter(R1==1) %>% left_join(datcsv1, by='id')
  data_phoneb <- data_emr2b[which(data_emr2b$R2==1),] %>% left_join(datcsv2, by='id')
  des2pha <- survey::svydesign(id = ~id, probs = ~probs, data = data_chartb)
  des3pha <- survey::svydesign(id = ~id, probs = ~probs2a, data = data_phoneb)
  cal2pha <- survey::calibrate(design=des2pha, calfun='raking', formula=~if1+if2+if3, population = totals1)
  cal3pha <- survey::calibrate(design=des3pha, calfun='raking', formula=~if1+if2+if3, population = totals2)
  
  pi1_new <- data.frame(id = cal2pha$cluster$id, pi1_new = cal2pha$prob) %>% group_by(id) %>% slice(1) %>% ungroup()
  pi2_new <- data.frame(id = cal3pha$cluster$id, pi2_new = cal3pha$prob) %>% group_by(id) %>% slice(1) %>% ungroup()
  datcsv <- data_phoneb %>% left_join(pi1_new, by='id') %>% left_join(pi2_new, by='id') %>% mutate(pprod_new = pi1_new*pi2_new)
  raking_new <- glm(y ~ X1 + X2, weights=1/pprod_new, family=binomial, data=datcsv)
  
  if_designR2 <- survey::svydesign(id = ~id, probs = ~ pprod_new, data = datcsv)
  resIPW_2ph  <- svyglm(y ~ X1 + X2, design=if_designR2, family=binomial)
  
  return(list(resIPW_2ph=resIPW_2ph))
}




fun_rak2ph_MI <- function(dd_emr, ifs2ph_step1, i=NULL) {
  if (is.null(i)){
    ifs2ph_step1_all <- Reduce('+', ifs2ph_step1)/length(ifs2ph_step1)
  } else{
    ifs2ph_step1_all <- ifs2ph_step1[[i]]
  }
  
  #  row_j <- matrix(NA, ncol=ncol(ifs2ph_step1[[1]]), nrow=nrow(ifs2ph_step1[[1]]))
  #  col_j <- matrix(NA, ncol=length(ifs2ph_step1), nrow=nrow(ifs2ph_step1[[1]]))
  #  for (jj in 2:ncol(ifs2ph_step1[[1]])){
  #    for (j in 1:length(ifs2ph_step1)){
  #      col_j[,j] <- c(ifs2ph_step1[[j]][,jj])
  #    }
  #    row_j[,jj] <- apply(col_j, 1, median)
  #  }
  #  ifs2ph_step1_all <- row_j
  #  ifs2ph_step1_all[,1] <- ifs2ph_step1[[1]][,1]
  #  ifs2ph_step1_all <- as.data.frame(ifs2ph_step1_all)
  #  
  colnames(ifs2ph_step1_all) <- c('id', 'if1', 'if2', 'if3')
  sum_ifs_2ph_MI <- c(length(unique(dd_emr$id)), colSums(ifs2ph_step1_all)[-1])
  names(sum_ifs_2ph_MI)[c(1)] <- c('(Intercept)')
  
  #  sum_ifs_2ph_MI[2:4] <- c(0,0,0)
  
  data_phoneb <- dd_emr[which(dd_emr$R2==1),] %>% left_join(ifs2ph_step1_all, by='id')
  desing_2phMI <- survey::svydesign(id = ~id, probs = ~ probs2, data = data_phoneb)
  if_cal2phMI  <- survey::calibrate(design=desing_2phMI, calfun='raking', formula=~if1+if2+if3, population = sum_ifs_2ph_MI)
  resCAL_2phMI <- svyglm(y ~ X1 + X2, design=if_cal2phMI, family=binomial)
  resCAL_2phMI
}


fun_rak3ph_MI <- function(dd_emr, ifs3ph_step1, ifs3ph_step2) {
  
  if_3ph_step1B_all <- Reduce('+', ifs3ph_step1)/length(ifs3ph_step1)
  #  row_j <- matrix(NA, ncol=ncol(ifs3ph_step1[[1]]), nrow=nrow(ifs3ph_step1[[1]]))
  #  col_j <- matrix(NA, ncol=length(ifs3ph_step1), nrow=nrow(ifs3ph_step1[[1]]))
  #  for (jj in 2:ncol(ifs3ph_step1[[1]])){
  #    for (j in 1:length(ifs3ph_step1)){
  #      col_j[,j] <- c(ifs3ph_step1[[j]][,jj])
  #    }
  #    row_j[,jj] <- apply(col_j, 1, median)
  #  }
  #  if_3ph_step1B_all <- as.data.frame(row_j)
  #  if_3ph_step1B_all[,1] <- ifs3ph_step1[[1]][,1]
  
  if_3ph_step2B_all <- Reduce('+', ifs3ph_step2)/length(ifs3ph_step2)
  #  row_j <- matrix(NA, ncol=ncol(ifs3ph_step2[[1]]), nrow=nrow(ifs3ph_step2[[1]]))
  #  col_j <- matrix(NA, ncol=length(ifs3ph_step2), nrow=nrow(ifs3ph_step2[[1]]))
  #  for (jj in 2:ncol(ifs3ph_step2[[1]])){
  #    for (j in 1:length(ifs3ph_step2)){
  #      col_j[,j] <- c(ifs3ph_step2[[j]][,jj])
  #    }
  #    row_j[,jj] <- apply(col_j, 1, median)
  #  }
  #  if_3ph_step2B_all <- as.data.frame(row_j)
  #  if_3ph_step2B_all[,1] <- ifs3ph_step2[[1]][,1]
  
  colnames(if_3ph_step1B_all) <- colnames(if_3ph_step2B_all) <- c('id', 'if1', 'if2', 'if3')
  ##
  totals_3ph2 <- cbind(length(unique(dd_emr$id)), matrix(colSums(if_3ph_step1B_all[,-1]), nrow=1))
  colnames(totals_3ph2) <- c('(Intercept)', colnames(if_3ph_step1B_all)[-1])
  totals_3ph1 <- cbind(length(unique(dd_emr$id[dd_emr$R1==1])), matrix(colSums(if_3ph_step2B_all[,-1]), nrow=1))
  colnames(totals_3ph1) <- c('(Intercept)', colnames(if_3ph_step2B_all)[-1])
  
  #  totals_3ph1[2:4] <- totals_3ph2[2:4] <- c(0,0,0)
  
  ## Raking, 3-phase
  data_chartb <- dd_emr %>% filter(R1==1) %>% left_join(if_3ph_step2B_all, by='id') #%>% filter(R1==1) %>% left_join(data_chartB, by='id')
  data_phoneb <- dd_emr[which(dd_emr$R2==1),] %>% left_join(if_3ph_step1B_all, by='id') #%>% filter(R2==1) %>% left_join(data_phoneB, by='id')
  des2pha <- survey::svydesign(id = ~id, probs = ~probs, data = data_chartb)
  des3pha <- survey::svydesign(id = ~id, probs = ~probs2a, data = data_phoneb)
  cal2pha <- survey::calibrate(design=des2pha, calfun='raking', formula=~if1+if2+if3, population = totals_3ph2[c(1:4)])
  cal3pha <- survey::calibrate(design=des3pha, calfun='raking', formula=~if1+if2+if3, population = totals_3ph1[c(1:4)])
  
  pi1_new <- data.frame(id = cal2pha$cluster$id, pi1_new = cal2pha$prob) %>% group_by(id) %>% slice(1)
  pi2_new <- data.frame(id = cal3pha$cluster$id, pi2_new = cal3pha$prob) %>% group_by(id) %>% slice(1)
  datcsv <- data_phoneb %>%
    left_join(pi1_new, by='id') %>%
    left_join(pi2_new, by='id') %>%
    mutate(pprod_new = pi1_new*pi2_new)
  raking_new <- glm(y ~ X1 + X2, weights=1/pprod_new, family=binomial, data=datcsv)
  
  raking_new
}

simulWeib <- function(N, lambda, rho, beta, rateC, nr=6, maxT=6, beta0_emr, beta0_chart, gamma_x1=FALSE)
{
  # covariate --> N Bernoulli trials
  x1 <- rnorm(N*nr*3, 0, 1)
  x2a <- rnorm(N, 0, 1)
  x2a <- rep(x2a, each=nr)
  x2b <- rnorm(N, 0, 1)
  x2b <- rep(x2b, each=nr)
  x2c <- rnorm(N, 0, 1)
  x2c <- rep(x2c, each=nr)
  x2 <- c(x2a, x2b, x2c)
  #### 
  rho <- .5
  #x1 <- x2*rho + sqrt(1-rho^2)*rnorm(length(x1),0,1)
  ####
  
  if (gamma_x1)
    x1 <- log(rgamma(length(x1), 10)) + .5*x2 + .1*x2^2
  
  # Weibull latent event times
  # v <- runif(n=N)
  # Tlat <- (- log(v) / (lambda * exp(x * beta)))^(1 / rho)
  
  # An alternative draw for event times
  lambda_wiki = lambda^(-1 / rho) #change definition of lambda to Wikipedia's
  lambda_prime = lambda_wiki / exp(c(x1 * beta[1] + x2 * beta[2]) / rho) #re-scale according to beta
  Tlat = rweibull(length(x1), shape=rho, scale=lambda_prime)
  
  # censoring times
  C <- rexp(n=length(x1), rate=rateC)
  C <- maxT
  
  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  
  # data set
  dat <- data.frame(id=rep(1:N, each=length(x1)/N),
                    time=time,
                    time_trunc=ifelse(time>10, 10, time),
                    n_followup=ceiling(ifelse(time>10, 10, time)),
                    #y=ifelse(time>10, 0, status),
                    y=status,
                    X1=x1,
                    X2=x2) %>%
    mutate(pr_chart = exp(beta0_chart[1] + beta0_chart[2]*y),
#           pr_chart = pr_chart/(1+pr_chart),
           pr_emr = exp(beta0_emr[1] + beta0_emr[2]*y),
#           pr_emr = pr_emr/(1+pr_emr),
           #           y_chart = rbinom(length(y), 1, pr_chart/(1+pr_chart)),
           #           y_emr = rbinom(length(y), 1, pr_emr/(1+pr_emr)),
           y_chart = rbinom(length(y), 1, ifelse(pr_chart>1, 1, pr_chart)),
           y_emr = rbinom(length(y), 1, ifelse(pr_emr>1, 1, pr_emr)),
           #####
           X1_emr = X1 + rnorm(length(X1), 0, sx_emr),
           X1_chart = X1 + rnorm(length(X1), 0, sx_chart),
           #####
           uid = paste0(id, x1, x2))
  
  dat
}

influence_func <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}

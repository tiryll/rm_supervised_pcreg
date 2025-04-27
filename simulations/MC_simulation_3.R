
rm(list=ls())

source(paste(getwd(),"/environment.R",sep=""))
source(paste(getwd(),"/functions/functions_pcr.R",sep=""))
source(paste(getwd(),"/functions/MC_simulation_functions.R",sep=""))
source(paste(getwd(),"/functions/functions_pls.R",sep=""))

############################
# Simulation 3
############################

set.seed(4519)

#model parameters
N <- 100
p <- 5000
p1 <- 50
p2 <- 0
X1_shift <- c(-1,1)
X2_shift <- c(-2,2)
sigma2_0 <- 1
sigma2_1 <- 1.5^2
beta<-2

Kpcr <- 2
Kspcr <- 2
Kpls <- 2

# simulation variables
Nsim <- 20
spcr_MSE <- numeric(Nsim)
pcr_MSE <- numeric(Nsim)
pls_MSE <- numeric(Nsim)
naive_MSE <- numeric(Nsim)
X_select <- numeric(Nsim)
X_p1_select <- numeric(Nsim)
coeff_spcr <- matrix(data = 0, nrow = Nsim, ncol = p)
coeff_pcr <- matrix(data = 0, nrow = Nsim, ncol = p)
coeff_pls <- matrix(data = 0, nrow = Nsim, ncol = p)

for (iter in 1:Nsim) {
  
  train_sim <- gene_simulation2(N,p,p1,p2,X1_shift,X2_shift,beta,sigma2_0,sigma2_1)
  X_train <- train_sim$X
  y_train <- train_sim$y
  
  spcr_mod <- spcreg_fit(y_train,X_train,Kspcr)
  pcr_mod <- pcreg_fit(y_train,X_train,Kpcr,scaling = "center")
  pls_mod <- plsreg_fit(y_train, X_train,Kpls, scaling = "center")
  
  coeff_pcr[iter,]<-pcr_mod$coefficients[2:length(pcr_mod$coefficients)]
  coeff_pls[iter,]<-pls_mod$coefficients[2:length(pls_mod$coefficients)]
  coeff_spcr[iter,]<-spcr_mod$coefficients[2:length(spcr_mod$coefficients)]
  
  X_select[iter] <- sum(spcr_mod$coefficients[2:length(spcr_mod$coefficients)] != 0)
  X_p1_select[iter] <- sum(spcr_mod$coefficients[2:(p1+1)] != 0)
  
  test_sim <- gene_simulation2(N,p,p1,p2,X1_shift,X2_shift,beta,sigma2_0,sigma2_1)
  X_test <- test_sim$X
  y_test <- test_sim$y
  
  spcr_pred <- spcreg_pred(spcr_mod,X_test)
  pcr_pred <- pcreg_pred(pcr_mod,X_test)
  pls_pred <- plsreg_pred(pls_mod, X_test)
  
  spcr_MSE[iter] <- 1/N*t(y_test-spcr_pred$prediction)%*%(y_test-spcr_pred$prediction)
  pcr_MSE[iter] <- 1/N*t(y_test-pcr_pred$prediction)%*%(y_test-pcr_pred$prediction)
  pls_MSE[iter] <- 1/N*t(y_test-pls_pred$prediction)%*%(y_test-pls_pred$prediction)
  naive_MSE[iter] <- 1/N*t(y_test - mean(y_train))%*%(y_test - mean(y_train))
  
  print(iter)
}

### Results Simulation 3

# Table 1 Simulation Parameters
sim3_param <-data.frame(t(c(Nsim,N,p,p1,p2,t(X1_shift),t(X2_shift),beta,sigma2_0,sqrt(sigma2_1))),row.names = NULL)

colnames(sim3_param) <- c('Nsims','Nobs','p','p$_1$','p$_2$','$\\mu_1$',
                          '$\\mu_2$','$\\mu_3$','$\\mu_4$','$\\beta$','$\\sigma_1$','$\\sigma_2$')


sim3_tab1 <- xtable(sim3_param,caption = "Simulation Parameters")

filename_sim3_tab1 <- paste(tables_dir,"/sim3_tab1.txt",sep="")

print(sim3_tab1,sanitize.colnames.function = identity, include.rownames=FALSE,file=filename_sim3_tab1)

# Table 2 MSE Results
sim3_summstat_MSE <- data.frame(rbind(compute_summary_stats(pcr_MSE),
                                      compute_summary_stats(pls_MSE),
                                      compute_summary_stats(spcr_MSE),
                                      compute_summary_stats(naive_MSE)))

rownames(sim3_summstat_MSE) <- c("PCR","PLS","SPCR","Naive")

sim3_tab2 <- xtable(sim3_summstat_MSE,caption = "Summary Statistics of Mean Squared Errors")

filename_sim3_tab2 <- paste(tables_dir,"/sim3_tab2.txt",sep="")

print(sim3_tab2,sanitize.colnames.function = identity, include.rownames=FALSE,file=filename_sim3_tab2)

# Table 3 Accuracy of SPCR

sim3_summstat_X <- data.frame(rbind(compute_summary_stats(X_select),
                                    compute_summary_stats(X_p1_select)))

rownames(sim3_summstat_X) <- c("Total","$X_\\{j<p_1\\}$")
colnames(sim3_summstat_X) <- c("Mean","Var","$25^{th}$","$50^{th}$","$75^{th}$")

sim3_tab3 <- xtable(sim3_summstat_X,caption = "Summary Statistics of the Number of Total Selected Variables and 
                    Correctly Selected Variables by the SPC method.")

filename_sim3_tab3 <- paste(tables_dir,"/sim3_tab3.txt",sep="")

print(sim3_tab3,sanitize.colnames.function = identity, sanitize.rownames.function = identity, 
      include.rownames=FALSE,file=filename_sim3_tab3)


### Export MSE and coefficients data for boxplots

filename_sim3_MSE_data <- paste(data_dir,"/sim3_MSE_data.txt",sep="")
filename_sim3_coeff_data <- paste(data_dir,"/sim3_coeff_data.rds",sep="")
filename_sim3_param <- paste(data_dir,"/sim3_param.rds",sep="")

write.table(data.frame(naive_MSE, pcr_MSE, spcr_MSE, pls_MSE),file=filename_sim3_MSE_data)

write_rds(list(coeff_spcr=coeff_spcr,coeff_pcr=coeff_pcr,coeff_pls=coeff_pls),
          file=filename_sim3_coeff_data,
          compress="none")

sim3_param <- list (N = N,
                    p = p,
                    p1 = p1,
                    p2 = p2,                     
                    X1_shift = X1_shift,
                    X2_shift = X2_shift,
                    sigma2_0 =sigma2_0,
                    sigma2_1 =sigma2_1,
                    beta=beta,
                    Kpcr = Kpcr,
                    Kspcr = Kspcr,
                    Kpls = Kpls,
                    Nsim =Nsim)

write_rds(sim3_param,file=filename_sim3_param,compress="none")
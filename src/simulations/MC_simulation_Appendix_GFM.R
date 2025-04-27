
rm(list=ls())

source(paste(getwd(),"/src/environment.R",sep=""))
source(paste(getwd(),"/src/functions_pcr.R",sep=""))
source(paste(getwd(),"/src/MC_simulation_functions.R",sep=""))
source(paste(getwd(),"/src/functions_pls.R",sep=""))

##################################
# Simulation Gaussian Factor Model
##################################

set.seed(4519)

#model parameters
N <- 100
p <- 5000
p1 <- 50
p2 <- 250
sigma2_U1 <- 1
sigma2_U2 <- 4
sigma2_0 <- 1
sigma2_1 <- 1.5^2
beta<-2

Kpcr <- 1
Kspcr <- 1
Kpls <- 1

# simulation variables
Nsim <- 20
spcr_MSE <- numeric(Nsim)
pcr_MSE <- numeric(Nsim)
pls_MSE <- numeric(Nsim)
naive_MSE <- numeric(Nsim)
X_select <- numeric(Nsim)
X_p1_select <- numeric(Nsim)

factor1 <- rnorm(N,0,sd=sqrt(sigma2_U1))
factor2 <- rnorm(N,0,sd=sqrt(sigma2_U2))

for (iter in 1:Nsim) {
  
  train_sim <- factor_simulation(N,p,p1,p2,factor1,factor2,beta,sigma2_0,sigma2_1)
  X_train <- train_sim$X
  y_train <- train_sim$y
  
  spcr_mod <- spcreg_fit(y_train,X_train,Kspcr)
  pcr_mod <- pcreg_fit(y_train,X_train,Kpcr,scaling = "center")
  pls_mod <- plsreg_fit(y_train, X_train,Kpls, scaling = "standardize")
  
  X_select[iter] <- sum(spcr_mod$coefficients[2:length(spcr_mod$coefficients)] != 0)
  X_p1_select[iter] <- sum(spcr_mod$coefficients[2:(p1+1)] != 0)
  
  test_sim <- factor_simulation(N,p,p1,p2,factor1,factor2,beta,sigma2_0,sigma2_1)
  X_test <- test_sim$X
  y_test <- test_sim$y
  
  spcr_pred <- spcreg_pred(spcr_mod,X_test)
  pcr_pred <- pcreg_pred(pcr_mod,X_test)
  pls_pred <- plsreg_pred(pls_mod, X_test)
  
  spcr_MSE[iter] <- 1/N*t(y_test-spcr_pred$prediction)%*%(y_test-spcr_pred$prediction)
  pcr_MSE[iter] <- 1/N*t(y_test-pcr_pred$prediction)%*%(y_test-pcr_pred$prediction)
  pls_MSE[iter] <- 1/N*t(y_test-pls_pred$prediction)%*%(y_test-pls_pred$prediction)
  naive_MSE[iter] <- 1/N*t(y_test - mean(y_train))%*%(y_test - mean(y_train))
  
  print(iter/Nsim)
}

# Results Gaussian Factor Model Simulation

# Table 1 Simulation Parameters
simGF_param <-data.frame(t(c(Nsim,N,p,p1,p2,beta,sqrt(sigma2_U1),sqrt(sigma2_U2),sqrt(sigma2_0),sqrt(sigma2_1))),row.names = NULL)

colnames(simGF_param) <- c('Nsims','Nobs','p','p$_1$','p$_2$','$\\beta$',
                           '$\\sigma_{U_1}$','$\\sigma_{U_2}$','$\\sigma_1$','$\\sigma_2$')


simGF_tab1 <- xtable(simGF_param,caption = "Simulation Parameters")

filename_simGF_tab1 <- paste(tables_dir,"/simGF_tab1.txt",sep="")

print(simGF_tab1,sanitize.colnames.function = identity, include.rownames=FALSE,file=filename_simGF_tab1)

# Table 2 MSE Results
simGF_summstat_MSE <- data.frame(rbind(compute_summary_stats(pcr_MSE),
                                      compute_summary_stats(pls_MSE),
                                      compute_summary_stats(spcr_MSE),
                                      compute_summary_stats(naive_MSE)))

rownames(simGF_summstat_MSE) <- c("PCR","PLS","SPCR","Naive")

simGF_tab2 <- xtable(simGF_summstat_MSE,caption = "Summary Statistics of Mean Squared Errors")

filename_simGF_tab2 <- paste(tables_dir,"/simGF_tab2.txt",sep="")

print(simGF_tab2,sanitize.colnames.function = identity, include.rownames=FALSE,file=filename_simGF_tab2)

# Table 3 Accuracy of SPCR

simGF_summstat_X <- data.frame(rbind(compute_summary_stats(X_select),
                                    compute_summary_stats(X_p1_select)))

rownames(simGF_summstat_X) <- c("Total","$X_\\{j<p_1\\}$")
colnames(simGF_summstat_X) <- c("Mean","Var","$25^{th}$","$50^{th}$","$75^{th}$")

simGF_tab3 <- xtable(simGF_summstat_X,caption = "Summary Statistics of the Number of Total Selected Variables and 
                    Correctly Selected Variables by the SPC method.")

filename_simGF_tab3 <- paste(tables_dir,"/simGF_tab3.txt",sep="")

print(simGF_tab3,sanitize.colnames.function = identity, sanitize.rownames.function = identity, 
      include.rownames=FALSE,file=filename_simGF_tab3)

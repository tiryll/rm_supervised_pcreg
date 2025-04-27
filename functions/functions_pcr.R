pcreg_fit <- function(y,X,M,scaling=NULL) {
  #' Fit Principal Component Regression
  #'
  #' Extracts M principal components of a scaled data matrix X by the 
  #' singular value decomposition and regresses y on these principal components
  #' by OLS. 
  #'
  #' @param y A numeric vector of dim n x 1 to be projected onto principal components 
  #' @param X A matrix of column vectors of dim n x p.  
  #' @param K number of principal components to be used in regression
  #' @param scaling optional parameter either "center" for mean centering the column vectors X
  #' or "standardized" for standardizing the column vectors X  
  #'
  #' @return Returns an list, containing the coefficients, fitted values, 
  #' the MSE for the training data, the scaling applied to X and also the 
  #' mean and standard Deviations of the column vectors in X.   
  #' 
  #' @examples
  #'  set.seed(4519)
  #'
  #' N <- 100
  #' p <- 50
  #'
  #' X_train <- matrix(rnorm(N*p),N,p)
  #' y_train <- 2*rowMeans(X_train[,1:5])+rnorm(1.5)
  #'
  #' pcr_train <- pcreg_fit(y_train,X_train,2,"standardize")
  #'   
  if (!is.null(scaling)) {
    if (scaling == "center") {
      X_sc <- scale(X,center=TRUE,scale=FALSE)   
      } else if (scaling == "standardize") {
        X_sc <- scale(X,center=TRUE,scale=TRUE)
        }
  } else {
    X_sc <- X
  }
  
  n=length(y)
  X_svd <- svd(X_sc, nu=0, nv=M+1)
  V <- X_svd$v
  
  Z <- X_sc%*%V[,1:M]
  
  beta_pcr <- V[,1:M]%*%solve(crossprod(Z))%*%crossprod(Z,y)
  
  pcr_coefs <- c(mean(y),beta_pcr)
  y_hat <- cbind(rep(1,n),X_sc)%*%pcr_coefs
  err_hat <-y-y_hat
  
  MSE_train <- 1/n*t(err_hat)%*%err_hat
    
  out <- list("coefficients"=pcr_coefs,"fittedvalues"=y_hat,"trainMSE"=MSE_train,
              "scaling"=scaling,"center"=colMeans(X),"StdDev"=apply(X,2,sd))
  return(out)
}


pcreg_pred <- function(pcr_mod,X){
  #' Prediction for Principal Component Regression
  #'
  #' uses the object returned by pcreg_fit() or pcreg_CV() to make a prediction for a
  #' test data set X. Note X is scaled according to the scaling performed in  
  #' the  pcreg_fit() or pcreg_CV() object by the column means and standard Deviations of the 
  #' training data used in pcreg_fit or pcreg_CV().
  #' 
  #' @param pcr_mod object returned by pcreg_fit() or pcreg_CV()
  #' @param X new test data of dim n_test x p
  #'
  #' @return predicted values
  #' 
  #' @examples
  #' 
  #'  set.seed(4519)
  #'
  #' N <- 100
  #' p <- 50
  #'
  #' X_train <- matrix(rnorm(N*p),N,p)
  #' y_train <- 2*rowMeans(X_train[,1:5])+rnorm(1.5)
  #'
  #' pcr_train <- pcreg_fit(y_train,X_train,2,"standardize")
  #'
  #' X_test <- matrix(rnorm(N*p),N,p)
  #' y_test <- 2*rowMeans(X_train[,1:5])+rnorm(1.5)
  #'
  #' pcr_test <- pcr_pred(pcr_train,X_test)
  #'
  #' y_hat <- pcr_test$prediction
  #' err_hat <- y_test - y_hat
  #'
  #' MSE_test <- 1/length(y_test)*t(err_hat)%*%err_hat
  #'
  
  if (pcr_mod$scaling == "center") {
    X_sc <- t(t(X) - pcr_mod$center)  
  } else if (pcr_mod$scaling == "standardize") {
    X_sc <- t((t(X) - pcr_mod$center)/pcr_mod$StdDev)
  } else {
    X_sc <- X
  }
  
  pcr_coefs <- pcr_mod$coefficients
  if (length(dim(X))>1){
    n <- length(X[,1])
    y_pred <- cbind(rep(1,n),X_sc)%*%pcr_coefs
    
  } else {
    y_pred <- t(c(1,X))%*%pcr_coefs 
  }
  
  out <- list("prediction"=y_pred)
  return(out)
}


pcreg_CV <- function(y,X,scaling,k){
  #' Fit Principal Component Regression via Cross-Validation
  #' 
  #' Extracts M principal components of a scaled data matrix X by the 
  #' singular value decomposition and regresses y on these principal components
  #' by OLS. Here the M is chosen by Cross-Validation using k partitions of X under
  #' MSE - loss.
  #' 
  #' @param y A numeric vector of dim n x 1 to be projected onto principal components 
  #' @param X A matrix of column vectors of dim n x p.  
  #' @param scaling optional parameter either "center" for mean centering the column vectors X
  #' or "standardized" for standardizing the column vectors X
  #' @param k A parameter indicating how many subgroups should be used for
  #' Cross-Validation. If k=n, this amounts to LOOCV.
  #'
  #' @return Returns an list, containing the coefficients, fitted values, 
  #' the MSE for the training data, the scaling applied to X, the 
  #' mean and standard Deviations of the column vectors in X and also the 
  #' Cross-Validation MSE and the number of components used in the final model.    
  #' 
  #' @examples
  #'  set.seed(4519)
  #'
  #' N <- 100
  #' p <- 50
  #'
  #' X_train <- matrix(rnorm(N*p),N,p)
  #' y_train <- 2*rowMeans(X_train[,1:5])+rnorm(1.5)
  #'
  #' pcr_mod10fCV <- pcreg_CV(y_train,X_train,"standardize",10)
  #'
  n <-length(y)
  
  if (n %% k > 0){
    return("k not a multiple of number of observations")
    break
  }
  
  n_k <- n/k
  p <- min(c(n-n_k-1,length(X[1,])))
  
  kp_grid <- expand.grid(1:k,1:p)
  kp_length <- length(kp_grid[,1])
  MSPE_kf <- numeric(kp_length)
  
  for (kp_i in 1:kp_length) {
    
    seq_nk <- seq((kp_grid[kp_i,1]-1)*n_k+1,n_k*kp_grid[kp_i,1],1)
    X_test <- X[seq_nk,]
    y_test <- y[seq_nk]
    X_train <- X[-seq_nk,]
    y_train <- y[-seq_nk]
    
    train_mod <- pcreg_fit(y_train,X_train,kp_grid[kp_i,2],scaling)
    yhat_test <- pcreg_pred(train_mod,X_test)$prediction
    
    MSPE_kf[kp_i]=1/n_k*(t(y_test-yhat_test)%*%(y_test-yhat_test)) 
  }
  
  MSPE_CV <- colMeans(matrix(MSPE_kf,k,p,byrow=TRUE))
  
  M_CV <- which.min(MSPE_CV)
  mod_CV <- pcreg_fit(y,X,M_CV,scaling)
  out <- append(mod_CV,list("MSPE_CV"=MSPE_CV[M_CV],"numbergroups"=k,"numbercomp"=M_CV))
  return(out)
  
}

spcreg_fit <- function(y,X,K){
  
  X_dm <- scale(X,center=TRUE,scale=FALSE) # demeaned
  p <- length(X[1,])
  N <- length(y)
  
  rhohat <- cor(X_dm,y)
  thresh <- seq(0,0.99,0.01)
  
  llikeratio_stat <- numeric(length(thresh))
  
  beta0hat <- mean(y)
  sigma2tilde <- 1/N*t(y-beta0hat)%*%(y-beta0hat)
  llike_H0 <- -N/2*(log(2*pi)+log(sigma2tilde))-N/2
  n_variables_last_Xthresh<-0
  
  for (t in 1:length(thresh)){
    
    if (sum(abs(rhohat) >= thresh[t])<=K){ #i.e no rho_j above threshold
      llikeratio_stat[t] <- 0
    } else{
      
      Xthresh <- X_dm[,which(abs(rhohat) >= thresh[t])]
      
      if (dim(Xthresh)[2]==n_variables_last_Xthresh){
        llikeratio_stat[t]<- llikeratio_stat[t-1]
      }else{
        Xthresh_svd <- svd(Xthresh, nu=K+1, nv=0)
        
        Uthresh <- Xthresh_svd$u[,1:K]
        
        beta <- t(Uthresh)%*%y
        sigma2hat <- 1/N*t(y-Uthresh%*%beta-beta0hat)%*%(y-Uthresh%*%beta-beta0hat)
        
        llikeratio_stat[t] <- 2*(-N/2*(log(2*pi)+log(sigma2hat))-N/2-llike_H0)
        n_variables_last_Xthresh <- dim(Xthresh)[2]
      }
      
    }
  }
  
  threshstar_ind <- which(llikeratio_stat==max(llikeratio_stat))
  threshstar <- min(thresh[threshstar_ind])
  Xthreshstar_select <- which(abs(rhohat) >= threshstar)
  
  Xthreshstar <- X_dm[,Xthreshstar_select]
  Xthreshstar_svd <- svd(Xthreshstar, nu=K+1, nv=K+1)
  
  Uthreshstar <- Xthreshstar_svd$u[,1:K]
  Dthreshstar <- Xthreshstar_svd$d[1:K]
  Vthreshstar <- Xthreshstar_svd$v[,1:K]
  
  gammastar <- t(Uthreshstar)%*%y
  Wthreshstar <- Vthreshstar%*%(diag(K)*c(1/Dthreshstar))
  betastarhat <- Wthreshstar%*%gammastar
  
  spcreg_coefs <- c(beta0hat,rep(0,p))
  
  for (i in seq_along(Xthreshstar_select)) {
    spcreg_coefs[Xthreshstar_select[i]+1] <- betastarhat[i]
  }
  
  y_hat <- cbind(rep(1,N),X_dm)%*%spcreg_coefs
  err_hat <-y-y_hat
  
  MSE_train <- 1/N*t(err_hat)%*%err_hat
  
  out <- list("coefficients"=spcreg_coefs,"fittedvalues"=y_hat,"trainMSE"=MSE_train,
              "threshold_index"=Xthreshstar_select,"center"=colMeans(X))
  return(out)
  
}

spcreg_pred <- function(spcreg_mod,X){
  
  X_dm <- t(t(X) - spcreg_mod$center)  
  spcreg_coefs <- spcreg_mod$coefficients
  
  if (length(dim(X))>1){
    n <- length(X[,1])
    y_pred <- cbind(rep(1,n),X_dm)%*%spcreg_coefs
    
  } else {
    y_pred <- t(c(1,X[spcreg_mod$threshold_index]))%*%pcr_coefs 
  }
  
  out <- list("prediction"=y_pred)
  return(out)
  
  
}
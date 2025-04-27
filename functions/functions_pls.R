plsreg_fit <- function(y,X,K,scaling="center") {
  #' Fit Partial Least Squares Regression
  #'
  #' Finds the Partial Least Squares Regression estimators for a vector of outcomes y and a matrix of predictors X.
  #' This is known as the PLS1 algorithm and only works for a regression model with one response variable y.
  #'  
  #'
  #' @param y A numeric vector of dim n x 1
  #' @param X A matrix of column vectors of dim n x p.  
  #' @param K number of partial least squares components to be used in regression
  #' @param scaling whether to standardize or center the variables. Default is center. Values can be "center" or "standardize".
  #'
  #' @return Returns an list, containing the coefficients, fitted values, 
  #' the MSE for the training data and also the 
  #' mean and standard Deviations of the column vectors in X.   
  #' 
  #' @examples
  #'  set.seed(4519)
  #'
  #' X_train <- rmvnorm(100,mean=rep(2,10),sigma=diag(10))
  #' y_train <- 2*rowMeans(X_train[,1:5])+rnorm(1.5)
  #'
  #' pls_train <- plsreg_fit(y_train,X_train,2,"standardize")
  #'
  if (scaling == "standardize") {
    X_pls <- scale(X,center=TRUE,scale=TRUE)   
  } else{
    X_pls <- scale(X,center=TRUE,scale=FALSE)
  }

  U<-scale(y,center = TRUE, scale = FALSE)
  
  N<-nrow(X)
  P<-ncol(X)
  
  y_pred<-vector(mode = "numeric",length = N)
  y_pred[]<-mean(y)
  W<-matrix(0,nrow = P,ncol = K)
  Z<-matrix(0, nrow = N, ncol = K)
  Q<-vector(mode = "numeric",length = K)
  Psi<-matrix(0,nrow = P,ncol = K)
  for (k in (1:K)) {
    W[,k]<- crossprod(X_pls,U) # calculate the univariate regression coefficient for P variables in X
    W[,k]<- W[,k]/ sqrt(as.numeric(crossprod(W[,k]))) # standardize the cross product of the weight coefficients to be 1 (following the maximazation criterion in Friedman et al.).
    Z[,k]<-X_pls%*%W[,k]
    Psi[,k]<-crossprod(X_pls,Z[,k])/as.numeric(crossprod(Z[,k]))
    Q[k]<- crossprod(U,Z[,k])/as.numeric(crossprod(Z[,k]))
    y_pred <- y_pred+Q[k]*Z[,k]
    X_pls<-X_pls-tcrossprod(Z[,k],Psi[,k])
    U<- U-Q[k]*Z[,k]
  }
  
  # Calculating beta coefficients for X directly
  if (K==1){
    R<-W
    Beta <- R*as.numeric(Q)
  }else{
    PsiW_inv<-solve(crossprod(Psi,W))
    R<-W%*%PsiW_inv
    Beta<-R%*% Q
  }
  
  
  pls_coefs<-c(mean(y),Beta)
  
  err_pred<-y-y_pred
  MSE_pls<-1/N*as.numeric(crossprod(err_pred))
  out<-list("coefficients"=pls_coefs,"fittedvalues"=y_pred,"residuals"=err_pred,
            "trainMSE"=MSE_pls, "center"=colMeans(X),"StdDev"=apply(X,2,sd),
            "scaling"=scaling)
  
  
  return(out)
}

plsreg_pred <- function(pls_mod,X){
  #' Prediction for Partial Least Squares Regression
  #'
  #' uses the object returned by plsreg_fit() to make a prediction for a
  #' test data set X.
  #' 
  #' @param pcr_mod object returned by plsreg_fit()
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
  #' pls_train <- plsreg_fit(y_train,X_train,2)
  #'
  #' X_test <- matrix(rnorm(N*p),N,p)
  #' y_test <- 2*rowMeans(X_train[,1:5])+rnorm(1.5)
  #'
  #' pls_test <- pls_pred(pls_train,X_test)
  #'
  #' y_hat <- pls_test$prediction
  #' err_hat <- y_test - y_hat
  #'
  #' MSE_test <- 1/length(y_test)*t(err_hat)%*%err_hat
  #'
  if (!is.vector(X)){
    X_sc <- sweep(X, MARGIN = 2, STATS = pls_mod$center, FUN = "-")
    if (pls_mod$scaling == "standardize") {
      X_sc <- sweep(X_sc, MARGIN = 2, STATS = pls_mod$StdDev, FUN = "/")
    }
  }else{
    X_sc <- (t(X) - pls_mod$center)
    if (pls_mod$scaling == "standardize") {
      X_sc <- t(t(X_sc)/pls_mod$StdDev)
    }
  }
  
  pls_coefs <- pls_mod$coefficients
  if (!is.vector(X)){
    n <- dim(X)[1]
    ones_vector<-rep(1,n)
    y_pred <- cbind(ones_vector,X_sc)%*%pls_coefs
    
  } else {
    y_pred <- t(c(1,X))%*%pls_coefs 
  }
  
  out <- list("prediction"=y_pred)
  return(out)
}
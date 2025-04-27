
gene_simulation <- function(N,p,p1,p2,X1_shift,X2_shift,beta,sigma2_0,sigma2_1){
  
  X <- matrix(rnorm(N*p,mean=0,sd=sqrt(sigma2_0)),N,p)
  
  factor1 <-c(rep(X1_shift[1],N/2),
              rep(X1_shift[2],N/2))
  
  factor2 <-c(rep(X2_shift[1],N/4),
              rep(X2_shift[2],N/4),
              rep(X2_shift[1],N/4),
              rep(X2_shift[2],N/4))
  
  X[,1:p1] <- X[,1:p1]+factor1
  X[,(p1+1):(p2+p1)] <- X[,(p1+1):(p2+p1)]+factor2
  
  y <- beta*factor1+rnorm(N,mean=0,sd=sqrt(sigma2_1))
  
  out <- list("y"=y,"X"=X)
  return(out)
}

factor_simulation <- function(N,p,p1,p2,factor1,factor2,beta,sigma2_0,sigma2_1) {
  
  X <- matrix(rnorm(N*p,mean=0,sd=sqrt(sigma2_0)),N,p)
  
  X[,1:p1] <- X[,1:p1]+factor1
  X[,(p1+1):(p2+p1)] <- X[,(p1+1):(p2+p1)]+factor2
  
  Y <- beta*factor1 + rnorm(N,mean=0,sd=sqrt(sigma2_1))
  
  out <- list("y"=Y,"X"=X)
  return(out)
}

gene_simulation2 <- function(N,p,p1,p2,X1_shift,X2_shift,beta,sigma2_0,sigma2_1){
  
  X <- matrix(rnorm(N*p,mean=0,sd=sqrt(sigma2_0)),N,p)
  
  factor1 <-c(rep(X1_shift[1],N/2),
              rep(X1_shift[2],N/2))
  factor2 <-c(rep(X2_shift[1],N/4),
              rep(X2_shift[2],N/4),
              rep(X2_shift[1],N/4),
              rep(X2_shift[2],N/4))
  
  X[,1:(p1/2)] <- X[,1:(p1/2)]+ factor1 + factor2
  X[,(p1/2+1):p1] <- X[,(p1/2+1):p1]+ factor1 - factor2
  
  y <- beta*factor1 + beta/8*factor2 + rnorm(N,mean=0,sd=sqrt(sigma2_1))
  
  out <- list("y"=y,"X"=X)
  return(out)
}

compute_summary_stats <- function(x) {
  
  name_of_X <- deparse(substitute(x))
  
  mean_value <- mean(x)
  variance_value <- var(x)
  median_value <- median(x)
  percentile_25 <- quantile(x, 0.25)
  percentile_75 <- quantile(x, 0.75)
  
  x_result_table <- data.frame(
    Mean = mean_value,
    Variance = variance_value,
    `25th Percentile` = percentile_25,
    Median = median_value,
    `75th Percentile` = percentile_75
  )
  
  row.names(x_result_table) <- NULL
  
  cat("Summary statistics for", name_of_X, ":\n")
  
  print(x_result_table)
  return(x_result_table)
  
}

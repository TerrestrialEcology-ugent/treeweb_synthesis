##########################

# helper functions for R script

# on TREEWEB synthesis

##############

#compute bayesian R2 for single function
bayes_R2 <- function(obs,pred,summary = TRUE){
  e <- -1 * sweep(pred,2,obs)
  var_ypred <- apply(pred,1,var)
  var_e <- apply(e,1,var)  
  if(summary){
    return(quantile(var_ypred / (var_ypred + var_e),probs=c(0.1,0.5,0.9)))
  }
  else{
    return(var_ypred / (var_ypred + var_e))    
  }
}
#apply this across the functions
#sapply(1:K,function(k) bayes_R2(dat_std[,(k+1)],ypred[,,k]))


#function to generate positive definite matrix (covariance)
Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}
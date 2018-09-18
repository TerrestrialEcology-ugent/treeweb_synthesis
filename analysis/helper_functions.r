##########################

# helper functions for R script

# on TREEWEB synthesis

##############

##### compute bayesian R2 for single function
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


#### function to generate positive definite matrix (covariance)
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

##### function to help initialization
init_fn <- function(){
  list(L_omega = chol(Posdef(K)))
}

######## function to plot ternary graphs
make_ind_gg <- function(predicted, fun = "Predation",fragm="Low", pal=viridis(10)){
  subs <- subset(predicted,Function==fun & Fragm==fragm)
  out <-   ggtern(subs,aes(fsyl,qrob,qrub))+
    theme_bw()+
    theme_nomask() +
    #facet_wrap(~Function) +
    
    geom_tri_tern(bins=4,fun=mean,aes(value=Med_real,fill=..stat..)) +
    stat_tri_tern(bins=4,fun=mean,geom="text",aes(value=Med_real,
                                                  label=sprintf("%.2f",..stat..)),
                  color="darkorange3",size=3,centroid=TRUE) +
    annotate(geom = "text",x=0.5,y=0.9,z=-0.3,label=fun) +
    scale_fill_gradientn(colours = pal) +
    theme(legend.position = "none")
  return(out)
}
################
# script to run the 10-fold
# cross-validation on the 2 models
# author: Lionel Hertzog
# date: XXXXX
##########


# load libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores())
library(rstudioapi) # for easy setwd

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# loading data
dat_std <- read.csv("../data/synthesis_responsedata_std.csv",sep=" ")
div_dd <- read.csv("../data/synthesis_expldata_raw.csv",sep=" ")

# load functions
source("kfold_functions.r")

# generate model input data
X_1 <- model.matrix(~ specrich_std * fragm_std,div_dd)
#generating data for stan fit
#dimension variables
N <- nrow(dat_std) #row number
K <- ncol(dat_std) - 1 #number of columns in Y
J <- ncol(X_1) #number of columns in X
Fr <- max(div_dd$id_frag)
F_id <- div_dd$id_frag
#define the number folds
n_fold <- 10
#create 10 blocks of data
hh <- sample(1:N,size = N,replace = FALSE)
holdout_10 <- matrix(0,nrow=N,ncol=n_fold)
id <- seq(1,50,by=5)
for(i in 1:n_fold){
  holdout_10[hh[id[i]:(id[i] + 4)],i] <- 1
}
#add the non-attributed points to the last fold (nasty debug)
holdout_10[which(apply(holdout_10,1,sum) == 0),10] <- 1
#turn into a list
holdout_10 <- split(holdout_10,rep(1:ncol(holdout_10),each=nrow(holdout_10)))

#define the data input
data_m <- list(N=N,K=K,J=J,X=X_1,F=Fr,F_id=F_id,y=dat_std[,-c(1)])
#create a list of data
data_l <- rep(list(data_m),10)
#add the holdout index to it
for(i in 1:10) data_l[[i]]$holdout <- holdout_10[[i]]

# run the model
ss <- stan_kfold(file="../model/multivariate_fragment_kfold.stan",data_l,chains=4,cores=3,init=init_fn) # this take some time
ee <- extract_log_lik_K(ss,holdout_10)
kk <- kfold(ee) 

# similar manipulations for the second model
div_dd$Eve <- with(div_dd,rba_fsyl * rba_qrob + rba_fsyl * rba_qrub + rba_qrob * rba_qrub)
div_dd$fragmCat <- cut(div_dd$fragm_std,c(-2,-0.78,0.70,2.5),labels = c("Low","Medium","High"))
X_3 <- model.matrix(~ -1 + total_ba_std + (rba_fsyl + rba_qrob + rba_qrub + Eve) * fragmCat,div_dd)

N <- nrow(dat_std)
K <- ncol(dat_std) - 1
J <- ncol(X_3) #number of columns in X
Fr <- max(div_dd$id_frag)
F_id <- div_dd$id_frag

# put the model input data together
data_m <- list(N=N,K=K,J=J,X=X_3,F=Fr,F_id=F_id,y=dat_std[,-c(1)])
#create a list of data
data_l <- rep(list(data_m),10)
#add the holdout index to it
for(i in 1:10) data_l[[i]]$holdout <- holdout_10[[i]]

ss2 <- stan_kfold(file="multivariate_fragment_kfold.stan",data_l,chains=4,cores=3,init=init_fn) # this take some time
ee2 <- extract_log_lik_K(ss2,holdout_10)

#compare the two models
k1 <- kfold(ee) #elpd -1932, se 106
k2 <- kfold(ee2) #elpd -2014, se 125
diff <- k1$pointwise - k2$pointwise
elpd_diff <- sum(diff)
se_diff <- sqrt(53) * sd(diff)
#diff: 82, SE: 55, so slight evidence for better second model but large uncertainties
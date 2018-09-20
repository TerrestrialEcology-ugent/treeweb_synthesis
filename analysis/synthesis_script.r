###########################

# R script for reproducing results

# of Hertzog et al ....

# author: Lionel Hertzog

# date: XXXXXX

###############

### loading libraries ###

library(MASS) # for mvrnorm
library(plyr)
library(reshape2)
library(tidyverse)
library(rstan)
library(corrr)
library(viridis)
library(gridExtra)
library(ggtern)
library(rstudioapi) # for easy setwd

### set working dorectory ###

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### loading data ###

dat_std <- read.csv("../data/synthesis_responsedata_std.csv",sep=" ")
div_dd <- read.csv("../data/synthesis_expldata_raw.csv",sep=" ")
impA <- read.csv("../data/synthesis_importance_scores.csv",sep=" ")

### loading helper functions ###

source("helper_functions.r")

### some helpful objects ###

nice_name <- c("C stock","Base sat.","pH","C:N","Phosphorus","Decomp.","Biomass","Veg. cover","Veg. div.",
               "Light trans.","LAI","Germination","Seed growth","Arth. div.","Herbivory","Predation","Spider fit.",
               "Spider diet","Bird SMI","Breed. succ.","Bird div.","Egg vol.","Egg bact.","Egg IgY")

### Model 1. Richness - fragmentation ###

# make the model matrix 
X_1 <- model.matrix(~ specrich_std * fragm_std,div_dd)

# make the matrix over which to derive predictions
X_1pred <- expand.grid(Int = 1, specrich = c(-0.96,0.45,1.88), fragm = c(-1.69, 0, 2.42))
X_1pred$interaction <- with(X_1pred, specrich * fragm)

# input data to the Stan model
N <- nrow(dat_std) # row number
N_pred <- nrow(X_1pred)
K <- ncol(dat_std) - 1 # number of columns in response matrix
J <- ncol(X_1) # number of columns in explonatory matrix
Fr <- max(div_dd$id_frag) # number of fragments
F_id <- div_dd$id_frag # fragment id for each plot

data_m <- list(N=N,N_pred = N_pred,K=K,J=J,X=X_1,X_pred=X_1pred,F=Fr,F_id=F_id,y=dat_std[,-c(1)],
               direction_manager = ifelse(impA$Direction=="maximize",1,-1),importance_manager=impA$Importance_weight)

# fit the model
m_rich <- stan(file="../model/multivariate_fragment.stan",data = data_m,init=init_fn) 

# model diagnostics
stan_diag(m_rich)
stan_rhat(m_rich, pars = "beta")
stan_ess(m_rich, pars = "beta")

# explore model
# library(shinystan)
# launch_shinystan(m_rich)

# extract posterior drws
post_m_rich <- rstan::extract(m_rich)
# save this for later use
# saveRDS(post_m_rich,"Stan_model_rich_posterior.rds")
# to reload the poseterior draws
# post_m_rich <- readRDS("post_m_rich.rds")


### Model 2. Diversity - interaction ###

# create 2 new variables
div_dd$Eve <- with(div_dd,rba_fsyl * rba_qrob + rba_fsyl * rba_qrub + rba_qrob * rba_qrub) # compute tree eveness
div_dd$fragmCat <- cut(div_dd$fragm_std,c(-2,-0.78,0.70,2.5),labels = c("Low","Medium","High")) # disretize the fragmentation intensity gradient

# the model matrices
X_3 <- model.matrix(~ -1 + total_ba_std + (rba_fsyl + rba_qrob + rba_qrub + Eve) * fragmCat,div_dd)
X_3pred <- expand.grid(total_ba=0,fsyl=seq(0,1,length=10),qrob=seq(0,1,length=10),qrub=seq(0,1,length=10))
# in the matrix to derive predicted values only keep sensible values (ie relative abundance that sum to 1)
tot <- rowSums(X_3pred)
X_3pred <- X_3pred[tot == 1,]

# add the fragmentation categories
X_3pred <- rbind(X_3pred,X_3pred,X_3pred)
X_3pred$fragmCatLow <- rep(c(1,0),times=c(46,46*2))
X_3pred$fragmCatMedium <- rep(c(0,1,0),each=46)
X_3pred$fragmCatHigh <- rep(c(0,1),times=c(46*2,46))
X_3pred$Eve <- with(X_3pred,fsyl*qrob+fsyl*qrub+qrob*qrub)
# add the interaction between fragmentation and the relative abundance
X_3pred <- cbind(X_3pred,with(X_3pred,fragmCatMedium*fsyl),with(X_3pred,fragmCatHigh*fsyl),
                 with(X_3pred,fragmCatMedium*qrob),with(X_3pred,qrob*fragmCatHigh),
                 with(X_3pred,qrub*fragmCatMedium),with(X_3pred,qrub*fragmCatHigh),
                 with(X_3pred,Eve*fragmCatMedium),with(X_3pred,Eve*fragmCatHigh))


# data to pass to the model
N <- nrow(dat_std)
N_pred <- nrow(X_3pred)
K <- ncol(dat_std) - 1
J <- ncol(X_3) #number of columns in X
Fr <- max(div_dd$id_frag)
F_id <- div_dd$id_frag

data_m_div <- list(N=N,N_pred = N_pred,K=K,J=J,X=X_3,X_pred=X_3pred,F=Fr,F_id=F_id,y=dat_std[,-c(1)],
               direction_manager = ifelse(impA$Direction == "maximize",1,-1),importance_manager=impA$Importance_weight)


#fit the model
m_div <- stan(file="../model/multivariate_fragment.stan",data=data_m_div,init=init_fn,control=list(max_treedepth=20)) #works and fit in 1min

# model diagnostics
stan_diag(m_div)
stan_rhat(m_div, pars = "beta")
stan_ess(m_div, pars = "beta")

# explore model
# library(shinystan)
# launch_shinystan(m_div)

# extract posterior drws
post_m_div <- rstan::extract(m_div)
# save this for later use
# saveRDS(post_m_div,"Stan_model_div_posterior.rds")
# to reload the saved posterior draws
# post_m_div <- readRDS("post_m_div.rds")


### Figures ###

# figure 1: slopes from richness model

bb <- post_m_rich$beta
bbm <- adply(bb,c(2,3),quantile,probs=c(0.05,0.5,0.95)) # get quantile per function and coefficient
bbm$Var <- names(dat_std)[-c(1)][bbm$X2] # add function names
bbm$Var <- factor(bbm$Var,levels = c(names(dat_std)[-1][1:15],names(dat_std)[-1][17:21],"Predation",names(dat_std)[-1][22:24])) # some re-ordering
bbm$Coef <- rep(c("Int","Diversity","Fragmentation","Interaction"),times=K)
names(bbm)[3:5] <- c("LCI","Med","UCI")

gg_fig1 <- ggplot(subset(bbm,Coef!="Int"),aes(x=Var,y=Med,ymin=LCI,ymax=UCI)) +
  geom_point(show.legend = FALSE) +
  geom_linerange(show.legend = FALSE) +
  facet_grid(~Coef) +
  coord_flip() +
  geom_hline(yintercept = 0,linetype="dashed") +
  annotate(geom = "rect",xmin=0,xmax=6.5,ymin=-Inf,ymax=Inf,fill="saddlebrown",alpha=0.3) +
  annotate(geom = "rect",xmin=6.5,xmax=13.5,ymin=-Inf,ymax=Inf,fill="springgreen2",alpha=0.3) +
  annotate(geom = "rect",xmin=13.5,xmax=17.5,ymin=-Inf,ymax=Inf,fill="darkgoldenrod2",alpha=0.3) +
  annotate(geom = "rect",xmin=17.5,xmax=25,ymin=-Inf,ymax=Inf,fill="tomato3",alpha=0.3) +
  theme(panel.background = element_rect(fill=NA),strip.background = element_rect(fill=NA,color="black")) +
  scale_x_discrete(labels=nice_name[c(1:15,17:21,16,22:24)]) +
  labs(x="",y="Estimated slopes with 90% credible intervals")

ggsave("01_slopes_diversity.png",gg_fig1,width=10,height=6,units="in",dpi=100)

# figure 2: desirability 

# gather desirability scores from the richness - fragmentation model
d_f <- post_m_rich$desirability_manager # extract the posterior draws
dd_f <- adply(d_f,2,quantile,probs=c(0.25,0.5,0.75)) # compute the quantile
# some data reshaping
names(dd_f)[2:4] <- c("LCI","Med","UCI")
dd_f$div <- X_1pred[dd_f$X1,"specrich"]
dd_f$fragm <- X_1pred[dd_f$X1,"fragm"]
dd_f$fragmF <- factor(dd_f$fragm,labels=c("Low","Medium","High"))
dd_f$divF <- factor(dd_f$div,labels=c(1,2,3))

ggd_div <- ggplot(dd_f,aes(x=div,y=Med,group=fragmF)) +
  geom_ribbon(aes(ymin=LCI,ymax=UCI,fill=fragmF),alpha=0.2) +
  geom_path(aes(color=fragmF)) +
  scale_x_continuous(breaks=c(-0.96,0.45,1.88),labels = c(1,2,3)) +
  labs(fill="Fragmentation level",color="Fragmentation level",
       x="Tree richness",y="Desirability score with 50% credible intervals") +
  theme_bw()

# gather desirability scores from the diversity - interaction model

ddf <- adply(post_m_div$desirability_manager,2,quantile,probs=0.5) #only get the median this time
ddf$fsyl <- X_3pred[ddf$X1,2]
ddf$qrob <- X_3pred[ddf$X1,3]
ddf$qrub <- X_3pred[ddf$X1,4]
ddf$Fragm <- rep(c("Low","Medium","High"),each=46)
names(ddf)[2] <- "Med"
ddf$Med100 <- ddf$Med * 100

gg_fl <- ggtern(subset(ddf,Fragm=="Low"),aes(fsyl,qrob,qrub))+
  theme_bw()+
  geom_tri_tern(bins=4,fun=mean,aes(value=Med,fill=..stat..)) +
  stat_tri_tern(bins=4,fun=mean,geom="text",aes(value=Med,
                                                label=sprintf("%.2f",..stat..)),
                color="lightgrey",centroid=TRUE) +
  scale_fill_viridis(option="C") +
  theme(legend.position = "none")

gg_fh <- ggtern(subset(ddf,Fragm=="High"),aes(fsyl,qrob,qrub))+
  theme_bw()+
  geom_tri_tern(bins=4,fun=mean,aes(value=Med,fill=..stat..)) +
  stat_tri_tern(bins=4,fun=mean,geom="text",aes(value=Med,
                                                label=sprintf("%.2f",..stat..)),
                color="lightgrey",centroid=TRUE) +
  scale_fill_viridis(option="C") +
  theme(legend.position = "none")

gg_d <- ggtern::arrangeGrob(grobs=list(gg_fl,gg_fh),ncol=2,
                            top="")
ggd <- ggtern::grid.arrange(ggd_div,gg_d)

#ggsave("02_desirability.png",ggd,width=20,height=20,units="cm",dpi=100) 

# figure 3: function - function correlations

# plot index for low and high fragmentation intensity
fragmCatLow <- which(div_dd$fragm_std < quantile(div_dd$fragm_std,0.25))
fragmCatHigh <- which(div_dd$fragm_std > quantile(div_dd$fragm_std,0.75))

# predicted correlation low fragmentation intensity for model 1
m1_corl <- apply(apply(post_m_rich$y_post[,fragmCatLow,],1,cor),1,median) # grab median correlation
# some data re-shaping
m1_ddl <- data.frame(Fun1 = rep(nice_name,each=24),Fun2 = rep(nice_name,times=24),Cor=m1_corl)
m1_ml <- dcast(m1_ddl, Fun1 ~ Fun2, value.var = "Cor")
# plot
n1 <- network_df(as_cordf(m1_ml[,-1]),min_cor = 0.3,legend=FALSE,repel=FALSE) + 
  labs(title="") + 
  theme(plot.margin = unit(c(0,0,0,0),"line"),plot.background = element_rect(fill="grey70"))

# same stuff for high fragmentation for model 1
m1_corh <- apply(apply(post_m_rich$y_post[,fragmCatHigh,],1,cor),1,median)
m1_ddh <- data.frame(Fun1 = rep(nice_name,each=24),Fun2 = rep(nice_name,times=24),Cor=m1_corh)
m1_mh <- dcast(m1_ddh, Fun1 ~ Fun2, value.var = "Cor")
n2 <- network_df(as_cordf(m1_mh[,-1]),min_cor = 0.3,legend=FALSE,repel=FALSE) + labs(title="") + theme(plot.margin = unit(c(0,0,0,0),"line"),
                                                                                                       plot.background = element_rect(fill="grey70"))
# same stuff for low fragmentation for model 2
m3_corl <- apply(apply(post_m_div$y_post[,fragmCatLow,],1,cor),1,median)
m3_ddl <- data.frame(Fun1 = rep(nice_name,each=24),Fun2 = rep(nice_name,times=24),Cor=m3_corl)
m3_ml <- dcast(m3_ddl, Fun1 ~ Fun2, value.var = "Cor")
n3 <- network_df(as_cordf(m3_ml[,-1]),min_cor = 0.3,legend=FALSE,repel=FALSE) + labs(title="") + theme(plot.margin = unit(c(0,0,0,0),"line"),
                                                                                                       plot.background = element_rect(fill="grey70"))
# same stuff for high fragmentation for model 2
m3_corh <- apply(apply(post_m_div$y_post[,fragmCatHigh,],1,cor),1,median)
m3_ddh <- data.frame(Fun1 = rep(nice_name,each=24),Fun2 = rep(nice_name,times=24),Cor=m3_corh)
m3_mh <- dcast(m3_ddh, Fun1 ~ Fun2, value.var = "Cor")
n4 <- network_df(as_cordf(m3_mh[,-1]),min_cor = 0.3,legend=FALSE,repel=FALSE) + labs(title="") + theme(plot.margin = unit(c(0,0,0,0),"line"),
                                                                                                       plot.background = element_rect(fill="grey70"))
# combine all panels
nn_all <- grid.arrange(n1,n2,n3,n4,ncol=2,top="Low fragmentation                                                       High fragmentation",
                       left="Composition model                                                  Richness model")
ggsave("../Output/Figures/03_trade_offs.png",nn_all,width=30,height=30,units="cm")
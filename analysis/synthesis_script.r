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


######## Figure ternary plots for all functions #########
predd <- post_m3$y_pred
predd <- adply(predd,c(2,3),quantile,probs=0.5)
predd$fsyl <- X_3pred[predd$X1,2]
predd$qrob <- X_3pred[predd$X1,3]
predd$qrub <- X_3pred[predd$X1,4]
predd$Function <- names(dat_std)[-1][predd$X2]
predd$Fragm <- rep(rep(c("Low","Medium","High"),each=46),24)
names(predd)[3] <- "Med"

#transform back the predicted values to their original scale
#get mean and sd for each function
dat %>%
  gather(Function,value,-id_plot) %>%
  group_by(Function) %>%
  summarise(M=mean(value,na.rm=TRUE),SD=sd(value,na.rm=TRUE)) -> mean_sd

predd %>%
  left_join(mean_sd,by="Function") %>%
  mutate(Med_real = (Med * SD) + M) %>%
  select(Fragm,Function,fsyl:qrub,Med_real) -> predd2

#the function that does the job

#go through vegetation functions
veg_fun <- c("Biomass","Veg_div","Cover","GLI","LAI","P_germ","Seed_biom")
fragm_lvl <- c("Low","High")
ss <- expand.grid(fragm = fragm_lvl, fun = veg_fun)
tmp <- list()
for(f in veg_fun){
  for(ff in fragm_lvl){
    tmp[[length(tmp)+1]] <- make_ind_gg(predd2,f,ff)
  }
}


gg_a <- ggtern::arrangeGrob(grobs=tmp,ncol=2,nrow=7,top = "Effect of varying tree composition on vegetation functions\nfor low fragmentation (left) vs high fragmentation (right) intensity")
ggf <- ggtern::grid.arrange(gg_a)
ggsave("../Output/Figures/17v2_ternary_veg.png",ggf,width=7,height=18,units="in")
#now for soil functions
veg_fun <- c("C_stock","CN","P","BS","pH","Decomp")
fragm_lvl <- c("Low","High")
ss <- expand.grid(fragm = fragm_lvl, fun = veg_fun)
tmp <- list()
for(f in veg_fun){
  for(ff in fragm_lvl){
    tmp[[length(tmp)+1]] <- make_ind_gg(predd2,f,ff)
  }
}


gg_a <- ggtern::arrangeGrob(grobs=tmp,ncol=2,nrow=6,top = "Effect of varying tree composition on soil functions\nfor low fragmentation (left) vs high fragmentation (right) intensity")
ggf <- ggtern::grid.arrange(gg_a)
ggsave("../Output/Figures/16v2_ternary_soil.png",ggf,width=7,height=16,units="in")
#now for arthropods functions
veg_fun <- c("Arth_div","Herbivory","Predation","Fit_spider","Spider_diet")
fragm_lvl <- c("Low","High")
ss <- expand.grid(fragm = fragm_lvl, fun = veg_fun)
tmp <- list()
for(f in veg_fun){
  for(ff in fragm_lvl){
    tmp[[length(tmp)+1]] <- make_ind_gg(predd2,f,ff)
  }
}


gg_a <- ggtern::arrangeGrob(grobs=tmp,ncol=2,nrow=5,
                            top = "Effect of varying tree composition on arthropod functions\nfor low fragmentation (left) vs high fragmentation (right) intensity")
ggf <- ggtern::grid.arrange(gg_a)
ggsave("../Output/Figures/14v2_ternary_arth.png",ggf,width=7,height=16,units="in")

#this figure rotated
gg_a <- ggtern::arrangeGrob(grobs=tmp,ncol=5,nrow=2,as.table = FALSE,
                            top = "Effect of varying tree composition on arthropod functions\nfor low fragmentation (top) vs high fragmentation (bottom) intensity")
ggf <- ggtern::grid.arrange(gg_a)
ggsave("../Output/Figures/ternary_arth_r.png",ggf,width=16,height=7,units="in")


#now for bird functions
veg_fun <- c("Bird_smi","Breed_succ","Bird_div","Egg_vol","Egg_bact","Egg_IgY")
fragm_lvl <- c("Low","High")
ss <- expand.grid(fragm = fragm_lvl, fun = veg_fun)
tmp <- list()
for(f in veg_fun){
  for(ff in fragm_lvl){
    tmp[[length(tmp)+1]] <- make_ind_gg(predd2,f,ff)
  }
}


gg_a <- ggtern::arrangeGrob(grobs=tmp,ncol=2,nrow=6,top = "Effect of varying tree composition on bird functions\nfor low fragmentation (left) vs high fragmentation (right) intensity")
ggf <- ggtern::grid.arrange(gg_a)
ggsave("../Output/Figures/15v2_ternary_bird.png",ggf,width=7,height=16,units="in")


########### figures desirability Kirwan #######
#look at the desirability
post_m3 <- readRDS("Stan_model_m3_posterior.rds")
ddf <- adply(post_m3$desirability_manager,2,quantile,probs=0.5)
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
  annotate(geom="text",x=100,y=0,z=100,label="a") +
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
ggsave("../Output/Figures/02_desirability.png",ggd,width=20,height=20,units="cm",dpi=100) 

ddn <- adply(post$desirability_natuur,2,quantile,probs=0.5)
ddn$fsyl <- X_4pred[ddn$X1,2]
ddn$qrob <- X_4pred[ddn$X1,3]
ddn$qrub <- X_4pred[ddn$X1,4]
ddn$Fragm <- rep(c("Low","Medium","High"),each=46)
names(ddn)[2] <- "Med"
ddn$Med100 <- ddn$Med * 100

gg_nl <- ggtern(subset(ddn,Fragm=="Low"),aes(fsyl,qrob,qrub))+
  theme_bw()+
  geom_tri_tern(bins=4,fun=mean,aes(value=Med,fill=..stat..)) +
  stat_tri_tern(bins=4,fun=mean,geom="text",aes(value=Med,
                                                label=sprintf("%.2f",..stat..)),
                color="lightgrey",centroid=TRUE) +
  scale_fill_viridis(option="C") +
  theme(legend.position = "none")

gg_nh <- ggtern(subset(ddn,Fragm=="High"),aes(fsyl,qrob,qrub))+
  theme_bw()+
  geom_tri_tern(bins=4,fun=mean,aes(value=Med,fill=..stat..)) +
  stat_tri_tern(bins=4,fun=mean,geom="text",aes(value=Med,
                                                label=sprintf("%.2f",..stat..)),
                color="lightgrey",centroid=TRUE) +
  scale_fill_viridis(option="C") +
  theme(legend.position = "none")

gg_d <- ggtern::arrangeGrob(grobs=list(gg_fl,gg_nl,gg_fh,gg_nh),ncol=2,
                            left="High fragmentation intensity          Low fragmentation intensity",
                            top="Desirability scores across variation in tree relative abundance\nfor foresters (left) or nature conservation (right)\nin the tree interaction model")
ggd <- ggtern::grid.arrange(gg_d)
ggsave("../Output/Figures/ternary_desirability_int_f.eps",ggd,width=20,height=20,units="cm",dpi=100) 


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

# plot desirability scores
d_f <- post_m1$desirability_manager
dd_f <- adply(d_f,2,quantile,probs=c(0.25,0.5,0.75))
names(dd_f)[2:4] <- c("LCI","Med","UCI")
dd_f$div <- X_1pred[dd_f$X1,"specrich"]
dd_f$fragm <- X_1pred[dd_f$X1,"fragm"]
dd_f$fragmF <- factor(dd_f$fragm,labels=c("High","Medium","Low"))
dd_f$divF <- factor(dd_f$div,labels=c(1,2,3))

ggd_div <- ggplot(dd_f,aes(x=div,y=Med,group=fragmF)) +
  geom_ribbon(aes(ymin=LCI,ymax=UCI,fill=fragmF),alpha=0.2) +
  geom_path(aes(color=fragmF)) +
  scale_x_continuous(breaks=c(-0.96,0.45,1.88),labels = c(1,2,3)) +
  labs(fill="Fragmentation level",color="Fragmentation level",
       x="Tree richness",y="Desirability score with 50% credible intervals") +
  theme_bw()

#save just that one
ggsave("../Output/Figures/02_model_forest_desirability.png",ggd_div,width=10,height=10,units="cm")


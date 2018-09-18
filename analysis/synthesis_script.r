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
div_dd <- read.csv("../Data/synthesis_expldata_raw.csv",sep=" ")

# loading helper functions

# first model
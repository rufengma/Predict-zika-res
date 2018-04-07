rm(list=ls())
setwd("C:/Users/Subho/Box Sync/Hunting Zika Virus with Machine Learning (Cary Institute)")
# setwd('/extdrive/Work/smajumd/Zika-codes')
library(plyr)
library(mi)
library(data.table)

prim = read.csv("Data/prim.csv", row.names = NULL)[,-1]
flavi = read.csv("Data/primates flavi.csv")[,-1]
mosquito = read.csv('Data/mosquito.csv')

# how much NA?
na.amount = apply(flavi, 2, FUN = function(x) length(which(is.na(x))))
na.df = data.frame(Variable=names(na.amount),
                   NAlength=as.numeric(na.amount))
na.df = na.df[order(na.df[,2], decreasing=T),]
row.names(na.df) = NULL

## take columns with more than 30 non-NA entries
which.in = with(na.df, Variable[which(NAlength < 310)])
flavi.small = flavi[,which(names(flavi) %in% which.in)]
# save(flavi.small, file='Data/flavi_small.Rda')

## number of NA entries for each sample
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))
flavi.small = flavi.small[-which(NA.amount==36),]

## impute
# step1: get data in shape
flavi.impute = missing_data.frame(flavi.small[-c(1:8, 45)])

# step2: change variable types
flavi.impute = change(flavi.impute, y='X6.1_DietBreadth', what='type', to='ordered-categorical')

set.seed(31052016)
# cl = makeCluster(detectCores()-1)
options(mc.cores=4)
system.time(flavi.imputed <- mi(flavi.impute, n.iter=2e3, parallel=TRUE))
save(flavi.imputed, file='../Zika-outputs/mi_model_nonlog_small.Rda')

# do same thing, but for continuous predictors only, and run chains longer
flavi.impute = missing_data.frame(flavi.small[,-c(1:9,13,15,26,27,45)])
set.seed(31052016)
options(mc.cores=4)
system.time(flavi.imputed <- mi(flavi.impute, n.iter=5e4, parallel=TRUE))
save(flavi.imputed, file=../Zika-outputs/mi_model_nonlog_big.Rda')


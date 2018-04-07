## Multi-response model for Mosquito-Flavivirus connections
## version 2.02: get out of sample scores
rm(list=ls())
setwd("D:/Study/My projects/Predict-zika-res/Code") ## change to your home directory
source("bmlml.R")

library(mice)
library(miceadds)
library(parallel)
library(caret)
library(pROC)

## load predictor data
mosquito = read.csv('mosquito.csv')

## process data, fill in NAs
na.X = which(is.na(mosquito), arr.ind=T)
mosquito[na.X] = 0
X1 = scale(as.matrix(mosquito[-c(1:3,24)]))

## load response data
virus = read.csv('FlaviData.csv')
mosq.matrix = read.csv('flavi_mos1234.csv')[1:174,]
mosq.matrix = mosq.matrix[,c(1,which(names(mosq.matrix) %in% paste(virus$code)))]

## process data, fill in NAs
Y1 = as.matrix(mosq.matrix[,-1])
na.Y = which(is.na(Y1), arr.ind=T)
Y1[na.Y] = 0
na.Y = which(Y1<4, arr.ind=T)
Y1[na.Y] = 0
one.Y = which(Y1==4, arr.ind=T)
Y1[one.Y] = 1

##############################################
# OOB scores
##############################################

N = nrow(X1)
D = ncol(X1)
L = ncol(Y1)

# full model with 1 latent factor
K = 1
r.hat = selectParams(Y1, X1, K=K, nU=1e2, maxit=5e2)
model.hat = bmlml.EM1(Y1, X1, K=K, r.tau=r.hat, maxit=5e2)
pred.mat = predict.bmlml(model.hat, X1)

# create stratified folds
set.seed(12102016)
pos.samples = as.numeric(which(rowSums(Y1)>0))
pos.folds = createFolds(as.numeric(pos.samples), k=10) # folds for positive samples
neg.samples = as.numeric(which(rowSums(Y1)==0))
neg.folds = createFolds(as.numeric(neg.samples), k=10) # folds for negative samples
idata.folds = lapply(1:10, function(i)
  c(pos.samples[pos.folds[[i]]], neg.samples[neg.folds[[i]]])) # plug together

# now train model and predict for each fold
for(i in 1:10){
  itest = idata.folds[[i]]
  
  r.hat.i = selectParams(Y1[-itest,], X1[-itest,], K=K, nU=1e2, maxit=5e2)
  model.hat.i = bmlml.EM1(Y1[-itest,], X1[-itest,], K=K, r.tau=r.hat.i, maxit=5e2)
  pred.mat[itest,] = predict.bmlml(model.hat.i, X1[itest,])
}

pred.df = data.frame(Species=mosquito$mosquitoSpecies,
                     NewWorld=with(mosquito, (northAmerica==1 | southAmerica==1)),
                     Detected=((1:N) %in% which(rowSums(Y1)>0)),
                     score=round(pred.mat[,L],3))

write.csv(pred.df, file='../Outputs/mosq scores_10fold.csv')

##############################################
# In-sample scores
##############################################

# full model with 1 latent factor
K = 1
r.hat = selectParams(Y1, scale(X1), K=K, nU=1e2, maxit=5e2)
model.hat = bmlml.EM1(Y1, scale(X1), K=K, r.tau=r.hat, maxit=5e2)
pred.mat = predict.bmlml(model.hat, scale(X1))

pred.df = data.frame(Species=mosquito$mosquitoSpecies,
                     NewWorld=with(mosquito, (northAmerica==1 | southAmerica==1)),
                     Detected=((1:N) %in% which(rowSums(Y1)>0)),
                     score=round(pred.mat[,L],3))

write.csv(pred.df, file='../Outputs/mosq scores.csv')

##############################################
# Validation
##############################################
in.score = read.csv('../Outputs/mosq scores.csv')
out.score = read.csv('../Outputs/mosq scores_10fold.csv')

vector.indicator = (rowSums(Y1)>0)

pdf('../Outputs/MultiModel_MosqFlavi_ROC.pdf', height=4, width=8)
par(mfrow=c(1,2))
plot.roc(vector.indicator, in.score$score, print.auc=T,
         main="In-sample")
plot.roc(vector.indicator, out.score$score, print.auc=T,
         main="10-fold CV")
par(mfrow=c(1,1))
dev.off()

##############################################
# Variable importance and table
##############################################
bar.data = abs(model.hat$W)
bar.names = names(mosquito[-c(1:3,24)])
bar.order = order(bar.data, decreasing=F)
defaultPar = par()

pdf('../Outputs/MultiModel_MosqFlavi_varimp.pdf', height=4, width=8)
par(las=2, mar=c(5,12,2,3),mfrow=c(1,1))
barplot(bar.data[bar.order], horiz=T,
        names.arg=bar.names[bar.order], cex.names=.7,
        main=paste("Average variable importance across factors"), cex.main=.9)
dev.off()

## mean of variables for low and high risk samples
pred.df = read.csv('../Outputs/mosq scores.csv')
score = pred.df$score
xl3 = which(score < quantile(score,.3))
xg7 = which(score > quantile(score,.7))
mean.mat = matrix(0, nrow=ncol(X1), ncol=2)
for(i in 1:ncol(X1)){
  ivar = as.numeric(paste(X1[,i]))
  mean.mat[i,] = c(mean(ivar[xl3]), mean(ivar[xg7]))
}

# put them all together
var.df = data.frame(Name=bar.names,
                    Importance=round(bar.data,2),
                    mean.low=round(mean.mat[,1],2),
                    mean.hi=round(mean.mat[,2],2))
var.df1 = var.df[order(bar.data, decreasing=T),]
write.csv(var.df1, file='../Outputs/variable means for low and high risk samples_mosq.csv', row.names=F)


############################################################
## End of file
############################################################
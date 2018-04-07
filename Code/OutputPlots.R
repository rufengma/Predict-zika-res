## Interact with outputs
rm(list=ls())
# setwd("C:/Users/IBM_ADMIN/Box Sync/Hunting Zika Virus with Machine Learning (Cary Institute)")
setwd("c:/Users/Subho/Box Sync/Hunting Zika Virus with Machine Learning (Cary Institute)")

load('Data/flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))
disease.matrix = read.csv('Data/flavi_prim01mid.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -1])

## histogram of mean scores
par(mfrow=c(1,2))

pred.df = read.csv("Validation/risk scores.csv")
hist(pred.df$mean.score, xlim=0:1, nclass=20,
     xlab="risk score", main="In-sample")
for(i in which(rowSums(Y1)>0)){
  abline(v=pred.df$mean.score[i], col=adjustcolor("red", alpha.f=.5))
}

pred.df = read.csv("Validation/risk scores_OOB.csv")
hist(pred.df$mean.score, xlim=0:1, nclass=20,
     xlab="risk score", main="Out-of-sample")
for(i in which(rowSums(Y1)>0)){
  abline(v=pred.df$mean.score[i], col=adjustcolor("red", alpha.f=.5))
}
par(mfrow=c(1,1))

## high risk species
pred.df1 = pred.df[order(pred.df$mean.score, decreasing=T),]
pred.df1[with(pred.df1, which(NewWorld & !Detected)),]

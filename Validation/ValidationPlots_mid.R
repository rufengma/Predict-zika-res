## SelectTuningParameter
rm(list=ls())
setwd("D:/Study/My projects/Predict-zika-res/Validation") ## change to your home directory
source('../Code/bmlml.R')

library(mice)
library(parallel)
library(pROC)

## load predictor data
load('../Outputs/mice_model_best.Rda')
load('../Code - Copy/flavi_small_trans.Rda')
load('../Code - Copy/flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))

# response data
disease.matrix = read.csv('../Data/flavi_prim01mid.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -c(1,4,6)])
Y1 = Y1[-175,] # set aside humans
flavi.small.trans = flavi.small.trans[-175,]
N = nrow(Y1)

# load predictions
pred.df = read.csv('../Outputs/risk scores_mid_NoHuPopDen.csv')
pred.df.cv = read.csv('../Outputs/risk scores_10fold_mid_NoHuPopDen.csv')
InPred = pred.df$mean.score
OutPred = pred.df.cv$mean.score

## In-sample prediction
load('../Outputs/pred_matrix_mid_NoHuPopDen.rda')
pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)

## save pdf plot, run if necessary
pdf('../Outputs/MultiModel_PrimFlavi_ROC.pdf', height=4, width=8)
par(mfrow=c(1,2))
plot.roc(flavi.small.trans$flavLabel, InPred, print.auc=T,
         main="In-sample")
for(i in 1:10){
  lines.roc(flavi.small.trans$flavLabel, pred.matrix1[,i],
            col=adjustcolor("black",alpha.f=.1))
}
lines.roc(flavi.small.trans$flavLabel, InPred, lwd=3)

## Out-of-sample prediction
load('../Outputs/pred_matrix_10fold_mid_NoHuPopDen.rda')
pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)
plot.roc(flavi.small.trans$flavLabel, OutPred, print.auc=T,
         main="10-fold CV")
for(i in 1:10){
  lines.roc(flavi.small.trans$flavLabel, pred.matrix1[,i],
            col=adjustcolor("black",alpha.f=.1))
}
lines.roc(flavi.small.trans$flavLabel, OutPred, lwd=3)
par(mfrow=c(1,1))
dev.off()

## Use this instead if leave-one-out results needed
# load('pred_matrix_OOB.rda')
# pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)
# plot.roc(flavi.small.trans$flavLabel, apply(pred.matrix1,1,mean), print.auc=T,
#          main="Out of sample")
# for(i in 1:10){
#   lines.roc(flavi.small.trans$flavLabel, pred.matrix1[,i],
#             col=adjustcolor("black",alpha.f=.1))
# }
# lines.roc(flavi.small.trans$flavLabel, apply(pred.matrix1,1,mean), lwd=3)
# par(mfrow=c(1,1))

## compare in-sample and validation predictions
# pdf("../Outputs/MultiModel_PrimFlavi_InVsOutPlot.pdf", width=5, height=5)
colvec = rep("black",nrow(Y1))
colvec[rowSums(Y1)>0] = "red"
plot(100*((ecdf(InPred))(InPred)),
     (100*(ecdf(OutPred))(OutPred)), pch=19,
     col=adjustcolor(colvec, alpha.f=.7),
     xlab="In-sample scores", ylab="10-fold CV scores",
     main="Training vs. validation percentile scores")
# legend(70, 30, c("Reservoirs","Unknown"),col=c("red","black"), pch=19)
legend("bottomright", c("Reservoirs","Unknown"),col=c("red","black"), pch=19)
abline(0,1, lwd=2)
# dev.off()

## top few species
pred.df$perc.score = round((100*(ecdf(InPred))(InPred)),1)
pred.df1 = pred.df[order(pred.df$perc.score, decreasing=T),]
# View(with(pred.df1, pred.df1[which(NewWorld & !Detected),]))
View(with(pred.df1, pred.df1[which(NewWorld & perc.score > 89),]))

## less out of sample score species
InOut = data.frame(cbind(paste(flavi.small.trans$Spp),
                         100*((ecdf(InPred))(InPred)),
                         (100*(ecdf(OutPred))(OutPred))))
View(InOut[which(pred.df$Detected),])

## top few species for each imputed dataset
load('../Outputs/pred_matrix_mid.rda')
pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)
for(i in 1:10){
  pred.df = data.frame(Species=paste(flavi.small.trans$Spp),
                       NewWorld=(flavi.small.trans$paleotropical==0),
                       Detected=((1:N) %in% which(rowSums(Y1)>0)),
                       mean.score=round(pred.matrix1[,i],3),
                       sd.score=round(apply(pred.matrix1, 1, sd),3),
                       NA.nums=rowSums(is.na(flavi.small.trans[,-(1:8)])))
  pred.df1 = pred.df[order(pred.df$mean.score, decreasing=T),]
  cat(paste("Imputed dataset",i,"-\n"))
  print(with(pred.df1, pred.df1[which(NewWorld & !Detected),])[1:5,1:4])
  cat("\n")
}

## compare BMLML and RF predictions 10-fold
BMLML.pred = read.csv('../Outputs/risk scores_mid.csv')[,-1]
RF.pred = read.csv('primates_RF_all_risk_scores.csv')
plot(rank(BMLML.pred$mean.score) ~ rank(RF.pred$mean_score), pch=19, col=colvec)
abline(0,1)

## save dataset-wise predictions
pred.frame = data.frame(pred.matrix1)
names(pred.frame) = paste0("Dataset.",1:10)
all.df = data.frame(Species=pred.df$Species)
all.df = data.frame(all.df, pred.frame, mean.score=pred.df$mean.score)
write.csv(all.df, file='pred_all.csv', row.names=F)

all.df1 = all.df
for(i in 2:11){
  all.df1[,i] = rank(all.df1[,i])/N
}
all.df1$mean.score = apply(all.df1[,2:11],1,mean)
write.csv(all.df1, file='pred_rank_all.csv', row.names=F)

all.df2 = all.df
for(i in 2:11){
  all.df2[,i] = all.df2[,i]/max(all.df2[,i])
}
all.df2$mean.score = apply(all.df2[,2:11],1,mean)
write.csv(all.df2, file='pred_scale_all.csv', row.names=F)

## partial dependency plots
flavi.imputed.data = read.csv('../Outputs/Imputed_data_1.csv')[,-1]
load('../Outputs/all_models_Tuned.Rda')
varInd = 20 # change this index to get PD curves for other variables
varnames = names(flavi.imputed.data[,-c(1:8,26)])
X = flavi.imputed.data[,-c(1:8,26)]
for(i in 1:5){
  X[,i] = as.numeric(paste(X[,i]))
}
X1 = scale(as.matrix(X))
N = nrow(X1)

## function for getting partial dependency plot for a variable
plotfun = function(varInd){
  defaultPar = par()
  
  ### plot 1: histogram
  par(fig=c(0,1,(2/3),1))
  hist(X[,varInd], freq=F, nclass=20, xlim=c(min(z[,1]),max(z[,1])),main="")
  
  # par(mfrow=c(2,1))
  
  ### plot 2
  par(fig=c(0,1,0,(2/3)))
  par(new=TRUE)
  # mean curve
  mv = mean(X[,varInd])
  sv = sd(X[,varInd])
  z = partial.bmlml(all.models[[1]], X1, 8, varInd, lims=c(-5,5), transform="none")
  z[,1] = mv + sv*z[,1]
  plot(z, type='l',ylim=c(0,1), lwd=2, col="blue", main=varnames[varInd], xlab="values", ylab="prob")
  
  # median curve
  z = partial.bmlml(all.models[[1]], X1, 8, varInd, lims=c(-5,5), median=T, transform="none")
  z[,1] = z[,1] = mv + sv*z[,1]
  lines(z, type='l',ylim=c(0,1), lwd=2, col="red")
  
  # individual curves
  for(n in 1:N){
    z = partial.bmlml(all.models[[1]], matrix(X1[n,], nrow=1), 8, varInd, lims=c(-5,5), transform="none")
    z[,1] = mv + sv*z[,1]
    lines(z, type='l', col=adjustcolor("black",alpha.f=.2))
  }
  
  par(defaultPar)
}

plotfun(20)

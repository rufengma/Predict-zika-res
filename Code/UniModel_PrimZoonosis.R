## Univariate model for Primate-Zoonosis connections
rm(list=ls())
setwd("D:/Study/My projects/Predict-zika-res/Code") ## change to your home directory

library(mice)
library(miceadds)
library(parallel)
library(randomForest)
library(pROC)

## load predictor data
load('../Outputs/mice_model_best.Rda')
load('flavi_small_trans.Rda')
load('flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))

# response data
disease.matrix = read.csv('../Data/flavi_prim01mid.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -c(1,4,6)])
ndata = 10
N = nrow(Y1) - 1

## get multiple predictions from random forest
pred.mat = matrix(0, nrow=N, ncol=ndata)
pred.cal.mat = pred.mat

for(i in 1:ndata){
  set.seed(31052016)
  # add one instance of imputed data to actual dataset
  flavi.imputed.data = read.csv(paste0('../Outputs/Imputed_data_',i,'.csv'))
  idata = flavi.imputed.data[-175,-c(1:9,22,34:36,47)] # add 26 if using bodymass_resid instead
  idata$flavLabel = as.numeric(rowSums(Y1[-175,])>0)
  
  ## make model
  rf.mod = randomForest(as.factor(flavLabel)~., data=idata, ntree=1e4, keep.inbag=T)
  rf.pred = predict(rf.mod, flavi.imputed.data, predict.all=T)
  
  ## get out of bag prediction for each sample
  ipred = rep(0, N)
  for(j in 1:N){
    which.oob = which(rf.mod$inbag[j,] == 0)
    ipred[j] = mean(as.numeric(rf.pred$individual[j, which.oob]))
  }
  
  ## calibration for PU problem
  true = idata$flavLabel
  e1 = mean(ipred[which(true==1)])
  e2 = sum(ipred[which(true==1)]) / sum(ipred)
  e3 = max(ipred)
  
  # method 1: direct calibration
  ipred.cal = ipred
  ipred.cal[which(true==0)] = apply(cbind(ipred[which(true==0)] / e3, 1), 1, min)
  
  # store
  pred.mat[,i] = ipred
  pred.cal.mat[,i] = ipred.cal
  cat(paste("Data",i,"done!\n"))
}

# average over all data and plot
pred = rowMeans(pred.mat)
pred.cal = rowMeans(pred.cal.mat)

# pdf(8,4,file="UniModel_preddensity.pdf")
plot(density(pred[which(true==0)]), lwd=2, col="darkgreen", xlim=c(0,1), ylim=c(0,25),
     main="Density plot of predicted probabilities", xlab="Probability")
lines(density(pred[which(true==1)]), lwd=2, col="darkred")
lines(density(pred.cal[which(true==0)]), lwd=2, col="darkgreen", lty=2)
legend('topright', c("FLAV+ primates", "Unknown", "Unknown calibrated"),
       col=c("darkred","darkgreen", "darkgreen"), lty=c(1,1,2), lwd=2)
# dev.off()

## ROC curves before and after calibration
# pdf(8,4,file="UniModel_roc.pdf")
par(mfrow=c(1,2))
plot.roc(true, pred, main="Before", print.auc=T)
plot.roc(true, pred.cal, main="After", print.auc=T)
par(mfrow=c(1,1))
# dev.off()

# interesting samples: top 10 New-world undetected primates
# which.ambi = which(pred.cal[which(true==0)] > .9*max(pred.cal[which(true==1)]))
pred.df = data.frame(Species=paste(flavi.imputed.data$Spp),
                     NewWorld=(flavi.imputed.data$paleotropical==0),
                     Detected=((1:length(true)) %in% which(flavi.imputed.data$ZLabel==1)),
                     Prob=round(pred.cal,2),
                     Prob.sd=round(apply(pred.cal.mat,1,sd),2),
                     NA.nums=rowSums(is.na(flavi.small.trans[,-(1:8)])))
pred.df1 = pred.df[order(pred.df$Prob, decreasing=TRUE),]
pred.df1[with(pred.df1, which(NewWorld & !Detected)),][1:10,-(2:3)]

# variable importance
varImpPlot(rf.mod)

# partial dependence
imp <- importance(rf.mod)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
par(mfrow=c(4,4))
for (i in seq_along(impvar)) {
  partialPlot(rf.mod, flavi.imputed.data[,-c(1,2,4:8)], impvar[i], xlab=impvar[i],
              # main=paste("Partial Dependence on \n", impvar[i]),
              main=paste(i),
              ylab="Probability")
}
par(mfrow=c(1,1))

############################################################
## End of code
############################################################

## Misc codes
# method 2: refitting with weighted data
# get augmented dataset
# w = rep(1, nrow(flavi.imputed.data))
# w[which(true==0)] = (1/e2-1) * 
#   predicted.cal[which(true==0)] / (1-predicted.cal[which(true==0)])
# newdata = rbind(flavi.imputed.data[which(true==1),],
#                 flavi.imputed.data[which(true==0),],
#                 flavi.imputed.data[which(true==0),])
# newdata$ZLabel = c(rep(1, nZ), rep(1, n-nZ), rep(0, n-nZ))
# wts = c(rep(1, nZ), w[which(true==0)], 1-w[which(true==0)])
# 
# # fit model
# glmod = glm(as.integer(ZLabel)~., family=binomial, data=newdata[,-c(1,2,4:13)], weights=wts)
# 
# predicted = predict(glmod, newdata=newdata, type="response")
# plot(density(predicted[which(true==0)]), lwd=2, col="darkgreen", xlim=c(0,1), ylim=c(0,10),
#      main="Density plot of predicted probabilities", xlab="Probability")
# lines(density(predicted[which(true==1)]), lwd=2, col="darkred")
# legend('topleft', c("Reservoirs", "Not known"), col=c("darkred","darkgreen"), lwd=2)

############################################################
## End of file
############################################################
## Multi-response model for Primate-Flavivirus connections
## Aggregate variable importance from model objects
setwd("C:/Users/Subho/Box Sync/Hunting Zika Virus with Machine Learning (Cary Institute)/Validation")
rm(list=ls())
library(mice)
library(miceadds)
library(parallel)

## load predictor data
load('../Outputs/mice_model_best.Rda')
load('../Code - Copy/flavi_small_trans.Rda')
load('../Code - Copy/flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))

# response data
disease.matrix = read.csv('../Data/flavi_prim01mid.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -c(1,4,6)])

# load imputed datasets
imputed.data = read.csv('../Outputs/Imputed_data_1.csv')[,-1]
X1 = imputed.data[,-c(1:8, 21, 33:35, 46)]

## set aside humans
Y1 = Y1[-175,]
X1 = X1[-175,]

## marginal interaction with risk scores
load('../Outputs/all_models_Tuned_mid.Rda')

# Variable importance
# importance of a variable in a factor is its average rank across runs and datasets
imp.mat = lapply(all.models, function(x) {abs(x$W)})
imp.mat = matrix(unlist(imp.mat), ncol=length(all.models), byrow=F)

bar.data = apply(imp.mat, 1, mean)
bar.names = names(X1)
bar.order = order(bar.data, decreasing=F)
defaultPar = par()

pdf('../Outputs/MultiModel_PrimFlavi_varimp.pdf', width=8, height=8)
par(las=2, mar=c(5,12,2,3),mfrow=c(1,1))
barplot(bar.data[bar.order], horiz=T,
        names.arg=bar.names[bar.order], cex.names=.7,
        main=paste("Average variable importance across factors"), cex.main=.9)
dev.off()

## mean of variables for low and high risk samples
pred.df = read.csv('../Outputs/risk scores.csv')
score = pred.df$mean.score
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
write.csv(var.df1, file='../Outputs/variable means for low and high risk samples.csv', row.names=F)

# make table with heading
# mytable = xtable::xtable(var.df1)
# names(mytable) = c(" ","mBIC2","RFGLS+BH", "q", paste0('t=',as.numeric(t.vec)/10))
# print(mytable, include.rownames = F)
# 
# par(defaultPar)

############################################################
## End of file
############################################################
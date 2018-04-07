## Multi-response model for Primate-Flavivirus connections
## Tuning parameter selection using Doss method,
## followed by validation by OOB and in-sample scores
## server version
setwd("C:/Users/Subho/Box Sync/Hunting Zika Virus with Machine Learning (Cary Institute)/Code - copy")
source('bmlml.R')

library(mice)
library(miceadds)
library(parallel)
library(caret)

## load predictor data
load('../Outputs/mice_model_best.Rda')
load('flavi_small_trans.Rda')
load('flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))

# response data
disease.matrix = read.csv('../Data/flavi_prim01mid.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -c(1,4,6)])

ndata = flavi.imputed.mice$m
N = nrow(flavi.small.trans) - 1
ndisease = ncol(Y1)
disease = ndisease
ncores = min(10, ndata)

## generate ndata instances of imputed dataset
imputed.frames = list()
for(i in 1:ndata){
  set.seed(i)
  idata = complete(flavi.imputed.mice, action=i)
  imputed.frames[[i]] = idata
}

# make list of columnwise min-max for thresholding
func = function(x){
  out = NULL
  if(class(x)=="numeric"){
    na.ind = which(is.na(x))
    out = c(min(x[-na.ind]), max(x[-na.ind]))
  }
  out
}
minmax.list = lapply(flavi.small.trans, func)

##############################################
# OOB scores
##############################################

loopfun = function(idata){
  set.seed(idata*1e3)
  
  # add one instance of imputed data to actual dataset
  flavi.imputed.data = flavi.small.trans
  flavi.imputed.data[,-c(1:8,25)] = imputed.frames[[idata]]
  
  ## threshold individual columns by min and max of non-NA elements
  for(j in 1:ncol(flavi.imputed.data)){
    if(class(flavi.imputed.data[,j])=="numeric"){
      imin = minmax.list[[j]][1]
      imax = minmax.list[[j]][2]
      flavi.imputed.data[which(flavi.imputed.data[,j]<imin),j] = imin
      flavi.imputed.data[which(flavi.imputed.data[,j]>imax),j] = imax
    }
  }
  
  ## create residual for log body mass
  # flavi.imputed.data$logBodyMass_g_Resid = residuals(lm(
  #   logBodyMass_g ~ X16.2_LittersPerYear_EXT+logAgeatFirstBirth_d+logGestationLen_d+
  #     logLitterSize+logNeonateBodyMass_g+logWeaningAge_d,
  #   data=flavi.imputed.data))
  
  ## save data
  # write.csv(flavi.imputed.data,
  #           file=paste0('../Outputs/Imputed_data_',idata,'.csv'))
  
  ## load response data
  X = flavi.imputed.data[,-c(1:8,21,33:35)] # add 26 if using bodymass_resid instead
  for(i in 1:5){
    X[,i] = as.numeric(paste(X[,i]))
  }
  X1 = scale(as.matrix(X))
  
  # set aside humans
  X1 = X1[-175,]
  Y1 = Y1[-175,]
  N = nrow(X1)
  D = ncol(X1)
  
  # full model with 1 latent factor
  K = 1
  r.hat = selectParams(Y1, X1, K=K, nU=1e2, maxit=5e2)
  model.hat = bmlml.EM1(Y1, X1, K=K, r.tau=r.hat, maxit=5e2)
  pred.mat = predict.bmlml(model.hat, X1, transform="scale")
  
  # create stratified folds
  set.seed(idata)
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
    pred.mat[itest,] = predict.bmlml(model.hat.i, X1[itest,], transform="scale")
  }
  
  # return score for Zika
  pred.mat[,disease]
}

# pred.matrix = mclapply(1:10, loopfun, mc.cores=10)
pred.matrix = lapply(1:10, loopfun)

save(pred.matrix, file='../Outputs/pred_matrix_10fold_mid_NoHuPopDen.rda')
pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)
pred.df = data.frame(Species=paste(flavi.small.trans$Spp[-175]),
                     NewWorld=(flavi.small.trans$paleotropical[-175]==0),
                     Detected=((1:N) %in% which(rowSums(Y1[-175,])>0)),
                     ZikaDetected=(Y1[-175,ncol(Y1)]==1),
                     mean.score=round(apply(pred.matrix1, 1, function(x) mean(x, na.rm=T)),3),
                     median.score=round(apply(pred.matrix1, 1, function(x) median(x, na.rm=T)),3),
                     sd.score=round(apply(pred.matrix1, 1, function(x) sd(x, na.rm=T)),3),
                     NA.nums=rowSums(is.na(flavi.small.trans[-175,-(1:8)])))

write.csv(pred.df, file='../Outputs/risk scores_10fold_mid_NoHuPopDen.csv')

##############################################
# In-sample scores
##############################################

loopfun = function(idata){
  set.seed(idata*1e3)
  flavi.imputed.data = flavi.small.trans
  flavi.imputed.data[,-c(1:8,25)] = imputed.frames[[idata]]
  
  ## threshold individual columns by min and max of non-NA elements
  for(j in 1:ncol(flavi.imputed.data)){
    if(class(flavi.imputed.data[,j])=="numeric"){
      imin = minmax.list[[j]][1]
      imax = minmax.list[[j]][2]
      flavi.imputed.data[which(flavi.imputed.data[,j]<imin),j] = imin
      flavi.imputed.data[which(flavi.imputed.data[,j]>imax),j] = imax
    }
  }
  
  ## create residual for log body mass
  # flavi.imputed.data$logBodyMass_g_Resid = residuals(lm(
  #   logBodyMass_g ~ X16.2_LittersPerYear_EXT+logAgeatFirstBirth_d+logGestationLen_d+
  #     logLitterSize+logNeonateBodyMass_g+logWeaningAge_d,
  #   data=flavi.imputed.data))
  
  ## save data
  # write.csv(flavi.imputed.data,
  #           file=paste0('../Outputs/Imputed_data_',idata,'.csv'))
  
  ## load response data
  X = flavi.imputed.data[,-c(1:8,21,33:35)] # add 26 if using bodymass_resid instead
  for(i in 1:5){
    X[,i] = as.numeric(paste(X[,i]))
  }
  X1 = as.matrix(X)
  
  # set aside humans
  X1 = X1[-175,]
  Y1 = Y1[-175,]
  N = nrow(X1)
  D = ncol(X1)

  # full model with 1 latent factor
  K = 1
  r.hat = selectParams(Y1, scale(X1), K=K, nU=1e2, maxit=5e2)
  model.hat = bmlml.EM1(Y1, scale(X1), K=K, r.tau=r.hat, maxit=5e2)
  pred.mat = predict.bmlml(model.hat, scale(X1), transform="scale")
  
  # return score for Zika
  pred.mat[,disease]
}

# pred.matrix = mclapply(1:10, loopfun, mc.cores=10)
pred.matrix = lapply(1:10, loopfun)

save(pred.matrix, file='../Outputs/pred_matrix_mid_NoHuPopDen.rda')

pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)
pred.df = data.frame(Species=paste(flavi.small.trans$Spp[-175]),
                     NewWorld=(flavi.small.trans$paleotropical[-175]==0),
                     Detected=((1:N) %in% which(rowSums(Y1[-175,])>0)),
                     ZikaDetected=(Y1[-175,ncol(Y1)]==1),
                     mean.score=round(apply(pred.matrix1, 1, function(x) mean(x, na.rm=T)),3),
                     median.score=round(apply(pred.matrix1, 1, function(x) median(x, na.rm=T)),3),
                     sd.score=round(apply(pred.matrix1, 1, function(x) sd(x, na.rm=T)),3),
                     NA.nums=rowSums(is.na(flavi.small.trans[-175,-(1:8)])))

write.csv(pred.df, file='../Outputs/risk scores_mid_NoHuPopDen.csv')

##############################################
# Store model objects
##############################################

loopfun = function(idata){
  set.seed(idata*1e3)
  flavi.imputed.data = flavi.small.trans
  flavi.imputed.data[,-c(1:8,25)] = imputed.frames[[idata]]
  
  ## threshold individual columns by min and max of non-NA elements
  for(j in 1:ncol(flavi.imputed.data)){
    if(class(flavi.imputed.data[,j])=="numeric"){
      imin = minmax.list[[j]][1]
      imax = minmax.list[[j]][2]
      flavi.imputed.data[which(flavi.imputed.data[,j]<imin),j] = imin
      flavi.imputed.data[which(flavi.imputed.data[,j]>imax),j] = imax
    }
  }
  
  ## create residual for log body mass
  # flavi.imputed.data$logBodyMass_g_Resid = residuals(lm(
  #   logBodyMass_g ~ X16.2_LittersPerYear_EXT+logAgeatFirstBirth_d+logGestationLen_d+
  #     logLitterSize+logNeonateBodyMass_g+logWeaningAge_d,
  #   data=flavi.imputed.data))
  
  ## save data
  # write.csv(flavi.imputed.data,
  #           file=paste0('../Outputs/Imputed_data_',idata,'.csv'))
  
  ## load response data
  X = flavi.imputed.data[,-c(1:8,21,33:35)] # add 26 if using bodymass_resid instead
  for(i in 1:5){
    X[,i] = as.numeric(paste(X[,i]))
  }
  X1 = as.matrix(X)
  
  # set aside humans
  X1 = X1[-175,]
  Y1 = Y1[-175,]
  N = nrow(X1)
  D = ncol(X1)

  # full model with 1 latent factor
  K = 1
  r.hat = selectParams(Y1, scale(X1), K=K, nU=1e2, maxit=5e2)
  model.hat = bmlml.EM1(Y1, scale(X1), K=K, r.tau=r.hat, maxit=5e2)

  # return model
  model.hat
}

# all.models = mclapply(1:10, loopfun, mc.cores=10)
all.models = lapply(1:10, loopfun)

save(all.models, file='../Outputs/all_models_Tuned_mid.rda')

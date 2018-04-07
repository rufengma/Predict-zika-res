## Multi-response model for Primate-Flavivirus connections
## Tuning parameter selection using Doss method,
## followed by validation by OOB and in-sample scores
## server version
# setwd("C:/Users/Subho/Box Sync/Hunting Zika Virus with Machine Learning (Cary Institute)/Code")
source('bmlml.R')

library(mice)
# library(miceadds)
library(parallel)

## load predictor data
load('../Outputs/mice_model_best.Rda')
load('flavi_small_trans.Rda')
load('flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))

# response data
disease.matrix = read.csv('flavi_prim01new.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -1])

ndata = flavi.imputed.mice$m
N = nrow(flavi.small.trans)
ndisease = 8
disease = 8
nrun = 1e2
ncores = min(10, ndata)
nfac = 8 # number of latent factors

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
  write.csv(flavi.imputed.data,
            file=paste0('../Outputs/Imputed_data_',idata,'.csv'))
  
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
  flavi.imputed.data$logBodyMass_g_Resid = residuals(lm(
    logBodyMass_g ~ X16.2_LittersPerYear_EXT+logAgeatFirstBirth_d+logGestationLen_d+
      logLitterSize+logNeonateBodyMass_g+logWeaningAge_d,
    data=flavi.imputed.data))
  
  ## load response data
  X = flavi.imputed.data[,-c(1:8,26)]
  for(i in 1:5){
    X[,i] = as.numeric(paste(X[,i]))
  }
  X1 = as.matrix(X)
  N = nrow(X1)
  D = ncol(X1)
  
  # full model with 1 latent factor
  K = 1
  r.hat = selectParams(Y1, scale(X1), K=K, nU=1e2, maxit=5e2)
  model.hat = bmlml.EM1(Y1, scale(X1), K=K, r.tau=r.hat, maxit=5e2)
  pred.mat = predict.bmlml(model.hat, scale(X1))
  
  # now do for each positive sample
  pos.ind = which(Y1>0, arr.ind=T)
  for(i in 1:nrow(pos.ind)){
    ix = pos.ind[i,1]
    iy = pos.ind[i,2]
    Y2 = Y1; Y2[ix,iy] = 0
    
    r.hat.i = selectParams(Y2, scale(X1), K=K, nU=1e2, maxit=5e2)
    model.hat.i = bmlml.EM1(Y2, scale(X1), K=K, r.tau=r.hat.i, maxit=5e2)
    pred.mat.i = predict.bmlml(model.hat.i, scale(X1))
    pred.mat[ix,] = pred.mat.i[ix,]
  }
  
  # return score for Zika
  pred.mat[,8]
}

pred.matrix = mclapply(1:10, loopfun, mc.cores=10)
save(pred.matrix, file='../Outputs/pred_matrix_OOB.rda')
pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)
pred.df = data.frame(Species=paste(flavi.small.trans$Spp),
                     NewWorld=(flavi.small.trans$paleotropical==0),
                     Detected=((1:N) %in% which(rowSums(Y1)>0)),
                     mean.score=round(apply(pred.matrix1, 1, mean),3),
                     sd.score=round(apply(pred.matrix1, 1, sd),3),
                     NA.nums=rowSums(is.na(flavi.small.trans[,-(1:8)])))

write.csv(pred.df, file='../Outputs/risk scores_OOB.csv')

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
  flavi.imputed.data$logBodyMass_g_Resid = residuals(lm(
    logBodyMass_g ~ X16.2_LittersPerYear_EXT+logAgeatFirstBirth_d+logGestationLen_d+
      logLitterSize+logNeonateBodyMass_g+logWeaningAge_d,
    data=flavi.imputed.data))
  
  ## load response data
  X = flavi.imputed.data[,-c(1:8,26)]
  for(i in 1:5){
    X[,i] = as.numeric(paste(X[,i]))
  }
  X1 = as.matrix(X)
  N = nrow(X1)
  D = ncol(X1)
  
  # full model with 1 latent factor
  K = 1
  r.hat = selectParams(Y1, scale(X1), K=K, nU=1e2, maxit=5e2)
  model.hat = bmlml.EM1(Y1, scale(X1), K=K, r.tau=r.hat, maxit=5e2)
  pred.mat = predict.bmlml(model.hat, scale(X1))
  
  # return score for Zika
  pred.mat[,8]
}

pred.matrix = mclapply(1:10, loopfun, mc.cores=10)
save(pred.matrix, file='../Outputs/pred_matrix.rda')

pred.matrix1 = matrix(unlist(pred.matrix), ncol=10, byrow=F)
pred.df = data.frame(Species=paste(flavi.small.trans$Spp),
                     NewWorld=(flavi.small.trans$paleotropical==0),
                     Detected=((1:N) %in% which(rowSums(Y1)>0)),
                     mean.score=round(apply(pred.matrix1, 1, mean),3),
                     sd.score=round(apply(pred.matrix1, 1, sd),3),
                     NA.nums=rowSums(is.na(flavi.small.trans[,-(1:8)])))

write.csv(pred.df, file='../Outputs/risk scores.csv')

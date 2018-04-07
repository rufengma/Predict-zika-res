## PD plots: aggregate version. PD plot of a variable is average of its PD plots over all datasets
setwd("D:/Study/My projects/Predict-zika-res/Validation") ## change to your home directory
source('../Code/bmlml.R')

## partial dependency plots
load('../Outputs/all_models_Tuned_mid.rda')
# varInd = 20 # change this index to get PD curves for other variables

## get response matrix
load('../Data/flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))
disease.matrix = read.csv('../Data/flavi_prim01mid.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -c(1,4,6)])
which.pos = as.numeric(which(rowSums(Y1)>0))

## function for getting partial dependency plot for a variable
get.PD = function(dataInd){
  
  # load data
  flavi.imputed.data = read.csv(paste0('../Outputs/Imputed_data_',dataInd,'.csv'))[,-1]
  varnames = names(flavi.imputed.data[,-c(1:8,21,33:35,46)])
  X = flavi.imputed.data[,-c(1:8,21,33:35,46)]
  for(i in 1:5){
    X[,i] = as.numeric(paste(X[,i]))
  }
  X1 = scale(as.matrix(X))
  
  # set aside humans
  X1 = X1[-175,]
  Y1 = Y1[-175,]
  N = nrow(X1)
  p = ncol(X1)
  
  y.list = list()
  x.mat = matrix(0, nrow=p, ncol=1e2)
  
  pb = txtProgressBar(0, p)
  for(varInd in 1:p){
    mv = mean(X[,varInd])
    sv = sd(X[,varInd])
    
    # now plot
    z = partial.bmlml(all.models[[dataInd]], matrix(X1[1,], nrow=1),
                      ncol(Y1), varInd, lims=c(-5,5), transform="none")
    z.mat = matrix(0, nrow=N, ncol=nrow(z))
    for(n in 1:N){
      z = partial.bmlml(all.models[[dataInd]], matrix(X1[n,], nrow=1),
                        ncol(Y1), varInd, lims=c(-5,5), transform="none")
      z.mat[n,] = z[,2]
    }
    
    # z = partial.bmlml(all.models[[dataInd]], X1, 8, varInd, lims=c(-5,5), transform="none")
    # z.mat = t(z)
    
    y.list[[varInd]] = z.mat
    x.mat[varInd,] = mv + sv*z[,1]
    setTxtProgressBar(pb, varInd)
  }
  close(pb)
  
  # return
  list(y.list=y.list, x.mat=x.mat)
}

# dataInd=1
PD.list = list()
for(dataInd in 1:10){
  PD.list[[dataInd]] = get.PD(dataInd=dataInd)
  
}
# save(PD.list, file="PDlist.Rda")

# now process
y.list1 = list()
y.list0 = y.list1
up.list0 = y.list0
up.list1 = y.list0
dn.list0 = y.list0
dn.list1 = y.list0

for(dataInd in 1:10){

  ## make summary measure, upper and lower bounds
  if(method=="median"){
    y.list0[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[-which.pos,,], 2:3, median)
    y.list1[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[which.pos,,], 2:3, median)
    
    up.list0[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[-which.pos,,], 2:3,
                                function(x) quantile(x, .75))
    up.list1[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[which.pos,,], 2:3,
                                function(x) quantile(x, .75))
    
    dn.list0[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[-which.pos,,], 2:3,
                                function(x) quantile(x, .25))
    dn.list1[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[which.pos,,], 2:3,
                                function(x) quantile(x, .25))
  } else{
    y.list0[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[-which.pos,,], 2:3, mean)
    y.list1[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[which.pos,,], 2:3, mean)
    
    up.list0[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[-which.pos,,], 2:3,
                                function(x) min(mean(x)+sd(x),1))
    up.list1[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[which.pos,,], 2:3,
                                function(x) min(mean(x)+sd(x),1))
    
    dn.list0[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[-which.pos,,], 2:3,
                                function(x) max(mean(x)-sd(x),0))
    dn.list1[[dataInd]] = apply(simplify2array(PD.list[[dataInd]]$y.list)[which.pos,,], 2:3,
                                function(x) max(mean(x)-sd(x),0))
  }
}
y.mat0 = apply(simplify2array(y.list0), 1:2, mean)
y.mat1 = apply(simplify2array(y.list1), 1:2, mean)
up.mat0 = apply(simplify2array(up.list0), 1:2, mean)
up.mat1 = apply(simplify2array(up.list1), 1:2, mean)
dn.mat0 = apply(simplify2array(dn.list0), 1:2, mean)
dn.mat1 = apply(simplify2array(dn.list1), 1:2, mean)

# load predicted scores to make quantile risks
pred.df = read.csv('../Outputs/risk scores_mid_NoHuPopDen.csv')
InPred = pred.df$mean.score
fPred = ecdf(InPred)

flavi.imputed.data = read.csv(paste0('../Outputs/Imputed_data_',dataInd,'.csv'))[,-1]
X = flavi.imputed.data[,-c(1:8,21,33:35,46)]

# par(mfrow=c(3,3))
plotInd = c(7,27,17,24,10,13,19,9,16,23) # indices to plot
indnames = c("Max latitude of geo range (dd)", "Log neonate body mass (g)",
             "Log body mass (g)", "log interbirth interval (d)",
             "Max longitude of geo range (dd)", "Mean precipitation (mm)",
             "Log age at first birth (d)", "Mid range latitude of geo range (dd)",
             "paleotropical","Log home range of individual (km2)")
xmat = matrix(c(-35, 35,
                NA, NA,
                NA, NA,
                NA, NA,
                -180, 180,
                NA, NA,
                NA, NA,
                -35, 35,
                NA, NA,
                NA, NA), ncol=2, byrow=T)
for(ii in 1:10){
  
  ## file save path
  i = plotInd[ii]
  folder = paste0("PD_Plots_mid_aggregate/")
  file = paste0(paste("PD", names(X)[i], method, sep="_"), ".pdf")
  path = paste0(folder,file)
  pdf(path, 5, 5)
  
  x.vec = PD.list[[1]]$x.mat[i,]
  y1.vec = 100*fPred(y.mat1[,i])
  y0.vec = 100*fPred(y.mat0[,i])
  up0.vec = 100*fPred(up.mat0[,i])
  up1.vec = 100*fPred(up.mat1[,i])
  dn0.vec = 100*fPred(dn.mat0[,i])
  dn1.vec = 100*fPred(dn.mat1[,i])
  
  plot(y1.vec~x.vec, type='l', main="",
       ylim=c(0, max(up1.vec)),
       xlab=indnames[ii], ylab="Percentile risk score")
  lines(y0.vec~x.vec, col="blue")
  
  polygon(c(rev(x.vec), x.vec), c(rev(dn1.vec), up1.vec),
          col = 'grey80', border = NA)
  lines(x.vec, y1.vec,  col = 'black')
  lines(x.vec, up1.vec, lty = 'dashed', col = 'black')
  lines(x.vec, dn1.vec, lty = 'dashed',ylim=c(0,1), lwd=2, col="black")
  
  polygon(c(rev(x.vec), x.vec), c(rev(dn0.vec), up0.vec),
          col = adjustcolor('skyblue',alpha.f=.2), border = NA)
  lines(x.vec, y0.vec,  col = 'blue')
  lines(x.vec, up0.vec, lty = 'dashed', col = 'blue')
  lines(x.vec, dn0.vec, lty = 'dashed', ylim=c(0,1), lwd=2, col="blue")
  
  # legend("topleft", c("Known","Unknown"), col=c("black","blue"), lwd=2)
  dev.off()
}

for(ii in c(1,5,8)){
  
  ## file save path
  i = plotInd[ii]
  folder = paste0("PD_Plots_mid_aggregate/")
  file = paste0(paste("PD", names(X)[i], method, sep="_"), ".pdf")
  path = paste0(folder,file)
  pdf(path, 5, 5)
  
  x.vec = PD.list[[1]]$x.mat[i,]
  y1.vec = 100*fPred(y.mat1[,i])
  y0.vec = 100*fPred(y.mat0[,i])
  up0.vec = 100*fPred(up.mat0[,i])
  up1.vec = 100*fPred(up.mat1[,i])
  dn0.vec = 100*fPred(dn.mat0[,i])
  dn1.vec = 100*fPred(dn.mat1[,i])
  
  plot(y1.vec~x.vec, type='l', main="",
       ylim=c(0, max(up1.vec)), xlim=xmat[ii,],
       xlab=indnames[ii], ylab="Percentile risk score")
  lines(y0.vec~x.vec, col="blue")
  
  polygon(c(rev(x.vec), x.vec), c(rev(dn1.vec), up1.vec),
          col = 'grey80', border = NA)
  lines(x.vec, y1.vec,  col = 'black')
  lines(x.vec, up1.vec, lty = 'dashed', col = 'black')
  lines(x.vec, dn1.vec, lty = 'dashed',ylim=c(0,1), lwd=2, col="black")
  
  polygon(c(rev(x.vec), x.vec), c(rev(dn0.vec), up0.vec),
          col = adjustcolor('skyblue',alpha.f=.2), border = NA)
  lines(x.vec, y0.vec,  col = 'blue')
  lines(x.vec, up0.vec, lty = 'dashed', col = 'blue')
  lines(x.vec, dn0.vec, lty = 'dashed', ylim=c(0,1), lwd=2, col="blue")
  
  # legend("topleft", c("Known","Unknown"), col=c("black","blue"), lwd=2)
  dev.off()
}
# par(mfrow=c(1,1))
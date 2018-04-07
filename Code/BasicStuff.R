rm(list=ls())
setwd("D:/Study/My projects/Predict-zika-res/Code") ## change to your home directory

library(plyr)
library(mi)
library(randomForest)
library(softImpute)
library(missMDA)

prim = read.csv("../Data/prim.csv", row.names = NULL)[,-1]
flavi = read.csv("../Data/primates flavi.csv")[,-1]
mosquito = read.csv('../Data/mosquito.csv')

# how much NA?
na.amount = apply(flavi, 2, FUN = function(x) length(which(is.na(x))))
na.df = data.frame(Variable=names(na.amount),
                   NAlength=as.numeric(na.amount))
na.df = na.df[order(na.df[,2], decreasing=T),]
row.names(na.df) = NULL

## plot NA values
defaultPar = par()
par(mar=c(12,4,2,1))
plot(na.df$NAlength, pch=19, xlab="", ylab="No. of NA samples")
axis(1, at=seq(1, nrow(na.df), by=1), labels = FALSE)
text(x = seq(0.5, nrow(na.df)-.5, by=1), par("usr")[3] - 0.2, labels = na.df$Variable,
     offset=7, srt = 90, pos = 1, xpd = TRUE, cex=.7)
abline(h=310, lty=2)
par(defaultPar)

## take columns with more than 30 non-NA entries
which.in = with(na.df, Variable[which(NAlength < 310)])
flavi.small = flavi[,which(names(flavi) %in% which.in)]
# save(flavi.small, file='Data/flavi_small.Rda')

## number of NA entries for each sample
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))
barplot(summary(as.factor(NA.amount)),
        main="Number of NA entries by sample",
        xlab="No. of NA predictors", ylab="No. of samples")
flavi.small = flavi.small[-which(NA.amount==36),]

## variable transformation.. using http://esapubs.org/archive/ecol/E090/184/metadata.htm
flavi.small.trans = mutate(flavi.small, logBodyMass_g = log(X5.1_AdultBodyMass_g),
                           logHeadBodyLen_mm = log(X13.1_AdultHeadBodyLen_mm),
                           logAgeatFirstBirth_d = log(X3.1_AgeatFirstBirth_d),
                           logGestationLen_d = log(X9.1_GestationLen_d),
                           logGR_Area_km2 = log(X26.1_GR_Area_km2),
                           logHomeRange_km2 = log(X22.1_HomeRange_km2),
                           logHomeRange_Indiv_km2 = log(X22.2_HomeRange_Indiv_km2),
                           logHuPopDen_Min_n.km2 = log(X27.1_HuPopDen_Min_n.km2+1),
                           logHuPopDen_Mean_n.km2 = log(X27.2_HuPopDen_Mean_n.km2),
                           logHuPopDen_5p_n.km2 = log(X27.3_HuPopDen_5p_n.km2+1),
                           logInterbirthInterval_d = log(X14.1_InterbirthInterval_d),
                           logLitterSize = log(X15.1_LitterSize),
                           logMaxLongevity_m = log(X17.1_MaxLongevity_m),
                           logNeonateBodyMass_g = log(X5.3_NeonateBodyMass_g),
                           logPopulationDensity_n.km2 = log(X21.1_PopulationDensity_n.km2),
                           logPopulationGrpSize = log(X10.1_PopulationGrpSize),
                           logSexualMaturityAge_d = log(X23.1_SexualMaturityAge_d),
                           logX10.2_SocialGrpSize = log(X10.2_SocialGrpSize),
                           logWeaningAge_d = log(X25.1_WeaningAge_d),
                           Temp_Mean_01degC4 = X28.2_Temp_Mean_01degC^4/1e8)

flavi.small.trans[,c('X5.1_AdultBodyMass_g', 'X13.1_AdultHeadBodyLen_mm',
                     'X3.1_AgeatFirstBirth_d', 'X9.1_GestationLen_d',
                     'X22.1_HomeRange_km2', 'X22.2_HomeRange_Indiv_km2',
                     'X14.1_InterbirthInterval_d', 'X15.1_LitterSize',
                     'X17.1_MaxLongevity_m', 'X26.1_GR_Area_km2',
                     'X5.3_NeonateBodyMass_g', 'X21.1_PopulationDensity_n.km2',
                     'X10.1_PopulationGrpSize', 'X10.2_SocialGrpSize',
                     'X23.1_SexualMaturityAge_d', 'X25.1_WeaningAge_d',
                     'X28.2_Temp_Mean_01degC','X27.1_HuPopDen_Min_n.km2',
                     'X27.2_HuPopDen_Mean_n.km2', 'X27.3_HuPopDen_5p_n.km2')] = c()
# save(flavi.small.trans, file="Data/flavi_small_trans.Rda")

## impute

# step1: get data in shape
flavi.impute = missing_data.frame(flavi.small.trans[-c(1:8, 25)])

# step2: change variable types
flavi.impute = change(flavi.impute, y='X6.1_DietBreadth', what='type', to='ordered-categorical')

# nothing for now
show(flavi.impute)
image(flavi.impute)

set.seed(31052016)
# cl = makeCluster(detectCores()-1)
options(mc.cores=2)
system.time(flavi.imputed <- mi(flavi.impute, n.iter=2000, parallel=TRUE))
# save(flavi.imputed, file='mi_model.Rda')
# stopCluster(cl)
flavi.imputed.data = flavi.small.trans
flavi.imputed.data[,-c(1:8, 25)] = mi::complete(flavi.imputed, m=1)[,1:36]
# save(flavi.imputed.data, file='Data/prim_flavi_imputed.Rda')

# check convergence
Rstat = Rhats(flavi.imputed)
df = data.frame(names=paste(lapply(names(Rstat)[1:36], function(x) substring(x, 6))),
                Rstat_mean=as.numeric(Rstat[1:36]),
                Rstat_sd=as.numeric(Rstat[37:72]))
View(df)

## random forest impute
flavi.small.rfimp = rfImpute(as.factor(ZLabel)~., flavi.small.trans[,-c(1,2,4:8)])
flavi.small.rfimp = cbind(flavi.small.trans[,c(1,2,4:8)], flavi.small.rfimp)
names(flavi.small.rfimp)[8] = 'ZLabel'
# save(flavi.small.rfimp, file='Data/prim_flavi_rfimputed.Rda')

## PCA impute
flavi.small.pcaimp = flavi.small.trans
flavi.small.pcaimp[,c(14:24,26:45)] = imputePCA(flavi.small.pcaimp[,c(14:24,26:45)])

## matrix completion
lam = lambda0(as.matrix(flavi.small.trans[,-c(1:13, 25)]))
soft.model = softImpute(as.matrix(flavi.small.trans[,-c(1:13, 25)]),
                              lambda=0, rank.max=30)
flavi.small.soft = flavi.small.trans
flavi.small.soft[,-c(1:13, 25)] = complete(flavi.small.trans[,-c(1:13, 25)], soft.model)

which.impute = (1:ncol(flavi.small.trans))[-c(1:13, 25)]
par(mfrow=c(3,3))
for(i in which.impute){
  plot(density(flavi.small.trans[,i], na.rm=T),
       main=paste(names(flavi.small.trans)[i]))
  lines(density(flavi.small.soft[,i], na.rm=T), lty=2)
}
par(mfrow=c(1,1))


#########################################################################
# test for normality
logt = function(x){
  minx = min(x, na.rm=T)
  minx = ifelse( minx < 0, minx-1, 0)
  log(x - minx)
}

sqt = function(x){
  minx = min(x, na.rm=T)
  minx = ifelse( minx < 0, minx, 0)
  sqrt(x - minx)
}

ntest = function(x){
  unt = shapiro.test(x)$p.value

  logt = shapiro.test(logt(x))$p.value
  sqt = shapiro.test(sqt(x))$p.value
  c(as.numeric(unt), as.numeric(logt), as.numeric(sqt))
}

box123 = function(x){
  minx = min(x, na.rm=T)
  minx = ifelse( minx < 0, minx-1, 0)
  
  par(mfrow=c(1,3))
  boxplot(x, main="no trans")
  boxplot(logt(x), main="log trans")
  boxplot(sqt(x), main="sqrt trans")
  par(mfrow=c(1,1))
}

which.nonnum = c(1:9, 15, 19, 23, 25, 26, 27, 37:40, 45)
pvals = apply(flavi.small[,-which.nonnum], 2, ntest)
pvals = matrix(unlist(pvals), ncol=3, byrow=T)
pvals = data.frame(Varname = names(flavi.small[-which.nonnum]),
                   pValue.unt = pvals[,1],
                   pValue.log = pvals[,2],
                   pValue.sqrt = pvals[,3])
pvals

# max.trans = rep(0, nrow(pvals))
flavi.small.trans = flavi.small
## apply transformation as per largest p-value
for(i in (1:ncol(flavi.small))[-which.nonnum]){
  var.ind = which(pvals$Varname == names(flavi.small)[i])
  max.trans = which.max(pvals[var.ind, 2:4])

  # apply required transformation
  if(max.trans==2){
    flavi.small.trans[,i] = logt(flavi.small.trans[,i])
  }
  else if(max.trans==3){
    flavi.small.trans[,i] = sqt(flavi.small.trans[,i])
  }
}

## OR figure out transformations individually.. change and check
box123(flavi.small[,46])
flavi.small.trans[,44] = (flavi.small[,44])

flavi.small.trans$X10.1_PopulationGrpSize = logt(flavi.small$X10.1_PopulationGrpSize)
flavi.small.trans$X10.2_SocialGrpSize = logt(flavi.small$X10.2_SocialGrpSize)


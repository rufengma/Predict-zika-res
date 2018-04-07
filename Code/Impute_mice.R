rm(list=ls())
#setwd("C:/Users/IBM_ADMIN/Box Sync/Hunting Zika Virus with Machine Learning (Cary Institute)")
setwd('/extdrive/Work/smajumd/Zika-codes')
library(plyr)
library(mi)
library(miceadds)
library(data.table)

prim = read.csv("prim.csv", row.names = NULL)[,-1]
flavi = read.csv("primates flavi.csv")[,-1]
mosquito = read.csv('mosquito.csv')

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

flavi.small.trans = data.table(flavi.small.trans)
flavi.small.trans[ ,c('X5.1_AdultBodyMass_g', 'X13.1_AdultHeadBodyLen_mm',
                     'X3.1_AgeatFirstBirth_d', 'X9.1_GestationLen_d',
                     'X22.1_HomeRange_km2', 'X22.2_HomeRange_Indiv_km2',
                     'X14.1_InterbirthInterval_d', 'X15.1_LitterSize',
                     'X17.1_MaxLongevity_m', 'X26.1_GR_Area_km2',
                     'X5.3_NeonateBodyMass_g', 'X21.1_PopulationDensity_n.km2',
                     'X10.1_PopulationGrpSize', 'X10.2_SocialGrpSize',
                     'X23.1_SexualMaturityAge_d', 'X25.1_WeaningAge_d',
                     'X28.2_Temp_Mean_01degC','X27.1_HuPopDen_Min_n.km2',
                     'X27.2_HuPopDen_Mean_n.km2', 'X27.3_HuPopDen_5p_n.km2'):=NULL]
flavi.small.trans = data.frame(flavi.small.trans)
# save(flavi.small.trans, file="flavi_small_trans.Rda")

## turn columns 9 to 13 to factors
for(i in 9:13){
  flavi.small.trans[,i] = as.factor(flavi.small.trans[,i])
}

## impute
system.time(flavi.imputed.mice <- mice(flavi.small.trans[,-c(1:8,25)],
                                       m=10, maxit=200, printFlag=T, 
                                       method=c(rep("polr",3),
                                                rep("logreg",2),
                                                rep("norm", 31))))
Rhat.mice(flavi.imputed.mice)
save(flavi.imputed.mice, file='../Zika-outputs/mice_model.Rda')


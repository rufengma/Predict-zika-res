rm(list=ls())
setwd("D:/Study/My projects/Predict-zika-res/Code") ## change to your home directory
library(plyr)
library(data.table)

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
na.df = na.df[-c(50:57),]
defaultPar = par()
par(mar=c(12,4,2,1))
# pdf('NumNabyPredictor.pdf',7,5)
plot(na.df$NAlength, pch=19, xlab="", ylab="No. of NA samples")
axis(1, at=seq(1, nrow(na.df), by=1), labels = FALSE)
text(x = seq(0.5, nrow(na.df)-.5, by=1), par("usr")[3] - 0.2, labels = na.df$Variable,
     offset=7, srt = 90, pos = 1, xpd = TRUE, cex=.7)
abline(h=310, lty=2)
# dev.off()
par(defaultPar)

## take columns with more than 30 non-NA entries
which.in = with(na.df, Variable[which(NAlength < 310)])
flavi.small = flavi[,which(names(flavi) %in% which.in)]
flavi.small = flavi.small[,-(29:32)] # set aside hupopden variables
# save(flavi.small, file='Data/flavi_small.Rda')

## number of NA entries for each sample
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))
barplot(summary(as.factor(NA.amount)),
        main="Number of NA entries by sample",
        xlab="No. of NA predictors", ylab="No. of samples")
# flavi.small = flavi.small[-which(NA.amount==36),]

## heatmap
flavi.small.mi = mi::missing_data.frame(flavi.small)
defaultPar = par()
par(mar=c(20,4,2,1))
image(flavi.small.mi)
par(defaultPar)

rm(list=ls())
setwd("D:/Study/My projects/Predict-zika-res/Code") ## change to your home directory

library(rgdal)
library(leaflet)

# read in output
load('../Data/flavi_small.Rda')
NA.amount = rowSums(is.na(flavi.small[,-(1:8)]))
disease.matrix = read.csv('../../Data/flavi_prim01new.csv')
Y1 = as.matrix(disease.matrix[-which(NA.amount==36), -1])

# read in risk scores
pred.df = read.csv('../Outputs/risk scores.csv')
scores = with(pred.df, mean.score[which(Y1[,8]==0)])
species = with(pred.df, gsub("_", " ", Species[which(Y1[,8]==0)]))
species = species[order(scores, decreasing=F)]
scores = scores[order(scores, decreasing=F)]
pal <- colorQuantile(c("white","#b80000"), NULL, n = 10) # change colors if needed
colors = pal(scores)

## overlay new world primate geographic range layers
# reads in shapefile from each folder in the directory and puts it in map object
setwd('new world/')
file.list = list.files()
file.list = species[which(species %in% file.list)]

## create map object and overlay OSM tiles
m = leaflet() %>% addTiles()

pb = txtProgressBar(0,length(file.list))
for(i in 1:length(file.list)){
  ipath = file.list[i]
  
  # get color
  icolor = colors[which(species == ipath)]
  if(length(icolor)==0){
    icolor = NA
  }
  
  ilist = list.files(ipath)
  if(length(ilist)==1){
    ipath = paste0(ipath, "/",ilist)
    iname = (list.files(ipath))[3]
  }
  else{
    iname = ilist[3]
  }
  
  iname = sub("\\..*$", "", iname)
  istring = paste("Species: ", file.list[i], ", Risk score ", scores[which(species == ipath)])
  imap = readOGR(dsn=ipath, layer=iname, verbose=F)
  m <- addPolygons(m, stroke = FALSE, smoothFactor = 0.1,
                   fillOpacity = ifelse(is.na(icolor), 0, .3),
                   color=icolor, data=imap, popup=istring)
  setTxtProgressBar(pb,i)
}
close(pb)

## overlay Africa and Madagascar geographic range layers
setwd('../Africa and Madagascar/')
file.list = list.files()
file.list = species[which(species %in% file.list)]

pb = txtProgressBar(0,length(file.list))
for(i in 1:length(file.list)){
  ipath = file.list[i]
  
  # get color
  icolor = colors[which(species == ipath)]
  if(length(icolor)==0){
    icolor = NA
  }
  
  ilist = list.files(ipath)
  if(length(ilist)==1){
    ipath = paste0(ipath, "/",ilist)
    iname = (list.files(ipath))[4]
  }
  else{
    iname = ilist[4]
  }
  
  iname = sub("\\..*$", "", iname)
  istring = paste("Species: ", file.list[i], ", Risk score ", scores[which(species == ipath)])
  imap = readOGR(dsn=ipath, layer=iname, verbose=F)
  m <- addPolygons(m, stroke = FALSE, smoothFactor = 0.1,
                   fillOpacity = ifelse(is.na(icolor), 0, .3),
                   color=icolor, data=imap, popup=istring)
  setTxtProgressBar(pb,i)
}
close(pb)

## overlay Asia primate geographic range layers
setwd('../Asia and SE Asia/')
file.list = list.files()
file.list = species[which(species %in% file.list)]

pb = txtProgressBar(0,length(file.list))
for(i in 1:length(file.list)){
  ipath = file.list[i]
  
  # get color
  icolor = colors[which(species == ipath)]
  if(length(icolor)==0){
    icolor = NA
  }
  
  ilist = list.files(ipath)
  if(length(ilist)==1){
    ipath = paste0(ipath, "/",ilist)
    iname = (list.files(ipath))[4]
  }
  else{
    iname = ilist[4]
  }
  
  iname = sub("\\..*$", "", iname)
  istring = paste("Species: ", file.list[i], ", Risk score ", scores[which(species == ipath)])
  imap = readOGR(dsn=ipath, layer=iname, verbose=F)
  m <- addPolygons(m, stroke = FALSE, smoothFactor = 0.1,
                   fillOpacity = ifelse(is.na(icolor), 0, .3),
                   color=icolor, data=imap, popup=istring)
  setTxtProgressBar(pb,i)
}
close(pb)

## country lines.. optional
download.file(file.path('http://www.naturalearthdata.com/http/',
                        'www.naturalearthdata.com/download/50m/cultural',
                        'ne_50m_admin_0_countries.zip'), 
              f <- tempfile())
unzip(f, exdir=tempdir())

world = readOGR(tempdir(), 'ne_50m_admin_0_countries', encoding='UTF-8')
m1 = addPolygons(m, stroke=TRUE, fill=FALSE, weight=1,
                  color=grey(.2), data=world) %>% 
    addProviderTiles("CartoDB.Positron") # changes tiles, optional
m1 %>% addLegend("bottomright", pal=pal, title="Zika risk score", values=scores, opacity=1)

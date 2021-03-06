#####################################################################
## This script is used to import and format psyllid and meteo data ##
#####################################################################
## source(here("src/importData.R"))
library(here)
library(readxl)
library(tibble)
library(dplyr)

#########################
## Identify data files ##
#########################
dataFiles = dir(here("data"))
xlxsFiles = dataFiles[grep("xlsx", dataFiles)]
xlxsFiles = xlxsFiles[grep("Data_2005_Comptage", xlxsFiles)]
treeNames = sub(",",".",sub(".xlsx","",sub("Data_2005_Comptage_", "", xlxsFiles)))
nTrees    = length(treeNames)
stages    = c("mature", "egg", "L1",  "L2", "L3", "L4", "L5", "imago")

###########################
## Import data to a list ##
###########################
rawDataList = vector("list", nTrees)
for (ii in 1:nTrees) {
  rawData       = readxl::read_excel(paste0("data/",xlxsFiles[ii]), col_types=c("guess",rep("numeric",15),"guess"))
  rawData$date  = as.POSIXct(rawData$date, format="%Y-%m-%d %H:%M:%S")
  rawData$stage = factor(rawData$stage, levels=c("mature","œuf","L1","L2","L3","L4","L5","imago"))
  levels(rawData$stage) = stages
  rawDataList[[ii]] = cbind(tree = treeNames[ii], rawData)
  rm(rawData)
}

###################
## Reformat data ##
###################
dataList = vector("list", nTrees)
for (ii in 1:nTrees) {
  x = rawDataList[[ii]]
  y = x %>% select("date") %>% distinct() %>% arrange()
  for (st in stages)
    eval(parse(text=paste0("y$",st,"=0")))
  for (iDate in 1:nrow(y)) {
    subX = x %>% filter(date == y$date[iDate]) %>% select(!date)
    y[iDate,"leaves"] = sum(apply( !is.na(subX %>% select(!stage & !tree)), 2, "any"))
    for (rr in 1:nrow(subX)){
      stage_rr = as.character(subX[rr,]$stage)
      y[iDate, stage_rr] =
        subX[rr,] %>% select(!stage & !tree) %>% rowSums(na.rm=TRUE)
    }
  }
  dataList[[ii]] = cbind(tree = treeNames[ii], y) %>% relocate(tree, .after=date)
  if (dataList[[ii]]$mature[1] == 0) {
    dataList[[ii]]$mature[1] = NA
  }
}

###########################################
## Find first and last dates in the list ##
###########################################
dateRange = range(dataList[[1]]$date)
for (ii in 2:nTrees) {
  dateRange = range(c(dateRange, dataList[[ii]]$date))
}
## "2005-03-16" "2005-07-13"


###################################
## Original "small-data" version ##
###################################
## Data_2005_Comptage_12.1     <- as.data.frame(readxl::read_excel("data/Data_2005_Comptage_12,1.xlsx"))[,-1]
## Data_2005_Comptage_12.1$f8  <- as.numeric(Data_2005_Comptage_12.1$f8)
## Data_2005_Comptage_12.1$f9  <- as.numeric(Data_2005_Comptage_12.1$f9)
## Data_2005_Comptage_12.1$f10 <- as.numeric(Data_2005_Comptage_12.1$f10)
## Data_2005_Comptage_12.1$f11 <- as.numeric(Data_2005_Comptage_12.1$f11)
## Data_2005_Comptage_12.1$f12 <- as.numeric(Data_2005_Comptage_12.1$f12)
## Data_2005_Comptage_12.1$f13 <- as.numeric(Data_2005_Comptage_12.1$f13)
## Data_2005_Comptage_12.1$f14 <- as.numeric(Data_2005_Comptage_12.1$f14)
## Data_2005_Comptage_12.1$f15 <- as.numeric(Data_2005_Comptage_12.1$f15)

## D_12.1 <- Data_2005_Comptage_12.1

## colnames(D_12.1) <- c("date","leaves","egg","L1","L2","L3","L4","L5","imago")
## D_12.1$egg <- as.numeric(D_12.1$egg)
## D_12.1$L1 <- as.numeric(D_12.1$L1)
## D_12.1$L2 <- as.numeric(D_12.1$L2)
## D_12.1$L3 <- as.numeric(D_12.1$L3)
## D_12.1$L4 <- as.numeric(D_12.1$L4)
## D_12.1$L5 <- as.numeric(D_12.1$L5)
## D_12.1$imago <- as.numeric(D_12.1$imago)

## for (i in D_12.1$date) {
##   D_12.1[D_12.1$date==i,3] <- sum(as.numeric(Data_2005_Comptage_12.1[which(Data_2005_Comptage_12.1$date==i&Data_2005_Comptage_12.1$stage=="oeuf"),]),na.rm = T)
##   D_12.1[D_12.1$date==i,4] <- sum(as.numeric(Data_2005_Comptage_12.1[which(Data_2005_Comptage_12.1$date==i&Data_2005_Comptage_12.1$stage=="L1"),]),na.rm = T)
##   D_12.1[D_12.1$date==i,5] <- sum(as.numeric(Data_2005_Comptage_12.1[which(Data_2005_Comptage_12.1$date==i&Data_2005_Comptage_12.1$stage=="L2"),]),na.rm = T)
##   D_12.1[D_12.1$date==i,6] <- sum(as.numeric(Data_2005_Comptage_12.1[which(Data_2005_Comptage_12.1$date==i&Data_2005_Comptage_12.1$stage=="L3"),]),na.rm = T)
##   D_12.1[D_12.1$date==i,7] <- sum(as.numeric(Data_2005_Comptage_12.1[which(Data_2005_Comptage_12.1$date==i&Data_2005_Comptage_12.1$stage=="L4"),]),na.rm = T)
##   D_12.1[D_12.1$date==i,8] <- sum(as.numeric(Data_2005_Comptage_12.1[which(Data_2005_Comptage_12.1$date==i&Data_2005_Comptage_12.1$stage=="L5"),]),na.rm = T)
##   D_12.1[D_12.1$date==i,9] <- sum(as.numeric(Data_2005_Comptage_12.1[which(Data_2005_Comptage_12.1$date==i&Data_2005_Comptage_12.1$stage=="imago"),]),na.rm = T)
##   D_12.1[D_12.1$date==i,2] <- length(Data_2005_Comptage_12.1[Data_2005_Comptage_12.1$date==i,2:16][!is.na(Data_2005_Comptage_12.1[Data_2005_Comptage_12.1$date==i,2:16])])
## }

## write.csv2(D_12.1,"D_12.1.csv")

################
## Meteo data ##
################
meteo = readxl::read_excel(paste0("data/",xlxsFiles[ii]), col_types=c("guess",rep("numeric",15),"guess"))
meteo = read.csv(file="data/Data_ALLMeteo_Montpellier.csv", sep=";", dec=",")[,-1] %>% as_tibble() %>% select(!ID)
meteo$date = as.POSIXct(meteo$date, format="%Y-%m-%d %H:%M:%S")
meteo$time = format(meteo$date,'%H:%M:%S')
meteo$temperature = round(meteo$temperature)


## Subset of meteo over study period only
subMeteo = meteo %>% filter(date >= dateRange[1] & date <= dateRange[2])

## A vector of the times of day in meteo
times = unique(meteo$time)

## A plot of temperature over the study period
## Three NAs need imputing
with(subMeteo, plot(date, temperature, typ="l", lwd=3))
abline(v = (subMeteo %>% filter(is.na(temperature)))$date, col=rgb(1,0,0,0.4))
(threeNAs = subMeteo %>% filter(is.na(temperature)))

## Another plot, zooming in on the period around the NAs
buffer = 7
plotRange = c(min(as.Date(threeNAs$date)) - buffer, max(as.Date(threeNAs$date)) + buffer)
par(mfrow = c(1,1))
with(subMeteo %>% filter(as.Date(date)>=plotRange[1] & as.Date(date)<=plotRange[2]), plot(date, temperature, typ="l", lwd=3))
abline(v=(subMeteo %>% filter(time=="00:00:00"))$date, col=rgb(1,0,0,0.5))


####################################################
## Possible evidence of non-3-hour time intervals ##
####################################################
table(diff(meteo$date))
which(diff(subMeteo$date)==2)
subMeteo$date[85:87] - subMeteo$date[84:86]

##################################
## Export data to an Rdata file ##
##################################
meteo    = subMeteo
psyllids = dataList
save(meteo, psyllids, treeNames, nTrees, file="data/data4nimble.Rdata" )

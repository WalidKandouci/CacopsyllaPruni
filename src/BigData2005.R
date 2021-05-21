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

###########################
## Import data to a list ##
###########################
rawDataList = vector("list", nTrees)
for (ii in 1:nTrees) {
  rawData       = readxl::read_excel(paste0("data/",xlxsFiles[ii]), col_types=c("guess","guess",rep("numeric",15),"guess"))[,-1]
  rawData$date  = as.Date(rawData$date, format="%Y-%m-%d %H:%M:%S")
  rawData$stage = factor(rawData$stage, levels=c("mature","oeuf","L1","L2","L3","L4","L5","imago"))
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
  stages = c("mature", "oeuf", "L1",  "L2", "L3", "L4", "L5", "imago")
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
  dataList[[ii]] = cbind(tree = treeNames[ii], y)
}


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

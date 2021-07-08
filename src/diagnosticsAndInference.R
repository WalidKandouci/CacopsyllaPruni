## rm(list=ls())
## source(here("src/diagnosticsAndInference.R"))
library(here)
library(dplyr)
library(nimbleTempDev)

########################
## Set some constants ##
########################
SDmodel = 1 # 2, 3, 4, 5 ## Identifywhich model to use for SD
nTemps  = 4 # 8 12 16 20 ## Number of temperatures in APT samplers
thin    = 10
setConstantsElsewhere = TRUE ## Prevents a redefinition in modelDefinition.R

## ###########################################
## Take arguments from script, if available ##
CA <- commandArgs(TRUE)
if (length(CA)==0) {
  UseScript <- FALSE
} else {
  UseScript <- TRUE
}

if (UseScript) { 
  print(CA)
  print(SDmodel <- as.integer(CA)[1])
  print(qsubID  <- as.integer(CA)[2])
  ## Ensure R can find nimble
  library(dplyr)
  library(here)
  Rlibs = here() %>% sub(pattern="work",replacement="save") %>% sub(pattern="CacopsyllaPruni",replacement="R")
  Rdirs = Rlibs %>% dir()
  x86dir = Rdirs[grep("x86",Rdirs)]
  libPath = paste0(Rlibs,"/",x86dir,"/4.0")
  .libPaths(new=libPath)
  ## .libPaths()
} else {
  qsubID <- 123
}

###########################
## Create rPsyllid model ##
###########################
source(here::here("src/modelDefinition.R"))

##########################
## Compile model to C++ ##
##########################
cPsyllid = compileNimble(rPsyllid)


###############################################
## Load APT output to analyse in this script ##
###############################################
samplesFileStem = ("Jun-19_23-08-57_2021_Temps4")
samples  = read.table(here(paste0("APT/",samplesFileStem,".txt")), header=TRUE)
samples2 = read.table(here(paste0("APT/",samplesFileStem,"_loglik.txt")), header=TRUE)

####################
## Remove burn-in ##
####################
burn = 1:4000
samples  = samples[-burn,]
samples2 = samples2[-burn,]

##############################
## Node lists to work with  ##
##############################
paraNames = gsub("\\..*","", colnames(samples)) %>% unique() # Names of parameter nodes
depNodes  = gsub("\\[.*","", cPsyllid$getDependencies(paraNames, self = FALSE, includeData = FALSE)) %>% unique()

#################
## Plot stages ##
#################

# The idea I have is to plot at first each stage on it's own, and maybe have the curv of the tree that we want to see in a difrent color ?
# I'm not sure of what and how to plot because "stage" is a big array and I don't know if we should plot everything ?
# and if we don't plot everything what should we really plot ?

# DP: Well first, you need to take one row of parameters from the MCMC output, plug it into the model & update and dependant nodes
#     This code chunk will probably need plugging into your loops somewhere

nMcmcSamples = 10 # 1000
pStage = array(NA, dim = c(nMcmcSamples, lMeteo, nTrees, nStagesTot)) # This will be a ragged array - i.e. nSteps[tree] is heterogeneous, so some elements of this array will remain as NAs

for (iMCMC in 1:nMcmcSamples) {
  print(iMCMC)
  iRow    = sample(nrow(samples),1)                          # Index for a random row. Later this can be replaced by a loop.
  for (ii in 1:length(paraNames)) {                            # A loop to fill each parameter node of the model that can handle alternative sets of parameter names
    myCommand = paste0("cPsyllid$", paraNames[ii], " = as.numeric(samples[iRow, grep(paraNames[ii],colnames(samples))])")
    eval(parse(text= myCommand))
    ## Here's a quick check that this works
    ## cPsyllid$shapeSD
    ## samples[iRow, grep("shapeSD",colnames(samples))]
  }
  simulate(cPsyllid, depNodes) # Updates "paras" "devKernel" "states" "pStage"
  for(iTree in 1:nTrees){
    for (iStage in 1:nStagesTot){
      #plot(cPsyllid$states[iTree, iStage, substage],type='l', col="grey", xlab="time", ylab="dev")
      #points(cPsyllid$states[iTree, iStage, substage],type='l', col="grey")
      ## 1) R won't find states. You need to use cPsyllid$states
      ## 2) states includes all sub-stages. We don't really care about substages, so you need to sum them (in the same way that the model does to obtain pStage
      iMeteo = min(iMeteoForObs[[iTree]]):max(iMeteoForObs[[iTree]])
      nSteps = length(iMeteo)
      if(iStage < nStagesTot){
        pStage[iMCMC,1:nSteps,iTree,iStage] = rowSums(cPsyllid$states[iTree, 1:nSteps, ((iStage-1)*res+1):(iStage*res)])
      } else{
        pStage[iMCMC,1:nSteps,iTree,iStage] = (cPsyllid$states[iTree, 1:nSteps, ((iStage-1)*res+1)])
      }
      ## So first construct pStage (the proportion of the population in each stage)
      ## Then plot it... but also think about how you will handle multiple lines of MCMC output (& will you use lines with transparency or polygons for showing CIs?)
    }
    #browser()
    pStage[iMCMC,1:nSteps,iTree,] = pStage[iMCMC,1:nSteps,iTree,]/ rowSums(pStage[iMCMC,1:nSteps,iTree,])
  }
}


## iMCMC = 1
## iStage = 1
## iTree = 1 # the tree
## iMeteo = min(iMeteoForObs[[iTree]]):max(iMeteoForObs[[iTree]])
## nSteps = length(iMeteo)
## plot(pStage[iMCMC,1:nSteps,iTree,iStage]~meteo$date[iMeteo], type='n',ylim = c(0,1),xlab = "time", ylab= "proportion")
## for(iStage in 1:nStagesTot){
##   lines(pStage[iMCMC,1:nSteps,iTree,iStage]~meteo$date[iMeteo])
##   pStage[iMCMC,1:nSteps,iTree,]
## }


## plot(pStage[iMCMC,1:nSteps,iTree,iStage]~meteo$date[iMeteo], type='n',ylim = c(0,1),xlab = "time", ylab= "proportion")
## for(iStage in 1:nStagesTot){
##   quant = apply(pStage[,1:nSteps,iTree,iStage], 2, "quantile", prob=c(0.025,0.5,0.975))
##   lines(quant[1,]~meteo$date[iMeteo], col="red")
##   lines(quant[2,]~meteo$date[iMeteo], col="blue")
##   lines(quant[3,]~meteo$date[iMeteo], col="red")
##   pStage[iMCMC,1:nSteps,iTree,]
## }

pdf(file = here(paste0("APT/",samplesFileStem, "_proportions.pdf")))
for(iTree in 1:nTrees) {
  iMeteo = min(iMeteoForObs[[iTree]]):max(iMeteoForObs[[iTree]])
  nSteps = length(iMeteo)
  plot(
    pStage[iMCMC, 1:nSteps, iTree, iStage] ~ meteo$date[iMeteo],
    type = 'n',
    ylim = c(0, 1),
    xlab = "time",
    ylab = "proportion",
    main = paste("Tree", psyllids[[iTree]][1, 2])
  )
  vecColor = c("red", "blue", "green", "grey", "orange", "purple", "pink")
  for (iStage in 1:nStagesTot) {
    quant = apply(pStage[, 1:nSteps, iTree, iStage], 2, "quantile", prob = c(0.025, 0.5, 0.975))
    #lines(quant[1,]~meteo$date[iMeteo], col="red")
    polygon(
      x = c(meteo$date[iMeteo], rev(meteo$date[iMeteo])),
      y = c(quant[1, ], rev(quant[3, ])),
      col = adjustcolor(vecColor[iStage], alpha.f = 0.50),
      border = NA
    )
    points(
      meteo$date[iMeteoForObs[[iTree]]][-1],
      (psyllidsArray[iTree, 1:nObs[iTree], ] / rowSums(psyllidsArray[iTree, 1:nObs[iTree], ]))[-1, iStage],
      pch = 4,
      lwd = 2,
      col = vecColor[iStage]
    )
    #lines(quant[2,]~meteo$date[iMeteo], col="blue")
    #lines(quant[3,]~meteo$date[iMeteo], col="red")
    #pStage[iMCMC,1:nSteps,iTree,]
  }
}
dev.off()

############################
## Plot the Briere curves ##
############################

if (FALSE) { # TRUE
  par(mfrow=n2mfrow(Const$nStagesDev))
  for (stage in 1:Const$nStagesDev){
    ttemp = -20:60
    quant = apply(paras[stage,iTemp,2],2, "quantile",prob = c(0.025, 0.5, 0.975))
    plot(ttemp,stBriere(T=-20:60, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]),
         type = 'n', xlab = "time",ylab = "proportion")
    polygon(
      x = c(ttemp, rev(ttemp)),
      y = c(quant[1, ], rev(quant[3, ])),
      col = adjustcolor("red", alpha.f = 0.25),
      border = NA
    )
    lines(ttemp, ttemp,stBriere(T=-20:60, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]), 
          col="red",lwd=1.5)
  }
}


################################################################################
################################################################################
# Fix model
SDmodel=1 #2,3,4,5,6

# Load data 

# Plot
# library(resample)
APTresults_mean <- colMeans(APTresults)
APTresults_sd   <- colStdevs(APTresults)



curve(stBriere(x,Tmin=APTresults$Tmax.1.[1],Tmax=APTresults$Tmin.1.[1],shape=APTresults$logit_shapeMean.1.[1],amplitude=exp(APTresults$logit_amplitudeMean.1.[1])), 
      -10,60, n=1001)


stBriere(-60:60,Tmin=APTresults$Tmax.1.[1],Tmax=APTresults$Tmin.1.[1],shape=APTresults$logit_shapeMean.1.[1],amplitude=exp(APTresults$logit_amplitudeMean.1.[1]))
         

polygon(
  x = c(ttemp, rev(ttemp)),
  y = c(stBriere(x,Tmin=(mean(APTresults$Tmax.1.)-sd(APTresults$Tmax.1.)),Tmax=APTresults$Tmin.1.[1],shape=APTresults$logit_shapeMean.1.[1],amplitude=exp(APTresults$logit_amplitudeMean.1.[1])), 
        rev(quant[3, ])),
  col = adjustcolor("red", alpha.f = 0.25),
  border = NA
)

#############
## To plot ##
#############
library(matrixStats)  # for the sdMeans function
APTresults <- read.csv("C:/Users/Walid/Desktop/model6_3636152_Jul-_3_1934_Temps4.txt", sep="")

SdModel_Oeuf <- matrix(NA, nrow=dim(APTresults)[1], ncol = length(-80:80))
SdModel_L1   <- matrix(NA, nrow=dim(APTresults)[1], ncol = length(-80:80))
SdModel_L2   <- matrix(NA, nrow=dim(APTresults)[1], ncol = length(-80:80))
SdModel_L3   <- matrix(NA, nrow=dim(APTresults)[1], ncol = length(-80:80))
SdModel_L4   <- matrix(NA, nrow=dim(APTresults)[1], ncol = length(-80:80))
SdModel_L5   <- matrix(NA, nrow=dim(APTresults)[1], ncol = length(-80:80))

for (i in 1: dim(APTresults)[1]){
  SdModel_Oeuf[i,] <- stBriere(-80:80,
                             Tmin=APTresults$Tmin.1.[i],
                             Tmax=APTresults$Tmax.1.[i],
                             shape=ilogit(APTresults$logit_shapeMean.1.[i]),
                             amplitude=ilogit(APTresults$logit_amplitudeMean.1.[i]))
  SdModel_L1[i,]   <- stBriere(-80:80,
                                  Tmin=APTresults$Tmin.2.[i],
                                  Tmax=APTresults$Tmax.2.[i],
                                  shape=ilogit(APTresults$logit_shapeMean.2.[i]),
                                  amplitude=ilogit(APTresults$logit_amplitudeMean.2.[i]))
  SdModel_L2[i,]   <- stBriere(-80:80,
                                  Tmin=APTresults$Tmin.3.[i],
                                  Tmax=APTresults$Tmax.3.[i],
                                  shape=ilogit(APTresults$logit_shapeMean.3.[i]),
                                  amplitude=ilogit(APTresults$logit_amplitudeMean.3.[i]))
  SdModel_L3[i,]   <- stBriere(-80:80,
                                  Tmin=APTresults$Tmin.4.[i],
                                  Tmax=APTresults$Tmax.4.[i],
                                  shape=ilogit(APTresults$logit_shapeMean.4.[i]),
                                  amplitude=ilogit(APTresults$logit_amplitudeMean.4.[i]))
  SdModel_L4[i,]   <- stBriere(-80:80,
                                  Tmin=APTresults$Tmin.5.[i],
                                  Tmax=APTresults$Tmax.5.[i],
                                  shape=ilogit(APTresults$logit_shapeMean.5.[i]),
                                  amplitude=ilogit(APTresults$logit_amplitudeMean.5.[i]))
  SdModel_L5[i,]   <- stBriere(-80:80,
                                  Tmin=APTresults$Tmin.6.[i],
                                  Tmax=APTresults$Tmax.6.[i],
                                  shape=ilogit(APTresults$logit_shapeMean.6.[i]),
                                  amplitude=ilogit(APTresults$logit_amplitudeMean.6.[i]))
}

meanOeuf <- colMeans(SdModel_Oeuf[1:1000,])
sdOeuf   <- colSds(as.matrix(SdModel_Oeuf[1:1000,]))

meanL1 <- colMeans(SdModel_L1[1:1000,])
sdL1   <- colSds(SdModel_L1[1:1000,])

meanL2 <- colMeans(SdModel_L2[1:1000,])
sdL2   <- colSds(SdModel_L2[1:1000,])

meanL3 <- colMeans(SdModel_L3[1:1000,])
sdL3   <- colSds(SdModel_L3[1:1000,])

meanL4 <- colMeans(SdModel_L4[1:1000,])
sdL4   <- colSds(SdModel_L4[1:1000,])

meanL5 <- colMeans(SdModel_L5[1:1000,])
sdL5   <- colSds(SdModel_L5[1:1000,])

##############
## SDmodel1 ##
##############

par(mfrow=c(3,2))
{
  # Oeuf
  plot(-80:80,meanOeuf,xlim=c(0,40),ylim=c(0,0.035),main="Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf-sdOeuf),rev(meanOeuf+sdOeuf)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,xlim=c(-10,60),ylim=c(0,0.3),meanL1,main="Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1-sdL1),rev(meanL1+sdL1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2,xlim=c(-10,30),ylim=c(0,0.07),main="Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2-sdL2),rev(meanL2+sdL2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,xlim=c(0,40),ylim=c(0,0.053),meanL3,main="Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3-sdL3),rev(meanL3+sdL3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,xlim=c(0,30),ylim=c(0,0.055),meanL4,main="Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4-sdL4),rev(meanL4+sdL4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5,xlim=c(-20,40),ylim=c(0,0.02),main="Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5-sdL5),rev(meanL5+sdL5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}

##############
## SDmodel2 ##
##############

SdModel2_Oeuf <- SdModel_Oeuf
SdModel2_L1   <- SdModel_L1
SdModel2_L2   <- SdModel_L2
SdModel2_L3   <- SdModel_L3
SdModel2_L4   <- SdModel_L4
SdModel2_L5   <- SdModel_L5

for (i in 1: dim(APTresults)[1]){
  SdModel2_Oeuf[i,] <- SdModel2_Oeuf[i,] * APTresults$beta0.1.[i]
  SdModel2_L1[i,] <- SdModel2_L1[i,]* APTresults$beta0.2.[i]
  SdModel2_L2[i,] <- SdModel2_L2[i,]* APTresults$beta0.3.[i]
  SdModel2_L3[i,] <- SdModel2_L3[i,]* APTresults$beta0.4.[i]
  SdModel2_L4[i,] <- SdModel2_L4[i,]* APTresults$beta0.5.[i]
  SdModel2_L5[i,] <- SdModel2_L5[i,]* APTresults$beta0.6.[i]
    }


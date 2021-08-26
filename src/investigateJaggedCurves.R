########################
## Set some constants ##
########################
SDmodel               = 5 # 2 # 1, 2, 3, 4, 5 ## Identifywhich model to use for SD
nTemps                = 4 # 8 12 16 20 ## Number of temperatures in APT samplers
thin                  = 10
nMcmcSamples          = 5 # 1000
setConstantsElsewhere = TRUE ## Prevents a redefinition in modelDefinition.R
## rm(list=ls())

###############
## Libraries ##
###############
library(here)
library(dplyr)
library(nimbleTempDev)
nimPrint("model = ", SDmodel)
source(here::here("src/functions.R"))

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
samplesFileStem = sub(aptOutputFile, pat=".txt",rep="")
samples  = read.table(here(paste0("APT/",samplesFileStem,".txt")), header=TRUE)
samples2 = read.table(here(paste0("APT/",samplesFileStem,"_loglik.txt")), header=TRUE)

## nrow(samples[!duplicated(samples),])
## sum(duplicated(samples))
## nrow(samples[!duplicated(samples),])

####################
## Remove burn-in ##
####################
burn = 1:1000
samples  = samples[-burn,]
samples2 = samples2[-burn,]

##############################
## Node lists to work with  ##
##############################
(paraNames = gsub("\\..*","", colnames(samples)) %>% unique()) # Names of parameter nodes
(depNodes  = gsub("\\[.*","", cPsyllid$getDependencies(paraNames, self = FALSE, includeData = FALSE)) %>% unique())

#####################################
## Plot population stage structure ##
## STEP 1: construct pStage array  ##
#####################################
nimPrint("#################")
nimPrint("LOOP ON MCMC ROWS")
nimPrint("#################")
pStage       = array(NA, dim = c(nMcmcSamples, lMeteo, nTrees, nStagesTot)) # This will be a ragged array - i.e. nSteps[tree] is heterogeneous, so some elements of this array will remain as NAs
devLogMeanSD = array(NA, dim = c(nMcmcSamples, nStagesDev, lTempVec, 2))
devMean      = array(NA, dim = c(nMcmcSamples, nStagesDev, lTempVec))
devStdev     = array(NA, dim = c(nMcmcSamples, nStagesDev, lTempVec))
devKernel    = array(NA, dim = c(nMcmcSamples, nStagesDev, lTempVec, res+1))
for (iMCMC in 1:nMcmcSamples) {
  nimPrint("iMCMC = ", iMCMC)
  iRow    = sample(nrow(samples),1)                          # Index for a random row. Later this can be replaced by a loop.
  ## Reparameterise model
  for (ii in 1:length(paraNames)) {                            # A loop to fill each parameter node of the model that can handle alternative sets of parameter names
    myCommand = paste0("cPsyllid$", paraNames[ii], " = as.numeric(samples[iRow, grep(paraNames[ii],colnames(samples))])")
    eval(parse(text= myCommand))
    ## Here's a quick check that this works
    ## cPsyllid$shapeSD
    ## samples[iRow, grep("shapeSD",colnames(samples))]
  }
  ## Update dependant nods
  simulate(cPsyllid, depNodes) # Updates "paras" "devKernel" "states" "pStage"
  #######################################################
  ## Loop on trees to extract pStage at each time step ##
  #######################################################
  for(iTree in 1:nTrees){
    for (iStage in 1:nStagesTot){
      iMeteo = min(iMeteoForObs[[iTree]]):max(iMeteoForObs[[iTree]])
      nSteps = length(iMeteo)
      if(iStage < nStagesTot){
        pStage[iMCMC,1:nSteps,iTree,iStage] = rowSums(cPsyllid$states[iTree, 1:nSteps, ((iStage-1)*res+1):(iStage*res)])
      } else {
        pStage[iMCMC,1:nSteps,iTree,iStage] = (cPsyllid$states[iTree, 1:nSteps, ((iStage-1)*res+1)])
      }
    }
    pStage[iMCMC,1:nSteps,iTree,] = pStage[iMCMC,1:nSteps,iTree,] / rowSums(pStage[iMCMC,1:nSteps,iTree,])
  }
  ####################################################################
  ## Loop on stages to extract mean development at each temperature ##
  ####################################################################
  for (iStage in 1:nStagesDev){
    devLogMeanSD[iMCMC,iStage,1:lTempVec,1:2] = t(apply(cPsyllid$paras[iStage,1:lTempVec,1:2], 1, "meanSd2logMeanSd"))
    devMean[iMCMC,iStage,1:lTempVec]  = cPsyllid$paras[iStage,1:lTempVec,1] ## mean  of development kernel
    devStdev[iMCMC,iStage,1:lTempVec] = cPsyllid$paras[iStage,1:lTempVec,2] ## stdev of development kernel
    devKernel[iMCMC,iStage,1:lTempVec, 1:(res+1)] = cPsyllid$devKernel[iStage,1:lTempVec,1:(res+1)]
  }
}

if (FALSE) {
  for (ii in 1:111) {
    nimPrint("iStage: ",iStage <- sample(nStagesDev, 1))
    nimPrint("iTemp: ",iTemp  <- 31) # sample(lTempVec, 1))
    nimPrint("Dup devMean: ", sum(duplicated(devMean[,iStage,iTemp])))
    nimPrint("Dup devStdev: ",sum(duplicated(devStdev[,iStage,iTemp])))
    nimPrint("Dup devLogMean: ", sum(duplicated(devLogMeanSD[,iStage,iTemp,1])))
    nimPrint("Dup devLogStdv: ", sum(duplicated(devLogMeanSD[,iStage,iTemp,2])))
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

#######################################
## Plot population stage structure   ##
## STEP 2: plot info in pStage array ##
#######################################
nimPrint("####################")
nimPrint("PLOTTING PROPORTIONS")
nimPrint("####################")
samplesFileStem <- samplesFileStem %>% sub(pat="-", rep="")
### pdf(file = here( paste0("figures/",samplesFileStem, "_nMCMC-",nMcmcSamples, "_proportions.pdf")))
for(iTree in 11) {
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
    #lines(quant[2,]~meteo$date[iMeteo], col="blue")
    #lines(quant[3,]~meteo$date[iMeteo], col="red")
    #pStage[iMCMC,1:nSteps,iTree,]
  }
}
abline(v=meteo$date[302])
which(iMeteo==302) # 222
iMeteo[222]


all(diff(rowSums(cPsyllid$states[iTree, 221:237,1:25])) < 0 )         ## TRUE
all(diff(rowSums(cPsyllid$states[iTree, 221:237,(1:25)+25*1])) > 0 )  ## TRUE
all(diff(rowSums(cPsyllid$states[iTree, 221:237,(1:25)+25*2])) < 0 )  ## TRUE

rowSums(cPsyllid$states[iTree, 221:237,(1:25)+25*2])

##for (time in 1:nSteps[tree])
# states[tree, time+1, 1:(nStagesDev*res+1)] <-
iTree = 11
iTime = 222
state222 = cPsyllid$states[iTree, iTime, 1:(nStagesDev*res+1)]
state223 = sparseTWstep(state222, cPsyllid$devKernel[1:nStagesDev, Const$iMeteoTemp[Const$iMeteoForObsMat[iTree,1] + iTime - 1], 1:(res+1)])


Const$iMeteoTemp[Const$iMeteoForObsMat[iTree,1] + iTime - 1]

tempVec[11]
cPsyllid$Tmin
cPsyllid$Tmax

cPsyllid$devKernel[1:nStagesDev, Const$iMeteoTemp[Const$iMeteoForObsMat[iTree,1] + iTime - 1], 1:(res+1)]



cbind(state222, state223, state223-state222)
state223
state223



####################################################
## This verifies that there is no population loss ##
## The total population remains very close to 1   ##
####################################################
for (iTree in 1:nTrees) {
  iMeteo = min(iMeteoForObs[[iTree]]):max(iMeteoForObs[[iTree]])
  nSteps = length(iMeteo)
  print(iTree)
  print(summary(1-rowSums(cPsyllid$states[iTree, 1:nSteps, ])))
}



##
iTree = 11
cPsyllid$states

iMeteo = min(iMeteoForObs[[iTree]]):max(iMeteoForObs[[iTree]])
nSteps = length(iMeteo)
for (iStage in 1:nStagesTot){
  if(iStage < nStagesTot){
    pStage[iMCMC,1:nSteps,iTree,iStage] = rowSums(cPsyllid$states[iTree, 1:nSteps, ((iStage-1)*res+1):(iStage*res)])
  } else {
    pStage[iMCMC,1:nSteps,iTree,iStage] = (cPsyllid$states[iTree, 1:nSteps, ((iStage-1)*res+1)])
  }
}
pStage[iMCMC,1:nSteps,iTree,] = pStage[iMCMC,1:nSteps,iTree,] / rowSums(pStage[iMCMC,1:nSteps,iTree,])

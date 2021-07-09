## rm(list=ls())
## source(here("src/diagnosticsAndInference.R"))
library(here)
library(dplyr)
library(nimbleTempDev)

########################
## Set some constants ##
########################
SDmodel               = 6 # 1, 2, 3, 4, 5 ## Identifywhich model to use for SD
nTemps                = 4 # 8 12 16 20 ## Number of temperatures in APT samplers
thin                  = 10
nMcmcSamples          = 10 # 1000
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
# calculate(cPsyllid)

###############################
## List available APT output ##
###############################
(aptOutputFiles = sort(dir(here("APT/"), pattern="*Temps4.txt")))
(aptOutputFile = aptOutputFiles[SDmodel])

###############################################
## Load APT output to analyse in this script ##
###############################################
samplesFileStem = sub(aptOutputFile, pat=".txt",rep="")
samples  = read.table(here(paste0("APT/",samplesFileStem,".txt")), header=TRUE)
samples2 = read.table(here(paste0("APT/",samplesFileStem,"_loglik.txt")), header=TRUE)

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
    devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2] = t(apply(cPsyllid$paras[iStage,1:lTempVec,1:2], 1, "meanSd2logMeanSd"))
    devMean[iMCMC,iStage,1:lTempVec]  = cPsyllid$paras[iStage,1:lTempVec,1] ## mean  of development kernel
    devStdev[iMCMC,iStage,1:lTempVec] = cPsyllid$paras[iStage,1:lTempVec,2] ## stdev of development kernel
    devKernel[iMCMC,iStage,1:lTempVec, 1:(res+1)] = cPsyllid$devKernel[iStage,1:lTempVec,1:(res+1)]
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
pdf(file = here(paste0("figures/",samplesFileStem, "_proportions.pdf")))
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



#####################################################
## Reduce dimension on uncertainty repressentation ##
## i.e. nMcmcSamples -> nQuantiles                 ##
#####################################################
pVec = c(0.01, 0.5, 0.99)
nQuantiles = length(pVec)
devKernelQuantiles = array(NA, dim = c(nQuantiles, nStagesDev, lTempVec, res+1))
devKernelMean      = array(NA, dim = c(nStagesDev, lTempVec, res+1))
devKernelStdv      = array(NA, dim = c(nStagesDev, lTempVec, res+1))
for (iStage in 1:nStagesDev) {
  nimPrint("iStage = ", iStage)
  for (iTemp in 1:lTempVec) {
    # Compute quantiles
    devKernelQuantiles[,iStage,iTemp,] = apply(devKernel[1:nMcmcSamples,iStage,iTemp,1:(res+1)], MARGIN=2, FUN=quantile, p=c(0.01, 0.5, 0.99))
    devKernelMean[iStage,iTemp,]   = apply(devKernel[1:nMcmcSamples,iStage,iTemp,1:(res+1)], MARGIN=2, FUN=mean)
    devKernelStdv[iStage,iTemp,]   = apply(devKernel[1:nMcmcSamples,iStage,iTemp,1:(res+1)], MARGIN=2, FUN=sd)
  }
}


#######################################################################
## Plot discrete representation of development kernel - not pretty!! ##
#######################################################################
library(ggplot2)

plotList = vector("list", 6)
for (iStage in 1:6) {
  meltDevKernMean = reshape::melt(devKernelMean[iStage,,])
  colnames(meltDevKernMean) = c("Temp","subStage","development")
  meltDevKernMean$Temp = tempVec[meltDevKernMean$Temp]
  head(meltDevKernMean)
  ggp <- ggplot(meltDevKernMean, aes(Temp, subStage)) + geom_tile(aes(fill = development))
  ## ggp
  plotList[[iStage]] = ggp + scale_fill_gradient(low = "green", high = "black", breaks=c(0, 0.001, 0.01, 0.1, 1))
}
multiplot(plotList[[1]],plotList[[2]],plotList[[3]],plotList[[4]],plotList[[5]],plotList[[6]])


##########################################################
## Plot continuous representation of development kernel ##
##########################################################
# devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2] = t(apply(cPsyllid$paras[iStage,1:lTempVec,1:2], 1, "meanSd2logMeanSd"))
# devMean[iMCMC,iStage,1:lTempVec]  = cPsyllid$paras[iStage,1:lTempVec,1] ## mean  of development kernel
# devStdev[iMCMC,iStage,1:lTempVec] = cPsyllid$paras[iStage,1:lTempVec,2] ## stdev of development kernel

## Quantiles for each line of MCMC
devQuantiles  = array(NA, dim = c(nMcmcSamples, nStagesDev, lTempVec, nQuantiles))
EdevQuantiles = array(NA, dim = c(nStagesDev, lTempVec, nQuantiles))
for (iStage in 1:nStagesDev)
  for (iMcmc in 1:nMcmcSamples)
    devQuantiles[iMcmc,iStage,,] = t(apply(devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2], 1, function(x) qlnorm(p=pVec, x[1], x[2])))

## Expected quantiles
for (iStage in 1:nStagesDev)
  for (iTemp in 1:nTemps)
    EdevQuantiles[iStage,iTemp,] = devQuantiles[,iStage,iTemp,] %>% colMeans()

iStage = 3
iQuant = 3
plot(tempVec, EdevQuantiles[iStage,,iQuant], typ="l", ylim=c(0, max(EdevQuantiles[iStage,,])))
lines(tempVec, devMean[iMCMC,iStage,1:lTempVec])

lines(tempVec, xyz[1,])
lines(tempVec, xyz[2,])

devMean[iMCMC,iStage,1:lTempVec]
qlnorm(p=c(0.01, 0.5, 0.99), -16.43441, 5.398823)
summary(rlnorm(n=1E4, -16.43441, 5.398823))

# pdf(file = here(paste0("figures/model",SDmodel,"_devKernelContinuous.pdf")))
par(mfrow=c(3,2))
for (iStage in 1:nStagesDev) {
  for (iQuant in 1:nQuantiles) {
    image(devKernelQuantiles[iQuant,iStage,,])

  apply(devKernel[1:nMcmcSamples,iStage,iTemp,1:(res+1)], MARGIN=2, FUN=quantile, p=pVec)


  devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2]

  devMean[iMCMC,iStage,1:lTempVec]
  devStdev[iMCMC,iStage,1:lTempVec]


  plot(-80:80,meanOeuf,main="Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf-sdOeuf),rev(meanOeuf+sdOeuf)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}
## dev.off()

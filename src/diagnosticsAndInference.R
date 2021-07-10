## rm(list=ls())
## source(here("src/diagnosticsAndInference.R"))
library(here)
library(dplyr)
library(nimbleTempDev)

########################
## Set some constants ##
########################
SDmodel               = 1 # 1, 2, 3, 4, 5 ## Identifywhich model to use for SD
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

#####################################################################
## Plot population stage structure: step 1, construct pStage array ##
#####################################################################
pStage       = array(NA, dim = c(nMcmcSamples, lMeteo, nTrees, nStagesTot)) # This will be a ragged array - i.e. nSteps[tree] is heterogeneous, so some elements of this array will remain as NAs
devLogMeanSD = array(NA, dim = c(nMcmcSamples, nStagesTot, lTempVec, 2))
devMean      = array(NA, dim = c(nMcmcSamples, nStagesTot, lTempVec))
devStdev     = array(NA, dim = c(nMcmcSamples, nStagesTot, lTempVec))
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
    #
    ## if (any(rowSums(pStage[iMCMC,1:nSteps,iTree,])!=1)) {
    ## print( max(1 - rowSums(pStage[iMCMC,1:nSteps,iTree,])) )
    pStage[iMCMC,1:nSteps,iTree,] = pStage[iMCMC,1:nSteps,iTree,] / rowSums(pStage[iMCMC,1:nSteps,iTree,])
    ## }
  }
  ####################################################################
  ## Loop on stages to extract mean development at each temperature ##
  ####################################################################
  for (iStage in 1:nStagesDev){
    devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2] = t(apply(cPsyllid$paras[iStage,1:lTempVec,1:2], 1, "meanSd2logMeanSd"))
    devMean[iMCMC,iStage,1:lTempVec]  = cPsyllid$paras[iStage,1:lTempVec,1] ## mean  of development kernel
    devStdev[iMCMC,iStage,1:lTempVec] = cPsyllid$paras[iStage,1:lTempVec,2] ## stdev of development kernel
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

########################################################################
## Plot population stage structure: step 2, plot info in pStage array ##
########################################################################
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




#######################################################################
## Plot continuous representation of development kernel as quantiles ##
#######################################################################

## devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2] =
##               t(apply(cPsyllid$paras[iStage,1:lTempVec,1:2], 1, "meanSd2logMeanSd"))
## devMean[iMCMC,iStage,1:lTempVec]  = cPsyllid$paras[iStage,1:lTempVec,1]
## devStdev[iMCMC,iStage,1:lTempVec] = cPsyllid$paras[iStage,1:lTempVec,2]


cbind(tempVec, devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2])

xyz = apply(devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2], 1, function(x) qlnorm(p=c(0.01, 0.5, 0.99), x[1], x[2]))
plot(tempVec, xyz[3,], typ="l", ylim=c(0, max(xyz)))
lines(tempVec, devMean[iMCMC,iStage,1:lTempVec])

lines(tempVec, xyz[1,])
lines(tempVec, xyz[2,])

devMean[iMCMC,iStage,1:lTempVec]
qlnorm(p=c(0.01, 0.5, 0.99), -16.43441, 5.398823)
summary(rlnorm(n=1E4, -16.43441, 5.398823))

# pdf(file = here(paste0("figures/model",SDmodel,"_devKernelContinuous.pdf")))
par(mfrow=c(3,2))
{
  # Oeuf

  d
  devLogMeanSD[iMCMC, iStage, 1:lTempVec, 1:2]

  devMean[iMCMC,iStage,1:lTempVec]
  devStdev[iMCMC,iStage,1:lTempVec]


  plot(-80:80,meanOeuf,main="Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf-sdOeuf),rev(meanOeuf+sdOeuf)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,meanL1,main="Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1-sdL1),rev(meanL1+sdL1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2,main="Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2-sdL2),rev(meanL2+sdL2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,meanL3,main="Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3-sdL3),rev(meanL3+sdL3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,meanL4,main="Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4-sdL4),rev(meanL4+sdL4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5,main="Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5-sdL5),rev(meanL5+sdL5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}
## dev.off()

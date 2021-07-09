## rm(list=ls())
## source(here("src/diagnosticsAndInference.R"))
library(here)
library(dplyr)
library(nimbleTempDev)

########################
## Set some constants ##
########################
SDmodel = 7 # 1, 2, 3, 4, 5 ## Identifywhich model to use for SD
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

#################################################
## Plot stages: step 1, construct pStage array ##
#################################################

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

####################################################
## Plot stages: step 2, plot info in pStage array ##
####################################################
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

curve(stBriere(x,Tmin=samples$Tmax.1.[1],Tmax=samples$Tmin.1.[1],shape=samples$logit_shapeMean.1.[1],amplitude=exp(samples$logit_amplitudeMean.1.[1])), -10,60, n=1001, ylab="Development", xlab="Temperature")

polygon(
  x = c(ttemp, rev(ttemp)),
  y = c(stBriere(x,
                 Tmin=(mean(samples$Tmax.1.)-sd(samples$Tmax.1.)), ## DP: that looks wrong
                 Tmax=samples$Tmin.1.[1],                          ## DP: This also looks wrong
                 shape=samples$logit_shapeMean.1.[1],
                 amplitude=exp(samples$logit_amplitudeMean.1.[1])), rev(quant[3, ])),
  col = adjustcolor("red", alpha.f = 0.25),
  border = NA
)
################################################################################

#############
## To plot ##
#############
library(matrixStats)  # for the sdMeans function
#### samples <- read.csv("C:/Users/Walid/Desktop/model6_3636152_Jul-_3_1934_Temps4.txt", sep="") ## DP: This is just repeating what is done on line 61, but is less elegent. Here you must copy and paste a file name, whereas line 61 will work conditionally on SDmodel (defined at start of script).

sdModel_Oeuf <- matrix(NA, nrow=dim(samples)[1], ncol = length(-80:80))
sdModel_L1   <- matrix(NA, nrow=dim(samples)[1], ncol = length(-80:80))
sdModel_L2   <- matrix(NA, nrow=dim(samples)[1], ncol = length(-80:80))
sdModel_L3   <- matrix(NA, nrow=dim(samples)[1], ncol = length(-80:80))
sdModel_L4   <- matrix(NA, nrow=dim(samples)[1], ncol = length(-80:80))
sdModel_L5   <- matrix(NA, nrow=dim(samples)[1], ncol = length(-80:80))

sdModel1_Oeuf = sdModel2_Oeuf = sdModel3_Oeuf = sdModel4_Oeuf = sdModel5_Oeuf = sdModel6_Oeuf = sdModel_Oeuf
sdModel1_L1 = sdModel2_L1 = sdModel3_L1 = sdModel4_L1 = sdModel5_L1 = sdModel6_L1 = sdModel_L1
sdModel1_L2 = sdModel2_L2 = sdModel3_L2 = sdModel4_L2 = sdModel5_L2 = sdModel6_L2 = sdModel_L2
sdModel1_L3 = sdModel2_L3 = sdModel3_L3 = sdModel4_L3 = sdModel5_L3 = sdModel6_L3 = sdModel_L3
sdModel1_L4 = sdModel2_L4 = sdModel3_L4 = sdModel4_L4 = sdModel5_L4 = sdModel6_L4 = sdModel_L4
sdModel1_L5 = sdModel2_L5 = sdModel3_L5 = sdModel4_L5 = sdModel5_L5 = sdModel6_L5 = sdModel_L5

##################################
## start with mean and mean+-SD ##
##################################

for (i in 1: dim(samples)[1]){
  sdModel_Oeuf[i,] <- stBriere(-80:80,
                             Tmin=samples$Tmin.1.[i],
                             Tmax=samples$Tmax.1.[i],
                             shape=ilogit(samples$logit_shapeMean.1.[i]),
                             amplitude=ilogit(samples$logit_amplitudeMean.1.[i]))
  sdModel_L1[i,]   <- stBriere(-80:80,
                                  Tmin=samples$Tmin.2.[i],
                                  Tmax=samples$Tmax.2.[i],
                                  shape=ilogit(samples$logit_shapeMean.2.[i]),
                                  amplitude=ilogit(samples$logit_amplitudeMean.2.[i]))
  sdModel_L2[i,]   <- stBriere(-80:80,
                                  Tmin=samples$Tmin.3.[i],
                                  Tmax=samples$Tmax.3.[i],
                                  shape=ilogit(samples$logit_shapeMean.3.[i]),
                                  amplitude=ilogit(samples$logit_amplitudeMean.3.[i]))
  sdModel_L3[i,]   <- stBriere(-80:80,
                                  Tmin=samples$Tmin.4.[i],
                                  Tmax=samples$Tmax.4.[i],
                                  shape=ilogit(samples$logit_shapeMean.4.[i]),
                                  amplitude=ilogit(samples$logit_amplitudeMean.4.[i]))
  sdModel_L4[i,]   <- stBriere(-80:80,
                                  Tmin=samples$Tmin.5.[i],
                                  Tmax=samples$Tmax.5.[i],
                                  shape=ilogit(samples$logit_shapeMean.5.[i]),
                                  amplitude=ilogit(samples$logit_amplitudeMean.5.[i]))
  sdModel_L5[i,]   <- stBriere(-80:80,
                                  Tmin=samples$Tmin.6.[i],
                                  Tmax=samples$Tmax.6.[i],
                                  shape=ilogit(samples$logit_shapeMean.6.[i]),
                                  amplitude=ilogit(samples$logit_amplitudeMean.6.[i]))
}

meanOeuf <- colMeans(sdModel_Oeuf)
sdOeuf   <- colSds(as.matrix(sdModel_Oeuf))

meanL1 <- colMeans(sdModel_L1)
sdL1   <- colSds(sdModel_L1)

meanL2 <- colMeans(sdModel_L2)
sdL2   <- colSds(sdModel_L2)

meanL3 <- colMeans(sdModel_L3)
sdL3   <- colSds(sdModel_L3)

meanL4 <- colMeans(sdModel_L4)
sdL4   <- colSds(sdModel_L4)

meanL5 <- colMeans(sdModel_L5)
sdL5   <- colSds(sdModel_L5)
pdf(file = "sdModelsResults.pdf")
par(mfrow=c(3,2))
{
  # Oeuf
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

########################
sdModel1_Oeuf = sdModel2_Oeuf = sdModel3_Oeuf = sdModel4_Oeuf = sdModel5_Oeuf = sdModel6_Oeuf= sdModel_Oeuf
sdModel1_L1 = sdModel2_L1 = sdModel3_L1 = sdModel4_L1 = sdModel5_L1 = sdModel6_L1 =  sdModel_L1
sdModel1_L2 = sdModel2_L2 = sdModel3_L2 = sdModel4_L2 = sdModel5_L2 = sdModel6_L2 = sdModel_L2
sdModel1_L3 = sdModel2_L3 = sdModel3_L3 = sdModel4_L3 = sdModel5_L3 = sdModel6_L3 = sdModel_L3
sdModel1_L4 = sdModel2_L4 = sdModel3_L4 = sdModel4_L4 = sdModel5_L4 = sdModel6_L4 = sdModel_L4
sdModel1_L5 = sdModel2_L5 = sdModel3_L5 = sdModel4_L5 = sdModel5_L5 = sdModel6_L5 = sdModel_L5
########################

##############
## sdModel1 ##
##############

for (i in 1: dim(samples)[1]){
  sdModel1_Oeuf[i,] <- sdModel1_Oeuf[i,] * exp(samples$beta0.1.[i])
  sdModel1_L1[i,]   <- sdModel1_L1[i,] * exp(samples$beta0.2.[i])
  sdModel1_L2[i,]   <- sdModel1_L2[i,] * exp(samples$beta0.3.[i])
  sdModel1_L3[i,]   <- sdModel1_L3[i,] * exp(samples$beta0.4.[i])
  sdModel1_L4[i,]   <- sdModel1_L4[i,] * exp(samples$beta0.5.[i])
  sdModel1_L5[i,]   <- sdModel1_L5[i,] * exp(samples$beta0.6.[i])
}

meanOeuf_sdModel1 <- colMeans(sdModel1_Oeuf)
sdOeuf_sdModel1   <- colSds(as.matrix(sdModel1_Oeuf))

meanL1_sdModel1 <- colMeans(sdModel1_L1)
sdL1_sdModel1   <- colSds(sdModel1_L1)

meanL2_sdModel1 <- colMeans(sdModel1_L2)
sdL2_sdModel1 <- colSds(sdModel1_L2)

meanL3_sdModel1 <- colMeans(sdModel1_L3)
sdL3_sdModel1 <- colSds(sdModel1_L3)

meanL4_sdModel1 <- colMeans(sdModel1_L4)
sdL4_sdModel1 <- colSds(sdModel1_L4)

meanL5_sdModel1 <- colMeans(sdModel1_L5)
sdL5_sdModel1 <- colSds(sdModel1_L5)

par(mfrow=c(3,2))
{
  # Oeuf
  plot(-80:80,meanOeuf_sdModel1,main="sdModel1 - Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf_sdModel1-sdOeuf_sdModel1),rev(meanOeuf_sdModel1+sdOeuf_sdModel1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf_sdModel1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,meanL1_sdModel1,main="sdModel1 - Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1_sdModel1-sdL1_sdModel1),rev(meanL1_sdModel1+sdL1_sdModel1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1_sdModel1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2_sdModel1,main="sdModel1 - Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2_sdModel1-sdL2_sdModel1),rev(meanL2_sdModel1+sdL2_sdModel1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2_sdModel1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,meanL3_sdModel1,main="sdModel1 - Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3_sdModel1-sdL3_sdModel1),rev(meanL3_sdModel1+sdL3_sdModel1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3_sdModel1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,meanL4_sdModel1,main="sdModel1 - Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4_sdModel1-sdL4_sdModel1),rev(meanL4_sdModel1+sdL4_sdModel1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4_sdModel1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5_sdModel1,main="sdModel1 - Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5_sdModel1-sdL5_sdModel1),rev(meanL5_sdModel1+sdL5_sdModel1)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5_sdModel1,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}

##############
## sdModel2 ##
##############

for (i in 1: dim(samples)[1]){
    if(sdModel2_Oeuf[i,]>0){
    sdModel2_Oeuf[i,] <- sdModel_Oeuf[i,] * exp(samples$beta0.1.[i])
    }
    if(sdModel2_L1[i,]>0){
    sdModel2_L1[i,]   <- sdModel_L1[i,] * exp(samples$beta0.2.[i])
    }
    if(sdModel2_L2[i,]>0){
    sdModel2_L2[i,]   <- sdModel_L2[i,] * exp(samples$beta0.3.[i])
    }
    if(sdModel2_L3[i,]>0){
    sdModel2_L3[i,]   <- sdModel_L3[i,] * exp(samples$beta0.4.[i])
    }
    if(sdModel2_L4[i,]>0){
    sdModel2_L4[i,]   <- sdmodel_l4[i,] * exp(samples$beta0.5.[i])
    }
    if(sdModel2_L5[i,]>0){
    sdModel2_L5[i,]   <- sdModel_L5[i,] * exp(samples$beta0.6.[i])
    }
}

meanOeuf_sdModel2 <- colMeans(sdModel2_Oeuf)
sdOeuf_sdModel2   <- colSds(as.matrix(sdModel2_Oeuf))

meanL1_sdModel2 <- colMeans(sdModel2_L1)
sdL1_sdModel2   <- colSds(sdModel2_L1)

meanL2_sdModel2 <- colMeans(sdModel2_L2)
sdL2_sdModel2 <- colSds(sdModel2_L2)

meanL3_sdModel2 <- colMeans(sdModel2_L3)
sdL3_sdModel2 <- colSds(sdModel2_L3)

meanL4_sdModel2 <- colMeans(sdModel2_L4)
sdL4_sdModel2 <- colSds(sdModel2_L4)

meanL5_sdModel2 <- colMeans(sdModel2_L5)
sdL5_sdModel2 <- colSds(sdModel2_L5)

par(mfrow=c(3,2))
{
  # Oeuf
  plot(-80:80,meanOeuf_sdModel2,main="sdModel2 - Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf_sdModel2-sdOeuf_sdModel2),rev(meanOeuf_sdModel2+sdOeuf_sdModel2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf_sdModel2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,meanL1_sdModel2,main="sdModel2 - Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1_sdModel2-sdL1_sdModel2),rev(meanL1_sdModel2+sdL1_sdModel2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1_sdModel2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2_sdModel2,main="sdModel2 - Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2_sdModel2-sdL2_sdModel2),rev(meanL2_sdModel2+sdL2_sdModel2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2_sdModel2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,meanL3_sdModel2,main="sdModel2 - Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3_sdModel2-sdL3_sdModel2),rev(meanL3_sdModel2+sdL3_sdModel2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3_sdModel2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,meanL4_sdModel2,main="sdModel2 - Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4_sdModel2-sdL4_sdModel2),rev(meanL4_sdModel2+sdL4_sdModel2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4_sdModel2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5_sdModel2,main="sdModel2 - Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5_sdModel2-sdL5_sdModel2),rev(meanL5_sdModel2+sdL5_sdModel2)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5_sdModel2,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}

##############
## sdModel3 ##
##############

for (i in 1: dim(samples)[1]){
  for (itemp in 1: lTempVec){
    sdModel3_Oeuf[i,] <- sdModel3_Oeuf[i,] * exp(samples$beta0.1.[i] + samples$beta1.1.[i] * tempVec[itemp])
    sdModel3_L1[i,]   <- sdModel3_L1[i,] * exp(samples$beta0.2.[i] + samples$beta1.2.[i] * tempVec[itemp])
    sdModel3_L2[i,]   <- sdModel3_L2[i,] * exp(samples$beta0.3.[i] + samples$beta1.3.[i] * tempVec[itemp])
    sdModel3_L3[i,]   <- sdModel3_L3[i,] * exp(samples$beta0.4.[i] + samples$beta1.4.[i] * tempVec[itemp])
    sdModel3_L4[i,]   <- sdModel3_L4[i,] * exp(samples$beta0.5.[i] + samples$beta1.5.[i] * tempVec[itemp])
    sdModel3_L5[i,]   <- sdModel3_L5[i,] * exp(samples$beta0.6.[i] + samples$beta1.6.[i] * tempVec[itemp])
  }
}

meanOeuf_sdModel3 <- colMeans(sdModel3_Oeuf)
sdOeuf_sdModel3   <- colSds(as.matrix(sdModel3_Oeuf))

meanL1_sdModel3 <- colMeans(sdModel3_L1)
sdL1_sdModel3   <- colSds(sdModel3_L1)

meanL2_sdModel3 <- colMeans(sdModel3_L2)
sdL2_sdModel3 <- colSds(sdModel3_L2)

meanL3_sdModel3 <- colMeans(sdModel3_L3)
sdL3_sdModel3 <- colSds(sdModel3_L3)

meanL4_sdModel3 <- colMeans(sdModel3_L4)
sdL4_sdModel3 <- colSds(sdModel3_L4)

meanL5_sdModel3 <- colMeans(sdModel3_L5)
sdL5_sdModel3 <- colSds(sdModel3_L5)

par(mfrow=c(3,2))
{
  # Oeuf
  plot(-80:80,meanOeuf_sdModel3,main="sdModel3 - Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf_sdModel3-sdOeuf_sdModel3),rev(meanOeuf_sdModel3+sdOeuf_sdModel3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf_sdModel3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,meanL1_sdModel3,main="sdModel3 - Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1_sdModel3-sdL1_sdModel3),rev(meanL1_sdModel3+sdL1_sdModel3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1_sdModel3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2_sdModel3,main="sdModel3 - Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2_sdModel3-sdL2_sdModel3),rev(meanL2_sdModel3+sdL2_sdModel3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2_sdModel3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,meanL3_sdModel3,main="sdModel3 - Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3_sdModel3-sdL3_sdModel3),rev(meanL3_sdModel3+sdL3_sdModel3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3_sdModel3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,meanL4_sdModel3,main="sdModel3 - Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4_sdModel3-sdL4_sdModel3),rev(meanL4_sdModel3+sdL4_sdModel3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4_sdModel3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5_sdModel3,main="sdModel3 - Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5_sdModel3-sdL5_sdModel3),rev(meanL5_sdModel3+sdL5_sdModel3)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5_sdModel3,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}

##############
## sdModel4 ##
##############

for (i in 1: dim(samples)[1]){
  for (itemp in 1: lTempVec){
    if(sdModel4_Oeuf[i,]>0){
      sdModel4_Oeuf[i,] <- sdModel_Oeuf[i,] * exp(samples$beta0.1.[i] + samples$beta1.1.[i] * tempVec[itemp])
      }
    if(sdModel4_L1[i,]>0){
      sdModel4_L1[i,]   <- sdModel_L1[i,] * exp(samples$beta0.2.[i] + samples$beta1.2.[i] * tempVec[itemp])
      }
    if(sdModel4_L2[i,]>0){
      sdModel4_L2[i,]   <- sdModel_L2[i,] * exp(samples$beta0.3.[i] + samples$beta1.3.[i] * tempVec[itemp])
      }
    if(sdModel4_L3[i,]>0){
      sdModel4_L3[i,]   <- sdModel_L3[i,] * exp(samples$beta0.4.[i] + samples$beta1.4.[i] * tempVec[itemp])
      }
    if(sdModel4_L4[i,]>0){
      sdModel4_L4[i,]   <- sdmodel_l4[i,] * exp(samples$beta0.5.[i] + samples$beta1.5.[i] * tempVec[itemp])
      }
    if(sdModel4_L5[i,]>0){
      sdModel4_L5[i,]   <- sdModel_L5[i,] * exp(samples$beta0.6.[i] + samples$beta1.6.[i] * tempVec[itemp])
      }
  }
}
    
meanOeuf_sdModel4 <- colMeans(sdModel4_Oeuf)
sdOeuf_sdModel4   <- colSds(as.matrix(sdModel4_Oeuf))

meanL1_sdModel4 <- colMeans(sdModel4_L1)
sdL1_sdModel4   <- colSds(sdModel4_L1)

meanL2_sdModel4 <- colMeans(sdModel4_L2)
sdL2_sdModel4 <- colSds(sdModel4_L2)

meanL3_sdModel4 <- colMeans(sdModel4_L3)
sdL3_sdModel4 <- colSds(sdModel4_L3)

meanL4_sdModel4 <- colMeans(sdModel4_L4)
sdL4_sdModel4 <- colSds(sdModel4_L4)

meanL5_sdModel4 <- colMeans(sdModel4_L5)
sdL5_sdModel4 <- colSds(sdModel4_L5)

par(mfrow=c(3,2))
{
  # Oeuf
  plot(-80:80,meanOeuf_sdModel4,main="sdModel4 - Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf_sdModel4-sdOeuf_sdModel4),rev(meanOeuf_sdModel4+sdOeuf_sdModel4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf_sdModel4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,meanL1_sdModel4,main="sdModel4 - Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1_sdModel4-sdL1_sdModel4),rev(meanL1_sdModel4+sdL1_sdModel4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1_sdModel4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2_sdModel4,main="sdModel4 - Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2_sdModel4-sdL2_sdModel4),rev(meanL2_sdModel4+sdL2_sdModel4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2_sdModel4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,meanL3_sdModel4,main="sdModel4 - Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3_sdModel4-sdL3_sdModel4),rev(meanL3_sdModel4+sdL3_sdModel4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3_sdModel4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,meanL4_sdModel4,main="sdModel4 - Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4_sdModel4-sdL4_sdModel4),rev(meanL4_sdModel4+sdL4_sdModel4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4_sdModel4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5_sdModel4,main="sdModel4 - Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5_sdModel4-sdL5_sdModel4),rev(meanL5_sdModel4+sdL5_sdModel4)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5_sdModel4,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}

##############
## sdModel5 ##
##############

for (i in 1: dim(samples)[1]){
  for (itemp in 1: lTempVec){
    sdModel5_Oeuf[i,] <- sdModel5_Oeuf[i,] * exp(samples$beta0.1.[i] + samples$beta1.1.[i] * tempVec[itemp] + samples$beta2.1.[i] * tempVec[itemp])
    sdModel5_L1[i,]   <- sdModel5_L1[i,] * exp(samples$beta0.2.[i] + samples$beta1.2.[i] * tempVec[itemp] + samples$beta2.2.[i] * tempVec[itemp])
    sdModel5_L2[i,]   <- sdModel5_L2[i,] * exp(samples$beta0.3.[i] + samples$beta1.3.[i] * tempVec[itemp] + samples$beta2.3.[i] * tempVec[itemp])
    sdModel5_L3[i,]   <- sdModel5_L3[i,] * exp(samples$beta0.4.[i] + samples$beta1.4.[i] * tempVec[itemp] + samples$beta2.4.[i] * tempVec[itemp])
    sdModel5_L4[i,]   <- sdModel5_L4[i,] * exp(samples$beta0.5.[i] + samples$beta1.5.[i] * tempVec[itemp] + samples$beta2.5.[i] * tempVec[itemp])
    sdModel5_L5[i,]   <- sdModel5_L5[i,] * exp(samples$beta0.6.[i] + samples$beta1.6.[i] * tempVec[itemp] + samples$beta2.6.[i] * tempVec[itemp])
  }
}

meanOeuf_sdModel5 <- colMeans(sdModel5_Oeuf)
sdOeuf_sdModel5   <- colSds(as.matrix(sdModel5_Oeuf))

meanL1_sdModel5 <- colMeans(sdModel5_L1)
sdL1_sdModel5   <- colSds(sdModel5_L1)

meanL2_sdModel5 <- colMeans(sdModel5_L2)
sdL2_sdModel5 <- colSds(sdModel5_L2)

meanL3_sdModel5 <- colMeans(sdModel5_L3)
sdL3_sdModel5 <- colSds(sdModel5_L3)

meanL4_sdModel5 <- colMeans(sdModel5_L4)
sdL4_sdModel5 <- colSds(sdModel5_L4)

meanL5_sdModel5 <- colMeans(sdModel5_L5)
sdL5_sdModel5 <- colSds(sdModel5_L5)

par(mfrow=c(3,2))
{
  # Oeuf
  plot(-80:80,meanOeuf_sdModel5,main="sdModel5 - Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf_sdModel5-sdOeuf_sdModel5),rev(meanOeuf_sdModel5+sdOeuf_sdModel5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf_sdModel5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,meanL1_sdModel5,main="sdModel5 - Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1_sdModel5-sdL1_sdModel5),rev(meanL1_sdModel5+sdL1_sdModel5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1_sdModel5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2_sdModel5,main="sdModel5 - Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2_sdModel5-sdL2_sdModel5),rev(meanL2_sdModel5+sdL2_sdModel5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2_sdModel5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,meanL3_sdModel5,main="sdModel5 - Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3_sdModel5-sdL3_sdModel5),rev(meanL3_sdModel5+sdL3_sdModel5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3_sdModel5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,meanL4_sdModel5,main="sdModel5 - Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4_sdModel5-sdL4_sdModel5),rev(meanL4_sdModel5+sdL4_sdModel5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4_sdModel5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5_sdModel5,main="sdModel5 - Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5_sdModel5-sdL5_sdModel5),rev(meanL5_sdModel5+sdL5_sdModel5)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5_sdModel5,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}

##############
## sdModel6 ##
##############

for (i in 1: dim(samples)[1]){
  for (itemp in 1: lTempVec){
    if(sdModel6_Oeuf[i,]>0){
      sdModel6_Oeuf[i,] <- sdModel_Oeuf[i,] * exp(samples$beta0.1.[i] + samples$beta1.1.[i] * tempVec[itemp] + samples$beta2.1.[i] * tempVec[itemp])
    }
    if(sdModel6_L1[i,]>0){
      sdModel6_L1[i,]   <- sdModel_L1[i,] * exp(samples$beta0.2.[i] + samples$beta1.2.[i] * tempVec[itemp] + samples$beta2.2.[i] * tempVec[itemp])
    }
    if(sdModel6_L2[i,]>0){
      sdModel6_L2[i,]   <- sdModel_L2[i,] * exp(samples$beta0.3.[i] + samples$beta1.3.[i] * tempVec[itemp] + samples$beta2.3.[i] * tempVec[itemp])
    }
    if(sdModel6_L3[i,]>0){
      sdModel6_L3[i,]   <- sdModel_L3[i,] * exp(samples$beta0.4.[i] + samples$beta1.4.[i] * tempVec[itemp] + samples$beta2.4.[i] * tempVec[itemp])
    }
    if(sdModel6_L4[i,]>0){
      sdModel6_L4[i,]   <- sdModel_L4[i,] * exp(samples$beta0.5.[i] + samples$beta1.5.[i] * tempVec[itemp] + samples$beta2.5.[i] * tempVec[itemp])
    }
    if(sdModel6_L5[i,]>0){
      sdModel6_L5[i,]   <- sdModel_L5[i,] * exp(samples$beta0.6.[i] + samples$beta1.6.[i] * tempVec[itemp] + samples$beta2.6.[i] * tempVec[itemp])
    }
  }
}

meanOeuf_sdModel6 <- colMeans(sdModel6_Oeuf)
sdOeuf_sdModel6   <- colSds(as.matrix(sdModel6_Oeuf))

meanL1_sdModel6 <- colMeans(sdModel6_L1)
sdL1_sdModel6   <- colSds(sdModel6_L1)

meanL2_sdModel6 <- colMeans(sdModel6_L2)
sdL2_sdModel6 <- colSds(sdModel6_L2)

meanL3_sdModel6 <- colMeans(sdModel6_L3)
sdL3_sdModel6 <- colSds(sdModel6_L3)

meanL4_sdModel6 <- colMeans(sdModel6_L4)
sdL4_sdModel6 <- colSds(sdModel6_L4)

meanL5_sdModel6 <- colMeans(sdModel6_L5)
sdL5_sdModel6 <- colSds(sdModel6_L5)

par(mfrow=c(3,2))
{
  # Oeuf
  plot(-80:80,meanOeuf_sdModel6,main="sdModel6 - Stade oeuf",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanOeuf_sdModel6-sdOeuf_sdModel6),rev(meanOeuf_sdModel6+sdOeuf_sdModel6)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanOeuf_sdModel6,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L1
  plot(-80:80,meanL1_sdModel6,main="sdModel6 - Stade L1",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL1_sdModel6-sdL1_sdModel6),rev(meanL1_sdModel6+sdL1_sdModel6)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL1_sdModel6,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L2
  plot(-80:80,meanL2_sdModel6,main="sdModel6 - Stade L2",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL2_sdModel6-sdL2_sdModel6),rev(meanL2_sdModel6+sdL2_sdModel6)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL2_sdModel6,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L3
  plot(-80:80,meanL3_sdModel6,main="sdModel6 - Stade L3",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL3_sdModel6-sdL3_sdModel6),rev(meanL3_sdModel6+sdL3_sdModel6)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL3_sdModel6,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L4
  plot(-80:80,meanL4_sdModel6,main="sdModel6 - Stade L4",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL4_sdModel6-sdL4_sdModel6),rev(meanL4_sdModel6+sdL4_sdModel6)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL4_sdModel6,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
  # L5
  plot(-80:80,meanL5_sdModel6,main="sdModel6 - Stade L5",xlab="température",ylab="devRate",type = "n",xaxs = "i",yaxs = "i")
  polygon(c(-80:80,rev(-80:80)),c((meanL5_sdModel6-sdL5_sdModel6),rev(meanL5_sdModel6+sdL5_sdModel6)),col = "springgreen", border = "springgreen",lwd=3)
  lines(-80:80,meanL5_sdModel6,lty="solid",col="black",lwd=1.5)
  axis(1, col = 'black')
  axis(2, col = 'black')
}


## DP comments
## The above is not too bad, although it could be improved by
## 1) not hard-wiring xlim and ylim - ideally the code can deduce good values based on the curves or polygons
## 2) You repeat everything for each stage... that's a lot of code to check and verify. It is a good exercise to put all that repetition inside a loop - it makes for more compact code.

#####################################################################
## This script provides the definition of the multi-stageIPM model ##
#####################################################################

## source(here::here("src/modelDefinition.R"))

## Transforming integer to date of class POSIXct
## as.POSIXct(as.integer(psyllids[[1]][1,1]), origin = "1970-01-01") == psyllids[[1]][1,1]

# rm(list=ls())
library(here)
library(tibble)
library(dplyr)
library(nimble)
library(nimbleAPT)
library(nimbleTempDev)
library(imputeTS)
library(lubridate)
library(coda)

baseDir = here()
setwd(baseDir)
load("data/data4nimble.Rdata")
source("src/functions.R")

#data4

if (!exists("setConstantsElsewhere")) { ## This flag permits other scripts (which source this script) to set the following constants
  SDmodel = 1 # 2, 3, 4, 5
}

# import imputed values & replace NA
imp <- na_kalman(meteo$temperature) # use the Kalman filter  to imput our missing values
#ggplot_na_distribution(meteo$temperature) # some nice plot
#ggplot_na_imputations(meteo$temperature, imp)
meteo$temperature <- round(imp) # NA free meteo dataset


################################
## BUGS code for nimble model ##
################################
psyllidCode <- nimbleCode ({
  #############################################################
  ## Development kernel at each temperature & for each stage ##
  #############################################################
  for (stage in 1:nStagesDev) { # iStag = index for {, L1, L2, L3, L4, L5, imago}
    ## Priors
    Tmin[stage]                 ~ dunif(-20, 20) # dnorm(0, sd=20)
    Tmax[stage]                 ~ dunif( 20, 60) # dnorm(40,sd=20)
    logit(amplitudeMean[stage]) ~ dLogitBeta(1,1)
    logit(shapeMean[stage])     ~ dLogitBeta(1,1)
    ## Briere functional response curves
    paras[stage,1:lTempVec,1] <- stBriere(t=tempVec[1:lTempVec], Tmin=Tmin[stage], Tmax=Tmax[stage], shape=shapeMean[stage], amplitude=amplitudeMean[stage]) # Mean of the development kernel
    if (SDmodel == 1) {
      logit(shapeSD[stage])       ~ dLogitBeta(1,1)
      logit(amplitudeSD[stage])   ~ dLogitBeta(1,1)
      paras[stage,1:lTempVec,2] <- stBriere(t=tempVec[1:lTempVec], Tmin=Tmin[stage], Tmax=Tmax[stage], shape=shapeSD[stage],   amplitude=amplitudeSD[stage])   # Standard deviation of the development kernel
    }
    ## etc
    for (iTemp in 1:lTempVec) { # iTemp = index for temperature
      ## Survival
      paras[stage,iTemp,3] <- 1
      ## Possibly add a parameter transformation step here ???
      devKernel[stage,iTemp,1:(res+1)] <- getKernel(paras=paras[stage,iTemp,1:3], res=res, devFunction = 1) ## Package currently has functions getM, setM and setMultiM... but we should write a function to just return the first column of getM and work with that (because the model matrix over many stages is very sparse).
    }
  }
  #####################
  ## Loop over trees ##
  #####################
  for (tree in 1:nTrees) { # Adding multiple trees means running the IPM seperately for each tree (due to different start dates)
    # IPM projections
    states[tree, 1, 1] <- 1
    for (substage in 2:(nStagesDev*res+1)){
      states[tree, 1, substage] <- 0
    }
    for (time in 1:nSteps[tree]) { # time = index for time-step
      states[tree, time+1, 1:(nStagesDev*res+1)] <- sparseTWstep(states[tree, time, 1:(nStagesDev*res+1)], devKernel[1:nStagesDev, iMeteoTemp[iMeteoForObsMat[tree,1] + time - 1], 1:(res+1)])
    }
    # Likelihood
    for (obs in 1:nObs[tree]) {
      for(stage in 1:nStagesDev){
        pStage[tree, obs, stage]  <- sum(states[tree, iMeteoForObsMat[tree,obs] - iMeteoForObsMat[tree,1] + 1, ((stage-1)*res+1):(stage*res)])
      }
      pStage[tree, obs, nStagesTot] <- states[tree, iMeteoForObsMat[tree,obs] - iMeteoForObsMat[tree,1] + 1, (nStagesDev*res+1)]
      psyllids[tree, obs, 1:nStagesTot] ~ dmultinom(prob = pStage[tree, obs, 1:nStagesTot], size = psyllidsTotal[tree, obs])
    }
  }
  ## A proxy node for returning logProbs
  ### sumLogProb ~ dnorm(0,1) ## dnorm(0,1) will not be used.  It just establishes sumLogProb as a stochastic node.
})

###############
## Constants ##
###############
nObs         = vector("numeric", length = length(psyllids))
nSteps       = vector("numeric", length = length(psyllids))
iMeteoForObs = vector("list",    length = length(psyllids))
for (tree in 1:length(iMeteoForObs)) {
  iMeteoForObs[[tree]] = sapply(psyllids[[tree]]$date, function(x) which(abs(meteo$date-x) == min(abs(meteo$date-x)))) %>% unlist()
  nObs[tree]           = length(iMeteoForObs[[tree]])
  nSteps[tree]         = length(min(iMeteoForObs[[tree]]):max(iMeteoForObs[[tree]]))
}
iMeteoForObsMat = matrix(NA,nrow = nTrees, ncol = max(nObs))
for (tree in 1:nTrees){
  iMeteoForObsMat[tree,1:nObs[tree]] = iMeteoForObs[[tree]]
}
stagesDev = c("egg", "L1", "L2", "L3", "L4", "L5")
stagesTot = c("egg", "L1", "L2", "L3", "L4", "L5", "imago")

Const             = list(
  SDmodel         = SDmodel,
  res             = (res        <- 25),                ## Resolution of within-stage development
  nTrees          = (nTrees     <- length(psyllids)),
  nStagesDev      = (nStagesDev <- length(stagesDev)), ## Number of developing stages (without imago)
  nStagesTot      = (nStagesTot <- length(stagesTot)), ## Total numer of stages (includes imago)
  tempMin         = (tempMin    <- min(meteo$temperature, na.rm = TRUE)),
  tempMax         = (tempMax    <- max(meteo$temperature, na.rm = TRUE)),
  tempVec         = (tempVec    <- tempMin:tempMax),
  lTempVec        = (lTempVec   <- length((tempVec <- tempMin:tempMax))),
  lMeteo          = (lMeteo     <- nrow(meteo)),
  meteoTemp       = meteo$temperature,
  iMeteoTemp      = sapply(meteo$temperature, function(x) which(x == tempVec)),
  iMeteoForObsMat = iMeteoForObsMat,
  nSteps          = nSteps
)

##############################################
## Initial (prior to MCMC) parameter values ##
##############################################
previousMCMCfile = here("APT/Jun-17_06-54-03_2021_Temps8.txt")
previous = read.table(previousMCMCfile, header=TRUE)
previous = previous[-round(1:nrow(previous)/2),]
##
previousSampScale = previous
previousSampScale[,-grep("T", colnames(previous))] = logit(previousSampScale[,-grep("T", colnames(previous))])
##
previous = previous %>% tail(1)

Inits = list(
  ## Last values retrned from a previous run
  Tmin                = previous %>% select(grep("Tmin",          colnames(previous))) %>% as.numeric(),
  Tmax                = previous %>% select(grep("Tmax",          colnames(previous))) %>% as.numeric(),
  logit_amplitudeMean = previous %>% select(grep("amplitudeMean", colnames(previous))) %>% logit() %>% as.numeric(),
  logit_shapeMean     = previous %>% select(grep("shapeMean",     colnames(previous))) %>% logit() %>% as.numeric(),
  logit_amplitudeSD   = previous %>% select(grep("amplitudeSD",   colnames(previous))) %>% logit() %>% as.numeric(),
  logit_shapeSD       = previous %>% select(grep("shapeSD",       colnames(previous))) %>% logit() %>% as.numeric()
  ### sumLogProb          = 0
)


################
## Model data ##
################
psyllidsArray = array(NA, dim = c(nTrees, max(nObs), nStagesTot))
psyllidsTotal = matrix(NA, nrow=nTrees, ncol=max(nObs))
for (tree in 1:nTrees) {
  psyllidsArray[tree, 1:nObs[tree], 1:nStagesTot] <- psyllids[[tree]] %>% select(all_of(stagesTot)) %>% as.matrix()
  psyllidsTotal[tree,1:nObs[tree]]                <- rowSums(psyllidsArray[tree, 1:nObs[tree], 1:nStagesTot])
}

Data = list(psyllids = psyllidsArray, psyllidsTotal = psyllidsTotal)

#####################################
## Build R version of nimble model ##
#####################################
rPsyllid = nimbleModel(psyllidCode, const=Const, init=Inits, data=Data, calculate = FALSE)

################
## Node lists ##
################
dataNodes    = rPsyllid$getNodeNames(dataOnly   = TRUE)                      # The data
stochNodes   = rPsyllid$getNodeNames(stochOnly  = TRUE, includeData = FALSE) # The parameters
detNodes     = rPsyllid$getNodeNames(determOnly = TRUE)                      # Deterministic nodes
monitorNodes = rPsyllid$getParents("paras", immediateOnly = TRUE)
## Filter out sumLogProb
## dataNodes  = dataNodes[which(dataNodes!="sumLogProb")]
## stochNodes = stochNodes[which(stochNodes!="sumLogProb")]
## detNodes   = detNodes[which(detNodes!="sumLogProb")]


#############################################################################
## See scripts 'cpuTimeExperiment.R' and 'fitModel.R' for using this model ##
#############################################################################

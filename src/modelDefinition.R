#####################################################################
## This script provides the definition of the multi-stageIPM model ##
#####################################################################

## Transforming integer to date of class POSIXct
## as.POSIXct(as.integer(psyllids[[1]][1,1]), origin = "1970-01-01") == psyllids[[1]][1,1]

# rm(list=ls())
library(here)
library(tibble)
library(dplyr)
library(nimble)
library(nimbleTempDev)
library(imputeTS)
library(lubridate)
library(coda)

baseDir = here()
setwd(baseDir)
load("data/data4nimble.Rdata")
source("src/functions.R")

#data4

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
    Tmin[stage]   ~ dunif(-20,20)
    Tmax[stage]   ~ dunif(20 ,60)
    aaMean[stage] ~ dexp(0.0001)
    aaSD[stage]   ~ dexp(0.0001)
    bbMean[stage] ~ dexp(2)
    bbSD[stage]   ~ dexp(2)
    ## Briere functional response curves
    paras[stage,1:lTempVec,1] <- briere(t=tempVec[1:lTempVec], Tmin=Tmin[stage], Tmax=Tmax[stage], aa=aaMean[stage], bb=bbMean[stage]) # Mean of the development kernel
    paras[stage,1:lTempVec,2] <- briere(t=tempVec[1:lTempVec], Tmin=Tmin[stage], Tmax=Tmax[stage], aa=aaSD[stage],   bb=bbSD[stage])   # Standard deviation of the development kernel
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
        pStage[tree, obs, stage]  <- sum(states[tree, iMeteoForObsMat[tree,obs] - iMeteoForObsMat[tree,1] + 1, ((stage-1)*res+1):(stage*res)]) ## BUG WITH INDEX FOR TIME
      }
      pStage[tree, obs, nStagesTot] <- states[tree, iMeteoForObsMat[tree,obs] - iMeteoForObsMat[tree,1] + 1, (nStagesDev*res+1)]
      psyllids[tree, obs, 1:nStagesTot] ~ dmultinom(prob = pStage[tree, obs, 1:nStagesTot], size = psyllidsTotal[tree, obs])
    }
  }
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
  res             = (res        <- 3),                                 ## Resolution of within-stage development
  nTrees          = (nTrees     <- length(psyllids)),
  nStagesDev      = (nStagesDev <- length(stagesDev)), ## Number of developing stages (without imago)
  nStagesTot      = (nStagesTot <- length(stagesTot)), ## Total numer of stages (includes imago)
  tempMin         = (tempMin    <- min(meteo$temperature, na.rm = TRUE)),
  tempMax         = (tempMax    <- max(meteo$temperature, na.rm = TRUE)),
  tempVec         = (tempVec    <- tempMin:tempMax),
  lTempVec        = (lTempVec   <- length(tempVec)),
  lMeteo          = (lMeteo     <- nrow(meteo)),
  meteoTemp       = meteo$temperature,
  iMeteoTemp      = sapply(meteo$temperature, function(x) which(x == tempVec)),
  iMeteoForObsMat = iMeteoForObsMat,
  nSteps          = nSteps
)

##############################################
## Initial (prior to MCMC) parameter values ##
##############################################
Inits     = list(
  aaMean = rep(0.0001, nStagesDev),
  aaSD   = rep(0.0001, nStagesDev),
  bbMean = rep(2,      nStagesDev),
  bbSD   = rep(2,      nStagesDev),
  Tmin   = rep(0,      nStagesDev),
  Tmax   = rep(40,     nStagesDev)
)

# Model data
# psyllids[tree, obs, 1:nStagesTot] ~ dmultinom(prob = pStage[tree, obs, 1:nStagesTot])#, psyllidsTotal = sum(psyllids[tree, obs, 1:nStagesTot]))
psyllidsArray = array(NA, dim = c(nTrees, max(nObs), nStagesTot))
psyllidsTotal = matrix(NA, nrow=nTrees, ncol=max(nObs))
for (tree in 1:nTrees) {
  psyllidsArray[tree, 1:nObs[tree], 1:nStagesTot] <- psyllids[[tree]] %>% select(all_of(stagesTot)) %>% as.matrix()
  psyllidsTotal[tree,1:nObs[tree]]                <- rowSums(psyllidsArray[tree, 1:nObs[tree], 1:nStagesTot])
}

Data = list(psyllids      = psyllidsArray,
            psyllidsTotal = psyllidsTotal)

# Build R version of nimble model
rPsyllid = nimbleModel(psyllidCode, const=Const, init=Inits, data=Data)

# Compile model to C++
cPsyllid = compileNimble(rPsyllid)

dataNodes  = cPsyllid$getNodeNames(dataOnly = TRUE)                       # The data
stochNodes = cPsyllid$getNodeNames(stochOnly = TRUE, includeData = FALSE) # The parameters
detNodes   = cPsyllid$getNodeNames(determOnly = TRUE)                     # Deterministic nodes

#simulate(rPsyllid, detNodes)
system.time(simulate(cPsyllid, detNodes))
##  user  system elapsed
## 0.667   0.000   0.667 ## This is great!!! I feared it could be much much slower.
## 1.71    0.03    1.73 for Walid




##########################
##########################
##########################
##########################
## DEBUGGING TO DO HERE ##
##########################
##########################
##########################
##########################
calculate(cPsyllid)                   ## NA - why?
calculate(cPsyllid, nodes=stochNodes) ## These work
calculate(cPsyllid, nodes=detNodes)   ## These work
cPsyllid$pStage                       ## tree, obs, stage
cPsyllid$states                       ##
cPsyllid$paras   ## Are these also NAs? Are their values too high for the devKernel ????
cPsyllid$devKernel
cPsyllid$states[1,2,]

cPsyllid$pStage[1,,]
cPsyllid$states[1,,]

## states[tree, time, subStage]
stateNodes = cPsyllid$getNodeNames()[grep("states", cPsyllid$getNodeNames())]
cPsyllid
simulate(cPsyllid, nodes=stateNodes)
calculate(cPsyllid, nodes=stateNodes)
getLogProb(cPsyllid, nodes=stateNodes)

tree = 1
time = 386
undebug(sparseTWstep)
(cPsyllid$states[tree, time+1, 1:(Const$nStagesDev*Const$res+1)] = sparseTWstep(cPsyllid$states[tree, time, 1:(Const$nStagesDev*Const$res+1)], cPsyllid$devKernel[1:Const$nStagesDev, Const$iMeteoTemp[Const$iMeteoForObsMat[tree,1] + time - 2], 1:(Const$res+1)])); time = time + 1

cPsyllid$states[tree,1:400,] # <- rPsyllid$states[tree,1:400,]
psyllids[[1]][,1]


###############################
## Next steps... set up MCMC ##
###############################

mcmcConf <- configureMCMC(model=rPsyllid, monitors=stochNodes) ## Sets a basic default MCMC configuration (I almost always end up adding block samplers to overcome problems/limitations with the dfault configguration)
mcmcConf$printSamplers() ## All univariate samplers. We'll probably have strong autocorrelation in the samples
mcmcConf$getMonitors()
mcmcConf$getMonitors2()
## Optionally ## mcmcConf$removeSamplers()
## Optionally ## mcmcConf$addSampler( Fill this gap)
## Build and compile the MCMC
Rmcmc = buildMCMC(mcmcConf)
Cmcmc = compileNimble(Rmcmc)

##################
## Run the MCMC ##
##################
nIter = 1000
STime <- run.time(Cmcmc$run(nIter, reset=FALSE))  ## 5.7 minutes for 1000 iterations -> we can do 100000 iterations over night, or 1E6 iterations in 5 days

samples <- coda::as.mcmc(as.matrix(Cmcmc$mvSamples))
summary(samples)

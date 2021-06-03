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
    Tmin[stage]        ~ dnorm( 0,20)
    Tmax[stage]        ~ dnorm(40,20)
    logit(amplitudeMean[stage]) ~ dLogitBeta(1,1)
    logit(amplitudeSD[stage])   ~ dLogitBeta(1,1)
    logit(shapeMean[stage])     ~ dLogitBeta(1,1)
    logit(shapeSD[stage])       ~ dLogitBeta(1,1)
    ## Briere functional response curves
    paras[stage,1:lTempVec,1] <- stBriere(t=tempVec[1:lTempVec], Tmin=Tmin[stage], Tmax=Tmax[stage], shape=shapeMean[stage], amplitude=amplitudeMean[stage]) # Mean of the development kernel
    paras[stage,1:lTempVec,2] <- stBriere(t=tempVec[1:lTempVec], Tmin=Tmin[stage], Tmax=Tmax[stage], shape=shapeSD[stage],   amplitude=amplitudeSD[stage])   # Standard deviation of the development kernel
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
  ## A proxy node for returning logProbs
  sumLogProb ~ dnorm(0,1) ## dnorm(0,1) will not be used.  It just establishes sumLogProb as a stochastic node.
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
Inits     = list(
  ## logit_amplitudeMean = logit(rep(0.08, nStagesDev)),
  ## logit_amplitudeSD   = logit(rep(0.1, nStagesDev)),
  ## logit_shapeMean     = logit(rep(0.8, nStagesDev)),
  ## logit_shapeSD       = logit(rep(0.8, nStagesDev)),
  ## Tmin                = rep(0,  nStagesDev),
  ## Tmax                = rep(40, nStagesDev),
  ## Values found via optim
  Tmin                = c(-0.840082505,1.890155937,0.964868738,-1.208933150,-1.894897093,1.074007783),
  Tmax                = c(44.145668347,42.015194827,41.383560659,42.012669226,41.961821721,41.866991422),
  logit_amplitudeMean = c(-2.538223744,-2.346815670,-2.148938115,-2.274139735,-2.109692812,-3.886697986),
  logit_shapeMean     = c(1.243808310,1.306726311,2.392200114,3.278973805,5.199365553,4.513076810),
  logit_amplitudeSD   = c(-0.379437151,0.835470536,1.410043012,2.437749389,2.159332219,-0.006822573),
  logit_shapeSD       = c(1.085060138,8.606434971,-1.124921897,0.517501802,0.302038131,7.553068749),
  sumLogProb          = 0
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
            psyllidsTotal = psyllidsTotal
            )

#####################################
## Build R version of nimble model ##
#####################################
rPsyllid = nimbleModel(psyllidCode, const=Const, init=Inits, data=Data, calculate = FALSE)

##########################
## Compile model to C++ ##
##########################
cPsyllid = compileNimble(rPsyllid)

################
## Node lists ##
################
dataNodes    = cPsyllid$getNodeNames(dataOnly   = TRUE)                      # The data
stochNodes   = cPsyllid$getNodeNames(stochOnly  = TRUE, includeData = FALSE) # The parameters
detNodes     = cPsyllid$getNodeNames(determOnly = TRUE)                      # Deterministic nodes
monitorNodes = cPsyllid$getParents("paras", immediateOnly = TRUE)
## Filter out sumLogProb
dataNodes  = dataNodes[which(dataNodes!="sumLogProb")]
stochNodes = stochNodes[which(stochNodes!="sumLogProb")]
detNodes   = detNodes[which(detNodes!="sumLogProb")]

####################################
## Initialise deterministic nodes ##
####################################
system.time(simulate(cPsyllid, detNodes))
##  user  system elapsed
## 0.667   0.000   0.667 ## This is great!!! I feared it could be much much slower.
## 1.71    0.03    1.73 for Walid

#########################################
## Test that the log-likelihood is finite
#########################################
# cPsyllid$tempVec = tempVec ## For some reason this vector gets set to silly values, so here we re-initialise
calculate(cPsyllid, c(stochNodes,dataNodes))
calculate(cPsyllid, "sumLogProb")
if (!is.finite(calculate(cPsyllid)))
  stop("Non-finite likelihood detected.")


##############################################
## Identify better initial parameter values ##
##############################################
if (FALSE) { # This step takes about 1/2 hour, the output has been pasted into the definition of 'Inits' above
  pVec = c(cPsyllid$Tmin, cPsyllid$Tmax, cPsyllid$logit_amplitudeMean, cPsyllid$logit_shapeMean, cPsyllid$logit_amplitudeSD, cPsyllid$logit_shapeSD)
  opt = optim(pVec, function(pVec){
    cPsyllid$Tmin                = pVec[1:6]
    cPsyllid$Tmax                = pVec[6 + 1:6]
    cPsyllid$logit_amplitudeMean = pVec[12 + 1:6]
    cPsyllid$logit_shapeMean     = pVec[18 + 1:6]
    cPsyllid$logit_amplitudeSD   = pVec[24 + 1:6]
    cPsyllid$logit_shapeSD       = pVec[30 + 1:6]
    calculate(cPsyllid, c(stochNodes,dataNodes))
  }, control = list(fnscale=-1), hessian = FALSE)
}



################################
## Plot the the Briere curves ##
################################
if (FALSE) {
  par(mfrow=n2mfrow(Const$nStagesDev*2))
  for (stage in 1:Const$nStagesDev) {
    curve(briere(t=x, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]), 0, 50, ylab="Mean", main=paste("stage",stage))
  }
  for (stage in 1:Const$nStagesDev) {
    curve(briere(t=x, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeSD[stage],   shape=cPsyllid$shapeSD[stage]), 0, 50, ylab="SD", main=paste("stage",stage))
  }
}

###############################
## Next steps... set up MCMC ##
###############################
thin = 10
mcmcConf <- configureMCMC(model=rPsyllid, monitors=monitorNodes, monitors2 = "sumLogProb", thin=thin, thin2=thin) ## Sets a basic default MCMC configuration (I almost always end up adding block samplers to overcome problems/limitations with the dfault configguration)
mcmcConf$getMonitors()
mcmcConf$printSamplers() ## All univariate samplers. We'll probably have strong autocorrelation in the samples
mcmcConf$removeSamplers()
configureStoreLogProb(mcmcConf, cPsyllid, 'sumLogProb') ## For tracking posterior log-likelihood
mcmcConf$addSampler(target=stochNodes, type="RW_block", control=list(scale=0.1))
mcmcConf

## Build and compile the MCMC
Rmcmc = buildMCMC(mcmcConf)
Cmcmc = compileNimble(Rmcmc)

# 1000 iterations in 2.5 hours on laptop with univariate samplers
# 1000 iterations in 5.4 minutes on workstation with 1 block sampler


##################
## Run the MCMC ##
##################
nIter =  1E4 # 1E5 ~ 9.7 hours on workstation
STime <- run.time(Cmcmc$run(nIter, reset=FALSE))  ## 5.7 minutes for 1000 iterations -> we can do 100000 iterations over night, or 1E6 iterations in 5 days

#####################################
## Extract samples and save tofile ##
#####################################
samples <- as.matrix(Cmcmc$mvSamples)
samples <- coda::as.mcmc(samples[!(is.na(samples[,1])),])

## summary(samples)
## plot(samples)

(fileName = paste0("MCMC/","mcmc",(date() %>% strsplit(" "))[[1]][c(2,4)] %>% paste(collapse=""),".txt"))
write.table(as.matrix(samples), file=fileName, row.names = FALSE)
## samples=read.table("MCMC/mcmcJun2.txt", header = TRUE)

if (FALSE) {
  ## These attempts to print to file are not working well... probably either too many samples, too many variables, or both.
  (fileName = paste0("MCMC/","mcmc",(date() %>% strsplit(" "))[[1]][c(2,4)] %>% paste(collapse=""),"_crosscorr.pdf"))
  pdf(file=fileName, width=25, height=25)
  crosscorr.plot(samples)
  dev.off()

  (fileName = paste0("MCMC/","mcmc",(date() %>% strsplit(" "))[[1]][c(2,4)] %>% paste(collapse=""),"_codaPlot.pdf"))
  pdf(file=fileName)
  plot(samples)
  dev.off()
}

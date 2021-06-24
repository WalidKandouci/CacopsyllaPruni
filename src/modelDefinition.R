#####################################################################
## This script provides the definition of the multi-stageIPM model ##
#####################################################################
## install.packages(c("here", "tibble", "dplyr", "nimble", "imputeTS", "lubridate", "coda")
## remotes::install_git(url="https://github.com/DRJP/nimbleAPT.git", subdir="nimbleAPT", build_vignettes = FALSE)

## source(here::here("src/modelDefinition.R"))

## Transforming integer to date of class POSIXct
## as.POSIXct(as.integer(psyllids[[1]][1,1]), origin = "1970-01-01") == psyllids[[1]][1,1]

# rm(list=ls())
library(here)          ## for locating a suitable base directory
library(tibble)        ## for tibbles
library(dplyr)         ## for piping %>%
library(nimble)        ## for working with BUGS type models in R
library(nimbleAPT)     ## for adaptive parallel tempering
library(nimbleTempDev) ## for temperature depandant IPM
library(lubridate)     ## for dates and times
library(coda)          ## for mcmc diagnostics

baseDir = here()
setwd(baseDir)
load("data/data4nimble.Rdata")
## source("src/functions.R")   ## Currently not using any of these funcitons

#data4

if (!exists("setConstantsElsewhere")) { ## This flag permits other scripts (which source this script) to set the following constants
  SDmodel = 1 # 2, 3, 4, 5
}


if (is.element("package:imputeTS", search())) {
    library(imputeTS)      ## for na_kalman
    ## Impute missing values
    imp <- na_kalman(meteo$temperature) # use the Kalman filter  to imput our missing values
                                        #ggplot_na_distribution(meteo$temperature) # some nice plot
                                        #ggplot_na_imputations(meteo$temperature, imp)
    meteo$temperature <- round(imp) # NA free meteo dataset
} else {
    ## A temporary hack because imputeTS won't install on migale
    meteo$temperature[is.na(meteo$temperature)] = c(10, 10, 9)
}


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
    beta0[stage]                ~ dnorm(0, tau=tauBeta0)
    if (SDmodel > 2) {
      beta1[stage] ~ ddexp(location=0, scale=scaleBeta1) ## Constrained version of Bhattacharya's "Dirichlet–Laplace Priors for Optimal Shrinkage"
    }
    if (SDmodel > 4) {
      beta2[stage] ~ ddexp(location=0, scale=scaleBeta2) ## Constrained version of Bhattacharya's "Dirichlet–Laplace Priors for Optimal Shrinkage"
    }
    ## Mean development as a function of temperature
    paras[stage,1:lTempVec,1] <- stBriere(T=tempVec[1:lTempVec], Tmin=Tmin[stage], Tmax=Tmax[stage], shape=shapeMean[stage], amplitude=amplitudeMean[stage]) #
    for (iTemp in 1:lTempVec) { # iTemp = index for temperature
      ## Within stage development at temperature tempVec[iTemp]
      devKernel[stage,iTemp,1:(res+1)] <- getKernel(paras=paras[stage,iTemp,1:3], res=res, devFunction = 1) ## Package currently has functions getM, setM and setMultiM... but we should write a function to just return the first column of getM and work with that (because the model matrix over many stages is very sparse).
      ## Survival
      paras[stage,iTemp,3] <- 1
      ## Standard deviation in development as a function of temperature
      if (SDmodel == 1) {
        paras[stage,iTemp,2] <- paras[stage,iTemp,1] * exp(beta0[stage])
      } else if (SDmodel == 2) {
        paras[stage,iTemp,2] <- (paras[stage,iTemp,1] > 0) * exp(beta0[stage])
      } else if (SDmodel == 3) {
        paras[stage,iTemp,2] <- paras[stage,iTemp,1] * exp(beta0[stage] + beta1[stage]*tempVec[iTemp])
      } else if (SDmodel == 4) {
        paras[stage,iTemp,2] <- (paras[stage,iTemp,1] > 0) * exp(beta0[stage] + beta1[stage]*tempVec[iTemp])
      } else if (SDmodel == 5) {
        paras[stage,iTemp,2] <- paras[stage,iTemp,1] * exp(beta0[stage] + beta1[stage]*tempVec[iTemp] + beta2[stage]*tempVec[iTemp]^2)
      } else if (SDmodel == 6) {
        paras[stage,iTemp,2] <- (paras[stage,iTemp,1] > 0) * exp(beta0[stage] + beta1[stage]*tempVec[iTemp] + beta2[stage]*tempVec[iTemp]^2)
      }
    }
  }
  ## Hyper-parameters for standard deviation of development models
  if (SDmodel > 2) {
    scaleBeta1 ~ dgamma(shape=1/nStagesDev, rate=1/2)
  }
  if (SDmodel > 4) {
    scaleBeta2 ~ dgamma(shape=1/nStagesDev, rate=1/2)
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

(res        <- 25)                ## Resolution of within-stage development
(nTrees     <- length(psyllids))
(nStagesDev <- length(stagesDev)) ## Number of developing stages (without imago)
(nStagesTot <- length(stagesTot)) ## Total numer of stages (includes imago)
(tempMin    <- min(meteo$temperature, na.rm = TRUE))
(tempMax    <- max(meteo$temperature, na.rm = TRUE))
(tempVec    <- tempMin:tempMax)
(lTempVec   <- length(tempVec))
(lMeteo     <- nrow(meteo))

Const             = list(
  SDmodel         = SDmodel,
  res             = res,
  nTrees          = nTrees,
  nStagesDev      = nStagesDev,
  nStagesTot      = nStagesTot,
  lTempVec        = lTempVec,
  lMeteo          = lMeteo,
  meteoTemp       = meteo$temperature,
  iMeteoTemp      = sapply(meteo$temperature, function(x) which(x == tempVec)),
  iMeteoForObsMat = iMeteoForObsMat,
  nSteps          = nSteps,
  tauBeta0        = 1E-11
)

##############################################
## Initial (prior to MCMC) parameter values ##
##############################################
fileStem         = paste0("model",SDmodel, "_")
previousMCMCfile = here("APT/Jun-18_19-38-37_2021_Temps4.txt")
previous         = read.table(previousMCMCfile, header=TRUE)
previous         = previous[-round(1:nrow(previous)/2),] # Remove 1/2 the samples as burn-in
##
previousSampScale = previous
previousSampScale[,-grep("T", colnames(previous))] = logit(previousSampScale[,-grep("T", colnames(previous))])
##
previous = previous %>% tail(1)

Inits = list(
  ## Last values retrned from a previous run
  Tmin                = previous %>% select(grep("Tmin",          colnames(previous))) %>% as.numeric(),
  Tmax                = previous %>% select(grep("Tmax",          colnames(previous))) %>% as.numeric(),
  ## For SDmodel == 1
  logit_amplitudeMean = previous %>% select(grep("amplitudeMean", colnames(previous))) %>% logit() %>% as.numeric(),
  logit_shapeMean     = previous %>% select(grep("shapeMean",     colnames(previous))) %>% logit() %>% as.numeric(),
  ## For SDmodel == 2
  beta0      = rep(1, nStagesDev),     # 1 for models 1 3 5
  beta1      = rep(1E-11, nStagesDev),
  beta2      = rep(1E-11, nStagesDev),
  scaleBeta1 = (1/nStagesDev) / (1/2), # Mean = shape / rate (gamma distribution)
  scaleBeta2 = (1/nStagesDev) / (1/2)  # Mean = shape / rate (gamma distribution)
)

if (is.element(SDmodel, c(2,4,6)))
    Inits$beta0 = rep(-1, nStagesDev) # -2 for models 2 & 6



################
## Model data ##
################
psyllidsArray = array(NA, dim = c(nTrees, max(nObs), nStagesTot))
psyllidsTotal = matrix(NA, nrow=nTrees, ncol=max(nObs))
for (tree in 1:nTrees) {
  psyllidsArray[tree, 1:nObs[tree], 1:nStagesTot] <- psyllids[[tree]] %>% select(all_of(stagesTot)) %>% as.matrix()
  psyllidsTotal[tree,1:nObs[tree]]                <- rowSums(psyllidsArray[tree, 1:nObs[tree], 1:nStagesTot])
}

Data = list(psyllids      = psyllidsArray,
            psyllidsTotal = psyllidsTotal,
            tempVec       = tempVec)

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

stochNodesUnique = unique(sub("\\[.*","",stochNodes))

# depNodes     = rPsyllid$getDependencies(stochNodes, self=FALSE, includeData = FALSE)
# detNodes %>% sub(pat="\\[.*", rep="") %>% unique() %>% sort()
# depNodes %>% sub(pat="\\[.*", rep="") %>% unique() %>% sort() # identical

#############################################################################
## See scripts 'cpuTimeExperiment.R' and 'fitModel.R' for using this model ##
#############################################################################

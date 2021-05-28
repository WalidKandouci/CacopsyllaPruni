#####################################################################
## This script provides the definition of the multi-stageIPM model ##
#####################################################################
library(here)
library(tibble)
library(dplyr)
library(nimble)
library(nimbleTempDev)
library(imputeTS)
library(lubridate)

baseDir = here()
setwd(baseDir)
load("data/data4nimble.Rdata")
# import imputed values & replace NA
imp <- na_kalman(meteo$temperature) # use the Kalman filter  to imput our missing values
#ggplot_na_distribution(meteo$temperature) # some nice plot
#ggplot_na_imputations(meteo$temperature, imp)
meteo$temperature <- round(imp) # NA free meteo dataset

# Construct stepForObs
#nplants <- length(psyllids)
#nobs <- 0
#for (ii in 1:nplants) {
# nobs[ii] <- length(psyllids[[ii]][,1])
#  }
#max(nobs)
#
#date_bidon <- ymd_hms("19940205120000")
#stepForObs <- rep(list(date_bidon),max(nobs))
#stepForObs <- data.frame(date_bidon,date_bidon)
#
#for (ii in 1:nplants) {
#  stepForObs[1:nobs[ii],ii] <-  psyllids[[ii]][,1]
#}
#
#class(stepForObs[1,1])
#colnames(stepForObs) <- paste0("plant", formatC(1:nplants, width = 2, flag = "0"))


# BUGS code for nimble model
# Just a very rough draft - there are many details which need adding or refining
psyllidCode <- nimbleCode ({
  #############################################################
  ## Development kernel at each temperature & for each stage ##
  #############################################################
  for (stage in 1:nStages) { # iStag = index for {, L1, L2, L3, L4, L5, imago}
    ## Priors
    Tmin[stage] ~ dunif(-20,20)
    Tmax[stage] ~ dunif(20 ,60)
    aaMean[stage] ~ dexp(0.0001)
    aaSD[stage]   ~ dexp(0.0001)
    bbMean[stage] ~ dexp(2)
    bbSD[stage]   ~ dexp(2)
    for (iTemp in 1:lTempVec) { # iTemp = index for temperature
      ## Briere
      paras[stage,iTemp,1] <- briere(t=tempVec[iTemp], Tmin=Tmin[stage], Tmax=Tmax[stage], aa=aaMean[stage],   bb=bbMean[stage])   # mean of the development kernel
      paras[stage,iTemp,2] <- briere(t=tempVec[iTemp], Tmin=Tmin[stage], Tmax=Tmax[stage], aa=aaSD[stage], bb=bbSD[stage]) # standard deviation for the development kernel
      paras[stage,iTemp,3] <- 1
      ## Possibly add a parameter transformation step here ???
      devKernel[stage,iTemp,1:res] <- getKernel(paras=paras[stage,iTemp,1:3], res=res, devFunction = 1) ## Package currently has functions getM, setM and setMultiM... but we should write a function to just return the first column of getM and work with that (because the model matrix over many stages is very sparse).
    }
  }
  #######################################################################
  ## Loop over trees 
  #######################################################################
  for (tree in 1:nTree) { # Adding multiple trees means running the IPM seperately for each tree (due to different start dates)
    # IPM projections
    for (tStep in 1:nTreeSteps[tree]) { # tStep = index for time-step
      iTemp <- iMeteoTemp[iMeteoForObsMat[tree,1] + tStep - 1]
      states[tree, tStep+1, 1:(nStages*res+1)] <- sparseTWstep(states[tree, tStep, 1:(nStages*res+1)],devKernel[iStage,iTemp,1:res])
    }
    # Likelihood
    for (obs in 1:nTreeDates[tree]) {
      for(stage in 1:nStages){
        pStage[tree, obs, stage]  <- sum(states[tree, iMeteoForObsMat[tree,obs], ((stage-1)*res+1):(stage*res)])
      }
      pStage[tree, obs, nStages1] <- states[tree, iMeteoForObsMat[tree,obs], (stage*res+1)]
      psyllids[tree, obs, 1:nStages1] ~ dmultinom(prob = pStage[tree, obs, 1:nStages1])#, size = sum(psyllids[tree, obs, 1:nStages1]))
      ## 1B is a magic number - we need to generalise some how - possibly via a ragged array -
      ## the prob vector will come from the IPM
    }
  }
})
## TO DO: add ragged array to importData.R

#########################
## Create nimble model ##
#########################
# Constants
tempMin  = min(meteo$temperature, na.rm = TRUE)
tempMax  = max(meteo$temperature, na.rm = TRUE)
tempVec  = tempMin:tempMax
lTempVec = length(tempVec)

iMeteoForObs = vector("list",length = length(psyllids))
nTreeDates =  vector("numeric", length = length(psyllids))
nTreeSteps = vector("numeric", length = length(psyllids))

for (tree in 1:length(iMeteoForObs)) {
  iMeteoForObs[[tree]] = sapply(psyllids[[tree]]$date, function(x) which(abs(meteo$date-x) == min(abs(meteo$date-x)))) %>% unlist()
  nTreeDates[tree] = length(iMeteoForObs[[tree]])
  nTreeSteps[tree] = length(min(iMeteoForObs[[tree]]):max(iMeteoForObs[[tree]]))
}

iMeteoForObsMat = matrix(NA,nrow = nTrees, ncol = max(nTreeDates))

for (tree in 1: nTrees){
  iMeteoForObsMat[tree,1:nTreeDates[tree]] = iMeteoForObs[[tree]] 
}

Const    = list(
  nTree      = length(psyllids),       ## Starting simple, just tree 1B to begin with
  res        = 3,                      ## We can increase the resolution once some rough code is working
  nStages    = 6,                      ## number of developing stages (without imago)
  nStages1   = 7,                      ## with imago
  tempMin    = tempMin,
  tempMax    = tempMax,
  tempVec    = tempVec,
  lTempVec   = lTempVec,
  lMeteo     = nrow(meteo),
  meteoTemp  = meteo$temperature,
  iMeteoTemp = sapply(meteo$temperature, function(x) which(x == tempVec)),
  iMeteoForObsMat = iMeteoForObsMat
  
  # Index of nearest meteo observation to each psyllid observation
  # Check for iDate_1B
  # data.frame(psyllid_date = psyllids[[which(treeNames=="1B")]]$date, nearest_meteo_date = meteo$date[Const$iDate_1B])
)

# Initial values for model parameters - tobe updated via MCMC
Inits = list(
  aaMean = rep(0.0001, Const$nStages),
  bbMean = rep(2, Const$nStages),
  aaMeanSD = rep(0.0001, Const$nStages),
  bbMeanSD = rep(2, Const$nStages),
  Tmin = rep(0,Const$nStages),
  Tmax = rep(40,Const$nStages)
)

# Model data
# psyllids[tree, obs, 1:nStages1] ~ dmultinom(prob = pStage[tree, obs, 1:nStages1])#, size = sum(psyllids[tree, obs, 1:nStages1]))
psyllidsArray = array(NA,dim = c(nTrees, max(sapply(psyllids, "nrow")), Const$nStages1))
Data = list(temperature = meteo$temperature,
            psyllids = psyllidsArray
            )

# Build R version of nimble model
rPsyllid = nimbleModel(psyllidCode, const=Const, init=Inits)#, data=Data)

# Compile model to C++
cPsyllid = compileNimble(rPsyllid)

seq(0,40,0.01)

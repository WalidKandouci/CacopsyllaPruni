#####################################################################
## This script provides the definition of the multi-stageIPM model ##
#####################################################################
library(here)
library(tibble)
library(dplyr)
library(nimble)
library(nimbleTempDev)

baseDir = here()
setwd(baseDir)
load("data/data4nimble.Rdata")

# BUGS code for nimble model
# Just a very rough draft - there are many details which need adding or refining
psyllidCode <- nimbleCode ({
  #############################################################
  ## Development kernel at each temperature & for each stage ##
  #############################################################
  for (iStage in 1:nStages) { # iStag = index for {egg, L1, L2, L3, L4, L5, imago}
    ## Priors
    aaMean[iStage]   ~ dWhatEver()
    aaMeanSD[iStage] ~ dWhatEver()
    bbMean[iStage]   ~ dWhatEver()
    bbMeanSD[iStage] ~ dWhatEver()
    for (iTemp in 1:lTempVec) { # iTemp = index for temperature
      parasEgg[iStage,iTemp,1] <- briere(t=tempVec[iTemp], Tmin=Tmin[iStage], Tmax=Tmax[iStage], aa=aaMean[iStage],   bb=bbMean[iStage])   # The mean of the development kernel
      parasEgg[iStage,iTemp,2] <- briere(t=tempVec[iTemp], Tmin=Tmin[iStage], Tmax=Tmax[iStage], aa=aaMeanSD[iStage], bb=bbMeanSD[iStage]) # The mean + 1 standard deviation for the development kernel
      ## Possibly add a parameter transformation step here ???
      devKernelEgg[iStage,iTemp,1:res] <- getMcol1(paras=parasEgg, res=res, devFunction = 1) ## Package currently has functions getM, setM and setMultiM... but we should write a function to just return the first column of getM and work with that (because the model matrix over many stages is very sparse).
    }
  }
  #######################################################################
  ## Loop over trees - this is tricky, so lets just start with tree 1B ##
  #######################################################################
  for (tree in 1:nTree) { # Adding multiple trees means running the IPM seperately for each tree (due to different start dates)
    # IPM projections
    for (tStep in 1:lTempVec) { # tStep = index for time-step
      IPMouput[iDate_1B, 1:nStages] <- newSparseTravellingWaveFunction()
    }
    # Likelihood
    for (iObs in 1:nObs[tree]) {
      psyllids1B[iObs, 1:nStages] ~ dmultinom(prob = IPMouput[iDate_1B, 1:nStages], )
      ## 1B is a magic number - we need to generalise some how
      ## the prob vector will come from the IPM
    }
  }
})

# Lists for nimble model: constants, initial values and data
tempMin  = min(meteo$temperature, na.rm = TRUE)
tempMax  = max(meteo$temperature, na.rm = TRUE)
tempVec  = tempMin:tempMax
lTempVec = length(tempVec)
Const    = list(
  nTree    = 1,       ## Starting simple, just tree 1B to begin with
  res      = 25,      ## We can increase the resolution once some rough code is working
  tempMin  = tempMin,
  tempMax  = tempMax,
  tempVec  = tempVec,
  lTempVec = lTempVec,
  lMeteo   = nrow(meteo),
  # Index of nearest meteo observation to each psyllid observation
  iDate_1B = sapply(psyllids[[which(treeNames=="1B")]]$date, function(x) which(abs(meteo$date-x) == min(abs(meteo$date-x)))) %>% unlist()
  # Check for iDate_1B
  # data.frame(psyllid_date = psyllids[[which(treeNames=="1B")]]$date, nearest_meteo_date = meteo$date[Const$iDate_1B])
)

Inits = list( )

Data = list(temperature = meteo$temperature,
            psyllids1B  = (psyllids[[which(treeNames=="1B")]][-1, c("œuf","L1","L2","L3","L4","L5","imago")])
            )

# Create model object
rPsyllid = nimbleModel(psyllidCode, const=Const, init=Inits, data=Data)

# Compile model
cPsyllid = compileNimble(rPsyllid)

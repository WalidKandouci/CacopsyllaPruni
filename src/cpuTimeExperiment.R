########################
## Set some constants ##
########################
SDmodel = 1 # 2, 3, 4, 5     ## Identifywhich model to use for SD
nTemps  = 12                 ## Number of temperatures in APT samplers
setConstantsElsewhere = TRUE ## Prevents a redefinition in modelDefinition.R

###########################
## Create rPsyllid model ##
###########################
source(here::here("src/modelDefinition.R"))

##########################
## Compile model to C++ ##
##########################
cPsyllid = compileNimble(rPsyllid)

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
## calculate(cPsyllid, "sumLogProb")
if (!is.finite(calculate(cPsyllid)))
  stop("Non-finite likelihood detected.")



#####################################
## CPU time experiment: res x time ##
#####################################
resVec <- seq(25,100,25) # with niter fix
dataResTime <- cbind(resVec,rep(0,length(resVec)))
colnames(dataResTime) <- c("res","time")
for (ii in 1: length(resVec)) {
  rPsyllidRes = nimbleModel(psyllidCode, const=list(
    SDmodel         = SDmodel,
    res             = (res        <- resVec[ii]),                ## Resolution of within-stage development
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
  ), init=Inits, data=Data, calculate = FALSE)
  cPsyllidRes = compileNimble(rPsyllidRes)
  ################
  ## Node lists ##
  ################
  dataNodesRes    = cPsyllidRes$getNodeNames(dataOnly   = TRUE)                      # The data
  stochNodesRes   = cPsyllidRes$getNodeNames(stochOnly  = TRUE, includeData = FALSE) # The parameters
  detNodesRes     = cPsyllidRes$getNodeNames(determOnly = TRUE)                      # Deterministic nodes
  monitorNodesRes = cPsyllidRes$getParents("paras", immediateOnly = TRUE)
  ## Filter out sumLogProb
  ## dataNodesRes  = dataNodesRes[which(dataNodesRes!="sumLogProb")]
  ## stochNodesRes = stochNodesRes[which(stochNodesRes!="sumLogProb")]
  ## detNodesRes   = detNodesRes[which(detNodesRes!="sumLogProb")]
  ####################################
  ## Initialise deterministic nodes ##
  ####################################
  dataResTime[ii,2] <- system.time(simulate(cPsyllidRes, detNodesRes))[3] # [3] is for the total time to calculate (User + CPU)
}

plot(dataResTime[,1], dataResTime[,2])

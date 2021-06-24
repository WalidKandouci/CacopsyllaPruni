## rm(list=ls())
## source(here::here("src/fitModel.R"))

########################
## Set some constants ##
########################
SDmodel = 4 # 2, 3, 4, 5 ## Identifywhich model to use for SD
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
SDmodel = 4
source(here::here("src/modelDefinition.R"))

##########################
## Compile model to C++ ##
##########################
cPsyllid = compileNimble(rPsyllid)

## cPsyllid$tempVec # Was becoming corrupted after simulation of detNodes due to a pointer bug in stBriere
## rPsyllid$tempVec #


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
calculate(cPsyllid)
calculate(cPsyllid, c(stochNodes,dataNodes))
## calculate(cPsyllid, "sumLogProb")
if (!is.finite(calculate(cPsyllid)))
  stop("Non-finite likelihood detected.")

##########################################################
## Optim for identifing better initial parameter values ##
##########################################################
nimPrint("Starting optim")
if (TRUE) { # This step takes about 1/2 hour, the output has been pasted into the definition of 'Inits' above
  for (iter in 1:3) {
    if (iter==1) {
      command = paste0("pVec=c(", paste0("cPsyllid$",stochNodesUnique,collapse=", "),")")
      eval(parse(text=command))
      names(pVec) = stochNodes
      print(pVec)
    }
    ## Run optim
    opt = optim(pVec, function(vec){
      names(vec) = stochNodes
      ## Initialise the model parameters using vec
      for (ii in 1:length(stochNodesUnique)) {
        command = paste0("cPsyllid$",stochNodesUnique[ii]," = ", paste0("vec[", paste0("c(",paste0(grep(stochNodesUnique[ii],names(vec)),collapse=","),")"), "]"))
        eval(parse(text=command))
      }
      ## Calculate posterior log likelihood
      calculate(cPsyllid)
    }, control = list(fnscale=-1, maxit=500), hessian = FALSE) # TRUE
    pVec = opt$par
    nimPrint("iter = ", iter)
    print(opt)
  }
}


cPsyllid$tempVec
rPsyllid$tempVec



## Ensure model is parameterised using optim optput
for (ii in 1:length(stochNodesUnique)) {
    command = paste0("cPsyllid$",stochNodesUnique[ii]," = ", paste0("pVec[", paste0("c(",paste0(grep(stochNodesUnique[ii],names(pVec)),collapse=","),")"), "]"))
    print(command)
    eval(parse(text=command))
}
cPsyllid$simulate(detNodes,includeData = FALSE)
nimPrint("logProbs following optim = ", calculate(cPsyllid, c(stochNodes,detNodes,dataNodes)))

## for (ii in 1:length(stochNodesUnique)) {
##   command = paste0("print(cPsyllid$",stochNodesUnique[ii],")")
##   print(stochNodesUnique[ii])
##   eval(parse(text=command))
##   command = paste0("print(cPsyllid$calculate(\"",stochNodesUnique[ii],"\"))")
##   eval(parse(text=command))
## }
## cPsyllid$simulate("paras")
## cPsyllid$paras
## cPsyllid$devKernel
## ## THE BUG - Mean and SD for stage 5 too close to zero. Why ???
## cPsyllid$paras[,,1:2]
## deleteme =
##   stage = 1
##   stBriere(T=tempVec[1:lTempVec], Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], shape=cPsyllid$shapeMean[stage], amplitude=cPsyllid$amplitudeMean[stage]) # Mean of the development kernel
## cPsyllid$simulate("paras[5,,1]")
## cPsyllid$paras[5,,1]


################################
## Plot the the Briere curves ##
################################
if (FALSE) { # TRUE
  par(mfrow=n2mfrow(Const$nStagesDev))
  for (stage in 1:Const$nStagesDev) {
    curve(stBriere(T=x, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]), -20, 60, ylab="E[Dev rate]", main=paste("stage",stage))
  }
}


if (TRUE) { # FALSE
  #####################################################################
  ## Standard MCMC. It tends to get stuck, so APT can perform better ##
  #####################################################################
  mcmcConf <- configureMCMC(model=cPsyllid,
                            monitors=stochNodes,    ## Needed for WAIC calculation
                            monitors2=monitorNodes, ## Prefered output
                            thin=thin, thin2=thin,
                            enableWAIC = TRUE)
  mcmcConf$getMonitors()
  mcmcConf$printSamplers() ## All univariate samplers. We'll probably have strong autocorrelation in the samples
  mcmcConf$removeSamplers()
  mcmcConf$addSampler(target=stochNodes, type="RW_block", control=list(scale=0.1)) # propCov=cov(previousSampScale),
  mcmcConf
  ##########################################
  ## Build and compile the MCMC algorithm ##
  ##########################################
  mcmcR <- buildMCMC(mcmcConf)
  mcmcC <- compileNimble(mcmcR)
  ##################
  ## Run the MCMC ##
  ##################
  nIter = 5E4 # 40
  RunTime <- run.time(mcmcC$run(nIter, thin = thin, thin2=thin, reset=TRUE)) ## 5.7 minutes for 1000 iterations -> we can do 100000 iterations over night, or 1E6 iterations in 5 days
  calculate(cPsyllid, c(stochNodes,dataNodes))
  ##
  samps <- tail(as.matrix(mcmcC$mvSamples), floor(nIter/thin)) ## Sampled parameters for T=1
  if (FALSE)
    plot(coda::as.mcmc(samps))
  covParas = cov(samps)
  if(any(eigen(covParas)$values <= 0)) {
    covParas = covParas + diag(diag(covParas))
    nimPrint("Adjusting covParas to ensure positive definite")
  }
  if(any(eigen(covParas)$values <= 0)) {
    covParas = "identity"
    nimPrint("Setting covParas to identity")
  }
}


########################################################################################################
## Configure an MCMC and add a block sampler from the nimbleAPT (adaptive parallel tempering) package ##
########################################################################################################
aptConf <- configureMCMC(model=rPsyllid,
                         monitors=stochNodes,     ## Needed for WAIC calculation
                         thin=thin, thin2=thin,
                         enableWAIC = TRUE)
aptConf$getMonitors()
aptConf$printSamplers() ## All univariate samplers. We'll probably have strong autocorrelation in the samples
aptConf$removeSamplers()
aptConf$addSampler(target=stochNodes, type="RW_block_tempered", control=list(scale=0.1, temperPriors=TRUE, propCov=covParas))
aptConf

####################################################
## Build and compile the APT based MCMC algorithm ##
####################################################
aptR <- buildAPT(aptConf, Temps = 1:nTemps, ULT = 1000, print= TRUE) # only 4 temperatures to avoid memory issues
aptC <- compileNimble(aptR)


###############################################################################
## Loop with short runs of APT, until mean loglik shows signs of convergence ##
###############################################################################
nIterDelta       <- 2E4 # 40  ## One iteration of the loop will take approx 1 day
TuneTemper       <- c(10, 1)  ## default value is c(10,1)
logliks          <- rnorm(nIter, cPsyllid$calculate(), 1)
logliks_previous <- logliks - rnorm(length(logliks),10, 1)
iter             <- 0
while( t.test(logliks_previous, logliks, alternative="less")$p.value < 0.05 ) {
# while( iter == 0 ) {
  iter <- iter+1
  nIter = nIterDelta * iter
  logliks_previous <- logliks
  print(paste0("iteration nb.", iter, "within while loop. meanL = ", mean(logliks_previous)))
  #################
  ## Short run 1 ##
  aptC$thinPrintTemps <- nIter / 10
  RunTime <- run.time(aptC$run(nIter,
                               reset          = TRUE,  ## Resets the adaptive MCMC. Let's proposal distributions change direction if required.
                               adaptTemps     = FALSE, ## Prevents temperature ladder adaptation (to avoid volatile behaviour when counter is reset)
                               resetTempering = TRUE,  ## Resets counter used in temperature ladder adaptation
                               printTemps     = TRUE,  ## Will print once only
                               tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2]))
  ## Update meanL
  nimPrint("While loop: 1st short run finished. RunTime (hours) = ", RunTime/60/60)
  logliks = tail(aptC$logProbs, floor(nIter/thin))
  meanL   = mean(logliks, na.rm=TRUE)
  nimPrint("meanL = ", meanL)
  #################
  ## Short run 2 ##
  aptC$thinPrintTemps <- round(nIter / 10) ## Ensures temps are only printed 10 times
  RunTime <- run.time(aptC$run(nIter,
                               reset          = FALSE, ## Do not reset the adaptive MCMC, let adaptation continue as it is
                               adaptTemps     = TRUE,  ## Allows temperature ladder to adjust
                               resetTempering = FALSE, ## Keeps the adjustments modest so avoids volatile behaviour
                               printTemps     = TRUE,  ## Prevents verbose printing of temperature ladder updates
                               tuneTemper1=TuneTemper[1], tuneTemper2=TuneTemper[2]))
  ## Update meanL
  nimPrint("While loop: 2nd short run finished. RunTime (hours) = ", RunTime/60/60)
  logliks = tail(aptC$logProbs, floor(nIter/thin))
  meanL   = mean(logliks, na.rm=TRUE)
  nimPrint("meanL = ", meanL)
  ## Calculate & print ESS
  samples   <- tail(as.matrix(aptC$mvSamples), floor(nIter/thin)) ## Sampled parameters for T=1
  mc        <- as.mcmc(samples)
  ESS       <- effectiveSize(mc)
  (ESS      <- ESS[order(ESS)])
  nimPrint("ESS: ", ESS)
  ## t-test for difference in means
  nimPrint("t-test for equivalence of means")
  print(t.test(logliks_previous, logliks, alternative="less")) # p is prob of obtaining result at least as extreme as the observed result, assuming H0 is true
  ## Generate output from aptC
  source(here::here("src/aptOutput.R"))
}
print(paste0("iteration nb.", iter, "within while loop. meanL = ", meanL))



#####################
## Long run of APT ##
#####################
nIterShort = nIter
nIter = 1E5 # 60
nimPrint("Estimated run-time (hours) = ", (RunTime/60/60) * nIter / nIterShort)
RunTime <- run.time(aptC$run(nIter, thin = 10, thin2=10, reset=FALSE, resetTempering=FALSE, adaptTemps=FALSE))
nimPrint("Run-time (hours) = ", RunTime / 60 / 60)
## SDmodel==2
#  4 temps - c(23.4)
#  8 temps - c(
# 12 temps -
# 16 temps -
# 20 temps -

## SDmodel==2
#  2 temps - 11.8 hours - 6E4 iterations

###############################
## Generate output from aptC ##
###############################
source(here::here("src/aptOutput.R"))

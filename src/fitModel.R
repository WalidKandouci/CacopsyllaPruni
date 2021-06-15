## source(here("src/fitModel.R"))

########################
## Set some constants ##
########################
SDmodel = 1 # 2, 3, 4, 5     ## Identifywhich model to use for SD
nTemps  = 4 # 8 12 16 20     ## Number of temperatures in APT samplers
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

########################################################################################################
## Configure an MCMC and add a block sampler from the nimbleAPT (adaptive parallel tempering) package ##
########################################################################################################
thin = 100
mcmcConf <- configureMCMC(model=rPsyllid, monitors=monitorNodes, ## monitors2 = "sumLogProb",
                          thin=thin, thin2=thin) ## Sets a basic default MCMC configuration (I almost always end up adding block samplers to overcome problems/limitations with the dfault configguration)
mcmcConf$getMonitors()
mcmcConf$printSamplers() ## All univariate samplers. We'll probably have strong autocorrelation in the samples
mcmcConf$removeSamplers()
mcmcConf$addSampler(target=stochNodes, type="RW_block_tempered", control=list(scale=0.1, propCov=cov(previousSampScale), temperPriors=TRUE))
mcmcConf

################################
## Build and compile the MCMC ##
################################
aptR <- buildAPT(mcmcConf, Temps = 1:nTemps, ULT = 1000, print= TRUE) # only 4 temperatures to avoid memory issues
aptC <- compileNimble(aptR)


##################
## Run the MCMC ##
##################
nIter =  6E4 # 28.8 hours with 4 temps
STime <- run.time(aptC$run(nIter, thin = 10, thin2=10, reset=TRUE)) ## 5.7 minutes for 1000 iterations -> we can do 100000 iterations over night, or 1E6 iterations in 5 days

#############################
## Extract log-likelihoods ##
#############################
STime / 60 / 60
logliks <- as.matrix(aptC$logProbs)
logliks <- coda::as.mcmc(logliks[!(is.na(logliks[,1])),])
if (FALSE)
  plot(logliks)

#####################################
## Extract samples and save tofile ##
#####################################
samples <- as.matrix(aptC$mvSamples)
samples <- coda::as.mcmc(samples[!(is.na(samples[,1])),])
summary(samples)
if (FALSE)
  plot(samples)


sort(effectiveSize(samples))
# sort(effectiveSize(samples[-(1:2500),]))

## crosscorr(samples)
## crosscorr.plot(samples)
## sort(apply(crosscorr(samples) - diag(1, 36), 1, function(x) max(abs(x))), dec=TRUE)

##########################
## Write output to file ##
##########################
(fileName = paste0("APT/", (date() %>% strsplit(" "))[[1]][c(2,4)] %>% paste0(collapse="") %>% stringr::str_replace_all(":","-"), "_",
                           (date() %>% strsplit(" "))[[1]][5] %>% substr(1,5) %>% stringr::str_replace(":",""),
                   "_Temps", nTemps,
                   ".txt"))
write.table(as.matrix(samples), file=fileName, row.names = FALSE)
## samples=read.table("APT/Jun10-20-00_2021_Temps8.txt", header = TRUE)

write.table(as.matrix(logliks), row.names = FALSE, file=sub("\\.","_loglik.",fileName))

if (TRUE) {
  ## Cross correlation plot
  pdf(file = sub("txt","pdf", sub("\\.","_crosscor.",fileName)), width=20, height=20)
  crosscorr.plot(samples)
  dev.off()
  ## Trajectories plot
  pdf(file = sub("txt","pdf", sub("\\.","_trajectories.",fileName)))
  plot(samples)
  dev.off()
  ## Log posterior likelihood trajectories
  pdf(file = sub("txt","pdf", sub("\\.","_trajecory-logliks.",fileName)))
  plot(logliks)
  dev.off()
}
## class(samples[-c(1:13000),])
## crosscorr.plot(as.mcmc(samples[-c(1:13000),]))

library(here)
library(dplyr)
library(nimbleTempDev)
source(here::here("src/modelDefinition.R"))
cPsyllid = compileNimble(rPsyllid)


###############################################
## Load APT output to analyse in this script ##
###############################################
samples  = read.table(here("APT/Jun-17_06-54-03_2021_Temps8.txt"), header=TRUE)
samples2 = read.table(here("APT/Jun-17_06-54-03_2021_Temps8_loglik.txt"), header=TRUE)

##############################
## Node lists to work with  ##
##############################
paraNames = gsub("\\..*","", colnames(samples)) %>% unique() # Names of parameter nodes
depNodes  = gsub("\\[.*","", cPsyllid$getDependencies(paraNames, self = FALSE, includeData = FALSE)) %>% unique()

#################
## Plot stages ##
#################

# The idea I have is to plot at first each stage on it's own, and maybe have the curv of the tree that we want to see in a difrent color ?
# I'm not sure of what and how to plot because "stage" is a big array and I don't know if we should plot everything ?
# and if we don't plot everything what should we really plot ?

# DP: Well first, you need to take one row of parameters from the MCMC output, plug it into the model & update and dependant nodes
#     This code chunk will probably need plugging into your loops somewhere
iMCMC     = sample(nrow(samples),1)                          # Index for a random row. Later this can be replaced by a loop.
for (ii in 1:length(paraNames)) {                            # A loop to fill each parameter node of the model that can handle alternative sets of parameter names
  eval(parse(text=paste0("cPsyllid$", paraNames[ii], " = as.numeric(samples[iMCMC, grep(paraNames[ii],colnames(samples))])")))
  ## Here's a quick check that this works
  ## cPsyllid$shapeSD
  ## samples[iMCMC, grep("shapeSD",colnames(samples))]
}
simulate(cPsyllid, depNodes) # Updates "paras" "devKernel" "states" "pStage"


nMcmcSamples = 11
pStage = array(NA, dim = c(nMcmcSamples, lMeteo, nTrees, nStagesTot)) # This will be a ragged array - i.e. nSteps[tree] is heterogeneous, so some elements of this array will remain as NAs

par(mfrow=c(4,4))
for (iMCMC in 1:nMcmcSamples) {
  for(iTree in 1:nTrees){
    for (iStage in 1:nStagesTot){
      ## plot(states[iTree, iStage, substage],type='l', col="grey", xlab="time", ylab="dev")
      ## points(states[iTree, iStage, substage],type='l', col="grey")
      ##
      ## 1) R won't find states. You need to use cPsyllid$states
      ## 2) states includes all sub-stages. We don't really care about substages, so you need to sum them (in the same way that the model does to obtain pStage
      pStage[iMCMC,,iTree,iStage] =
      ## So first construct pStage (the proportion of the population in each stage)
      ## Then plot it... but also think about how you will handle multiple lines of MCMC output (& will you use lines with transparency or polygons for showing CIs?)
    }
  }
}

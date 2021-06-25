##################################################################################
## This script can be called from fitModel.R to generate output files from aptC ##
##################################################################################


#############################
## Extract log-likelihoods ##
#############################
logliks = tail(aptC$logProbs, floor(nIter/thin))
logliks = logliks[!is.na(logliks)] ## It's unlikely, but just in case the first row gives NA
logliks = coda::as.mcmc(logliks)
#####################################
## Extract samples and save tofile ##
#####################################
samples = as.matrix(aptC$mvSamples)
samples = tail(samples, floor(nIter/thin))
samples = coda::as.mcmc(samples)

##########################
## Write output to file ##
##########################
(fileName = paste0(
  ## "MCMC/model",
  "APT/model",
  SDmodel, "_", qsubID, "_",
  (date() %>% strsplit(" "))[[1]][c(2,3)] %>% paste0(collapse="-") %>% paste0("_"),
  (date() %>% strsplit(" "))[[1]][c(4)] %>% stringr::str_replace_all(":","-"), "_",
  (date() %>% strsplit(" "))[[1]][5] %>% substr(1,5) %>% stringr::str_replace(":",""),
  "_Temps", nTemps, ".txt"))
write.table(as.matrix(samples), file=fileName, row.names = FALSE)
write.table(as.matrix(logliks), row.names = FALSE, file=sub("\\.","_loglik.",fileName))
if (TRUE) {
  ## Cross correlation plot
  try({
    pdf(file = sub("txt","pdf", sub("\\.","_crosscor.",fileName)), width=20, height=20)
    crosscorr.plot(samples)
    dev.off()
  })
  ## Trajectories plot
  try({
    pdf(file = sub("txt","pdf", sub("\\.","_trajectories.",fileName)))
    plot(samples)
    dev.off()
  })
  ## Log posterior likelihood trajectories
  try({
    pdf(file = sub("txt","pdf", sub("\\.","_trajecory-logliks.",fileName)))
    plot(logliks)
    dev.off()
  })
  ## Temperature trajectories
  try({
    pdf(file = sub("txt","pdf", sub("\\.","_APTtemp-trajecories.",fileName)))
    plotTempTraj(aptC)
    dev.off()
  })
}
## class(samples[-c(1:13000),])
## crosscorr.plot(as.mcmc(samples[-c(1:13000),]))

#############################################
## Widely Applicable Information Criterion ##
#############################################
waic = aptC$calculateWAIC()
write.table(waic, file=sub("\\.","_waic.",fileName), row.names = FALSE, col.names="waic")

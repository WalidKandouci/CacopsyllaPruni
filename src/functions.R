## A four stage model
lmean   = log(c(0.06, 0.04, 0.03, 0.02))
lsd     = c(0.6, 0.4, 0.8, 0.7)
surv    = c(0.97, 0.98, 0.99, 0.95)
res     = rep(10, 4)
devFunc = rep(1, 4)
paraMat = cbind(lmean, lsd, surv, res, devFunc)
paraMat

M1      = getM(paras=paraMat[1,1:3], res=paraMat[1,4], devFunction = 1)
M2      = getM(paras=paraMat[2,1:3], res=paraMat[2,4], devFunction = 1)
M3      = getM(paras=paraMat[3,1:3], res=paraMat[3,4], devFunction = 1)
Mlist = list(M1, M2, M3)
nstage = length(Mlist)

par(mfrow=n2mfrow(nstage))
for (ii in 1:nstage) {
  plot(seq(0-1/(2*res[ii]), 1+1/(2*res[ii]), l=res[ii]+1), Mlist[[ii]][,1], typ="s", col="darkgreen",
       xlab="Proportion of stage completed", main=paste("Stage", ii), ylab="Proportion of population")
}

M <- setMultiM(paras=paraMat[,1:3], res=paraMat[,4], devFunction=paraMat[,5], femfec=0)
options(width=10000)
options(max.print=999999)
M
V <- M[,1] * 0
V[1] <- 1
M %*% V

sparseTW1step = nimbleFunction(
  run = function(vecOld=double(1),
                 devMat=double(2)
                 ){
    nstage = nimDim(devMat)[1]
    res = nimDim(devMat)[2]-1
    vecLength = res*nstage+1
    vecNew = nimNumeric(length=vecLength)
    for (stage in 1:nstage){
      iWB = nimNumeric(length = res+1) # iWB: index within-between stage
      iWB[1:res] = (stage -1) * res + 1:res
      iWB[res+1] = max(iWB[1:res]) + 1
      for (ii in 1:res) {
        vecNew[iWB[ii:res]] = vecNew[iWB[ii:res]] + vecOld[iWB[ii]] * devMat[stage,1:(res-ii+1)]
        vecNew[iWB[res+1]]  = vecNew[iWB[res+1]]  + vecOld[iWB[ii]] * sum(devMat[stage,(res+1):(res-ii+2)])
      }
    }
    vecNew[vecLength] = vecNew[vecLength] + vecOld[vecLength]
    returnType(double(1))
    return(vecNew)
  }
)

############################
## Plot the Briere curves ##
############################

if (FALSE) { # TRUE
  par(mfrow=n2mfrow(Const$nStagesDev))
  for (stage in 1:Const$nStagesDev){
    ttemp = -20:60
    quant = apply(paras[stage,iTemp,2],2, "quantile",prob = c(0.025, 0.5, 0.975))
    plot(ttemp,stBriere(T=-20:60, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]),
      type = 'n', xlab = "time",ylab = "proportion")
    polygon(
      x = c(ttemp, rev(ttemp)),
      y = c(quant[1, ], rev(quant[3, ])),
      col = adjustcolor("red", alpha.f = 0.25),
      border = NA
    )
    lines(ttemp, ttemp,stBriere(T=-20:60, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]), 
          col="red",lwd=1.5)
  }
}
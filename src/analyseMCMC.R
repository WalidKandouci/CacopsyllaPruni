################################
## Plot the the Briere curves ##
################################

if (FALSE) { # TRUE
  par(mfrow=n2mfrow(Const$nStagesDev))
  for (stage in 1:Const$nStagesDev) {
    quant = apply(as.array(stBriere(T=-20:60, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage])),
                  1, "quantile",prob = c(0.025, 0.5, 0.975))
    plot(
      (-20:60),
      stBriere(T=-20:60, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]),
      type = 'l',
      ylim = c(0, 1),
      xlab = "time",
      ylab = "proportion",
      main = paste("Tree", psyllids[[iTree]][1, 2])
    )
    points((-20:60),quant[1,],col='red',type = 'l')
    points((-20:60),quant[3,],col='red',type = 'l')
    #polygon(
    #  x = c(-20:60, rev(-20:60)),
    #  y = c(quant[1, ], rev(quant[3, ])),
    #  col = adjustcolor("red", alpha.f = 0.50),
    #  border = NA
    #)
    #curve(stBriere(T=x, Tmin=cPsyllid$Tmin[stage], Tmax=cPsyllid$Tmax[stage], amplitude=cPsyllid$amplitudeMean[stage], shape=cPsyllid$shapeMean[stage]), -20, 60, ylab="E[Dev rate]", main=paste("stage",stage))
  }
}
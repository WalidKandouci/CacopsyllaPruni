#setwd( "C:/Users/Walid/Documents/STAGE/R/2005 File/Plots and Models")
################################################################################
# Import Data:
Data_Oeuf_L1 <- readRDS(file = "Data_Oeuf_L1.rds")
Data_Oeuf_L1[33,c(3,4)] <- rep(0.04042830,2)
Data_Oeuf_L1[34,c(3,4)] <- rep(0.07498052,2)

Data_L1_L2 <- readRDS(file = "Data_L1_L2.rds")
Data_L1_L2[29,c(3,4)] <- rep(0.07276033,2)
Data_L1_L2[30,c(3,4)] <- rep(0.14019129,2)

Data_L2_L3 <- readRDS(file = "Data_L2_L3.rds")
Data_L2_L3[21,] <- c(10,0.07838828,0.07838828,0.07838828)
Data_L2_L3[22,] <- c(17,0.23107143,0.23107143,0.23107143)

Data_L3_L4 <- readRDS(file = "Data_L3_L4.rds")
Data_L3_L4[19,c(3,4)] <- rep(0.2284921,2)

Data_L4_L5 <- readRDS(file = "Data_L4_L5.rds")
Data_L4_L5[16,c(3,4)] <- rep(0.1669841,2)

Data_L5_imago <- readRDS(file = "Data_L5_imago.rds")
Data_L5_imago[15,c(3,4)] <- rep(0.17523810,2)
################################################################################
# Les modeles ##################################################################
################################################################################
par(mfrow=c(3,2))

# Modele Briere
Modele_Briere <- function(tt, PP){
  Tmin <- PP[1]
  Tmax <- PP[2]
  aa <- PP[3]
  bb <- PP[4]
  if (tt<=Tmin || tt>=Tmax) {
    return(0)
    }
  else{
    return( aa * (10^-5) * (tt-Tmin)^2 * ( (Tmax-tt) ^ (1/bb) ) )
    }
  }

R2_Briere <- function(PP,data){
  xobs <- data[,1]
  yobs <- data[,2]
  yth <- sapply(xobs,Modele_Briere,PP)
  return ( sum((yobs-yth)^2) )
}
################################################################################
# Modele Kontodimas
Modele_Kontodimas <- function ( tt, PP ){
  Tmin <- PP[1]
  Tmax <- PP[2]
  aa <- PP[3]
  if (tt<=Tmin || tt>= Tmax) {
    return(0)
  }
  else{
    return(aa * (tt-Tmin) ^ 2 * (Tmax-tt))
  }
}

R2_Kontodimas <- function(PP,data){
  xobs <- data[,1]
  yobs <- data[,2]
  yth <- sapply(xobs,Modele_Kontodimas,PP)
  return ( sum((yobs-yth)^2) )
}
################################################################################
# Modele Taylor

Modele_Taylor <- function(tt, PP){
  Rm <- PP[1]
  Tm <- PP[2]
  To <- PP[3]
  
  if (tt<=Tmin || tt>= Tmax){
    return(0)
  }
  else{
    return(Rm * exp(-1/2 * ((tt - Tm)/To)^2))
  }
}


R2_Taylor <- function(PP,data){
  xobs <- data[,1]
  yobs <- data[,2]
  yth <- sapply(xobs,Modele_Taylor,PP)
  return ( sum((yobs-yth)^2) )
}

################################################################################
# Plot: Oeuf - L1 ##############################################################
################################################################################

# Briere

################################################################################
{
#plot(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19,
#     xlim=c(0,25),ylim=c(0,0.15),
#     xlab = "temperature",ylab = "devRate",main="Stade oeuf-L1: Modéle Briére")
#points(Data_Oeuf_L1[33:34,c(1,2)], pch=19, col="red")
#
tt <- seq(from=-40, to=40, by=1)
Tmin <- -3.471607
Tmax <- 23.910450
aa <- 18.5
#aa_lower <-16.5
#aa_uper <-20.5
bb <- 8.614587
#lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa_mean,bb))}), 
#     col="red",lwd=2,lty=1)
#lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa_lower,bb))}), 
#       col="red",lwd=2,lty=2)
#lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa_uper,bb))}), 
#       col="red",lwd=2,lty=2)

nlmo_Briere_mean <-  optim ( data=Data_Oeuf_L1[,c(1,2)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
nlmo_Briere_lower <- optim ( data=Data_Oeuf_L1[,c(1,3)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
nlmo_Briere_uper <-  optim ( data=Data_Oeuf_L1[,c(1,4)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )

Y0_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_mean$par[1], nlmo_Briere_mean$par[2], nlmo_Briere_mean$par[3], nlmo_Briere_mean$par[4]))})
Y1_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_lower$par[1],nlmo_Briere_lower$par[2], nlmo_Briere_lower$par[3], nlmo_Briere_lower$par[4]))})
Y2_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_uper$par[1], nlmo_Briere_uper$par[2], nlmo_Briere_uper$par[3], nlmo_Briere_uper$par[4]))})

#polygon(c(tt,rev(tt)),c(Y1_Briere,rev(Y2_Briere)),col ="#FFCDD2")
#lines ( tt, Y1_Briere,col ="#FFCDD2")
#lines ( tt, Y2_Briere,col ="#FFCDD2")
#lines ( tt, Y0_Briere,col="red",lwd=2)
#points(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19)
#points(Data_Oeuf_L1[33:34,c(1,2)], pch=19, col="red")
}
################################################################################

# Kontodimas

################################################################################
{
#plot(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19,xlim=c(0,25),ylim=c(0,0.15),
#     xlab = "temperature",ylab = "devRate",main="Stade oeuf-L1: Modéle de Kontodimas")
#points(Data_Oeuf_L1[33:34,c(1,2)], pch=19, col="red")
#
tt <- seq(from=-40, to=40, by=1)
Tmin <- -2.350632e-01
Tmax <- 4.179383e+01
aa <- 1.212191e-05


nlmo_Kontodimas_mean <- optim (data=Data_Oeuf_L1[,c(1,2)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
nlmo_Kontodimas_lower <- optim (data=Data_Oeuf_L1[,c(1,3)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
nlmo_Kontodimas_uper <- optim (data=Data_Oeuf_L1[,c(1,4)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)

Y0_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_mean$par[1], nlmo_Kontodimas_mean$par[2], nlmo_Kontodimas_mean$par[3]))})
Y1_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_lower$par[1], nlmo_Kontodimas_lower$par[2], nlmo_Kontodimas_lower$par[3]))})
Y2_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_uper$par[1], nlmo_Kontodimas_uper$par[2], nlmo_Kontodimas_uper$par[3]))})

#polygon(c(tt,rev(tt)),c(Y1_Kontodimas,rev(Y2_Kontodimas)),col ="#FFCDD2")
#lines ( tt, Y1_Kontodimas,col ="#FFCDD2",lwd=2)
#lines ( tt, Y2_Kontodimas,col ="#FFCDD2",lwd=2)
#lines ( tt, Y0_Kontodimas,col="red",lwd=2)
#points(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19)
#points(Data_Oeuf_L1[33:34,c(1,2)], pch=19, col="red")
}
################################################################################

# Taylor

################################################################################
{
#plot(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19,xlim=c(0,25),ylim=c(0,0.15),
#     xlab = "temperature",ylab = "devRate",main="Stade oeuf-L1: Modéle de Taylor")
#  
#points(Data_Oeuf_L1[33:34,c(1,2)], pch=19, col="red")
  
tt <- seq(from=-40, to=40, by=1)
Rm <- 0.1271879 
Tm <- 26.0212581 
To <- 10.7338761

nlmo_Taylor_mean <- optim (data=Data_Oeuf_L1[,c(1,2)], p=c(Rm,Tm,To), f=R2_Taylor)
nlmo_Taylor_lower <- optim (data=Data_Oeuf_L1[,c(1,3)], p=c(Rm,Tm,To), f=R2_Taylor)
nlmo_Taylor_uper <- optim (data=Data_Oeuf_L1[,c(1,4)], p=c(Rm,Tm,Tm), f=R2_Taylor)
  
Y0_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_mean$par[1], nlmo_Taylor_mean$par[2], nlmo_Taylor_mean$par[3]))})
Y1_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_lower$par[1], nlmo_Taylor_lower$par[2], nlmo_Taylor_lower$par[3]))})
Y2_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_uper$par[1], nlmo_Taylor_uper$par[2], nlmo_Taylor_uper$par[3]))})
  
#polygon(c(tt,rev(tt)),c(Y1_Taylor,rev(Y2_Taylor)),col ="#FFCDD2")
#lines ( tt, Y1_Taylor,col ="#FFCDD2",lwd=2)
##lines ( tt, Y2_Taylor,col ="#FFCDD2",lwd=2)
##lines ( tt, Y0_Taylor,col="red",lwd=2)
#points(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19)
#points(Data_Oeuf_L1[33:34,c(1,2)], pch=19, col="red")
}
################################################################################

# Model VS Model VS Model

{
  plot(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19,xlim=c(0,30),ylim=c(0,0.13),type='n',
       xlab = "temperature",ylab = "devRate",main="Stade oeuf-L1")
  lines ( tt, Y0_Kontodimas,col="navy",lwd=3)
  lines ( tt, Y0_Taylor,col="springgreen",lwd=3)
  lines ( tt, Y0_Briere,col="red",lwd=3)
  points(Data_Oeuf_L1$temp_moyenne,Data_Oeuf_L1$devRate_Moyen,pch=19,lwd=0.5)
  points(Data_Oeuf_L1[33:34,c(1,2)], pch=19, col="goldenrod",lwd=0.5)
legend("topleft", 
       c("Briere","Kontodimas","Taylor"),
       col=c("red","navy","springgreen"),
       lty = c(1,1,1),
       lwd = 3,
       cex=0.5,
       box.lty=1
)
}
################################################################################
# Plot: L1 - Oeuf ##############################################################
################################################################################

# Briere

################################################################################
#par(mfrow=c(3,1))
{
#  plot(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19,
#       xlim=c(5,23),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L1-L2: Modéle Briére")
#  points(Data_L1_L2[29:30,c(1,2)], pch=19, col="red")
#  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- -9.510325
  Tmax <- 22.845427
  aa <- 17.303670
  bb <- 13.912289
  #lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa,bb))}), 
  #    col="red",lwd=2,lty=1)
  
  nlmo_Briere_mean <-  optim ( data=Data_L1_L2[,c(1,2)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_lower <- optim ( data=Data_L1_L2[,c(1,3)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_uper <-  optim ( data=Data_L1_L2[,c(1,4)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  
  Y0_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_mean$par[1], nlmo_Briere_mean$par[2], nlmo_Briere_mean$par[3], nlmo_Briere_mean$par[4]))})
  Y1_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_lower$par[1],nlmo_Briere_lower$par[2], nlmo_Briere_lower$par[3], nlmo_Briere_lower$par[4]))})
  Y2_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_uper$par[1], nlmo_Briere_uper$par[2], nlmo_Briere_uper$par[3], nlmo_Briere_uper$par[4]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Briere,rev(Y2_Briere)),col ="#FFCDD2")
#  lines ( tt, Y1_Briere,col ="#FFCDD2")
#  lines ( tt, Y2_Briere,col ="#FFCDD2")
#  lines ( tt, Y0_Briere,col="red",lwd=2)
#  points(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19)
#  points(Data_L1_L2[29:30,c(1,2)], pch=19, col="red")
}
################################################################################

# Kontodimas

################################################################################
{
#  plot(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19,xlim=c(5,23),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L1-L2: Modéle de Kontodimas")
#  points(Data_L1_L2[29:30,c(1,2)], pch=19, col="red")
#  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- -4.277287e+00
  Tmax <- 4.389937e+01
  aa <- 1.131127e-05
  
#  lines( tt,sapply(tt,function(x){Modele_Kontodimas(x,c(Tmin,Tmax,aa))}), 
#      col="red",lwd=2,lty=1)
  
  nlmo_Kontodimas_mean <- optim (data=Data_L1_L2[,c(1,2)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_lower <- optim (data=Data_L1_L2[,c(1,3)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_uper <- optim (data=Data_L1_L2[,c(1,4)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  
  Y0_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_mean$par[1], nlmo_Kontodimas_mean$par[2], nlmo_Kontodimas_mean$par[3]))})
  Y1_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_lower$par[1], nlmo_Kontodimas_lower$par[2], nlmo_Kontodimas_lower$par[3]))})
  Y2_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_uper$par[1], nlmo_Kontodimas_uper$par[2], nlmo_Kontodimas_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Kontodimas,rev(Y2_Kontodimas)),col ="#FFCDD2")
#  lines ( tt, Y1_Kontodimas,col ="#FFCDD2",lwd=2)
#  lines ( tt, Y2_Kontodimas,col ="#FFCDD2",lwd=2)
#  lines ( tt, Y0_Kontodimas,col="red",lwd=2)
#  points(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19)
#  points(Data_L1_L2[29:30,c(1,2)], pch=19, col="red")
}
################################################################################

# Taylor

################################################################################
{
#  plot(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19,xlim=c(5,23),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L1-L2: Modéle de Taylor")
#  
#  points(Data_L1_L2[29:30,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Rm <- 0.182333
  Tm <- 26.609012 
  To <- 12.841648
  
#  lines( tt,sapply(tt,function(x){Modele_Taylor(x,c(Rm,Tm,To))}), 
#        col="red",lwd=2,lty=1)
#  
  nlmo_Taylor_mean <- optim (data=Data_L1_L2[,c(1,2)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_lower <- optim (data=Data_L1_L2[,c(1,3)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_uper <- optim (data=Data_L1_L2[,c(1,4)], p=c(Rm,Tm,Tm), f=R2_Taylor)
  
  Y0_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_mean$par[1], nlmo_Taylor_mean$par[2], nlmo_Taylor_mean$par[3]))})
  Y1_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_lower$par[1], nlmo_Taylor_lower$par[2], nlmo_Taylor_lower$par[3]))})
  Y2_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_uper$par[1], nlmo_Taylor_uper$par[2], nlmo_Taylor_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Taylor,rev(Y2_Taylor)),col ="#FFCDD2")
#  lines ( tt, Y1_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y2_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Taylor,col="red",lwd=2)
#  points(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19)
#  points(Data_L1_L2[29:30,c(1,2)], pch=19, col="red")
}
################################################################################

# Model VS Model VS Model

################################################################################
{
  plot(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19,xlim=c(-10,32),ylim=c(0,0.21),type="n",
       xlab = "temperature",ylab = "devRate",main="Stade L1-L2")
  lines ( tt, Y0_Kontodimas,col="navy",lwd=3)
  lines ( tt, Y0_Taylor,col="springgreen",lwd=3)
  lines ( tt, Y0_Briere,col="red",lwd=3)
  points(Data_L1_L2$temp_moyenne,Data_L1_L2$devRate_Moyen,pch=19)
  points(Data_L1_L2[29:30,c(1,2)], pch=19, col="goldenrod")
legend("topleft", 
       c("Briere","Kontodimas","Taylor"),
       col=c("red","navy","springgreen"),
       lty = c(1,1,1),
       lwd = 3,
       cex=0.5,
       box.lty=1)
}
################################################################################
# Plot: L2 - L3 ################################################################
################################################################################

# Briere

################################################################################
#par(mfrow=c(3,1))
{
#  plot(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19,
#       xlim=c(5,28),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L2-L3: Modéle Briére")
#  points(Data_L2_L3[21:22,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- -4.888451
  Tmax <- 28.378576
  aa <- 12.652170
  bb <- 2.898462
  #lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa,bb))}), 
      #col="red",lwd=2,lty=1)
  
  nlmo_Briere_mean <-  optim ( data=Data_L2_L3[,c(1,2)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_lower <- optim ( data=Data_L2_L3[,c(1,3)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_uper <-  optim ( data=Data_L2_L3[,c(1,4)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  
  Y0_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_mean$par[1], nlmo_Briere_mean$par[2], nlmo_Briere_mean$par[3], nlmo_Briere_mean$par[4]))})
  Y1_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_lower$par[1],nlmo_Briere_lower$par[2], nlmo_Briere_lower$par[3], nlmo_Briere_lower$par[4]))})
  Y2_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_uper$par[1], nlmo_Briere_uper$par[2], nlmo_Briere_uper$par[3], nlmo_Briere_uper$par[4]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Briere,rev(Y2_Briere)),col ="#FFCDD2")
#  lines ( tt, Y1_Briere,col ="#FFCDD2")
##  lines ( tt, Y2_Briere,col ="#FFCDD2")
##  lines ( tt, Y0_Briere,col="red",lwd=2)
#  points(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19)
#  points(Data_L2_L3[21:22,c(1,2)], pch=19, col="red")
}
################################################################################

# Kontodimas

################################################################################
{
#  plot(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19,xlim=c(5,28),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L2-L3: Modéle de Kontodimas")
#  points(Data_L2_L3[21:22,c(1,2)], pch=19, col="red")
#  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <--1.529734e+00
  Tmax <- 3.707186e+01
  aa <-2.048031e-05
  
    #lines( tt,sapply(tt,function(x){Modele_Kontodimas(x,c(Tmin,Tmax,aa))}), 
     #   col="red",lwd=2,lty=1)
  
  nlmo_Kontodimas_mean <- optim (data=Data_L2_L3[,c(1,2)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_lower <- optim (data=Data_L2_L3[,c(1,3)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_uper <- optim (data=Data_L2_L3[,c(1,4)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  
  Y0_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_mean$par[1], nlmo_Kontodimas_mean$par[2], nlmo_Kontodimas_mean$par[3]))})
  Y1_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_lower$par[1], nlmo_Kontodimas_lower$par[2], nlmo_Kontodimas_lower$par[3]))})
  Y2_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_uper$par[1], nlmo_Kontodimas_uper$par[2], nlmo_Kontodimas_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Kontodimas,rev(Y2_Kontodimas)),col ="#FFCDD2")
#  lines ( tt, Y1_Kontodimas,col ="#FFCDD2",lwd=2)
#  lines ( tt, Y2_Kontodimas,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Kontodimas,col="red",lwd=2)
#  points(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19)
#  points(Data_L2_L3[21:22,c(1,2)], pch=19, col="red")
}
################################################################################

# Taylor

################################################################################
{
#  plot(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19,xlim=c(5,28),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L2-L3: Modéle de Taylor")
#  
#  points(Data_L2_L3[21:22,c(1,2)], pch=19, col="red")
#  
  tt <- seq(from=-40, to=40, by=1)
  Rm <- 0.1740779
  Tm <- 23.9774241 
  To <- 10.7637230
  
  #lines( tt,sapply(tt,function(x){Modele_Taylor(x,c(Rm,Tm,To))}), 
   #      col="red",lwd=2,lty=1)
  
  nlmo_Taylor_mean <- optim (data=Data_L2_L3[,c(1,2)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_lower <- optim (data=Data_L2_L3[,c(1,3)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_uper <- optim (data=Data_L2_L3[,c(1,4)], p=c(Rm,Tm,Tm), f=R2_Taylor)
  
  Y0_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_mean$par[1], nlmo_Taylor_mean$par[2], nlmo_Taylor_mean$par[3]))})
  Y1_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_lower$par[1], nlmo_Taylor_lower$par[2], nlmo_Taylor_lower$par[3]))})
  Y2_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_uper$par[1], nlmo_Taylor_uper$par[2], nlmo_Taylor_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Taylor,rev(Y2_Taylor)),col ="#FFCDD2")
##  lines ( tt, Y1_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y2_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Taylor,col="red",lwd=2)
##  points(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19)
#  points(Data_L2_L3[21:22,c(1,2)], pch=19, col="red")
}
################################################################################

# Model VS Model VS Model

################################################################################
{
  plot(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19,xlim=c(-3,35),ylim=c(0,0.23),type="n",
       xlab = "temperature",ylab = "devRate",main="Stade L2-L3")
  lines ( tt, Y0_Kontodimas,col="navy",lwd=3)
  lines ( tt, Y0_Taylor,col="springgreen",lwd=3)
  lines ( tt, Y0_Briere,col="red",lwd=3)
  points(Data_L2_L3$temp_moyenne,Data_L2_L3$devRate_Moyen,pch=19)
  points(Data_L2_L3[21:22,c(1,2)], pch=19, col="goldenrod")
  legend("topleft", 
         c("Briere","Kontodimas","Taylor"),
         col=c("red","navy","springgreen"),
         lty = c(1,1,1),
         lwd = 3,
         cex=0.5,
         box.lty=1)
}###########################################################################
# Plot: L3 - L4 ################################################################
################################################################################

# Briere

################################################################################
#par(mfrow=c(3,1))
{
#  plot(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19,
#       xlim=c(13,25),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L3-L4: Modéle Briére")
#  points(Data_L3_L4[19,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- -0.5914309
  Tmax <- 30.0150014  
  aa <- 7.5764591
  bb <- 1.3427041
  #lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa,bb))}), 
  #col="red",lwd=2,lty=1)
  
  nlmo_Briere_mean <-  optim ( data=Data_L3_L4[,c(1,2)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_lower <- optim ( data=Data_L3_L4[,c(1,3)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_uper <-  optim ( data=Data_L3_L4[,c(1,4)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  
  Y0_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_mean$par[1], nlmo_Briere_mean$par[2], nlmo_Briere_mean$par[3], nlmo_Briere_mean$par[4]))})
  Y1_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_lower$par[1],nlmo_Briere_lower$par[2], nlmo_Briere_lower$par[3], nlmo_Briere_lower$par[4]))})
  Y2_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_uper$par[1], nlmo_Briere_uper$par[2], nlmo_Briere_uper$par[3], nlmo_Briere_uper$par[4]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Briere,rev(Y2_Briere)),col ="#FFCDD2")
#  lines ( tt, Y1_Briere,col ="#FFCDD2")
##  lines ( tt, Y2_Briere,col ="#FFCDD2")
##  lines ( tt, Y0_Briere,col="red",lwd=2)
#  points(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19)
#  points(Data_L3_L4[19,c(1,2)], pch=19, col="red")
}
################################################################################

# Kontodimas

################################################################################
{
#  plot(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19,
#       xlim=c(13,25),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L3-L4: Modéle de Kontodimas")
#  points(Data_L3_L4[19,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- 2.448643e+00
  Tmax <- 3.097268e+01
  aa <- 5.360027e-05
  
  #lines( tt,sapply(tt,function(x){Modele_Kontodimas(x,c(Tmin,Tmax,aa))}), 
  #   col="red",lwd=2,lty=1)
  
  nlmo_Kontodimas_mean <- optim (data=Data_L3_L4[,c(1,2)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_lower <- optim (data=Data_L3_L4[,c(1,3)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_uper <- optim (data=Data_L3_L4[,c(1,4)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  
  Y0_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_mean$par[1], nlmo_Kontodimas_mean$par[2], nlmo_Kontodimas_mean$par[3]))})
  Y1_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_lower$par[1], nlmo_Kontodimas_lower$par[2], nlmo_Kontodimas_lower$par[3]))})
  Y2_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_uper$par[1], nlmo_Kontodimas_uper$par[2], nlmo_Kontodimas_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Kontodimas,rev(Y2_Kontodimas)),col ="#FFCDD2")
#  lines ( tt, Y1_Kontodimas,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y2_Kontodimas,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Kontodimas,col="red",lwd=2)
#  points(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19)
#  points(Data_L3_L4[19,c(1,2)], pch=19, col="red")
}
################################################################################

# Taylor

################################################################################
{
#  plot(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19,xlim=c(5,28),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L3-L4: Modéle de Taylor")
#  
#  points(Data_L3_L4[19,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Rm <-  0.1830892
  Tm <- 21.1899669
  To <- 7.8732500
  
  #lines( tt,sapply(tt,function(x){Modele_Taylor(x,c(Rm,Tm,To))}), 
  #      col="red",lwd=2,lty=1)
  
  nlmo_Taylor_mean <- optim (data=Data_L3_L4[,c(1,2)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_lower <- optim (data=Data_L3_L4[,c(1,3)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_uper <- optim (data=Data_L3_L4[,c(1,4)], p=c(Rm,Tm,Tm), f=R2_Taylor)
  
  Y0_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_mean$par[1], nlmo_Taylor_mean$par[2], nlmo_Taylor_mean$par[3]))})
  Y1_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_lower$par[1], nlmo_Taylor_lower$par[2], nlmo_Taylor_lower$par[3]))})
  Y2_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_uper$par[1], nlmo_Taylor_uper$par[2], nlmo_Taylor_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Taylor,rev(Y2_Taylor)),col ="#FFCDD2")
#  lines ( tt, Y1_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y2_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Taylor,col="red",lwd=2)
#  points(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19)
#  points(Data_L3_L4[19,c(1,2)], pch=19, col="red")
}
################################################################################

# Model VS Model VS Model

################################################################################

{
  plot(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19,xlim=c(1,32),ylim=c(0,0.23),type="n",
       xlab = "temperature",ylab = "devRate",main="Stade L3-L4")
  lines ( tt, Y0_Kontodimas,col="navy",lwd=3)
  lines ( tt, Y0_Taylor,col="springgreen",lwd=3)
  lines ( tt, Y0_Briere,col="red",lwd=3)
  points(Data_L3_L4$temp_moyenne,Data_L3_L4$devRate_Moyen,pch=19)
  points(Data_L3_L4[19,c(1,2)], pch=19, col="goldenrod")
  legend("topleft", 
         c("Briere","Kontodimas","Taylor"),
         col=c("red","navy","springgreen"),
         lty = c(1,1,1),
         lwd = 3,
         cex=0.5,
         box.lty=1)
}
################################################################################
# Plot: L4 - L5 ################################################################
################################################################################

# Briere

################################################################################
#par(mfrow=c(3,1))
{
#  plot(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19,
#       xlim=c(13,28),ylim=c(0,0.25),
#       xlab = "temperature",ylab = "devRate",main="Stade L4-L5: Modéle Briére")
#  points(Data_L4_L5[16,c(1,2)], pch=19, col="red")
#  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- -6.054626
  Tmax <- 27.065688 
  aa <- 15.680779
  bb <- 4.138661
  #lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa,bb))}), 
  #col="red",lwd=2,lty=1)
  
  nlmo_Briere_mean <-  optim ( data=Data_L4_L5[,c(1,2)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_lower <- optim ( data=Data_L4_L5[,c(1,3)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_uper <-  optim ( data=Data_L4_L5[,c(1,4)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  
  Y0_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_mean$par[1], nlmo_Briere_mean$par[2], nlmo_Briere_mean$par[3], nlmo_Briere_mean$par[4]))})
  Y1_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_lower$par[1],nlmo_Briere_lower$par[2], nlmo_Briere_lower$par[3], nlmo_Briere_lower$par[4]))})
  Y2_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_uper$par[1], nlmo_Briere_uper$par[2], nlmo_Briere_uper$par[3], nlmo_Briere_uper$par[4]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Briere,rev(Y2_Briere)),col ="#FFCDD2")
##  lines ( tt, Y1_Briere,col ="#FFCDD2")
##  lines ( tt, Y2_Briere,col ="#FFCDD2")
##  lines ( tt, Y0_Briere,col="red",lwd=2)
##  points(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19)
#  points(Data_L4_L5[16,c(1,2)], pch=19, col="red")
}
################################################################################

# Kontodimas

################################################################################
{
#  plot(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19,
#       xlim=c(13,28),ylim=c(0,0.25),
#       xlab = "temperature",ylab = "devRate",main="Stade L4-L5: Modéle de Kontodimas")
#  points(Data_L4_L5[16,c(1,2)], pch=19, col="red")
#  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- 2.665918e+00
  Tmax <- 3.311551e+01
  aa <- 4.386207e-05
  
  #lines( tt,sapply(tt,function(x){Modele_Kontodimas(x,c(Tmin,Tmax,aa))}), 
  #   col="red",lwd=2,lty=1)
  
  nlmo_Kontodimas_mean <- optim (data=Data_L4_L5[,c(1,2)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_lower <- optim (data=Data_L4_L5[,c(1,3)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_uper <- optim (data=Data_L4_L5[,c(1,4)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  
  Y0_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_mean$par[1], nlmo_Kontodimas_mean$par[2], nlmo_Kontodimas_mean$par[3]))})
  Y1_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_lower$par[1], nlmo_Kontodimas_lower$par[2], nlmo_Kontodimas_lower$par[3]))})
  Y2_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_uper$par[1], nlmo_Kontodimas_uper$par[2], nlmo_Kontodimas_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Kontodimas,rev(Y2_Kontodimas)),col ="#FFCDD2")
##  lines ( tt, Y1_Kontodimas,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y2_Kontodimas,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Kontodimas,col="red",lwd=2)
##  points(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19)
#  points(Data_L4_L5[16,c(1,2)], pch=19, col="red")
}
################################################################################

# Taylor

################################################################################
{
#  plot(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19,
#       xlim=c(13,28),ylim=c(0,0.25),
#       xlab = "temperature",ylab = "devRate",main="Stade L4-L5: Modéle de Taylor")
#  
#  points(Data_L4_L5[16,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Rm <- 0.1820486
  Tm <- 23.0149564
  To <- 8.9655981
  
  #lines( tt,sapply(tt,function(x){Modele_Taylor(x,c(Rm,Tm,To))}), 
  #      col="red",lwd=2,lty=1)
  
  nlmo_Taylor_mean <- optim (data=Data_L4_L5[,c(1,2)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_lower <- optim (data=Data_L4_L5[,c(1,3)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_uper <- optim (data=Data_L4_L5[,c(1,4)], p=c(Rm,Tm,Tm), f=R2_Taylor)
  
  Y0_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_mean$par[1], nlmo_Taylor_mean$par[2], nlmo_Taylor_mean$par[3]))})
  Y1_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_lower$par[1], nlmo_Taylor_lower$par[2], nlmo_Taylor_lower$par[3]))})
  Y2_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_uper$par[1], nlmo_Taylor_uper$par[2], nlmo_Taylor_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Taylor,rev(Y2_Taylor)),col ="#FFCDD2")
##  lines ( tt, Y1_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y2_Taylor,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Taylor,col="red",lwd=2)
##  points(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19)
#  points(Data_L4_L5[16,c(1,2)], pch=19, col="red")
}
################################################################################

# Model VS Model VS Model

################################################################################
{
  plot(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19,type = "n",
       xlim=c(1,32),ylim=c(0,0.25),
       xlab = "temperature",ylab = "devRate",main="Stade L4-L5")
  lines ( tt, Y0_Kontodimas,col="navy",lwd=3)
  lines ( tt, Y0_Taylor,col="springgreen",lwd=3)
  lines ( tt, Y0_Briere,col="red",lwd=3)
  points(Data_L4_L5$temp_moyenne,Data_L4_L5$devRate_Moyen,pch=19)
  points(Data_L4_L5[16,c(1,2)], pch=19, col="goldenrod")
  legend("topleft", 
         c("Briere","Kontodimas","Taylor"),
         col=c("red","navy","springgreen"),
         lty = c(1,1,1),
         lwd = 3,
         cex=0.5,
         box.lty=1)
}
################################################################################
# Plot: L5 - imago ################################################################
################################################################################

# Briere

################################################################################
{
#  plot(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19,
#       xlim=c(13,25),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L5-imago: Modéle Briére")
#  points(Data_L5_imago[15,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- -3.144059
  Tmax <- 27.469951 
  aa <- 17.065835
  bb <- 3.170311
  #lines( tt,sapply(tt,function(x){Modele_Briere(x,c(Tmin,Tmax,aa,bb))}), 
  #col="red",lwd=2,lty=1)
  
  nlmo_Briere_mean <-  optim ( data=Data_L5_imago[,c(1,2)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_lower <- optim ( data=Data_L5_imago[,c(1,3)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  nlmo_Briere_uper <-  optim ( data=Data_L5_imago[,c(1,4)], p=c(Tmin,Tmax,aa,bb) ,f=R2_Briere )
  
  Y0_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_mean$par[1], nlmo_Briere_mean$par[2], nlmo_Briere_mean$par[3], nlmo_Briere_mean$par[4]))})
  Y1_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_lower$par[1],nlmo_Briere_lower$par[2], nlmo_Briere_lower$par[3], nlmo_Briere_lower$par[4]))})
  Y2_Briere <- sapply(tt,function(x){Modele_Briere (x,c(nlmo_Briere_uper$par[1], nlmo_Briere_uper$par[2], nlmo_Briere_uper$par[3], nlmo_Briere_uper$par[4]))})
#  
#  polygon(c(tt,rev(tt)),c(Y1_Briere,rev(Y2_Briere)),col ="#FFCDD2")
#  lines ( tt, Y1_Briere,col ="#FFCDD2")
#  lines ( tt, Y2_Briere,col ="#FFCDD2")
#  lines ( tt, Y0_Briere,col="red",lwd=2)
#  points(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19)
#  points(Data_L5_imago[15,c(1,2)], pch=19, col="red")
}
################################################################################

# Kontodimas

################################################################################
{
#  plot(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19,
#       xlim=c(13,25),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L5-imago: Modéle de Kontodimas")
#  points(Data_L5_imago[15,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Tmin <- 2.667722e+00
  Tmax <- 3.311468e+01
  aa <- 4.387384e-05
  
  #lines( tt,sapply(tt,function(x){Modele_Kontodimas(x,c(Tmin,Tmax,aa))}), 
  #   col="red",lwd=2,lty=1)
  
  nlmo_Kontodimas_mean <- optim (data=Data_L5_imago[,c(1,2)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_lower <- optim (data=Data_L5_imago[,c(1,3)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  nlmo_Kontodimas_uper <- optim (data=Data_L5_imago[,c(1,4)], p=c(Tmin,Tmax,aa), f=R2_Kontodimas)
  
  Y0_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_mean$par[1], nlmo_Kontodimas_mean$par[2], nlmo_Kontodimas_mean$par[3]))})
  Y1_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_lower$par[1], nlmo_Kontodimas_lower$par[2], nlmo_Kontodimas_lower$par[3]))})
  Y2_Kontodimas <- sapply(tt,function(x){Modele_Kontodimas (x,c(nlmo_Kontodimas_uper$par[1], nlmo_Kontodimas_uper$par[2], nlmo_Kontodimas_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Kontodimas,rev(Y2_Kontodimas)),col ="#FFCDD2")
#  lines ( tt, Y1_Kontodimas,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y2_Kontodimas,col ="#FFCDD2",lwd=2)
##  lines ( tt, Y0_Kontodimas,col="red",lwd=2)
#  points(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19)
#  points(Data_L5_imago[15,c(1,2)], pch=19, col="red")
}
################################################################################

# Taylor

################################################################################
{
#  plot(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19,xlim=c(10,29),ylim=c(0,0.3),
#       xlab = "temperature",ylab = "devRate",main="Stade L5-imago: Modéle de Taylor")
#  
#  points(Data_L5_imago[15,c(1,2)], pch=19, col="red")
  
  tt <- seq(from=-40, to=40, by=1)
  Rm <-  0.1820486
  Tm <- 23.0149564
  To <- 8.9655981
  
  #lines( tt,sapply(tt,function(x){Modele_Taylor(x,c(Rm,Tm,To))}), 
  #      col="red",lwd=2,lty=1)
  
  nlmo_Taylor_mean <- optim (data=Data_L5_imago[,c(1,2)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_lower <- optim (data=Data_L5_imago[,c(1,3)], p=c(Rm,Tm,To), f=R2_Taylor)
  nlmo_Taylor_uper <- optim (data=Data_L5_imago[,c(1,4)], p=c(Rm,Tm,Tm), f=R2_Taylor)
  
  Y0_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_mean$par[1], nlmo_Taylor_mean$par[2], nlmo_Taylor_mean$par[3]))})
  Y1_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_lower$par[1], nlmo_Taylor_lower$par[2], nlmo_Taylor_lower$par[3]))})
  Y2_Taylor <- sapply(tt,function(x){Modele_Taylor (x,c(nlmo_Taylor_uper$par[1], nlmo_Taylor_uper$par[2], nlmo_Taylor_uper$par[3]))})
  
#  polygon(c(tt,rev(tt)),c(Y1_Taylor,rev(Y2_Taylor)),col ="#FFCDD2")
##  lines ( tt, Y1_Taylor,col ="#FFCDD2",lwd=2)
###  lines ( tt, Y2_Taylor,col ="#FFCDD2",lwd=2)
###  lines ( tt, Y0_Taylor,col="red",lwd=2)
##  points(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19)
#  points(Data_L5_imago[15,c(1,2)], pch=19, col="red")
}
################################################################################

# Model VS Model VS Model

################################################################################

{
  plot(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19,xlim=c(-37,37),ylim=c(0,0.2),type="n",
       xlab = "temperature",ylab = "devRate",main="Stade L5-imago")
  lines ( tt, Y0_Kontodimas,col="navy",lwd=3)
  lines ( tt, Y0_Taylor,col="springgreen",lwd=3)
  lines ( tt, Y0_Briere,col="red",lwd=3)
  points(Data_L5_imago$temp_moyenne,Data_L5_imago$devRate_Moyen,pch=19)
  points(Data_L5_imago[15,c(1,2)], pch=19, col="goldenrod")
  legend("topleft", 
         c("Briere","Kontodimas","Taylor"),
         col=c("red","navy","springgreen"),
         lty = c(1,1,1),
         lwd = 3,
         cex=0.5,
         box.lty=1)
}


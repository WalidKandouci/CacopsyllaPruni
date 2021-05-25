# Time series
library(readxl)
library(forecast)
library(ggplot2)
library(tidyverse)
library(rio)
# 2005 Meteo Data 
DataMeteo2005 <- read_excel("DataMeteo2005.xltx")
DataMeteo2005$temperature <- as.numeric(DataMeteo2005$temperature)
DataMeteo2005$humidity <- as.numeric(DataMeteo2005$humidity)
DataMeteo2005$date <- as.Date(DataMeteo2005$date)
Meteo2005 <- DataMeteo2005[,c(1,2)]
# Option 1: linear type with the function "na.interp"

plot(missingMeteo[1:110],
     type = "l",
     main = "",
     xlab = "jours",
     ylab = "températures",
     lwd=3,
     col="red")
lines(Meteo2005$temperature[1:110],
      type = "l",
      lwd=3)

missingMeteo <- na.interp(Meteo2005$temperature)


length(missingMeteo) == length(Meteo2005$temperature)
which(is.na(Meteo2005$temperature))
missingMeteo[which(is.na(Meteo2005$temperature))]
# Option 2: package "imputeTS"
library(imputeTS)
ggplot_na_distribution(Meteo2005$temperature)
#ggplot_na_gapsize(Meteo2005$temperature)
imp <- na_kalman(Meteo2005$temperature)
ggplot_na_imputations(Meteo2005$temperature, imp)

imp[which(is.na(Meteo2005$temperature))]

# Option 3: auto.arima()
y <- Meteo2005$temperature
fit <- auto.arima(Meteo2005$temperature)
kr <- KalmanRun(Meteo2005$temperature, fit$model)
id.na <- which(is.na(Meteo2005$temperature))
for (i in id.na){
  y[i] <- fit$model$Z %*% kr$states[i,]
}
which(is.na(y))
y[id.na]



#################
# PLOT
TS <- ts(y, start = c(1,1), frequency = 7)
autoplot(TS)
#################
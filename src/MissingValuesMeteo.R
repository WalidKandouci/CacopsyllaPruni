# Packages needed
library(readxl)
library(forecast)
library(ggplot2)
library(tidyverse)
library(imputeTS)
#  Import our 2005 Meteo Data 
DataMeteo2005 <- read_excel("DataMeteo2005.xltx")
DataMeteo2005$temperature <- as.numeric(DataMeteo2005$temperature)
DataMeteo2005$humidity <- as.numeric(DataMeteo2005$humidity)
DataMeteo2005$date <- as.Date(DataMeteo2005$date)
# Option 1: linear type with the function "na.interp"
missingMeteo <- na.interp(DataMeteo2005$temperature)

plot(missingMeteo[1:110],
     type = "l",
     main = "",
     xlab = "jours",
     ylab = "températures",
     lwd=3,
     col="red")
lines(DataMeteo2005$temperature[1:110],
      type = "l",
      lwd=3)

length(missingMeteo) == length(DataMeteo2005$temperature)
which(is.na(DataMeteo2005$temperature))
missingMeteo[which(is.na(DataMeteo2005$temperature))]
# Option 2: package "imputeTS"
ggplot_na_distribution(DataMeteo2005$temperature)
imp <- na_kalman(DataMeteo2005$temperature)
ggplot_na_imputations(DataMeteo2005$temperature, imp)

imp[which(is.na(DataMeteo2005$temperature))]

# Option 3: auto.arima()
y <- DataMeteo2005$temperature
fit <- auto.arima(DataMeteo2005$temperature)
kr <- KalmanRun(DataMeteo2005$temperature, fit$model)
id.na <- which(is.na(DataMeteo2005$temperature))
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
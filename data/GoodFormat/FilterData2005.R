library(readr)
Data_ALLMeteo_Montpellier <- as.data.frame(read_csv2("Data_ALLMeteo_Montpellier.csv"))
Data_ALLMeteo_Montpellier <- Data_ALLMeteo_Montpellier[,-1]
Data_ALLMeteo_Montpellier$temperature <- round(Data_ALLMeteo_Montpellier$temperature)
DataMeteo_Montpellier2005 <- Data_ALLMeteo_Montpellier[Data_ALLMeteo_Montpellier$ID==2005,]
DataMeteo_Montpellier2005 <- DataMeteo_Montpellier2005[,-4]
rownames(DataMeteo_Montpellier2005) <- NULL
DataMeteo_Montpellier2005_Psyll <- DataMeteo_Montpellier2005[697:1552,]

DataMeteo_Montpellier2005_Psyll[,4] <- 1:856
colnames(DataMeteo_Montpellier2005_Psyll) <- c("date","temperature","humidity","index")

write.csv2(DataMeteo_Montpellier2005_Psyll,"DataMeteo2005.csv")
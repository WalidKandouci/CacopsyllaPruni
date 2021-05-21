################################################################################
Data_2005_Comptage_13.3 <- as.data.frame(read_excel("Data_2005_Comptage_13,3.xlsx"))
Data_2005_Comptage_13.3 <- Data_2005_Comptage_13.3[,-1]
Data_2005_Comptage_13.3$f8 <- as.numeric(Data_2005_Comptage_13.3$f8)
Data_2005_Comptage_13.3$f9 <- as.numeric(Data_2005_Comptage_13.3$f9)
Data_2005_Comptage_13.3$f10 <- as.numeric(Data_2005_Comptage_13.3$f10)
Data_2005_Comptage_13.3$f11 <- as.numeric(Data_2005_Comptage_13.3$f11)
Data_2005_Comptage_13.3$f12 <- as.numeric(Data_2005_Comptage_13.3$f12)
Data_2005_Comptage_13.3$f13 <- as.numeric(Data_2005_Comptage_13.3$f13)
Data_2005_Comptage_13.3$f14 <- as.numeric(Data_2005_Comptage_13.3$f14)
Data_2005_Comptage_13.3$f15 <- as.numeric(Data_2005_Comptage_13.3$f15)

D_13.3 <- Data_2005_Comptage_13.3

colnames(D_13.3) <- c("date","leaves","egg","L1","L2","L3","L4","L5","imago")
D_13.3$egg <- as.numeric(D_13.3$egg)
D_13.3$L1 <- as.numeric(D_13.3$L1)
D_13.3$L2 <- as.numeric(D_13.3$L2)
D_13.3$L3 <- as.numeric(D_13.3$L3)
D_13.3$L4 <- as.numeric(D_13.3$L4)
D_13.3$L5 <- as.numeric(D_13.3$L5)
D_13.3$imago <- as.numeric(D_13.3$imago)

for (i in D_13.3$date) {
  D_13.3[D_13.3$date==i,3] <- sum(as.numeric(Data_2005_Comptage_13.3[which(Data_2005_Comptage_13.3$date==i&Data_2005_Comptage_13.3$stage=="Oeuf"),]),na.rm = T)
  D_13.3[D_13.3$date==i,4] <- sum(as.numeric(Data_2005_Comptage_13.3[which(Data_2005_Comptage_13.3$date==i&Data_2005_Comptage_13.3$stage=="L1"),]),na.rm = T)
  D_13.3[D_13.3$date==i,5] <- sum(as.numeric(Data_2005_Comptage_13.3[which(Data_2005_Comptage_13.3$date==i&Data_2005_Comptage_13.3$stage=="L2"),]),na.rm = T)
  D_13.3[D_13.3$date==i,6] <- sum(as.numeric(Data_2005_Comptage_13.3[which(Data_2005_Comptage_13.3$date==i&Data_2005_Comptage_13.3$stage=="L3"),]),na.rm = T)
  D_13.3[D_13.3$date==i,7] <- sum(as.numeric(Data_2005_Comptage_13.3[which(Data_2005_Comptage_13.3$date==i&Data_2005_Comptage_13.3$stage=="L4"),]),na.rm = T)
  D_13.3[D_13.3$date==i,8] <- sum(as.numeric(Data_2005_Comptage_13.3[which(Data_2005_Comptage_13.3$date==i&Data_2005_Comptage_13.3$stage=="L5"),]),na.rm = T)
  D_13.3[D_13.3$date==i,9] <- sum(as.numeric(Data_2005_Comptage_13.3[which(Data_2005_Comptage_13.3$date==i&Data_2005_Comptage_13.3$stage=="imago"),]),na.rm = T)
  D_13.3[D_13.3$date==i,2] <- length(Data_2005_Comptage_13.3[Data_2005_Comptage_13.3$date==i,2:16][!is.na(Data_2005_Comptage_13.3[Data_2005_Comptage_13.3$date==i,2:16])])
}

write.csv2(D_13.3,"D_13.3.csv")
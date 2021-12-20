print("%redu")
rrun_dfvx <- as.data.frame(rrun_dfvx)
lastcol <- (((1) * (combn)) + 1) * steps
y1 <- seq(1, lastcol, steps)

pcreduI <- matrix(0, 16, length(vacnames))
pcreduI35 <- matrix(0, 16, length(vacnames))
pcreduM <- matrix(0, 9, length(vacnames))
pcreduM35 <- matrix(0, 9, length(vacnames))

pcredrow <- c("All", "All0-14years", "All15plus", "HAll", "H014years", "H15plus", "NAll", "N014years", "N15plus", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")
agecols <- c("black", "red", "blue", "purple", "green")

pcreduI[1, ] <- 100 * ((rrun_dfvx$TBIalltot[151] - rrun_dfvx$TBIalltot[(y1[-1] + 150)])) / (rrun_dfvx$TBIalltot[151])
pcreduI[2, ] <- 100 * ((rrun_dfvx$"TBIall0-14"[151] - rrun_dfvx$"TBIall0-14"[(y1[-1] + 150)])) / (rrun_dfvx$"TBIall0-14"[151])
pcreduI[3, ] <- 100 * ((rrun_dfvx$TBIall15p[151] - rrun_dfvx$TBIall15p[(y1[-1] + 150)])) / (rrun_dfvx$TBIall15p[151])
pcreduI[4, ] <- 100 * ((rrun_dfvx$TBIHtot[151] - rrun_dfvx$TBIHtot[(y1[-1] + 150)])) / (rrun_dfvx$TBIHtot[151])
pcreduI[5, ] <- 100 * ((rrun_dfvx$"TBIH0-14"[151] - rrun_dfvx$"TBIH0-14"[(y1[-1] + 150)])) / (rrun_dfvx$"TBIH0-14"[151])
pcreduI[6, ] <- 100 * ((rrun_dfvx$TBIH15p[151] - rrun_dfvx$TBIH15p[(y1[-1] + 150)])) / (rrun_dfvx$TBIH15p[151])
pcreduI[7, ] <- 100 * ((rrun_dfvx$TBItot[151] - rrun_dfvx$TBItot[(y1[-1] + 150)])) / (rrun_dfvx$TBItot[151])
pcreduI[8, ] <- 100 * ((rrun_dfvx$"TBI0-14"[151] - rrun_dfvx$"TBI0-14"[(y1[-1] + 150)])) / (rrun_dfvx$"TBI0-14"[151])
pcreduI[9, ] <- 100 * ((rrun_dfvx$"TBI15+"[151] - rrun_dfvx$"TBI15+"[(y1[-1] + 150)])) / (rrun_dfvx$"TBI15+"[151])
pcreduI[10:16, ] <- t(rrun_dfvx[y1[-1], c("type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")])

pcreduI35[1, ] <- 100 * ((rrun_dfvx$TBIalltot[136] - rrun_dfvx$TBIalltot[(y1[-1] + 135)])) / (rrun_dfvx$TBIalltot[136])
pcreduI35[2, ] <- 100 * ((rrun_dfvx$"TBIall0-14"[136] - rrun_dfvx$"TBIall0-14"[(y1[-1] + 135)])) / (rrun_dfvx$"TBIall0-14"[136])
pcreduI35[3, ] <- 100 * ((rrun_dfvx$TBIall15p[136] - rrun_dfvx$TBIall15p[(y1[-1] + 135)])) / (rrun_dfvx$TBIall15p[136])
pcreduI35[4, ] <- 100 * ((rrun_dfvx$TBIHtot[136] - rrun_dfvx$TBIHtot[(y1[-1] + 135)])) / (rrun_dfvx$TBIHtot[136])
pcreduI35[5, ] <- 100 * ((rrun_dfvx$"TBIH0-14"[136] - rrun_dfvx$"TBIH0-14"[(y1[-1] + 135)])) / (rrun_dfvx$"TBIH0-14"[136])
pcreduI35[6, ] <- 100 * ((rrun_dfvx$TBIH15p[136] - rrun_dfvx$TBIH15p[(y1[-1] + 135)])) / (rrun_dfvx$TBIH15p[136])
pcreduI35[7, ] <- 100 * ((rrun_dfvx$TBItot[136] - rrun_dfvx$TBItot[(y1[-1] + 135)])) / (rrun_dfvx$TBItot[136])
pcreduI35[8, ] <- 100 * ((rrun_dfvx$"TBI0-14"[136] - rrun_dfvx$"TBI0-14"[(y1[-1] + 135)])) / (rrun_dfvx$"TBI0-14"[136])
pcreduI35[9, ] <- 100 * ((rrun_dfvx$"TBI15+"[136] - rrun_dfvx$"TBI15+"[(y1[-1] + 135)])) / (rrun_dfvx$"TBI15+"[136])
pcreduI35[10:16, ] <- t(rrun_dfvx[y1[-1], c("type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")])

colnames(pcreduI) <- vacnames
rownames(pcreduI) <- pcredrow
colnames(pcreduI35) <- vacnames
rownames(pcreduI35) <- pcredrow

pcredrowM <- c("All", "HAll", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")

pcreduM[1, ] <- 100 * ((rrun_dfvx$TBMtot[151] - rrun_dfvx$TBMtot[(y1[-1] + 150)])) / (rrun_dfvx$TBMtot[151])
pcreduM[2, ] <- 100 * ((rrun_dfvx$TBMHdatot[151] - rrun_dfvx$TBMHdatot[(y1[-1] + 150)])) / (rrun_dfvx$TBMHdatot[151])
pcreduM[3:9, ] <- t(rrun_dfvx[y1[-1], c("type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")])

pcreduM35[1, ] <- 100 * ((rrun_dfvx$TBMtot[136] - rrun_dfvx$TBMtot[(y1[-1] + 135)])) / (rrun_dfvx$TBMtot[136])
pcreduM35[2, ] <- 100 * ((rrun_dfvx$TBMHdatot[136] - rrun_dfvx$TBMHdatot[(y1[-1] + 135)])) / (rrun_dfvx$TBMHdatot[136])
pcreduM35[3:9, ] <- t(rrun_dfvx[y1[-1], c("type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")])

colnames(pcreduM) <- vacnames
rownames(pcreduM) <- pcredrowM
colnames(pcreduM35) <- vacnames
rownames(pcreduM35) <- pcredrowM

setwd("Vxoutput_M72_refit")
setwd(vaxfol)
write.table(pcreduI, paste("2050_reduction_incidence_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(pcreduI35, paste("2035_reduction_incidence_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(pcreduM, paste("2050_reduction_mortality_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(pcreduM35, paste("2035_reduction_mortality_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
setwd(home)

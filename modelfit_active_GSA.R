library(plyr)
library(grid)
library(ggplot2)
library(gtable)
library(grid)

C <- 1
if (C == 0) {
  home <- "/Users/lsh355020/R/Gates_sa"
}
if (C == 1) {
  home <- "/home/lsh355020/SA_Gates/"
}

if (C == 0) {
  input <- "/Users/lsh355020/R/Gates_sa/Vxoutput"
}
if (C == 1) {
  input <- "/home/lsh355020/SA_Gates/Vxoutput"
}

if (C == 0) {
  output <- "/Users/lsh355020/R/Gates_sa/Epi"
}
if (C == 1) {
  output <- "/home/lsh355020/SA_Gates/Epi"
}

setwd(home)
setwd("Data")
fitdata <- as.data.frame(read.csv("fitdata_sa.csv", header = TRUE, check.names = F))

setwd(home)
setwd(input)

setwd(home)

rrun <- 20

hivoff <- 0
artoff <- 0
art100 <- 0
artred100 <- 0

source("#DataGrabGSA.R")
setwd(home)
source("CFunctions_GSA_uHfix.R")

setwd(home)
para <- as.matrix(drop.levels(read.csv("fitcheckinterim.csv", header = TRUE, check.names = F)))

nm <- c(pararange[, 1], "p0")

para <- as.data.frame(para)

cntry <- "South Africa"

setwd(home)

yrTBItot <- matrix(0, 51, rrun)
yrTBI_014 <- matrix(0, 51, rrun)
yrTBI_15 <- matrix(0, 51, rrun)

yrTBIHtot <- matrix(0, 51, rrun)

yrTBMtot <- matrix(0, 51, rrun)
yrTBMHtot <- matrix(0, 51, rrun)

yrTBNtot <- matrix(0, 51, rrun)
yrTBN_014 <- matrix(0, 51, rrun)
yrTBN_15 <- matrix(0, 51, rrun)

yrTBNHtot <- matrix(0, 51, rrun)

yrTBPtot <- matrix(0, 51, rrun)

yrHIVpc <- matrix(0, 51, rrun)

yrpoptot <- matrix(0, 51, rrun)
yrpop014 <- matrix(0, 51, rrun)
yrpop1564 <- matrix(0, 51, rrun)
yrpop65 <- matrix(0, 51, rrun)

Pcactv2000 <- matrix(0, rrun, 5)
Pcactv2025 <- matrix(0, rrun, 5)
Pcactv2050 <- matrix(0, rrun, 5)

yrTBRa <- matrix(0, 51, rrun)
yrTBRi <- matrix(0, 51, rrun)
yrTBRa_014 <- matrix(0, 51, rrun)
yrTBRi_014 <- matrix(0, 51, rrun)
yrTBRa_1564 <- matrix(0, 51, rrun)
yrTBRi_1564 <- matrix(0, 51, rrun)
yrTBRa_65 <- matrix(0, 51, rrun)
yrTBRi_65 <- matrix(0, 51, rrun)
yrTBRaR <- matrix(0, 51, rrun)
yrTBRaL <- matrix(0, 51, rrun)

yrTBPI <- matrix(0, 51, rrun)
yrTBPI014 <- matrix(0, 51, rrun)
yrTBPI15p <- matrix(0, 51, rrun)
TBPI2015 <- matrix(0, rrun, 6)

xout <- c()
eee <- c()
for (kkk in 1:rrun) {
  print(kkk)
  for (i in 1:length(nm)) {
    assign(nm[i], as.numeric(para[kkk, i]))
  }
  neta2 <- neta

  system.time(Xn <- FitGo(cntry, 1, c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, (1 / 2), c(0.02, 0.02, 0.8, 0.07)), c(1900, 2050), 0, 0))

  Xn <- as.data.frame(Xn)

  yrTBItot[, kkk] <- Xn$TBIalltot[101:151]
  yrTBI_014[, kkk] <- Xn$"TBIall0-14"[101:151]
  yrTBI_15[, kkk] <- Xn$"TBIall15p"[101:151]

  yrTBIHtot[, kkk] <- Xn$"TBIHdatot"[101:151]

  yrTBMtot[, kkk] <- Xn$TBMtot[101:151]
  yrTBMHtot[, kkk] <- Xn$TBMHdatot[101:151]

  yrTBNtot[, kkk] <- Xn$TBNalltot[101:151]
  yrTBN_014[, kkk] <- Xn$"TBNall0-14"[101:151]
  yrTBN_15[, kkk] <- Xn$TBNall15p[101:151]

  yrTBNHtot[, kkk] <- Xn$TBNHdatot[101:151]

  yrTBPtot[, kkk] <- Xn$TBPalltot[101:151]

  yrHIVpc[, kkk] <- Xn$TBpcHIV[101:151]
  yrpoptot[, kkk] <- Xn$PSIZEALL[101:151]
  yrpop014[, kkk] <- Xn$"YearPsize0-14"[101:151] + Xn$"YearPsize0-14H"[101:151]
  yrpop1564[, kkk] <- (Xn$YearPsize15plus[101:151] - Xn$"YearPsize65+"[101:151]) + (Xn$"YearPsize15+H"[101:151] - Xn$"YearPsize65+H"[101:151])
  yrpop65[, kkk] <- Xn$"YearPsize65+"[101:151] + Xn$"YearPsize65+H"[101:151]

  Xn <- as.data.frame(Xn)

  yrTBRa[, kkk] <- Xn$"TBRatot"[101:151]
  yrTBRi[, kkk] <- Xn$"TBRitot"[101:151]

  yrTBRaR[, kkk] <- Xn$"TBRatot"[101:151]
  yrTBRaL[, kkk] <- Xn$"TBRitot"[101:151]

  yrTBRa_014[, kkk] <- Xn$"TBRa0-14"[101:151]
  yrTBRi_014[, kkk] <- Xn$"TBRi0-14"[101:151]

  yrTBRa_1564[, kkk] <- (Xn$"TBRa15-64"[101:151])
  yrTBRi_1564[, kkk] <- (Xn$"TBRi15-64"[101:151])

  yrTBRa_65[, kkk] <- Xn$"TBRa65+"[101:151]
  yrTBRi_65[, kkk] <- Xn$"TBRi65+"[101:151]

  TBPI <- as.data.frame(TBPI)
  yrTBPI[, kkk] <- TBPI$"All ages"[101:151]
  yrTBPI014[, kkk] <- TBPI$"0-14"[101:151]
  yrTBPI15p[, kkk] <- (TBPI$"15-54"[101:151] + TBPI$"55+"[101:151]) / 2

  adultLTBI <- (TBPI$"15-54"[116] + TBPI$"55+"[116]) / 2
  print(adultLTBI)
  TBPI2015[kkk, 1:5] <- t(TBPI[116, 1:5])

  TBPI2015[kkk, 6] <- adultLTBI
}

setwd(output)
write.table(yrTBPI, "LTBI_annual.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBPI014, "LTBI_annual_014y.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBPI15p, "LTBI_annual_15plus.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(TBPI2015, "LTBI_2015.csv", sep = ",", row.names = FALSE, col.names = FALSE)

colnames(Pcactv2000) <- c("0-14", "15-19", "20-64", "65+", "run")
colnames(Pcactv2025) <- c("0-14", "15-19", "20-64", "65+", "run")
colnames(Pcactv2050) <- c("0-14", "15-19", "20-64", "65+", "run")

write.table(Pcactv2000, "Percent actv per agegrp_2000.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(Pcactv2025, "Percent actv per agegrp_2025.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(Pcactv2050, "Percent actv per agegrp_2050.csv", sep = ",", row.names = FALSE, col.names = FALSE)

write.table(yrTBRa, "Percent_reactivation_tot.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBRa_014, "Percent_reactivation_014.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBRa_1564, "Percent_reactivation_1564.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBRa_65, "Percent_reactivation_65.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBRi, "Percent_reinfection_tot.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBRi_014, "Percent_reinfection_014.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBRi_1564, "Percent_reinfection_1564.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBRi_65, "Percent_reinfection_65.csv", sep = ",", row.names = FALSE, col.names = FALSE)

med_TBN <- matrix(0, 51, 9)
med_TBNH <- matrix(0, 51, 3)
med_TBM <- matrix(0, 51, 3)
med_TBMH <- matrix(0, 51, 3)
med_TBP <- matrix(0, 51, 3)
med_TBI <- matrix(0, 51, 9)
med_TBIH <- matrix(0, 51, 3)
med_pcH <- matrix(0, 51, 3)
med_pop <- matrix(0, 51, 12)

Pcactv2000 <- as.matrix(t(Pcactv2000))
Pcactv2025 <- as.matrix(t(Pcactv2025))
Pcactv2050 <- as.matrix(t(Pcactv2050))

med_actv <- matrix(0, 12, 3)

for (jj in 1:4) {
  med_actv[jj, ] <- c(median(Pcactv2000[jj, ]), min(Pcactv2000[jj, ]), max(Pcactv2000[jj, ]))
  med_actv[(4 + jj), ] <- c(median(Pcactv2025[jj, ]), min(Pcactv2025[jj, ]), max(Pcactv2025[jj, ]))
  med_actv[(8 + jj), ] <- c(median(Pcactv2050[jj, ]), min(Pcactv2050[jj, ]), max(Pcactv2050[jj, ]))
}

med_actv <- as.data.frame(cbind((rep(c("0-14", "15-19", "20-64", "65+"), 3)), rep(c("2000", "2025", "2050"), each = 4), med_actv))
colnames(med_actv) <- c("Age", "Year", "Median", "Min", "Max")

write.table(med_actv, "Active_pc_median.csv", sep = ",", row.names = FALSE)

for (jj in 1:51) {
  med_TBI[jj, c(1:3)] <- c(median(yrTBItot[jj, ]), min(yrTBItot[jj, ]), max(yrTBItot[jj, ]))
  med_TBI[jj, c(4:6)] <- c(median(yrTBI_014[jj, ]), min(yrTBI_014[jj, ]), max(yrTBI_014[jj, ]))
  med_TBI[jj, c(7:9)] <- c(median(yrTBI_15[jj, ]), min(yrTBI_15[jj, ]), max(yrTBI_15[jj, ]))

  med_TBIH[jj, c(1:3)] <- c(median(yrTBIHtot[jj, ]), min(yrTBIHtot[jj, ]), max(yrTBIHtot[jj, ]))

  med_TBM[jj, c(1:3)] <- c(median(yrTBMtot[jj, ]), min(yrTBMtot[jj, ]), max(yrTBMtot[jj, ]))
  med_TBMH[jj, c(1:3)] <- c(median(yrTBMHtot[jj, ]), min(yrTBMHtot[jj, ]), max(yrTBMHtot[jj, ]))

  med_TBN[jj, c(1:3)] <- c(median(yrTBNtot[jj, ]), min(yrTBNtot[jj, ]), max(yrTBNtot[jj, ]))
  med_TBN[jj, c(4:6)] <- c(median(yrTBN_014[jj, ]), min(yrTBN_014[jj, ]), max(yrTBN_014[jj, ]))
  med_TBN[jj, c(7:9)] <- c(median(yrTBN_15[jj, ]), min(yrTBN_15[jj, ]), max(yrTBN_15[jj, ]))

  med_TBNH[jj, c(1:3)] <- c(median(yrTBNHtot[jj, ]), min(yrTBNHtot[jj, ]), max(yrTBNHtot[jj, ]))

  med_TBP[jj, c(1:3)] <- c(median(yrTBPtot[jj, ]), min(yrTBPtot[jj, ]), max(yrTBPtot[jj, ]))

  med_pcH[jj, c(1:3)] <- c(median(yrHIVpc[jj, ]), min(yrHIVpc[jj, ]), max(yrHIVpc[jj, ]))

  med_pop[jj, c(1:3)] <- c(median(yrpoptot[jj, ]), min(yrpoptot[jj, ]), max(yrpoptot[jj, ]))
  med_pop[jj, c(4:6)] <- c(median(yrpop014[jj, ]), min(yrpop014[jj, ]), max(yrpop014[jj, ]))
  med_pop[jj, c(7:9)] <- c(median(yrpop1564[jj, ]), min(yrpop1564[jj, ]), max(yrpop1564[jj, ]))
  med_pop[jj, c(10:12)] <- c(median(yrpop65[jj, ]), min(yrpop65[jj, ]), max(yrpop65[jj, ]))
}

setwd(output)

Year <- seq(2000, 2050, 1)
med_TBN <- as.data.frame(cbind(Year, med_TBN))
med_TBNH <- as.data.frame(cbind(Year, med_TBNH))
med_TBM <- as.data.frame(cbind(Year, med_TBM))
med_TBMH <- as.data.frame(cbind(Year, med_TBMH))
med_TBP <- as.data.frame(cbind(Year, med_TBP))
med_TBI <- as.data.frame(cbind(Year, med_TBI))
med_TBIH <- as.data.frame(cbind(Year, med_TBIH))
med_pcH <- as.data.frame(cbind(Year, med_pcH))
med_pop <- as.data.frame(cbind(Year, med_pop))

write.table(med_TBN, "med_TBN.csv", sep = ",", row.names = F)
write.table(med_TBNH, "med_TBNH.csv", sep = ",", row.names = F)
write.table(med_TBM, "med_TBM.csv", sep = ",", row.names = F)
write.table(med_TBMH, "med_TBMH.csv", sep = ",", row.names = F)
write.table(med_TBP, "med_TBP.csv", sep = ",", row.names = F)
write.table(med_TBI, "med_TBI.csv", sep = ",", row.names = F)
write.table(med_TBIH, "med_TBIH.csv", sep = ",", row.names = F)
write.table(med_pcH, "med_pcH.csv", sep = ",", row.names = F)
write.table(med_pop, "med_pop.csv", sep = ",", row.names = F)

med_TBRa <- matrix(0, 51, 12)
med_TBRi <- matrix(0, 51, 12)

for (jj in 1:51) {
  med_TBRa[jj, c(1:3)] <- c(median(yrTBRa[jj, ]), min(yrTBRa[jj, ]), max(yrTBRa[jj, ]))
  med_TBRa[jj, c(4:6)] <- c(median(yrTBRa_014[jj, ]), min(yrTBRa_014[jj, ]), max(yrTBRa_014[jj, ]))
  med_TBRa[jj, c(7:9)] <- c(median(yrTBRa_1564[jj, ]), min(yrTBRa_1564[jj, ]), max(yrTBRa_1564[jj, ]))
  med_TBRa[jj, c(10:12)] <- c(median(yrTBRa_65[jj, ]), min(yrTBRa_65[jj, ]), max(yrTBRa_65[jj, ]))

  med_TBRi[jj, c(1:3)] <- c(median(yrTBRi[jj, ]), min(yrTBRi[jj, ]), max(yrTBRi[jj, ]))
  med_TBRi[jj, c(4:6)] <- c(median(yrTBRi_014[jj, ]), min(yrTBRi_014[jj, ]), max(yrTBRi_014[jj, ]))
  med_TBRi[jj, c(7:9)] <- c(median(yrTBRi_1564[jj, ]), min(yrTBRi_1564[jj, ]), max(yrTBRi_1564[jj, ]))
  med_TBRi[jj, c(10:12)] <- c(median(yrTBRi_65[jj, ]), min(yrTBRi_65[jj, ]), max(yrTBRi_65[jj, ]))
}

colnames(med_TBRa) <- rep(c("tot", "0-14", "15-64", "65+"), each = 3)
colnames(med_TBRi) <- rep(c("tot", "0-14", "15-64", "65+"), each = 3)
Year <- seq(2000, 2050, 1)
med_TBRa2 <- as.data.frame(cbind(Year, med_TBRa))
med_TBRi2 <- as.data.frame(cbind(Year, med_TBRi))

write.table(med_TBRa2, "TBRa_median.csv", sep = ",", row.names = FALSE)
write.table(med_TBRi2, "TBRi_median.csv", sep = ",", row.names = FALSE)

med_TBRa <- as.data.frame(med_TBRa)
med_TBRi <- as.data.frame(med_TBRi)

ggTBN <- ggplot(med_TBN, aes(x = Year, y = med_TBN[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBN[, 3], ymax = med_TBN[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_TBN[, 2]), colour = "black", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_ribbon(aes(ymin = med_TBNH[, 3], ymax = med_TBNH[, 4], width = 3), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = med_TBNH[, 2]), colour = "orange", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1000), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "All-age notifcation rate (/100,000 pop)") +
  geom_point(aes(x = 2000, y = fitdata[11, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2000, ymin = fitdata[11, 4], ymax = fitdata[11, 5]), size = 0.6) +
  geom_point(aes(x = 2015, y = fitdata[12, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[12, 4], ymax = fitdata[12, 5]), size = 0.6) +
  geom_point(aes(x = 2015, y = fitdata[15, 3]), colour = "orange", size = 2) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[15, 4], ymax = fitdata[15, 5]), colour = "orange", size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBN

ggTBN014 <- ggplot(med_TBN, aes(x = Year, y = med_TBN[, 5])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBN[, 6], ymax = med_TBN[, 7], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_TBN[, 5]), colour = "red", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1000), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "TB notification rate \n 0-14 years (/100,000 pop)") +
  geom_point(aes(x = 2015, y = fitdata[13, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[13, 4], ymax = fitdata[13, 5]), size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBN014

ggTBN15 <- ggplot(med_TBN, aes(x = Year, y = med_TBN[, 8])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBN[, 9], ymax = med_TBN[, 10], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_TBN[, 8]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1000), xlim = c(2000, 2051)) +
  labs(x = "Year", y = "TB notification rate \n 15+ years (/100,000 pop)") +
  geom_point(aes(x = 2015, y = fitdata[14, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[14, 4], ymax = fitdata[14, 5]), size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBN15

ggTBNall <- grid.arrange(ggTBN, ggTBN014, ggTBN15, ncol = 1)

ggsave("ggTBN.pdf", plot = ggTBN, width = 16, height = 9, dpi = 120)
ggsave("ggTBNall.pdf", plot = ggTBNall, width = 7, height = 10, dpi = 120)
ggsave("ggTBN014.pdf", plot = ggTBN014, width = 16, height = 9, dpi = 120)
ggsave("ggTBN15.pdf", plot = ggTBN1554, width = 16, height = 9, dpi = 120)

ggTBM <- ggplot(med_TBM, aes(x = Year, y = med_TBM[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBM[, 3], ymax = med_TBM[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_TBM[, 2], width = 5), colour = "black", size = 1.1) +
  geom_ribbon(aes(ymin = med_TBMH[, 3], ymax = med_TBMH[, 4], width = 3), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = med_TBMH[, 2], width = 5), colour = "orange", size = 1.1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 500), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "All-age mortality rate (/100,000population)") +
  geom_point(aes(x = 2000, y = fitdata[7, 3], size = 8), colour = "black", show.legend = FALSE) +
  geom_errorbar(aes(x = 2000, ymin = fitdata[7, 4], ymax = fitdata[7, 5]), size = 1.1) +
  geom_point(aes(x = 2015, y = fitdata[8, 3], size = 8), colour = "black", show.legend = FALSE) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[8, 4], ymax = fitdata[8, 5]), size = 1.1) +
  geom_point(aes(x = 2000, y = fitdata[9, 3], size = 8), colour = "orange", show.legend = FALSE) +
  geom_errorbar(aes(x = 2000, ymin = fitdata[9, 4], ymax = fitdata[9, 5]), colour = "orange", size = 0.5) +
  geom_point(aes(x = 2015, y = fitdata[10, 3], size = 8), colour = "orange", show.legend = FALSE) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[10, 4], ymax = fitdata[10, 5]), colour = "orange", size = 0.5) +
  theme(axis.line.x = element_line(size = 0.7, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.7, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBM

ggsave("ggTBM.pdf", plot = ggTBM, width = 16, height = 9, dpi = 120)

ggTBP <- ggplot(med_TBP, aes(x = Year, y = med_TBP[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBP[, 3], ymax = med_TBP[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_TBP[, 2]), colour = "black", size = 1.1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1600), xlim = c(2000, 2051)) +
  labs(x = "Year", y = "All-age prevalence rate  (/100,000 population)") +
  geom_point(aes(x = 2015, y = fitdata[16, 3], size = 10), colour = "black", show.legend = FALSE) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[16, 4], ymax = fitdata[16, 5]), size = 0.7) +
  theme(axis.line.x = element_line(size = 0.7, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.7, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBP

ggsave("ggTBP.pdf", plot = ggTBP, width = 16, height = 9, dpi = 120)

ggTBI <- ggplot(med_TBI, aes(x = Year, y = med_TBI[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBI[, 3], ymax = med_TBI[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_TBI[, 2]), colour = "black", size = 1) +
  geom_ribbon(aes(ymin = med_TBIH[, 3], ymax = med_TBIH[, 4], width = 3), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = med_TBIH[, 2]), colour = "orange", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1400), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "All-age incidence rate \n (/100,000 pop)") +
  geom_point(aes(x = 2000, y = fitdata[1, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2000, ymin = fitdata[1, 4], ymax = fitdata[1, 5]), size = 1.1) +
  geom_point(aes(x = 2016, y = fitdata[2, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2016, ymin = fitdata[2, 4], ymax = fitdata[2, 5]), size = 1.1) +
  geom_point(aes(x = 2000, y = fitdata[5, 3]), colour = "orange", size = 2) +
  geom_errorbar(aes(x = 2000, ymin = fitdata[5, 4], ymax = fitdata[5, 5]), colour = "orange", size = 0.5) +
  geom_point(aes(x = 2016, y = fitdata[6, 3]), colour = "orange", size = 2) +
  geom_errorbar(aes(x = 2016, ymin = fitdata[6, 4], ymax = fitdata[6, 5]), colour = "orange", size = 0.5) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBI

ggsave("ggTBI.pdf", plot = ggTBI, width = 16, height = 9, dpi = 120)

ggTBI014 <- ggplot(med_TBI, aes(x = Year, y = med_TBI[, 5])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBI[, 6], ymax = med_TBI[, 7], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_TBI[, 5]), colour = "red", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1400), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "Incidence rate 0-14 years \n (/100,000 pop)") +
  geom_point(aes(x = 2016, y = fitdata[3, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2016, ymin = fitdata[3, 4], ymax = fitdata[3, 5]), size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBI014

ggTBI15 <- ggplot(med_TBI, aes(x = Year, y = med_TBI[, 8])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBI[, 9], ymax = med_TBI[, 10], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_TBI[, 8]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1600), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "Incidence rate 15+ years \n (/100,000 pop)") +
  geom_point(aes(x = 2016, y = fitdata[4, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2016, ymin = fitdata[4, 4], ymax = fitdata[4, 5]), size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBI15
ggTBIall <- grid.arrange(ggTBI, ggTBI014, ggTBI15, ncol = 1)

ggsave("ggTBI.pdf", plot = ggTBI, width = 16, height = 9, dpi = 120)
ggsave("ggTBIall.pdf", plot = ggTBIall, width = 7, height = 10, dpi = 120)
ggsave("ggTBI014.pdf", plot = ggTBI014, width = 16, height = 9, dpi = 120)
ggsave("ggTBI15.pdf", plot = ggTBI15, width = 16, height = 9, dpi = 120)

ggpcH <- ggplot(med_pcH, aes(x = Year, y = med_pcH[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_pcH[, 3], ymax = med_pcH[, 4], width = 3), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = med_pcH[, 2]), colour = "orange", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "All-age percent of TB cases in HIV positive populations") +
  geom_point(aes(x = 2000, y = fitdata[17, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2000, ymin = fitdata[17, 4], ymax = fitdata[17, 5]), size = 0.6) +
  geom_point(aes(x = 2016, y = fitdata[18, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2016, ymin = fitdata[18, 4], ymax = fitdata[18, 5]), size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 20), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 1, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggpcH

ggsave("ggpcHIV.pdf", plot = ggpcH, width = 16, height = 9, dpi = 120)

ggTBR <- ggplot(med_TBRa2, aes(x = Year, y = med_TBRa2[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBRa2[, 3], ymax = med_TBRa2[, 4], width = 3), fill = "grey70") +
  geom_ribbon(aes(ymin = med_TBRi2[, 3], ymax = med_TBRi2[, 4], width = 3), fill = "red", alpha = 0.4) +
  geom_line(aes(y = med_TBRa2[, 2]), colour = "black") +
  geom_line(aes(y = med_TBRi2[, 2]), colour = "red") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  labs(x = "Year", y = "Proportion of All New Cases (%)")
ggTBR

ggTBR014 <- ggplot(med_TBRa2, aes(x = Year, y = med_TBRa2[, 5])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBRa2[, 6], ymax = med_TBRa2[, 7]), fill = "grey70") +
  geom_ribbon(aes(ymin = med_TBRi2[, 6], ymax = med_TBRi2[, 7]), fill = "red", alpha = 0.4) +
  geom_line(aes(y = med_TBRa2[, 5])) +
  geom_line(aes(y = med_TBRi2[, 5], colour = "red")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  labs(x = "Year", y = "Proportion of New Cases aged 0-14years (%)")
ggTBR014

ggTBR1564 <- ggplot(med_TBRa2, aes(x = Year, y = med_TBRa2[, 8])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBRa2[, 9], ymax = med_TBRa2[, 10]), fill = "grey70") +
  geom_ribbon(aes(ymin = med_TBRi2[, 9], ymax = med_TBRi2[, 10]), fill = "red", alpha = 0.4) +
  geom_line(aes(y = med_TBRa2[, 8])) +
  geom_line(aes(y = med_TBRi2[, 8], colour = "red")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  labs(x = "Year", y = "Proportion of New Cases aged 15-64years (%)")
ggTBR1564

ggTBR65 <- ggplot(med_TBRa2, aes(x = Year, y = med_TBRa2[, 11])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBRa2[, 12], ymax = med_TBRa2[, 13]), fill = "grey70") +
  geom_ribbon(aes(ymin = med_TBRi2[, 12], ymax = med_TBRi2[, 13]), fill = "red", alpha = 0.4) +
  geom_line(aes(y = med_TBRa2[, 11])) +
  geom_line(aes(y = med_TBRi2[, 11], colour = "red")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  labs(x = "Year", y = "Proportion of New Cases aged 65+years (%)")
ggTBR65

ggsave("ggTBR.pdf", plot = ggTBR)
ggsave("ggTBR014.pdf", plot = ggTBR014)
ggsave("ggTBR1564.pdf", plot = ggTBR1564)
ggsave("ggTBR65.pdf", plot = ggTBR65)

med_TBPI <- matrix(0, 51, 9)

for (jj in 1:51) {
  med_TBPI[jj, 1:3] <- c(median(yrTBPI[jj, ]), min(yrTBPI[jj, ]), max(yrTBPI[jj, ]))
  med_TBPI[jj, 4:6] <- c(median(yrTBPI014[jj, ]), min(yrTBPI014[jj, ]), max(yrTBPI014[jj, ]))
  med_TBPI[jj, 7:9] <- c(median(yrTBPI15p[jj, ]), min(yrTBPI15p[jj, ]), max(yrTBPI15p[jj, ]))
}

med_TBPI2015 <- matrix(0, 3, 3)
med_TBPI2015[1, 1:3] <- c(median(TBPI2015[, 1]), min(TBPI2015[, 1]), max(TBPI2015[, 1]))
med_TBPI2015[2, 1:3] <- c(median(TBPI2015[, 2]), min(TBPI2015[, 2]), max(TBPI2015[, 2]))
med_TBPI2015[3, 1:3] <- c(median(TBPI2015[, 6]), min(TBPI2015[, 6]), max(TBPI2015[, 6]))

write.table(med_TBPI, "LTBI_median_annual.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(med_TBPI2015, "LTBI_median_2015.csv", sep = ",", row.names = FALSE, col.names = FALSE)

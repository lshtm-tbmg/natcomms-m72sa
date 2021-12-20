library(plyr)
library(grid)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(here)
home <- here()
input <- here("Data")
outputv <- here("Vxoutput_M72_refit")
outpute <- here("Vxoutput_M72_refit/0.epi")
C <- 0

setwd(home)
setwd("Data")
fitdata <- as.data.frame(read.csv("fitdata_sa.csv", header = TRUE, check.names = F))

setwd(home)
setwd(input)

setwd(home)

rrun <- 1000

hivoff <- 0
artoff <- 0
art100 <- 0
artred100 <- 0

source("#DataGrabGSA_v_M72.R")
setwd(home)
source("CFunctions_GSA_uHfix_vratio_M72.R")

setwd(home)
setwd(input)

para <- as.matrix(drop.levels(read.csv("fit1000finalSA_191005_1615seeds.csv", header = TRUE, check.names = F)))

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

Pcactv2000 <- matrix(0, 5, rrun)
Pcactv2025 <- matrix(0, 5, rrun)
Pcactv2050 <- matrix(0, 5, rrun)

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

source("#vxSA_M72.R")

vaxfol <- paste(typen, effDis, durs, vage, cover, sep = "_")
vaxfol_path <- here(paste0("Vxoutput_M72_refit/", vaxfol, "/1.heatmap"))
if (!dir.exists(vaxfol_path)) dir.create(vaxfol_path, recursive = TRUE)
if (!dir.exists(here("Vxoutput_M72_refit/0.epi"))) dir.create(here("Vxoutput_M72_refit/0.epi"), recursive = TRUE)

vcumDALY <- matrix(0, rrun, 7)
cumDALY <- matrix(0, rrun, 11)
annualDALY <- c()
cumDALY_list <- list()
vcumDALY_list <- list()
inc1550 <- c()
incnum <- c()
treatnum <- c()
numart <- c()
CAV2535all <- c()
CAV2550all <- c()
DAV2550all <- c()
NVax2550all <- c()

xout <- c()
eee <- c()
for (kkk in 1:rrun) {
  print(kkk)

  setwd(home)

  cumulvx <- c()
  vaxgive <- c()
  vaxgiveyr <- c()
  cumulvxyrM <- c()
  cumulvxyrI <- c()
  NumV <- c()
  inc2050 <- matrix(0, (combn + 1), 11)
  mort2050 <- c()
  inc2035 <- c()
  mort2035 <- c()
  rrun_dfvx <- c()
  vacnames <- c()

  inc2050 <- matrix(0, (combn + 1), 10)
  mort2050 <- matrix(0, (combn + 1), 7)
  inc2035 <- matrix(0, (combn + 1), 10)
  mort2035 <- matrix(0, (combn + 1), 7)

  LTBI2015 <- rep(0, 2)

  dt <- (1 / 2)
  year1 <- 1900
  yearend <- 2050

  for (i in 1:length(nm)) {
    assign(nm[i], as.numeric(para[kkk, i]))
  }
  neta2 <- neta

  Xn <- FitGo(cntry, 1, c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, dt, c(0.02, 0.02, 0.8, 0.07)), c(year1, yearend), 0, C)

  Xn <- as.data.frame(Xn)

  treatnum <- newtreat[125:151, ]
  numart <- artnum[126:151, ]

  setwd(home)
  setwd("Vxoutput_M72_refit")
  setwd(vaxfol)
  write.table(Xn, paste("Xnbaseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBP_age, paste("TBPneg_age_baseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBPH_age, paste("TBPHIV_age_baseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBMo_age, paste("TBdeaths_age_baseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(treatnum, paste("TBtreat_baseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(numart, paste("ART_baseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBPIage, paste("TBPI_age_baseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

  setwd(home)

  new_active <- cbind(TBAc, 0, 0)
  new_ac_age <- cbind(TBAc_age, 0, 0)
  new_mort <- cbind(TBMo, 0, 0)

  inc1550 <- rbind(inc1550, TBIall[116:151, 1])
  incnum <- rbind(incnum, t(newinc[116:151]))

  inc2050[1, ] <- c(TBIall[151, ], TBIHda[151, ], TBI[151, c(1, 2, 10)], "baseline")
  mort2050[1, ] <- c(TBM[151, 1:5], TBMHda[151, ], "baseline")
  inc2035[1, ] <- c(TBIall[136, ], TBIHda[136, ], TBI[136, c(1, 2, 10)], "baseline")
  mort2035[1, ] <- c(TBM[136, 1:5], TBMHda[136, ], "baseline")
  LTBI2015[1] <- LTBIc
  LTBI2015[2] <- LTBIa
  write.csv(LTBI2015, "LTBI_2015.csv")

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

  casesall <- sum(new_actv[201:202, ] + new_actvH[201:202, ])

  Pcactv2000[1, kkk] <- sum(new_actv[201:202, 1:15] + new_actvH[201:202, 1:15]) / casesall * 100
  Pcactv2000[2, kkk] <- sum(new_actv[201:202, 16:20] + new_actvH[201:202, 16:20]) / casesall * 100
  Pcactv2000[3, kkk] <- sum(new_actv[201:202, 21:65] + new_actvH[201:202, 21:65]) / casesall * 100
  Pcactv2000[4, kkk] <- sum(new_actv[201:202, 66:Mnage] + new_actvH[201:202, 66:Mnage]) / casesall * 100

  Pcactv2000[5, kkk] <- kkk

  Pcactv2025[1, kkk] <- sum(new_actv[251:252, 1:15] + new_actvH[251:252, 1:15]) / casesall * 100
  Pcactv2025[2, kkk] <- sum(new_actv[251:252, 16:20] + new_actvH[251:252, 16:20]) / casesall * 100
  Pcactv2025[3, kkk] <- sum(new_actv[251:252, 21:65] + new_actvH[251:252, 21:65]) / casesall * 100
  Pcactv2025[4, kkk] <- sum(new_actv[251:252, 66:Mnage] + new_actvH[251:252, 66:Mnage]) / casesall * 100

  Pcactv2025[5, kkk] <- kkk
  casesall <- sum(new_actv[301:302, ] + new_actvH[301:302, ])

  Pcactv2050[1, kkk] <- sum(new_actv[301:302, 1:15] + new_actvH[301:302, 1:15]) / casesall * 100
  Pcactv2050[2, kkk] <- sum(new_actv[301:302, 16:20] + new_actvH[301:302, 16:20]) / casesall * 100
  Pcactv2050[3, kkk] <- sum(new_actv[301:302, 21:65] + new_actvH[301:302, 21:65]) / casesall * 100
  Pcactv2050[4, kkk] <- sum(new_actv[301:302, 66:Mnage] + new_actvH[301:302, 66:Mnage]) / casesall * 100

  Pcactv2050[5, kkk] <- kkk

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
  TBPI2015[kkk, 1:5] <- t(TBPI[116, 1:5])

  TBPI2015[kkk, 6] <- adultLTBI

  eee <- cbind(Xn, 0, 0, 0, 0, 0, 0, 0)
  colnames(eee) <- c(colnames(Xn), "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")

  rrun_dfvx <- eee

  for (nn in typen:typen) {
    count <- 0
    coms <- matrix(0, combn, 5)
    for (vv in 1:length(effInf)) {
      for (zz in 1:length(effDis)) {
        for (xx in 1:length(durs)) {
          for (ss in 1:length(Hsafe)) {
            for (hh in 1:length(Hvax)) {
              count <- count + 1
              coms[count, ] <- c(effInf[vv], effDis[zz], durs[xx], Hsafe[ss], Hvax[hh])

              ticI <- effInf[vv]
              ticD <- effDis[zz]
              toc <- durs[xx]
              covH <- Hsafe[ss]
              VEH <- Hvax[hh]
              print(c(nn, ticI, ticD, toc, covH, VEH))

              X <- FitGo(cntry, c(nn, cover, coverM, ticI, ticD, toc, fms, vage, covH, VEH), c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, dt, c(0.02, 0.02, 0.8, 0.07)), c(year1, yearend), 0, C)

              if (nn == 1) {
                vtp <- "PPI"
              } else if (nn == 2) {
                vtp <- "PRI"
              } else if (nn == 3) {
                vtp <- "PSI_LR"
              }

              vxtyp <- paste(vtp, "_", ticI, "_", ticD, "_", toc, "_", covH, "_", VEH)
              vacnames <- c(vacnames, vxtyp)

              eee <- cbind(X, nn, ticI, ticD, toc, covH, VEH, count)
              colnames(eee) <- c(colnames(X), "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")

              rrun_dfvx <- rbind(rrun_dfvx, eee)

              new_active <- cbind(new_active, TBAc, nn, count)
              new_ac_age <- cbind(new_ac_age, TBAc_age, nn, count)
              new_mort <- cbind(new_mort, TBMo, nn, count)
              NumV <- cbind(NumV, NV, nn, count)

              inc2050[(count + 1), ] <- c(TBIall[151, ], TBIHda[151, ], TBI[151, c(1, 2, 10)], vxtyp)
              mort2050[(count + 1), ] <- c(TBM[151, 1:5], TBMHda[151, ], vxtyp)
              inc2035[(count + 1), ] <- c(TBIall[136, ], TBIHda[136, ], TBI[136, c(1, 2, 10)], vxtyp)
              mort2035[(count + 1), ] <- c(TBM[136, 1:5], TBMHda[136, ], vxtyp)

              colnames(inc2050) <- c(colnames(TBI[, c(1, 2, 10)]), "HAllage", "H0-14", "H15plus", "NAllage", "N0-14", "N15plus", "vx")
              colnames(inc2035) <- c(colnames(TBI[, c(1, 2, 10)]), "HAllage", "H0-14", "H15plus", "NAllage", "N0-14", "N15plus", "vx")
              colnames(mort2035) <- c(colnames(TBM[, 1:5]), "TBMHda_all", "vx")
              colnames(mort2050) <- c(colnames(TBM[, 1:5]), "TBMHda_all", "vx")
              inc1550 <- rbind(inc1550, TBIall[116:151, 1])
              incnum <- rbind(incnum, t(newinc[116:151]))
              treatnum <- newtreat[125:151, ]
              numart <- artnum[126:151, ]
            }
          }
        }
      }
    }
  }

  assign("new_active", new_active, envir = .GlobalEnv)
  assign("new_ac_age", new_ac_age, envir = .GlobalEnv)
  assign("new_mort", new_mort, envir = .GlobalEnv)
  assign("NumV", NumV, envir = .GlobalEnv)
  assign("inc2050", inc2050, envir = .GlobalEnv)
  assign("mort2050", mort2050, envir = .GlobalEnv)
  assign("inc2035", inc2035, envir = .GlobalEnv)
  assign("mort2035", mort2035, envir = .GlobalEnv)
  assign("rrun_dfvx", rrun_dfvx, envir = .GlobalEnv)
  assign("vacnames", vacnames, envir = .GlobalEnv)
  assign("kkk", kkk, envir = .GlobalEnv)
  assign("inc1550", inc1550, envir = .GlobalEnv)
  assign("incnum", incnum, envir = .GlobalEnv)
  assign("treatnum", treatnum, envir = .GlobalEnv)
  assign("numart", numart, envir = .GlobalEnv)

  setwd(home)
  setwd("Vxoutput_M72_refit")
  setwd(vaxfol)
  write.table(new_active, paste("new_active_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(new_mort, paste("new_mort_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(NumV, paste("number_vaccinated_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(inc2050, paste("inc_rates_2050_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(mort2050, paste("mort_rates_2050_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(inc2035, paste("inc_rates_2035", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(mort2035, paste("mort_rates_2035_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

  write.table(inc1550, paste("inc_rates_2015to50_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(incnum, paste("number_incident", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBP_age, paste("TBPneg_age", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBPH_age, paste("TBPHIV_age", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBMo_age, paste("TBdeaths_age", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(TBPIage, paste("TBPI_age", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

  write.table(treatnum, paste("TBtreat", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(numart, paste("ART", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

  print(vacnames)
  setwd(home)
  source("#%reduction_clusterGSA_M72.R")
  setwd(home)
  source("#DALY_gen_SA.R")
}

setwd(home)
setwd("Vxoutput_M72_refit")
setwd(vaxfol)

cumDALY <- rbindlist(cumDALY_list)
vcumDALY <- rbindlist(vcumDALY_list)

write.table(vcumDALY, "1_cumulativeVaccinecostDALY.csv", sep = ",", row.names = FALSE)

write.table(cumDALY, "1_cumulativecostDALY.csv", sep = ",", row.names = FALSE)

med_cD2_list <- list()

med_cD2_list[[1]] <- cumDALY[, lapply(.SD, median), .SDcol = c("costall", "d_costall", "costall_nd", "d_costall_nd", "costall_na", "d_costall_na", "costall_nad", "d_costall_nad", "DALYdiff", "d_DALYdiff", "cumcostDALY", "cumcostDALYdis", "cumcostDALY_nd", "cumcostDALYdis_nd", "cumcostDALY_na", "cumcostDALYdis_na", "cumcostDALY_nad", "cumcostDALYdis_nad")]
med_cD2_list[[2]] <- cumDALY[, lapply(.SD, min), .SDcol = c("costall", "d_costall", "costall_nd", "d_costall_nd", "costall_na", "d_costall_na", "costall_nad", "d_costall_nad", "DALYdiff", "d_DALYdiff", "cumcostDALY", "cumcostDALYdis", "cumcostDALY_nd", "cumcostDALYdis_nd", "cumcostDALY_na", "cumcostDALYdis_na", "cumcostDALY_nad", "cumcostDALYdis_nad")]
med_cD2_list[[3]] <- cumDALY[, lapply(.SD, max), .SDcol = c("costall", "d_costall", "costall_nd", "d_costall_nd", "costall_na", "d_costall_na", "costall_nad", "d_costall_nad", "DALYdiff", "d_DALYdiff", "cumcostDALY", "cumcostDALYdis", "cumcostDALY_nd", "cumcostDALYdis_nd", "cumcostDALY_na", "cumcostDALYdis_na", "cumcostDALY_nad", "cumcostDALYdis_nad")]

med_cD2 <- rbindlist(med_cD2_list)
med_cD2 <- t(med_cD2)
colnames(med_cD2) <- c("median", "min", "max")

write.table(med_cD2, "1_median_CcD.csv", sep = ",", row.names = TRUE, col.names = TRUE)

numyr <- nrow(ledf[, 1])
med_annual <- matrix(0, numyr, (3 * (ncol(annualDALY) - 1)))

annualDALY <- as.data.frame(annualDALY)

for (ii in 2:(ncol(annualDALY))) {
  for (jj in 1:numyr) {
    seq1 <- seq(1, (rrun * numyr), numyr)
    med_annual[jj, (3 * (ii - 1) - 2)] <- median(annualDALY[seq1 + jj - 1, ii])
    med_annual[jj, (3 * (ii - 1) - 1)] <- min(annualDALY[seq1 + jj - 1, ii])
    med_annual[jj, (3 * (ii - 1))] <- max(annualDALY[seq1 + jj - 1, ii])
  }
}

med_annual <- cbind(ledf[, 1], med_annual)

colnames(med_annual) <- c("Year", rep((colnames(annualDALY[, -1])), times = 1, each = 3))

write.table(med_annual, paste("1_median_CcD_annual.csv", sep = ""), sep = ",", row.names = FALSE)

setwd(home)
source("#DALYplots_SA.R")

setwd(home)
setwd("Vxoutput_M72_refit")
setwd("0.epi")
write.table(yrTBPI, "LTBI_annual.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBPI014, "LTBI_annual_014y.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(yrTBPI15p, "LTBI_annual_15plus.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(TBPI2015, "LTBI_2015.csv", sep = ",", row.names = FALSE, col.names = FALSE)

Pcactv2000 <- as.matrix(t(Pcactv2000))
Pcactv2025 <- as.matrix(t(Pcactv2025))
Pcactv2050 <- as.matrix(t(Pcactv2050))
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
  coord_cartesian(ylim = c(0, 1250), xlim = c(1999, 2051)) +
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
  coord_cartesian(ylim = c(0, 1250), xlim = c(1999, 2051)) +
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
  coord_cartesian(ylim = c(0, 1250), xlim = c(2000, 2051)) +
  labs(x = "Year", y = "TB notification rate \n 15+ years (/100,000 pop)") +
  geom_point(aes(x = 2015, y = fitdata[14, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[14, 4], ymax = fitdata[14, 5]), size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBN15

ggTBNall <- grid.arrange(ggTBN, ggTBN014, ggTBN15, ncol = 1)

ggsave("ggTBN.png", plot = ggTBN, width = 16, height = 9, dpi = 120)
ggsave("ggTBNall.png", plot = ggTBNall, width = 7, height = 10, dpi = 120)
ggsave("ggTBN014.png", plot = ggTBN014, width = 16, height = 9, dpi = 120)
ggsave("ggTBN15.png", plot = ggTBN15, width = 16, height = 9, dpi = 120)

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
    axis.text = element_text(size = 22), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBM

ggsave("ggTBM.png", plot = ggTBM, width = 16, height = 9, dpi = 120)

ggTBP <- ggplot(med_TBP, aes(x = Year, y = med_TBP[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBP[, 3], ymax = med_TBP[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_TBP[, 2]), colour = "black", size = 1.1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 2000), xlim = c(2000, 2051)) +
  labs(x = "Year", y = "All-age prevalence rate  (/100,000 population)") +
  geom_point(aes(x = 2015, y = fitdata[16, 3], size = 10), colour = "black", show.legend = FALSE) +
  geom_errorbar(aes(x = 2015, ymin = fitdata[16, 4], ymax = fitdata[16, 5]), size = 0.7) +
  theme(axis.line.x = element_line(size = 0.7, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.7, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 22), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBP

ggsave("ggTBP.png", plot = ggTBP, width = 16, height = 9, dpi = 120)

ggTBI <- ggplot(med_TBI, aes(x = Year, y = med_TBI[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_TBI[, 3], ymax = med_TBI[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_TBI[, 2]), colour = "black", size = 1) +
  geom_ribbon(aes(ymin = med_TBIH[, 3], ymax = med_TBIH[, 4], width = 3), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = med_TBIH[, 2]), colour = "orange", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 2000), xlim = c(1999, 2051)) +
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
  coord_cartesian(ylim = c(0, 2000), xlim = c(1999, 2051)) +
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
  coord_cartesian(ylim = c(0, 2000), xlim = c(1999, 2051)) +
  labs(x = "Year", y = "Incidence rate 15+ years \n (/100,000 pop)") +
  geom_point(aes(x = 2016, y = fitdata[4, 3]), colour = "black", size = 2) +
  geom_errorbar(aes(x = 2016, ymin = fitdata[4, 4], ymax = fitdata[4, 5]), size = 0.6) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBI15
ggTBIall <- grid.arrange(ggTBI, ggTBI014, ggTBI15, ncol = 1)

ggsave("ggTBI.png", plot = ggTBI, width = 16, height = 9, dpi = 120)
ggsave("ggTBIall.png", plot = ggTBIall, width = 7, height = 10, dpi = 120)
ggsave("ggTBI014.png", plot = ggTBI014, width = 16, height = 9, dpi = 120)
ggsave("ggTBI15.png", plot = ggTBI15, width = 16, height = 9, dpi = 120)

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

ggsave("ggpcHIV.png", plot = ggpcH, width = 16, height = 9, dpi = 120)

med_TBRa2b <- med_TBRa2
colnames(med_TBRa2) <- NULL

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
  labs(x = "Year", y = "Proportion of All New Cases (%)") +
  theme(
    axis.text = element_text(size = 22), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  )
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
  labs(x = "Year", y = "Proportion of New Cases aged 0-14years (%)") +
  theme(
    axis.text = element_text(size = 22), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  )
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
  labs(x = "Year", y = "Proportion of New Cases aged 15-64years (%)") +
  theme(
    axis.text = element_text(size = 22), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  )
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
  labs(x = "Year", y = "Proportion of New Cases aged 65+years (%)") +
  theme(
    axis.text = element_text(size = 22), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  )
ggTBR65

ggsave("ggTBR.png", plot = ggTBR, width = 9, height = 8, dpi = 120)
ggsave("ggTBR014.png", plot = ggTBR014, width = 9, height = 8, dpi = 120)
ggsave("ggTBR1564.png", plot = ggTBR1564, width = 9, height = 8, dpi = 120)
ggsave("ggTBR65.png", plot = ggTBR65, width = 9, height = 8, dpi = 120)

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

setwd(home)
setwd("Vxoutput_M72_refit")
setwd(vaxfol)

prevLTBI <- c()
for (i in 1:rrun) {
  holder <- fread(paste("TBPI_age_baseline", i, ".csv", sep = ""))

  prevLTBI <- rbind(prevLTBI, holder[127, ])
}

prevLTBI <- t(data.matrix(prevLTBI))
med_TBPIage <- matrix(0, 101, 3)

for (jj in 1:101) {
  med_TBPIage[jj, 1:3] <- c(median(prevLTBI[jj, ]), min(prevLTBI[jj, ]), max(prevLTBI[jj, ]))
}

write.table(prevLTBI, "1.prevLTBI.csv", sep = ",", row.names = FALSE, col.names = FALSE)

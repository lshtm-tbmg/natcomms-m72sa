library(plyr)

C <- 1
if (C == 0) {
  home <- "/Users/lsh355020/R/Gates_sa"
}
if (C == 1) {
  home <- "/home/lsh355020/SA_Gates/"
}

if (C == 0) {
  output <- "/Users/lsh355020/R/Gates_sa/output"
}
if (C == 1) {
  output <- "/home/lsh355020/SA_Gates/output"
}

setwd(home)

cntry <- "South Africa"

source("#DataGrabGSA.R")
setwd(home)

source("CFunctions_GSA_uHfix.R")

nm <- c(pararange[, 1], "p0")

setwd(home)
para <- read.csv("fit1000finalSA.csv", sep = ",")
n_p <- nrow(para)
setwd(home)

source("#vxSA_s15.R")

setwd(home)
cumulvx <- c()
vaxgive <- c()
vaxgiveyr <- c()
cumulvxyrM <- c()
cumulvxyrI <- c()
NumV <- c()
inc2050 <- matrix(0, (combn + 1), 10)
mort2050 <- matrix(0, (combn + 1), 7)
inc2035 <- matrix(0, (combn + 1), 10)
mort2035 <- matrix(0, (combn + 1), 7)
rrun_dfvx <- c()
vacnames <- c()
LTBI2015 <- rep(0, 2)

if (C == 0) {
  kkk <- 1
}
if (C == 1) {
  kkk <- as.numeric(Sys.getenv("SGE_TASK_ID"))
}

print(kkk)

dt <- (1 / 2)
year1 <- 1900
yearend <- 2050

year <- c(seq(year1, yearend, 1), rep(0, ((1 / dt) * (yearend - year1 + 1) - (yearend - year1 + 1))))

for (i in 1:length(nm)) {
  assign(nm[i], as.numeric(para[kkk, i]))
}
neta2 <- neta

Xn <- FitGo(cntry, 1, c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, dt, c(0.02, 0.02, 0.8, 0.07)), c(year1, yearend), 0, C)

setwd(home)

new_active <- cbind(TBAc, 0, 0)
new_ac_age <- cbind(TBAc_age, 0, 0)
new_mort <- cbind(TBMo, 0, 0)
inc2050[1, ] <- c(TBIall[151, ], TBIHda[151, ], TBI[151, c(1, 2, 10)], "baseline")
mort2050[1, ] <- c(TBM[151, 1:5], TBMHda[151, ], "baseline")
inc2035[1, ] <- c(TBIall[136, ], TBIHda[136, ], TBI[136, c(1, 2, 10)], "baseline")
mort2035[1, ] <- c(TBM[136, 1:5], TBMHda[136, ], "baseline")
LTBI2015[1] <- LTBIc
LTBI2015[2] <- LTBIa
write.csv(LTBI2015, "LTBI_2015.csv")

eee <- cbind(Xn, 0, 0, 0, 0, 0, 0, 0)
colnames(eee) <- c(colnames(Xn), "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count")
rrun_dfvx <- eee

for (nn in 2:typen) {
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
          }
        }
      }
    }
  }
}

assign("new_active", new_active, envir = .GlobalEnv)
assign("new_ac_age", new_ac_age, envir = .GlobalEnv)
assign("new_mort", new_mort, envir = .GlobalEnv)
assign("inc2050", inc2050, envir = .GlobalEnv)
assign("mort2050", mort2050, envir = .GlobalEnv)
assign("inc2035", inc2035, envir = .GlobalEnv)
assign("mort2035", mort2035, envir = .GlobalEnv)
assign("rrun_dfvx", rrun_dfvx, envir = .GlobalEnv)
assign("vacnames", vacnames, envir = .GlobalEnv)
assign("kkk", kkk, envir = .GlobalEnv)

setwd(home)
setwd("Vxoutputs15")
write.table(new_active, paste("new_active_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(new_mort, paste("new_mort_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(inc2050, paste("inc_rates_2050_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(mort2050, paste("mort_rates_2050_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(inc2035, paste("inc_rates_2035", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(mort2035, paste("mort_rates_2035_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

print(vacnames)
setwd(home)
source("#%reduction_clusterGSAs15.R")
setwd(home)

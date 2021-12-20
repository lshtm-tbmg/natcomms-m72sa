library(plyr)

C <- 1

if (C == 0) {
  home <- "/Users/Rebecca/GatesVx_india"
}
if (C == 1) {
  home <- "/home/lsh355020/India_Gates/"
}

if (C == 0) {
  output <- "/Users/Rebecca/GatesVx_india/output"
}
if (C == 1) {
  output <- "/home/lsh355020/India_Gates/output"
}

if (C == 0) {
  input <- "/Users/Rebecca/GatesVx_india/Data"
}
if (C == 1) {
  input <- "/home/lsh355020/India_Gates/Data"
}

setwd(home)

source("#DataGrabGIn.R")
setwd(home)

source("CFunctions_GIn.R")

setwd(home)
setwd(input)

para <- as.matrix(drop.levels(read.csv("Afithit_13hits_para_171123.csv", header = TRUE, check.names = F)))

setwd(home)

nm <- c(pararange[, 1], "p0")
cntry <- "India"
p1950 <- 376325

source("#vxIn.R")

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

for (i in 1:length(nm)) {
  assign(nm[i], as.numeric(para[kkk, i]))
}
neta2 <- neta

Xn <- FitGo(cntry, 1, c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, dt, c(0.02, 0.02, 0.8, 0.07)), c(year1, yearend), 0, C)

setwd("Vxoutput")
write.table(Xn, paste("Xnbaseline", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
setwd(home)

new_active <- cbind(TBAc, 0, 0)
new_ac_age <- cbind(TBAc_age, 0, 0)
new_mort <- cbind(TBMo, 0, 0)
inc2050[1, ] <- c(TBI[151, ], "baseline")
mort2050 <- rbind(TBM[151, ])
inc2035 <- rbind(TBI[136, ])
mort2035 <- rbind(TBM[136, ])

eee <- cbind(Xn, 0, 0, 0, 0, 0)
colnames(eee) <- c(colnames(Xn), "type", "VE_I", "VE_D", "dur", "count")
rrun_dfvx <- eee

for (nn in 1:typen) {
  count <- 0
  coms <- matrix(0, combn, 3)
  for (vv in 1:length(effInf)) {
    for (zz in 1:length(effDis)) {
      for (xx in 1:length(durs)) {
        count <- count + 1
        coms[count, ] <- c(effInf[vv], effDis[zz], durs[xx])

        ticI <- effInf[vv]
        ticD <- effDis[zz]
        toc <- durs[xx]
        print(c(nn, ticI, ticD, toc))

        X <- FitGo(cntry, c(nn, cover, coverM, ticI, ticD, toc, fms, vage), c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, dt, c(0.02, 0.02, 0.8, 0.07)), c(year1, yearend), 0, C)

        if (nn == 1) {
          vtp <- "PPI"
        } else if (nn == 2) {
          vtp <- "PRI"
        } else if (nn == 3) {
          vtp <- "PSI_LR"
        }

        vxtyp <- paste(vtp, "_", ticI, "_", ticD, "_", toc)
        vacnames <- c(vacnames, vxtyp)

        eee <- cbind(X, nn, ticI, ticD, toc, count)
        colnames(eee) <- c(colnames(X), "type", "VE_I", "VE_D", "dur", "count")
        rrun_dfvx <- rbind(rrun_dfvx, eee)

        new_active <- cbind(new_active, TBAc, nn, count)
        new_ac_age <- cbind(new_ac_age, TBAc_age, nn, count)
        new_mort <- cbind(new_mort, TBMo, nn, count)
        NumV <- cbind(NumV, NV, nn, count)

        inc2050[(count + 1), ] <- c(TBI[151, ], vxtyp)
        colnames(inc2050) <- c(colnames(TBI), "vx")

        mort2050 <- rbind(mort2050, TBM[151, ])

        inc2035 <- rbind(inc2035, TBI[136, ])
        mort2035 <- rbind(mort2035, TBM[136, ])
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

setwd(home)
setwd("Vxoutput")
write.table(new_active, paste("new_active_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(new_mort, paste("new_mort_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(NumV, paste("number_vaccinated_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(inc2050, paste("inc_rates_2050_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(mort2050, paste("mort_rates_2050_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(inc2035, paste("inc_rates_2035", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(mort2035, paste("mort_rates_2035_", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

print(vacnames)
setwd(home)
source("#%reduction_clusterGIn.R")
setwd(home)

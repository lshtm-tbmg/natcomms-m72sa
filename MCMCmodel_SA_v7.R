MCMCmodel <- function(para) {
  print("testMCMC")
  para <- as.data.frame(para)
  para
  print(para)

  para[34, 1] <- para[4, 1] * (para[27, 1] / para[16, 1])
  print("fHchild")
  print(para[34, 1])
  print("pHchild")
  para[30, 1] <- para[8, 1] * (para[29, 1] / para[10, 1])
  print(para[30, 1])
  if (para[30, 1] > 1) {
    para[30, 1] <- 1
  }
  print(para[30, 1])

  para[35, 1] <- para[5, 1] * (para[28, 1])
  print("vHchild")
  print(para[35, 1])
  if (para[35, 1] > 1) {
    para[35, 1] <- 1
  }
  print(para[35, 1])

  para[36, 1] <- para[6, 1] * (para[31, 1] / para[14, 1])
  print("rHchild")
  print(para[36, 1])
  if (para[36, 1] > 1) {
    para[36, 1] <- 1
  }
  print(para[36, 1])

  para[21, 1] <- para[11, 1]
  para[13, 1] <- para[14, 1]
  para[22, 1] <- para[18, 1]
  para[26, 1] <- para[19, 1]
  para[9, 1] <- para[10, 1]
  para[17, 1] <- para[16, 1]
  para[24, 1] <- para[23, 1]

  if ((para[12, 1] >= para[23, 1]) & (para[24, 1] >= para[23, 1]) & (para[16, 1] >= para[27, 1])) {
    paracheck <- 1
  } else {
    paracheck <- 0
  }
  print(paracheck)

  if (paracheck == 0) {
    AR <- c(0, 0)
    hittrack <- c(hittrack, 99)
  } else {
    C <- 1
    if (C == 0) {
      home <- "/Users/lsh355020/R/Gates_sa"
    }
    if (C == 1) {
      home <- "/home/lsh355020/SA_Gates"
    }

    setwd(home)

    source("#DataGrabGSA_nopara_v.R")

    setwd(home)
    source("CFunctions_GSA_uHfix_vratio.R")
    cntry <- "South Africa"
    nm <- c(pararange[, 1], "p0")

    dt <- (1 / 2)
    year1 <- 1900
    yearend <- 2050

    typen <- 0

    for (i in 1:(length(nm) + 1)) {
      assign(nm[i], as.numeric(para[i, ]), envir = .GlobalEnv)
    }

    vHadult <- vHratio * vadult
    assign("vHadult", vHadult, envir = .GlobalEnv)

    print("vHadult")
    print(vHratio)
    print(vadult)
    print(vHadult)
    print("rHadult")
    print(radult)
    print(rHadult)

    neta2 <- neta
    Xn <- FitGo(cntry, 1, c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, dt, c(0.02, 0.02, 0.8, 0.07)), c(year1, yearend), 0, C)
    head(Xn)

    year <- c(seq(year1, yearend, 1), rep(0, ((1 / dt) * (yearend - year1 + 1) - (yearend - year1 + 1))))
    xout <- cbind(Xn, year)

    assign("xout", xout, envir = .GlobalEnv)

    source("fithit_SA.R")

    fithit
    hittrack <- c(hittrack, fithit)
    hitvar <- rbind(hitvar, inbound)

    AR <- c(0, 0)
    print(fithit)
    print(hittrack)
    print(length(hittrack))

    if (fithit >= hits) AR <- c(1, 0)
    if (fithit >= hits) accvar <- rbind(accvar, inbound)
    if (fithit >= hits) para16 <- rbind(para16, para)

    if (fithit >= hits) write.table(para16, paste("ABCacceptinterim_v7b_", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)
  }

  print("AR")
  print(AR)
  assign("hittrack", hittrack, envir = .GlobalEnv)
  assign("hitvar", hitvar, envir = .GlobalEnv)
  assign("accvar", accvar, envir = .GlobalEnv)
  assign("para16", para16, envir = .GlobalEnv)

  return(AR)
}

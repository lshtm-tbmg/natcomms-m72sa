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
  output <- "/home/lsh355020/SA_Gates/output_vHratio"
}

setwd(home)

cntry <- "South Africa"

source("#DataGrabGSA_v.R")
setwd(home)

source("CFunctions_GSA_uHfix_vratio.R")

nm <- c(pararange[, 1], "p0")

n_p <- 1000

if (C == 0) {
  job <- 0
}
if (C == 1) {
  job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
}

randparam <- mat.or.vec(n_p, (length(nm) + 1))
colnames(randparam) <- c("set", nm)
randparam[1:n_p, 1] <- seq(1, n_p, 1)

for (dd in 1:(length(nm) - 1)) {
  randparam[1:n_p, dd + 1] <- runif(n_p, as.numeric(pararange[dd, 2]), as.numeric(pararange[dd, 3]))
}

for (dd in 1:n_p) {
  randparam[dd, "uiscaleC"] <- runif(1, randparam[dd, "uiscaleA"], as.numeric(pararange[12, 3]))

  if (randparam[dd, "fadult"] < randparam[dd, "fHadult"]) {
    randparam[dd, "fHadult"] <- runif(1, as.numeric(pararange[27, 2], randparam[dd, "fadult"]))
  }

  randparam[dd, "rmortTBH"] <- rnorm(1, randparam[dd, "rmortTB"], 0.25)
  if (randparam[dd, "rmortTBH"] >= 1) {
    randparam[dd, "rmortTBH"] <- 0.99
  }
  if (randparam[dd, "rmortTBH"] <= -1) {
    randparam[dd, "rmortTBH"] <- -0.99
  }

  randparam[dd, "fHchild"] <- randparam[dd, "fchild"] * (randparam[dd, "fHadult"] / randparam[dd, "fadult"])
  randparam[dd, "pHchild"] <- randparam[dd, "pchild"] * (randparam[dd, "pHadult"] / randparam[dd, "padult"])
  randparam[dd, "vHchild"] <- randparam[dd, "vchild"] * randparam[dd, "vHratio"]
  randparam[dd, "rHchild"] <- randparam[dd, "rchild"] * (randparam[dd, "rHadult"] / randparam[dd, "radult"])

  randparam[dd, "nelderly"] <- randparam[dd, "n"]
  randparam[dd, "relderly"] <- randparam[dd, "radult"]
  randparam[dd, "velderly"] <- randparam[dd, "vadult"]
  randparam[dd, "CDRscaleE"] <- randparam[dd, "CDRscale"]
  randparam[dd, "pelderly"] <- randparam[dd, "padult"]
  randparam[dd, "felderly"] <- randparam[dd, "fadult"]
  randparam[dd, "uiscaleE"] <- randparam[dd, "uiscaleA"]
}

head(randparam)

randparam[1:n_p, (length(nm) + 1)] <- rep(psz1900, n_p)
print(randparam)

setwd(home)
setwd("parameters_vHratio")
write.table(randparam, paste("paraout_", cntry, "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)

para <- read.csv(paste("paraout_", cntry, "_", job, ".csv", sep = ""))[-1]

setwd(home)

dt <- (1 / 2)
year1 <- 1900
yearend <- 2050

typen <- 0
year <- c(seq(year1, yearend, 1), rep(0, ((1 / dt) * (yearend - year1 + 1) - (yearend - year1 + 1))))

count <- 1
nn <- 1
nmbr <- count * nn

xout <- mat.or.vec(((yearend - year1 + 1) * (1 / dt) * n_p * (typen + 1)), (173 + 5))

LTBI2015 <- matrix(0, n_p, 2)

for (kkk in 1:n_p)
{
  print(kkk)
  for (i in 1:(length(nm) + 1)) {
    assign(nm[i], as.numeric(para[kkk, i]))
  }
  neta2 <- neta
  Xn <- FitGo(cntry, 1, c(p0, rmort, neta2, rmortTB, CDRscale, CDRscaleE, alpha), c(2, dt, c(0.02, 0.02, 0.8, 0.07)), c(year1, yearend), 0, C)

  xout[((((yearend - year1 + 1) * (1 / dt) * kkk * nmbr) - ((1 / dt) * (yearend - year1 + 1) - 1)):((1 / dt) * (yearend - year1 + 1) * kkk * nmbr)), ] <- cbind(Xn, times, year, (rep(nn, length(times))), (rep(count, length(times))), (rep(kkk, length(times))))

  LTBI2015[kkk, 1] <- LTBIc
  LTBI2015[kkk, 2] <- LTBIa
}

JOB <- rep(job, n_p)
JOB <- as.data.frame(JOB)
xout <- as.data.frame(xout)
xout <- cbind(xout, JOB)
colnames(xout) <- c(colnames(Xn), "timestep", "year", "type", "vxint", "fit", "job")

setwd(home)
source("#likelihood_16_SA.R")

xoutsub <- cbind(
  xout[xout[, "year"] %in% 2000:2050, "TBIalltot"], xout[xout[, "year"] %in% 2000:2050, "TBIall0-14"],
  xout[xout[, "year"] %in% 2000:2050, "TBIall15p"], xout[xout[, "year"] %in% 2000:2050, "TBIHdatot"],
  xout[xout[, "year"] %in% 2000:2050, "TBpcHIV"], xout[xout[, "year"] %in% 2000:2050, "TBMtot"],
  xout[xout[, "year"] %in% 2000:2050, "TBMHdatot"], xout[xout[, "year"] %in% 2000:2050, "TBNalltot"],
  xout[xout[, "year"] %in% 2000:2050, "TBNall0-14"], xout[xout[, "year"] %in% 2000:2050, "TBNall15p"],
  xout[xout[, "year"] %in% 2000:2050, "TBNHdatot"], xout[xout[, "year"] %in% 2000:2050, "TBPalltot"],
  xout[xout[, "year"] %in% 2000:2050, "year"],
  xout[xout[, "year"] %in% 2000:2050, "fit"], xout[xout[, "year"] %in% 2000:2050, "job"]
)

colnames(xoutsub) <- c("TBIalltot", "TBIall014", "TBIall15", "TBIHdatot", "TBpcHIV", "TBMtot", "TBMHdatot", "TBNalltot", "TBNall0-14", "TBNall15p", "TBNHdatot", "TBPalltot", "year", "fit", "job")

startyr <- 1900
xoutsubH <- cbind(
  xout[xout[, "year"] %in% startyr:2050, "HIVPtot"], xout[xout[, "year"] %in% startyr:2050, "HIVP014"],
  xout[xout[, "year"] %in% startyr:2050, "HIVP15plus"], xout[xout[, "year"] %in% startyr:2050, "HIVP1549"],
  xout[xout[, "year"] %in% startyr:2050, "year"],
  xout[xout[, "year"] %in% startyr:2050, "fit"], xout[xout[, "year"] %in% startyr:2050, "job"]
)
colnames(xoutsubH) <- c("HIVPtot", "HIVP014", "HIVP15plus", "HIVP1549", "year", "fit", "job")

setwd(home)
setwd(output)
getwd()

if (C == 0) {
  write.table(xout, paste("xout", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)
}

write.table(xoutsub, paste("xoutsub", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(xoutsubH, paste("xoutsubH", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.csv(LTBI2015, paste("LTBI_2015_", job, ".csv", sep = ""))

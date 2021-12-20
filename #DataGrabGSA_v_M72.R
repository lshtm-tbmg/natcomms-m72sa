C <- 0
input <- here("Data")

dt <- 0.5
Mnage <- 101

setwd(input)
cntry <- "South Africa"

suppressWarnings(library("gdata"))

countries <- as.matrix(read.csv("CountryList_22HBC.csv", check.names = F))

myneta <- suppressWarnings(as.matrix(read.csv("myneta_sa.csv", header = FALSE, check.names = F)[, -1]))

births <- read.csv("BirthR_17.csv", header = TRUE, check.names = F, stringsAsFactors = FALSE)
bb <- births[, cntry] / 1000

mortage <- read.csv("Lifetables_SA.csv", header = TRUE, check.names = F)
mort <- as.matrix(mortage)

mortageH <- read.csv("AIDSmort_age.csv", header = TRUE, check.names = F)
mortH <- as.matrix(mortageH)

pstruc <- suppressWarnings(read.csv("AgeStruc1950.csv", header = TRUE, check.names = F))
ps <- pstruc[1:101, cntry]

suctv <- suppressWarnings(read.csv("SucT_3co.csv", header = TRUE, check.names = F))[1:22, ]
para <- as.matrix(drop.levels(read.csv("para_SA_v.csv", header = TRUE, check.names = F)))[1:29, ]
pararange <- as.matrix(drop.levels(read.csv("pararanges_v3.csv", header = TRUE, check.names = F)))

yHIVi <- as.matrix(read.csv("HIVinc.csv", header = FALSE, check.names = F))
yHIVin <- yHIVi[, 2]
dim(yHIVin) <- c(19, 61)
dim(yHIVin)
yHIVi <- yHIVin[2:18, ]
yHIVi <- as.matrix(yHIVi)
yHIVi <- as.numeric(yHIVi)
dim(yHIVi) <- c(17, 61)

sHIVi <- t(yHIVi * dt / 1000)
sHIViA <- matrix(0, 61, Mnage)
sHIViA[, 1] <- sHIVi[, 1]
sHIViA[, 2:5] <- 0
sHIViA[, 6:10] <- sHIVi[, 2]
sHIViA[, 11:15] <- sHIVi[, 3]
sHIViA[, 16:20] <- sHIVi[, 4]
sHIViA[, 21:25] <- sHIVi[, 5]
sHIViA[, 26:30] <- sHIVi[, 6]
sHIViA[, 31:35] <- sHIVi[, 7]
sHIViA[, 36:40] <- sHIVi[, 8]
sHIViA[, 41:45] <- sHIVi[, 9]
sHIViA[, 46:50] <- sHIVi[, 10]
sHIViA[, 51:55] <- sHIVi[, 11]
sHIViA[, 56:60] <- sHIVi[, 12]
sHIViA[, 61:65] <- sHIVi[, 13]
sHIViA[, 66:70] <- sHIVi[, 14]
sHIViA[, 71:75] <- sHIVi[, 15]
sHIViA[, 76:80] <- sHIVi[, 16]
sHIViA[, 81:Mnage] <- sHIVi[, 17]
hivall <- rbind((rep(0, Mnage)), sHIViA)

H04rate <- hivall[, 1] / dt

paramsG <- suppressWarnings(as.matrix(read.table("CDRfit.csv", sep = ",", header = TRUE)))

paramsG <- as.data.frame(paramsG)

cdrm <- matrix(0, 1, (2050 - 1990 + 1))

K <- paramsG$K
A <- paramsG$A
Qv <- paramsG$Qv
slope <- paramsG$slope
inflect <- paramsG$inflect

cdrm[1, (1:(2050 - 1990 + 1))] <- A + ((K - A) / ((1 + (exp(-slope * (((1990:2050) - inflect)))))^(1 / Qv)))

cdrm <- cdrm / 100

write.table(cdrm, "cdrm.csv", sep = ",")

cdr <- matrix(0, length(cdrm), Mnage)
cdr[, 1:Mnage] <- cdrm
cdrH <- cdr
suctm <- matrix(0, length(countries), (2050 - 1994 + 1))
rownames(suctm) <- countries

suctm[cntry, (1:(2015 - 1994 + 1))] <- t(suctv[(1:22), cntry])
suctm[cntry, (2016 - 1994 + 1):(2050 - 1994 + 1)] <- suctm[cntry, (2015 - 1994 + 1)]

suctt <- suctm[cntry, ]
sucttH <- suctt

art <- as.matrix(drop.levels(read.csv("ARTcovHIVall.csv", header = F, check.names = F)))[, -1]

probt <- 0.1

widthage <- 1
mm <- 101
Mnage <- ceiling(mm / widthage)

for (i in 1:length(para[, 1])) {
  assign(para[i, 1], as.numeric(para[i, 2]))
}

setwd(input)

rm(list = ls())

library(EasyABC)
library(coda)

C <- 1
if (C == 0) {
  home <- "/Users/lsh355020/R/Gates_sa"
}
if (C == 1) {
  home <- "/home/lsh355020/SA_Gates"
}

setwd(home)

print("testABC")

source("MCMCmodel_SA_v7.R")

if (C == 0) {
  input <- "/Users/lsh355020/R/Gates_sa/Data"
}
if (C == 1) {
  input <- "/home/lsh355020/SA_Gates/Data"
}
setwd(input)

para2 <- as.matrix(read.csv("para_SA_v.csv", header = TRUE, check.names = F))[1:29, ]
para2 <- as.numeric(para2[, 2])
assign("uni", para2[14], envir = .GlobalEnv)
assign("ui", para2[20], envir = .GlobalEnv)
assign("uniH", para2[24], envir = .GlobalEnv)
assign("uiH", para2[23], envir = .GlobalEnv)
assign("z", para2[8], envir = .GlobalEnv)

if (C == 0) {
  job <- 1
}
if (C == 1) {
  job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
}

pararange <- as.matrix(read.csv("pararanges_v2.csv", header = TRUE, check.names = F))
print(pararange)
seeds <- as.matrix(read.csv("Afithit_16hits_180322_rmortfix_v7.csv", header = TRUE, check.names = F))
print(seeds)
setwd(home)

pararange <- cbind(as.numeric(as.character(pararange[, 2])), as.numeric(as.character(pararange[, 3])))
print(pararange)
paraprior <- list(
  c("unif", pararange[1, 1], pararange[1, 2]),
  c("unif", pararange[2, 1], pararange[2, 2]),
  c("unif", pararange[3, 1], pararange[3, 2]),
  c("unif", pararange[4, 1], pararange[4, 2]),
  c("unif", pararange[5, 1], pararange[5, 2]),
  c("unif", pararange[6, 1], pararange[6, 2]),
  c("unif", pararange[7, 1], pararange[7, 2]),
  c("unif", pararange[8, 1], pararange[8, 2]),
  c("unif", pararange[9, 1], pararange[9, 2]),
  c("unif", pararange[10, 1], pararange[10, 2]),
  c("unif", pararange[11, 1], pararange[11, 2]),
  c("unif", pararange[12, 1], pararange[12, 2]),
  c("unif", pararange[13, 1], pararange[13, 2]),
  c("unif", pararange[14, 1], pararange[14, 2]),
  c("unif", pararange[15, 1], pararange[15, 2]),
  c("unif", pararange[16, 1], pararange[16, 2]),
  c("unif", pararange[17, 1], pararange[17, 2]),
  c("unif", pararange[18, 1], pararange[18, 2]),
  c("unif", pararange[19, 1], pararange[19, 2]),
  c("unif", pararange[20, 1], pararange[20, 2]),
  c("unif", pararange[21, 1], pararange[21, 2]),
  c("unif", pararange[22, 1], pararange[22, 2]),
  c("unif", pararange[23, 1], pararange[23, 2]),
  c("unif", pararange[24, 1], pararange[24, 2]),
  c("unif", pararange[25, 1], pararange[25, 2]),
  c("unif", pararange[26, 1], pararange[26, 2]),
  c("unif", pararange[27, 1], pararange[27, 2]),
  c("unif", pararange[28, 1], pararange[28, 2]),
  c("unif", pararange[29, 1], pararange[29, 2]),
  c("unif", pararange[30, 1], pararange[30, 2]),
  c("unif", pararange[31, 1], pararange[31, 2]),
  c("unif", pararange[32, 1], pararange[32, 2]),
  c("unif", pararange[33, 1], pararange[33, 2]),
  c("unif", pararange[34, 1], pararange[34, 2]),
  c("unif", pararange[35, 1], pararange[35, 2]),
  c("unif", pararange[36, 1], pararange[36, 2]),
  c("unif", pararange[37, 1], pararange[37, 2]),
  c("unif", 13628.428, 13628.428)
)

sum_stat_fits <- c(1, 0)

propfrac <- 1 / 100
propfracB <- 1 / 100
prop_ran <- c(
  0, ((pararange[2, 2] - pararange[2, 1]) * propfrac), ((pararange[3, 2] - pararange[3, 1]) * propfrac),
  ((pararange[4, 2] - pararange[4, 1]) * propfrac), ((pararange[5, 2] - pararange[5, 1]) * propfrac),
  ((pararange[6, 2] - pararange[6, 1]) * propfrac), ((pararange[7, 2] - pararange[7, 1]) * propfrac),
  ((pararange[8, 2] - pararange[8, 1]) * propfrac), ((pararange[9, 2] - pararange[9, 1]) * propfrac),
  ((pararange[10, 2] - pararange[10, 1]) * propfrac), ((pararange[11, 2] - pararange[11, 1]) * propfrac),
  ((pararange[12, 2] - pararange[12, 1]) * propfrac), ((pararange[13, 2] - pararange[13, 1]) * propfrac),
  ((pararange[14, 2] - pararange[14, 1]) * propfrac), ((pararange[15, 2] - pararange[15, 1]) * propfrac),
  ((pararange[16, 2] - pararange[16, 1]) * propfrac), ((pararange[17, 2] - pararange[17, 1]) * propfrac),
  ((pararange[18, 2] - pararange[18, 1]) * propfrac), ((pararange[19, 2] - pararange[19, 1]) * propfrac),
  ((pararange[20, 2] - pararange[20, 1]) * propfrac), ((pararange[21, 2] - pararange[21, 1]) * propfrac),
  ((pararange[22, 2] - pararange[22, 1]) * propfrac), ((pararange[23, 2] - pararange[23, 1]) * propfrac),
  ((pararange[24, 2] - pararange[24, 1]) * propfrac), ((pararange[25, 2] - pararange[25, 1]) * propfrac),
  ((pararange[26, 2] - pararange[26, 1]) * propfrac), ((pararange[27, 2] - pararange[27, 1]) * propfrac),
  ((pararange[28, 2] - pararange[28, 1]) * propfracB), ((pararange[29, 2] - pararange[29, 1]) * propfrac),
  ((pararange[30, 2] - pararange[30, 1]) * propfrac), ((pararange[31, 2] - pararange[31, 1]) * propfrac),
  ((pararange[32, 2] - pararange[32, 1]) * propfrac), ((pararange[33, 2] - pararange[33, 1]) * propfrac),
  ((pararange[34, 2] - pararange[34, 1]) * propfrac), ((pararange[35, 2] - pararange[35, 1]) * propfrac),
  ((pararange[36, 2] - pararange[36, 1]) * propfrac), ((pararange[37, 2] - pararange[37, 1]) * propfrac), 0
)

print(prop_ran)

print(job)

inipa <- seeds[job, ]

print("inipa")
print(inipa)

para <- c()
hittrack <- c()
hitvar <- c()
accvar <- c()
samppoint <- 20000
hits <- 16
para16 <- c()

ABCmcmc <- ABC_mcmc(
  method = "Marjoram_original", model = MCMCmodel, prior = paraprior, summary_stat_target = sum_stat_fits, n_rec = samppoint, n_between_sampling = 1, proposal_range = prop_ran,
  init_param = inipa, acceptance = TRUE, rejection = TRUE, verbose = TRUE
)

ABCmcmc
hittrack
hitvar
accvar

para16_t <- para16
para16_t <- as.matrix(para16_t)
dim(para16_t) <- c(38, nrow(accvar))
para16_t2 <- t(para16_t)

hitcount <- colSums(hitvar)
hitcount <- hitcount / samppoint * 100
hitcount

count99 <- count(hittrack == 99)
critpc <- count99[2, 2] / samppoint * 100
critpc

countacc <- count(hittrack >= hits & hittrack != 99)
acceptpc <- countacc[2, 2] / samppoint * 100
acceptpc

acccritpc <- cbind(critpc, acceptpc)

colnames(ABCmcmc$param) <- nm
colnames(para16_t2) <- nm

setwd(home)
setwd("output_v35")
getwd()
write.table(para16_t2, paste("ABCacceptfinal_v7b_", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(ABCmcmc$param, paste("ABCmcmc_v7b_", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(acccritpc, paste("acccritpc_v7b_", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(hitcount, paste("hitcount_v7b_", "_", job, ".csv", sep = ""), sep = ",", row.names = FALSE)

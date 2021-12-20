if (grepl(x = Sys.info()["nodename"], pattern = "comp")) {
  SGE_TASK_ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else {
  SGE_TASK_ID <- 1
}

library(data.table)
data.table::setDTthreads(1)

vax_prof_table <- data.table(
  typen = c(1, 1, 1, 1, 3, 3),
  effInf = c(0, 0, 0, 0, 0, 0),
  effDis = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  durs = c(5, 5, 10, 10, 5, 10),
  cover = c(0.8, 0.5, 0.8, 0.5, 0.8, 0.8),
  coverM = c(0, 0, 0, 0, 0, 0),
  vage = c(15, 18, 15, 18, 10, 10),
  fms = c(100, 100, 100, 100, 100, 100)
)

list2env(vax_prof_table[SGE_TASK_ID, ], envir = .GlobalEnv)

Hvax <- c(0.8)
Hsafe <- c(1)

combn <- length(effInf) * length(effDis) * length(durs) * length(Hvax) * length(Hsafe) * length(coverM)

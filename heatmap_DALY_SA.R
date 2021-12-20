library(plyr)
library(grid)
library(ggplot2)

C <- 0
if (C == 0) {
  vaccin <- "/Users/lsh355020/R/M72/M72Ind/Vxoutput_M72_WHO"
}
if (C == 1) {
  vaccin <- "/home/lsh355020/Gates_China/Vxoutput"
}

if (C == 0) {
  vaccout <- "/Users/lsh355020/R/M72/M72Ind/Vxoutput_M72_WHO/1.heatmap"
}
if (C == 1) {
  vaccout <- "/home/lsh355020/Gates_China/Vxoutput/heatmap"
}

if (C == 0) {
  home <- "/Users/lsh355020/R/M72/M72Ind/"
}
if (C == 1) {
  home <- "/home/lsh355020/Gates_China/"
}

setwd(vaccin)

rrun <- 1000

setwd(home)
source("#vxIn_M72.R")
numvx <- combn * typen

setwd(vaccin)
setwd(vaxfol)

RImatrix <- c()
RImatrix35 <- c()
RImatrixM <- c()

for (xx in 1:rrun) {
  print(xx)
  RIcall <- t(read.csv(paste("2050_reduction_incidence_", xx, ".csv", sep = "")))
  RIcall <- cbind(RIcall, rep(xx, nrow(RIcall)))
  RImatrix <- rbind(RImatrix, RIcall[, c(1, 10:17)])

  RIcall35 <- t(read.csv(paste("2035_reduction_incidence_", xx, ".csv", sep = "")))
  RIcall35 <- cbind(RIcall35, rep(xx, nrow(RIcall35)))
  RImatrix35 <- rbind(RImatrix35, RIcall35[, c(1, 10:17)])

  RIcallM <- t(read.csv(paste("2050_reduction_mortality_", xx, ".csv", sep = "")))
  RIcallM <- cbind(RIcallM, rep(xx, nrow(RIcallM)))
  RImatrixM <- rbind(RImatrixM, RIcallM[, c(1, 10:17)])
}

RImatrix <- round(RImatrix, 2)
RImatrix35 <- round(RImatrix35, 2)
RImatrixM <- round(RImatrixM, 2)

colnames(RImatrix) <- c("redu", "type", "VE_I", "VE_D", "dur", "covM", "mmax", "count", "run")
colnames(RImatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covM", "mmax", "count", "run")
colnames(RImatrixM) <- c("redu", "type", "VE_I", "VE_D", "dur", "covM", "mmax", "count", "run")

setwd("1.heatmap")
write.table(RImatrix, "redu_inc_alldata.csv", sep = ",", row.names = F)
write.table(RImatrixM, "redu_mort_alldata.csv", sep = ",", row.names = F)

RImatrix <- as.matrix(RImatrix)
RImatrix35 <- as.matrix(RImatrix35)
RImatrixM <- as.matrix(RImatrixM)

med_RI <- matrix(0, numvx, 3)

for (jj in 1:numvx) {
  seq1 <- seq(1, (rrun * numvx), numvx)
  med_RI[jj, 1] <- median(RImatrix[seq1 + jj - 1, 1])
  med_RI[jj, 2] <- min(RImatrix[seq1 + jj - 1, 1])
  med_RI[jj, 3] <- max(RImatrix[seq1 + jj - 1, 1])
}

colnames(med_RI) <- c("median", "min", "max")
med_RI <- cbind(med_RI, RImatrix[1:numvx, 2:7])

write.table(med_RI, paste("1.Ind_redu_inc_median.csv", sep = ""), sep = ",", row.names = F)

setwd(home)
source("#char4impact.R")

library(ggplot2)
library(RColorBrewer)

my_pal <- brewer.pal(11, "Spectral")
use_pal <- colorRampPalette(my_pal)

groups <- c("PPI", "PRI", "PSI")

for (ig in 1:typen) {
  med_RI <- as.data.frame(med_RI)
  if (typen >= 1) {
    RI50_PPI <- subset(med_RI, type == 1)
  }
  if (typen >= 2) {
    RI50_PRI <- subset(med_RI, type == 2)
  }
  if (typen >= 3) {
    RI50_PSI <- subset(med_RI, type == 3)
  }

  y <- effInf * 100
  x <- effDis * 100

  for (hh in 1:length(durs)) {
    if (ig == 1) assign(paste0("RI50_PPI_D", durs[hh]), subset(RI50_PPI, dur == durs[hh]))
    if (ig == 2) assign(paste0("RI50_PRI_D", durs[hh]), subset(RI50_PRI, dur == durs[hh]))
    if (ig == 3) assign(paste0("RI50_PSI_D", durs[hh]), subset(RI50_PSI, dur == durs[hh]))

    png(paste(vaccout, "/A.plot/", groups[ig], "_", durs[hh], "_heatmap.png", sep = ""), width = 6, height = 6, units = "in", res = 600)

    z <- get(paste0("RI50_", groups[ig], "_D", durs[hh]))
    z <- matrix(z[, 1], nrow = length(effDis), ncol = length(effInf), byrow = TRUE)

    ctlns <- contourLines(x, y, z, levels = c(0))

    print(filled.contour(x, y, z, xlab = "Vaccine efficacy against infection (%)", ylab = "Vaccine efficacy against disease (%)", levels = seq(0, 100, by = 10), color = use_pal, main = paste("Incidence rate reduction in 2050 compared to no new vaccine\nbaseline for a ", groups[ig], " vaccine providing ", durs[hh], " years protection", sep = ""), cex.main = 0.9, key.title = title(main = "IRR\n(%)", cex.main = 0.9)))
    dev.off()
  }
}

assign("RI50_UA", subset(med_RI, (dur == 5 | dur == 10) & (VE_I == 0.2 | VE_I == 1) & (VE_D == 0.2 | VE_D == 1)))

dodge <- position_dodge(width = 0.9)

ggRI_UA <- ggplot(RI50_UA, aes(x = factor(count), y = median, fill = factor(type))) +
  geom_bar(stat = "identity", position = dodge) +
  scale_x_discrete("Vaccine characteristics", labels = c("219" = "VE_I 20%\nVE_D 20%\nDur 5yr", "221" = "VE_I 20%\nVE_D 20%\nDur 10yr", "291" = "VE_I 20%\nVE_D 100%\nDur 5yr", "293" = "VE_I 20%\nVE_D 100%\nDur 10yr", "1011" = "VE_I 100%\nVE_D 20%\nDur 5yr", "1013" = "VE_I 100%\nVE_D 20%\nDur 10yr", "1083" = "VE_I 100%\nVE_D 100%\nDur 5yr", "1085" = "VE_I 100%\nVE_D 100%\nDur 10yr")) +
  geom_errorbar(aes(ymin = RI50_UA$min, ymax = RI50_UA$max, width = 0.25), position = dodge, width = 0.25) +
  theme_classic() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(vjust = 0.5, size = 8), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) +
  scale_fill_discrete(labels = c("PPI", "PRI", "PSI")) +
  labs(x = "Vaccine", y = "TB incidence rate reduction in 2050 compared to no new vaccine scenario (%)")

g3 <- ggplot_gtable(ggplot_build(ggRI_UA))
g3$layout$clip[g3$layout$name == "panel"] <- "off"
grid.draw(g3)

if (C == 0) {
  ggsave("ggRI_UA.pdf", plot = g3)
}

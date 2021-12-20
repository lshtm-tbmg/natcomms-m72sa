library(plyr)
library(grid)
library(ggplot2)

C <- 0
vaccin <- here("Vxoutput_M72")

vaccout <- here("Vxoutput_M72")

home <- here()

setwd(home)
setwd(vaccin)

rrun <- 1000

setwd(home)
source("#vxSA_M72.R")
numvx <- combn * typen

setwd(vaccin)
setwd(vaxfol)

RImatrix <- c()
RIHmatrix <- c()
RINmatrix <- c()
RImatrix35 <- c()
RIHmatrix35 <- c()
RINmatrix35 <- c()
RImatrixM <- c()
RImatrixMH <- c()

for (xx in 1:rrun) {
  print(xx)
  RIcall <- t(read.csv(paste("2050_reduction_incidence_", xx, ".csv", sep = "")))
  RIcall <- cbind(RIcall, rep(xx, nrow(RIcall)))
  RImatrix <- rbind(RImatrix, RIcall[, c(1, 10:19)])
  RIHmatrix <- rbind(RIHmatrix, RIcall[, c(4, 10:19)])
  RINmatrix <- rbind(RINmatrix, RIcall[, c(7, 10:19)])

  RIcall35 <- t(read.csv(paste("2035_reduction_incidence_", xx, ".csv", sep = "")))
  RIcall35 <- cbind(RIcall35, rep(xx, nrow(RIcall35)))
  RImatrix35 <- rbind(RImatrix35, RIcall35[, c(1, 10:19)])
  RIHmatrix35 <- rbind(RIHmatrix35, RIcall35[, c(4, 10:19)])
  RINmatrix35 <- rbind(RINmatrix35, RIcall35[, c(7, 10:19)])

  RIcallM <- t(read.csv(paste("2050_reduction_mortality_", xx, ".csv", sep = "")))
  RIcallM <- cbind(RIcallM, rep(xx, nrow(RIcallM)))
  RImatrixM <- rbind(RImatrixM, RIcallM[, c(1, 3:12)])
  RImatrixMH <- rbind(RImatrixMH, RIcallM[, c(2, 3:12)])

  RIcallactive <- t(read.csv(paste("new_active_", xx, ".csv", sep = "")))
  RIcallactive <- cbind(RIcallactive, rep(xx, nrow(RIcallM)))
  RImatrixM <- rbind(RImatrixM, RIcallM[, c(1, 3:12)])
  RImatrixMH <- rbind(RImatrixMH, RIcallM[, c(2, 3:12)])
}

RImatrix <- round(RImatrix, 2)
RIHmatrix <- round(RIHmatrix, 2)
RINmatrix <- round(RINmatrix, 2)
RImatrix35 <- round(RImatrix35, 2)
RIHmatrix35 <- round(RIHmatrix35, 2)
RINmatrix35 <- round(RINmatrix35, 2)
RImatrixM <- round(RImatrixM, 2)
RImatrixMH <- round(RImatrixMH, 2)

colnames(RImatrix) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RIHmatrix) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RINmatrix) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RImatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RImatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RIHmatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RINmatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RImatrixM) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")
colnames(RImatrixMH) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "count", "run")

colnames(RImatrix) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RIHmatrix) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RINmatrix) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RImatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RImatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RIHmatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RINmatrix35) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RImatrixM) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")
colnames(RImatrixMH) <- c("redu", "type", "VE_I", "VE_D", "dur", "covH", "VEH", "covM", "mmax", "count", "run")

setwd("1.heatmap")

write.table(RImatrix, "redu_inc_alldata.csv", sep = ",", row.names = F)
write.table(RIHmatrix, "redu_incH_alldata.csv", sep = ",", row.names = F)
write.table(RINmatrix, "redu_incN_alldata.csv", sep = ",", row.names = F)
write.table(RImatrix35, "redu_inc35_alldata.csv", sep = ",", row.names = F)
write.table(RIHmatrix35, "redu_inc35H_alldata.csv", sep = ",", row.names = F)
write.table(RINmatrix35, "redu_inc35N_alldata.csv", sep = ",", row.names = F)
write.table(RImatrixM, "redu_mort_alldata.csv", sep = ",", row.names = F)
write.table(RImatrixMH, "redu_mortH_alldata.csv", sep = ",", row.names = F)

RImatrix <- as.matrix(RImatrix)
RImatrix35 <- as.matrix(RImatrix35)
RImatrixM <- as.matrix(RImatrixM)

med_RI <- matrix(0, numvx, 3)
med_RIH <- matrix(0, numvx, 3)
med_RIN <- matrix(0, numvx, 3)
med_RI35 <- matrix(0, numvx, 3)
med_RI35H <- matrix(0, numvx, 3)
med_RI35N <- matrix(0, numvx, 3)
med_RIM <- matrix(0, numvx, 3)
med_RIMH <- matrix(0, numvx, 3)

for (jj in 1:numvx) {
  seq1 <- seq(1, ((rrun) * numvx), numvx)
  med_RI[jj, 1] <- median(RImatrix[seq1 + jj - 1, 1])
  med_RI[jj, 2] <- min(RImatrix[seq1 + jj - 1, 1])
  med_RI[jj, 3] <- max(RImatrix[seq1 + jj - 1, 1])

  med_RIH[jj, 1] <- median(RIHmatrix[seq1 + jj - 1, 1])
  med_RIH[jj, 2] <- min(RIHmatrix[seq1 + jj - 1, 1])
  med_RIH[jj, 3] <- max(RIHmatrix[seq1 + jj - 1, 1])

  med_RIN[jj, 1] <- median(RINmatrix[seq1 + jj - 1, 1])
  med_RIN[jj, 2] <- min(RINmatrix[seq1 + jj - 1, 1])
  med_RIN[jj, 3] <- max(RINmatrix[seq1 + jj - 1, 1])
}

colnames(med_RI) <- c("median", "min", "max")
med_RI <- cbind(med_RI, RImatrix[1:numvx, 2:9])
colnames(med_RIH) <- c("median", "min", "max")
med_RIH <- cbind(med_RIH, RIHmatrix[1:numvx, 2:9])
colnames(med_RIN) <- c("median", "min", "max")
med_RIN <- cbind(med_RIN, RINmatrix[1:numvx, 2:9])

colnames(med_RIM) <- c("median", "min", "max")
med_RIM <- cbind(med_RIM, RImatrixM[1:numvx, 2:9])
colnames(med_RIMH) <- c("median", "min", "max")
med_RIMH <- cbind(med_RIMH, RImatrixMH[1:numvx, 2:9])

write.table(med_RI, paste("1.SA_redu_inc_median", vaxfol, ".csv", sep = ""), sep = ",", row.names = F)
write.table(med_RIH, paste("1.SA_redu_incH_median", vaxfol, ".csv", sep = ""), sep = ",", row.names = F)
write.table(med_RIN, paste("1.SA_redu_incN_median", vaxfol, ".csv", sep = ""), sep = ",", row.names = F)
write.table(med_RI, paste("1.SA_redu_inc_median.csv", sep = ""), sep = ",", row.names = F)
write.table(med_RIH, paste("1.SA_redu_incH_median.csv", sep = ""), sep = ",", row.names = F)
write.table(med_RIN, paste("1.SA_redu_incN_median.csv", sep = ""), sep = ",", row.names = F)
write.table(med_RI35, "redu_inc2035_median.csv", sep = ",", row.names = F)
write.table(med_RI35H, "redu_inc2035H_median.csv", sep = ",", row.names = F)
write.table(med_RI35N, "redu_inc2035N_median.csv", sep = ",", row.names = F)
write.table(med_RIM, "redu_mort_median.csv", sep = ",", row.names = F)
write.table(med_RIMH, "redu_mortH_median.csv", sep = ",", row.names = F)

setwd(vaccout)
med_RI <- read.csv("redu_inc_median.csv", sep = ",")
med_RI35 <- read.csv("redu_inc2035_median.csv", sep = ",")
med_RIM <- read.csv("redu_mort_median.csv", sep = ",")

setwd(home)

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)

my_pal <- brewer.pal(11, "Spectral")
use_pal <- colorRampPalette(my_pal)

groups <- c("PPI", "PRI", "PSI")

for (ig in 1:1) {
  med_RI <- as.data.frame(med_RI)
  med_RIM <- as.data.frame(med_RIM)
  med_RI35 <- as.data.frame(med_RI35)
  if (typen >= 1) {
    RI50_PPI <- subset(med_RI, type == 1)
  }
  if (typen >= 1) {
    RM50_PPI <- subset(med_RIM, type == 1)
  }
  if (typen >= 1) {
    RI35_PPI <- subset(med_RI35, type == 1)
  }
  if (typen >= 2) {
    RI50_PRI <- subset(med_RI, type == 2)
  }
  if (typen >= 2) {
    RM50_PRI <- subset(med_RIM, type == 2)
  }
  if (typen >= 2) {
    RI35_PRI <- subset(med_RI35, type == 2)
  }
  if (typen >= 3) {
    RI50_PSI <- subset(med_RI, type == 3)
  }
  if (typen >= 3) {
    RM50_PSI <- subset(med_RIM, type == 3)
  }
  if (typen >= 3) {
    RI35_PSI <- subset(med_RI35, type == 3)
  }

  png(paste(vaccout, "/A.plot/", groups[ig], "_heatmap.png", sep = ""), width = 20, height = 20, units = "in", res = 600)
  par(mfrow = c(3, 3))

  y <- effInf * 100
  x <- effDis * 100

  for (hh in 1:length(durs)) {
    if (ig == 1) assign(paste0("RI50_PPI_D", durs[hh]), subset(RI50_PPI, dur == durs[hh]))
    if (ig == 2) assign(paste0("RI50_PRI_D", durs[hh]), subset(RI50_PRI, dur == durs[hh]))
    if (ig == 2) assign(paste0("RM50_PRI_D", durs[hh]), subset(RM50_PRI, dur == durs[hh]))
    if (ig == 2) assign(paste0("RI35_PRI_D", durs[hh]), subset(RI35_PRI, dur == durs[hh]))
    if (ig == 3) assign(paste0("RI50_PSI_D", durs[hh]), subset(RI50_PSI, dur == durs[hh]))
    if (ig == 3) assign(paste0("RM50_PSI_D", durs[hh]), subset(RM50_PSI, dur == durs[hh]))
    if (ig == 3) assign(paste0("RI35_PSI_D", durs[hh]), subset(RI35_PSI, dur == durs[hh]))

    png(paste(vaccout, "/A.plot/", groups[ig], "_", durs[hh], "_heatmap.png", sep = ""), width = 6, height = 6, units = "in", res = 600)

    z <- get(paste0("RI50_", groups[ig], "_D", durs[hh]))
    z <- matrix(z[, 1], nrow = length(effDis), ncol = length(effInf), byrow = TRUE)

    ctlns <- contourLines(x, y, z, levels = c(0))

    print(filled.contour(x, y, z, xlab = "Vaccine efficacy against infection (%)", ylab = "Vaccine efficacy against disease (%)", levels = seq(0, 100, by = 10), color = use_pal, main = paste("Incidence rate reduction in 2050 compared to no new vaccine\nbaseline for a ", groups[ig], " vaccine providing ", durs[hh], " years protection", sep = ""), cex.main = 0.9, key.title = title(main = "IRR\n(%)", cex.main = 0.9)))

    dev.off()
  }
}

setwd(home)
source("heatplot_2SA.R")
setwd(home)

dev.off()

setwd(vaccout)

assign("RI50_UA", subset(med_RI, (dur == 5 | dur == 10) & (VE_I == 0.2 | VE_I == 1) & (VE_D == 0.2 | VE_D == 1)))
assign("RI50_UA7", subset(med_RI, (dur == 5 | dur == 10) & (VE_I == 0.2 | VE_I == 0.7) & (VE_D == 0.2 | VE_D == 0.7)))

dodge <- position_dodge(width = 0.9)

png(paste(vaccout, "/A.plot/RI_UR.png", sep = ""), width = 6, height = 6, units = "in", res = 600)

ggRI_UA2 <- ggplot(RI50_UA, aes(x = factor(count), y = median, ymin = min, ymax = max, fill = factor(type))) +
  geom_bar(stat = "identity", position = dodge) +
  scale_fill_manual(values = use_pal) +
  scale_x_discrete("Vaccine characteristics", labels = c("1" = "VE-POI 20%\nVE-POD 20%\nDur 5yr", "2" = "VE-POI 20%\nVE-POD 20%\nDur 10yr", "3" = "VE-POI 20%\nVE-POD 100%\nDur 5yr", "4" = "VE-POI 20%\nVE-POD 100%\nDur 10yr", "5" = "VE-POI 100%\nVE-POD 20%\nDur 5yr", "6" = "VE-POI 100%\nVE-POD 20%\nDur 10yr", "7" = "VE-POI 100%\nVE-POD 100%\nDur 5yr", "8" = "VE-POI 100%\nVE-POD 100%\nDur 10yr")) +
  geom_errorbar(aes(width = 0.25), position = dodge, width = 0.25) +
  facet_wrap(c("VE_I", "VE_D"), nrow = 2, scales = "free_x", as.table = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(vjust = 0.5, size = 8), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_fill_discrete(labels = c("PPI", "PRI", "PSI")) +
  labs(x = "Vaccine", y = "TB incidence rate reduction in 2050 compared to no new vaccine (%)") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(size = 10.2)) +
  theme(axis.title.x = element_text(size = 11))

g4 <- ggplot_gtable(ggplot_build(ggRI_UA2))
g4$layout$clip[g4$layout$name == "panel"] <- "off"
grid.draw(g4)

ggsave("ggRI_UA2.pdf")

png(paste(vaccout, "/A.plot/RI_UR7.png", sep = ""), width = 6, height = 6, units = "in", res = 600)

ggRI_UA27 <- ggplot(RI50_UA7, aes(x = factor(count), y = median, ymin = min, ymax = max, fill = factor(type))) +
  geom_bar(stat = "identity", position = dodge) +
  scale_fill_manual(values = use_pal) +
  scale_x_discrete("Vaccine characteristics", labels = c("219" = "VE-POI 20%\nVE-POD 20%\nDur 5yr", "221" = "VE-POI 20%\nVE-POD 20%\nDur 10yr", "291" = "VE-POI 20%\nVE-POD 70%\nDur 5yr", "293" = "VE-POI 20%\nVE-POD 70%\nDur 10yr", "1011" = "VE-POI 70%\nVE-POD 20%\nDur 5yr", "1013" = "VE-POI 70%\nVE-POD 20%\nDur 10yr", "1083" = "VE-POI 70%\nVE-POD 70%\nDur 5yr", "1085" = "VE-POI 70%\nVE-POD 70%\nDur 10yr")) +
  geom_errorbar(aes(width = 0.25), position = dodge, width = 0.25) +
  facet_wrap(c("VE_I", "VE_D"), nrow = 2, scales = "free_x", as.table = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(vjust = 0.5, size = 8), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_fill_discrete(labels = c("PPI", "PRI", "PSI")) +
  labs(x = "Vaccine", y = "TB incidence rate reduction in 2050 compared to no new vaccine (%)") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(size = 10.2)) +
  theme(axis.title.x = element_text(size = 11))

g5 <- ggplot_gtable(ggplot_build(ggRI_UA27))
g5$layout$clip[g5$layout$name == "panel"] <- "off"
grid.draw(g5)

ggsave("ggRI_UA27.pdf")

library(plyr)
library(grid)
library(ggplot2)
library(leaflet)
library(gtable)
library(grid)
library(gridExtra)

C <- 0

if (C == 0) {
  vacc <- "/Users/lsh355020/R/M72/M72SA/Vxoutput_M72_WHO"
}
if (C == 0) {
  vaccout <- "/Users/lsh355020/R/M72/M72SA/Vxoutput_M72_WHO/1.heatmap"
}

setwd(home)

rrun <- 1000

setwd(vacc)

callactive <- c()
CAV2535all <- c()
CAV2550all <- c()
NVax2550all <- c()

activecols2 <- seq(6, (((typen) * count * 5) + 1), 5)

for (xx in 1:rrun) {
  print(xx)
  CAV <- matrix(0, (yearend - year1 + 1), ((typen) * count))
  callactive <- read.csv(paste("new_active_", xx, ".csv", sep = ""))
  CAV <- callactive[, 1] - callactive[, (activecols2)]

  write.table(CAV, paste("annual_cases_averted_", xx, ".csv", sep = ""), sep = ",", row.names = FALSE)

  CAV2535 <- sum(callactive[126:136, 1]) - colSums(callactive[126:136, (activecols2)])

  CAV2550 <- colSums(CAV)

  write.table(CAV2535, paste("CA_2535_", xx, ".csv", sep = ""), sep = ",", row.names = FALSE)
  write.table(CAV2550, paste("CA_2550_", xx, ".csv", sep = ""), sep = ",", row.names = FALSE)

  NVax2550 <- read.csv(paste("Nvax_2550_", xx, ".csv", sep = ""))

  CAV2535all <- rbind(CAV2535all, t(as.matrix(CAV2535)))
  CAV2550all <- rbind(CAV2550all, t(as.matrix(CAV2550)))
  NVax2550all <- rbind(NVax2550all, t(as.matrix(NVax2550)))

  NNVc <- NVax2550all / CAV2550all
}

write.table(NNVc, paste("NNVc_2550_", ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(CAV2535all, paste("CAV2535all", ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(CAV2550all, paste("CAV2550all", ".csv", sep = ""), sep = ",", row.names = FALSE)
write.table(NVax2550all, paste("NVax2550all", ".csv", sep = ""), sep = ",", row.names = FALSE)

vacc <- rep(seq(1, 12, 1), 8)

Dur <- rep(c("10", "20"), 48)
VE <- rep(c("40", "60", "80"), each = 2, 16)
Cov <- rep(c("30", "70"), each = 6, 8)
Type <- rep(c("PRI", "PSI-LR", "PPI", "PRI", "PSI-LR", "PPI", "PSI-L", "PSI-L"), each = 12)
Age <- c(rep(c("Older Adult", "Adolescent and Adult"), each = 36), rep(c("Older Adult", "Adolescent and Adult"), each = 12))

Type <- rep(c("Pre & Post Infection", "Post Infection"), each = 12)
VE <- rep(c("50", "80"), each = 6, 2)

setwd(vacc)
period_CA50_range <- as.data.frame(read.csv("period_CA_range_2050.csv"))
period_DA50_range <- as.data.frame(read.csv("period_DA_range_2050.csv"))
period_CA35_range <- as.data.frame(read.csv("period_CA_range_2035.csv"))
period_NVax50_range <- as.data.frame(read.csv("period_NVax50_range.csv"))
period_NNVCA50_range <- as.data.frame(read.csv("period_NNVCA_range_2050.csv"))
period_NNVDA50_range <- as.data.frame(read.csv("period_NNVDA_range_2050.csv"))

dodge <- position_dodge(width = 0.9)

med_RI <- as.data.frame(med_RI)
assign("RI_PPIPSI", subset(med_RI, (type == 1 | type == 3)))

med_RIb <- cbind(RI_PPIPSI, Type, VE)
med_RI2 <- med_RIb[order(med_RIb$Type, med_RIb$VE, med_RIb$covM, med_RIb$mmax), ]

medCA <- cbind(period_CA50_range[c(1:12, 25:36), 1:3], med_RIb[, 5:13])
med_CA2 <- medCA[order(medCA$Type, medCA$VE, medCA$covM, medCA$mmax), ]

medDA <- cbind(period_DA50_range[c(1:12, 25:36), 1:3], med_RIb[, 5:13])
med_DA2 <- medDA[order(medDA$Type, medDA$VE, medDA$covM, medDA$mmax), ]

medCA35 <- cbind(period_CA35_range[c(1:12, 25:36), 1:3], med_RIb[, 5:13])
med_CA235 <- medCA35[order(medCA35$Type, medCA35$VE, medCA35$covM, medCA35$mmax), ]

medNvax <- cbind(period_NVax50_range[c(1:12, 25:36), 1:3], med_RIb[, 5:13])
med_Nvax2 <- medNvax[order(medNvax$Type, medNvax$VE, medNvax$covM, medNvax$mmax), ]

medNNVCA50 <- cbind(period_NNVCA50_range[c(1:12, 25:36), 1:3], med_RIb[, 5:13])
medNNVCA502 <- medNNVCA50[order(medNNVCA50$Type, medNNVCA50$VE, medNNVCA50$covM, medNNVCA50$mmax), ]

medNNVDA50 <- cbind(period_NNVDA50_range[c(1:12, 25:36), 1:3], med_RIb[, 5:13])
medNNVDA502 <- medNNVDA50[order(medNNVDA50$Type, medNNVDA50$VE, medNNVDA50$covM, medNNVDA50$mmax), ]

cov_VE2 <- paste("Mass\n ", (med_RI2$covM) * 100, "% cov\n17-", med_RI2$mmax, "yrs", sep = "")

med_RI3 <- cbind(med_RI2, cov_VE2)
med_CA3 <- cbind(med_CA2, cov_VE2)
med_DA3 <- cbind(med_DA2, cov_VE2)
med_CA335 <- cbind(med_CA235, cov_VE2)
med_Nvax3 <- cbind(med_Nvax2, cov_VE2)
medNNVCA503 <- cbind(medNNVCA502, cov_VE2)
medNNVDA503 <- cbind(medNNVDA502, cov_VE2)

med_RI3$cov_VE2 <- as.character(med_RI3$cov_VE2)
med_CA3$cov_VE2 <- as.character(med_CA3$cov_VE2)
med_DA3$cov_VE2 <- as.character(med_DA3$cov_VE2)
med_CA335$cov_VE2 <- as.character(med_CA335$cov_VE2)
medNNVCA503$cov_VE2 <- as.character(medNNVCA503$cov_VE2)
medNNVDA503$cov_VE2 <- as.character(medNNVDA503$cov_VE2)

med_RI3$cov_VE2 <- factor(med_RI3$cov_VE2, levels = unique(med_RI3$cov_VE2))
med_CA3$cov_VE2 <- factor(med_CA3$cov_VE2, levels = unique(med_CA3$cov_VE2))
med_DA3$cov_VE2 <- factor(med_DA3$cov_VE2, levels = unique(med_DA3$cov_VE2))
med_CA335$cov_VE2 <- factor(med_CA335$cov_VE2, levels = unique(med_CA335$cov_VE2))
med_Nvax3$cov_VE2 <- factor(med_Nvax3$cov_VE2, levels = unique(med_Nvax3$cov_VE2))
medNNVCA503$cov_VE2 <- factor(medNNVCA503$cov_VE2, levels = unique(medNNVCA503$cov_VE2))
medNNVDA503$cov_VE2 <- factor(medNNVDA503$cov_VE2, levels = unique(medNNVDA503$cov_VE2))

setwd(vaccout)
write.table(med_RI3, "ReductInc2050_median.csv", sep = ",", row.names = F)
write.table(med_RIM3, "ReductMort2050_median.csv", sep = ",", row.names = F)

ggRIall <- ggplot(med_RI3, aes(x = cov_VE2, y = median, group = factor(VE), fill = VE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = med_RI3$min, ymax = med_RI3$max, width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics and implementation", y = "TB incidence rate reduction in 2050 \n(vaccine vs. no vaccine, %)") +
  ggtitle("SA") +
  facet_wrap("Type", nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 16, hjust = 0.4, face = "bold")) +
  coord_cartesian(ylim = c(-0.5, 100)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggRIall

ggsave("ggRIall_SA.pdf", plot = ggRIall, width = 9, height = 5, dpi = 120)

ggCAall <- ggplot(med_CA3, aes(x = cov_VE2, y = V1 / 1000, group = factor(VE), fill = VE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = med_CA3$V2 / 1000, ymax = med_CA3$V3 / 1000, width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics and implementation", y = "Cumulative TB cases averted\nby vaccination 2028-2050 (million)") +
  ggtitle("SA") +
  facet_wrap("Type", nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 16, hjust = 0.4, face = "bold")) +
  coord_cartesian(ylim = c(-0.05, 4)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggCAall

ggsave("ggCA_2850_SA.pdf", plot = ggCAall, width = 9, height = 5, dpi = 120)

ggDAall <- ggplot(med_DA3, aes(x = cov_VE2, y = (V1 / 1000), group = factor(VE), fill = VE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (med_DA3$V2 / 1000), ymax = (med_DA3$V3 / 1000), width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics and implementation", y = "Cumulative TB deaths averted\nby vaccination 2028-2050 (million)") +
  ggtitle("SA") +
  facet_wrap("Type", nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 16, hjust = 0.4, face = "bold")) +
  coord_cartesian(ylim = c(-0.05, 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggDAall

ggsave("ggDA_2850_SA.pdf", plot = ggDAall, width = 9, height = 5, dpi = 120)

ggCA35 <- ggplot(med_CA335, aes(x = cov_VE2, y = (V1 / 1000), group = factor(VE), fill = VE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (med_CA335$V2 / 1000), ymax = (med_CA335$V3 / 1000), width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics and implementation", y = "Cumulative TB cases averted\nby vaccination 2028-2035 (million)") +
  ggtitle("SA") +
  facet_wrap("Type", nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 16, hjust = 0.4, face = "bold")) +
  coord_cartesian(ylim = c(-0.05, 4)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggCA35

ggsave("ggCA_2835_SA.pdf", plot = ggCA35, width = 9, height = 5, dpi = 120)

ggNvax <- ggplot(med_Nvax3, aes(x = cov_VE2, y = (V1 / 1000), group = factor(VE), fill = VE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (med_Nvax3$V2 / 1000), ymax = (med_Nvax3$V3 / 1000), width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics and implementation", y = "Cumulative number vaccinated 2028-2050 (million)") +
  ggtitle("SA") +
  facet_wrap("Type", nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 16, hjust = 0.4, face = "bold")) +
  coord_cartesian(ylim = c(-5, 150)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggNvax

ggsave("ggNVax_2850_SA.pdf", plot = ggNvax, width = 9, height = 5, dpi = 120)

ggNNVc <- ggplot(medNNVCA503, aes(x = cov_VE2, y = V1, group = factor(VE), fill = VE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (medNNVCA503$V2), ymax = (medNNVCA503$V3), width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics and implementation", y = "Cumulative number needed to vaccinate\nper case averted 2028-2050") +
  ggtitle("SA") +
  facet_wrap("Type", nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 16, hjust = 0.4, face = "bold")) +
  coord_cartesian(ylim = c(-5, 700)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggNNVc

ggsave("ggNNVc_2850_SA.pdf", plot = ggNNVc, width = 9, height = 5, dpi = 120)

ggNNVd <- ggplot(medNNVDA503, aes(x = cov_VE2, y = V1, group = factor(VE), fill = VE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (medNNVDA503$V2), ymax = (medNNVDA503$V3), width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics and implementation", y = "Cumulative number needed to vaccinate\nper death averted 2028-2050") +
  ggtitle("SA") +
  facet_wrap("Type", nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 16, hjust = 0.4, face = "bold")) +
  coord_cartesian(ylim = c(-5, 3500)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggNNVd

ggsave("ggNNVd_2850_SA.pdf", plot = ggNNVd, width = 9, height = 5, dpi = 120)

PPIplot <- med_RI3[med_RI3$type == "1", ]
ggRIppi <- ggplot(PPIplot, aes(x = cov_VE2, y = V1, fill = Dur)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = PPIplot$V10, ymax = PPIplot$V19, width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics", y = "TB incidence reduction (%) in 2050 - PPI vaccine compared to no vaccine (%)") +
  facet_wrap(~Age) +
  theme(axis.text = element_text(size = 8)) +
  coord_cartesian(ylim = c(0, 40)) +
  scale_y_continuous(expand = c(0, 0))
ggRIppi
ggsave("ggRIppi.pdf", plot = ggRIppi)

PSIplot <- med_RI3[med_RI3$Type == "PSI-LR", ]
ggRIpsi <- ggplot(PSIplot, aes(x = cov_VE2, y = V1, fill = Dur)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = PSIplot$V10, ymax = PSIplot$V19, width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics", y = "TB incidence reduction (%) in 2050 - PSI-LR vaccine compared to no vaccine (%)") +
  facet_wrap(~Age) +
  theme(axis.text = element_text(size = 8)) +
  coord_cartesian(ylim = c(0, 40)) +
  scale_y_continuous(expand = c(0, 0))

ggRIpsi
ggsave("ggRIpsi.pdf", plot = ggRIpsi)
gg7 <- ggplotGrob(ggRIall)

gg7 <- gtable_add_rows(gg7, gg7$height[5], pos = 2)
gg7 <- gtable_add_rows(gg7, unit(3 / 10, "line"), 3)
gg7 <- gtable_add_rows(gg7, gg7$height[3], pos = 8)
gg7 <- gtable_add_rows(gg7, unit(3 / 10, "line"), 9)

gg8 <- gtable_add_grob(gg7,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Pre- and Post-Infection (PPI)", gp = gpar(cex = 1, fontface = "bold", col = "black", just = "top"))
  ),
  t = 3, l = 4, b = 3, r = 7, name = c("a", "b")
)

gg9 <- gtable_add_grob(gg8,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Pre-Infection (PRI)", gp = gpar(cex = 1, fontface = "bold", col = "black"))
  ),
  t = 3, l = 10, b = 3, r = 13, name = c("a", "b")
)

gg10 <- gtable_add_grob(gg9,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Post-Infection in Latency and Recovered (PSI-LR)", gp = gpar(cex = 1, fontface = "bold", col = "black"))
  ),
  t = 9, l = 10, b = 9, r = 13, name = c("a", "b")
)

gg11 <- gtable_add_grob(gg10,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Post-Infection in Latency (PSI-L)", gp = gpar(cex = 1, fontface = "bold", col = "black"))
  ),
  t = 9, l = 4, b = 9, r = 7, name = c("a", "b")
)

gg12 <- gtable_add_grob(gg11,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 4, b = 5, r = 4, name = c("a", "b")
)

gg13 <- gtable_add_grob(gg12,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 7, b = 5, r = 7, name = c("a", "b")
)

gg14 <- gtable_add_grob(gg13,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 10, b = 5, r = 10, name = c("a", "b")
)

gg15 <- gtable_add_grob(gg14,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 13, b = 5, r = 13, name = c("a", "b")
)

gg16 <- gtable_add_grob(gg15,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 4, b = 11, r = 4, name = c("a", "b")
)

gg17 <- gtable_add_grob(gg16,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 7, b = 11, r = 7, name = c("a", "b")
)

gg18 <- gtable_add_grob(gg17,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 10, b = 11, r = 10, name = c("a", "b")
)

gg19 <- gtable_add_grob(gg18,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 13, b = 11, r = 13, name = c("a", "b")
)

grid.newpage()
grid.draw(gg19)

ggsave("gg19.pdf", plot = gg19, width = 16, height = 9, dpi = 120)

med_RI2 <- as.data.frame(med_RI2)
med_RI_vacctype <- med_RI2[med_RI2$vacc == 9, ]

med_RI_vacctype <- med_RI_vacctype[order(Type, Age), ]

ggRI_vt <- ggplot(med_RI_vacctype, aes(x = Type, y = V1, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = med_RI_vacctype$V10, ymax = med_RI_vacctype$V19, width = 0.25), position = position_dodge(width = 0.9)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  labs(x = "Vaccine Mechanism", y = "TB incidence reduction in 2050 (vaccine vs. no vaccine, %)")

ggsave("ggRI_vt.pdf", plot = ggRI_vt)

seq1 <- seq(1, (rrun * 9), 9)

for (kk in 1:12) {
  pcRLvL[, kk] <- (RImatrix[seq1, (kk + 12)] - RImatrix[seq1, (kk + 72)]) / RImatrix[seq1, (kk + 12)] * 100
}
rangeRLvL <- matrix(0, 2, 12)

for (ii in 1:12) {
  rangeRLvL[1, ii] <- min(pcRLvL[, ii])
  rangeRLvL[2, ii] <- max(pcRLvL[, ii])
}

ggRIMall <- ggplot(med_RIM3, aes(x = cov_VE2, y = V1, fill = Dur)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = med_RIM3$V10, ymax = med_RIM3$V19, width = 0.25), position = position_dodge(width = 0.9)) +
  labs(x = "Vaccine characteristics", y = "TB mortality reduction in 2050 (vaccine vs. no vaccine, %)") +
  facet_wrap(c("Type", "Age"), nrow = 2, scales = "free_x") +
  theme(axis.text = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12)) +
  theme(
    panel.margin = unit(c(0.5, 1, 0.5, 1), "lines"),
    axis.line = element_line(colour = "black", size = 1)
  ) +
  coord_cartesian(ylim = c(-0.5, 40)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()

ggRIMall

ggM7 <- ggplotGrob(ggRIMall)

ggM7 <- gtable_add_rows(ggM7, ggM7$height[5], pos = 2)
ggM7 <- gtable_add_rows(ggM7, unit(3 / 10, "line"), 3)
ggM7 <- gtable_add_rows(ggM7, ggM7$height[3], pos = 8)
ggM7 <- gtable_add_rows(ggM7, unit(3 / 10, "line"), 9)

ggM8 <- gtable_add_grob(ggM7,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Pre- and Post-Infection (PPI)", gp = gpar(cex = 1, fontface = "bold", col = "black", just = "top"))
  ),
  t = 3, l = 4, b = 3, r = 7, name = c("a", "b")
)

ggM9 <- gtable_add_grob(ggM8,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Pre-Infection (PRI)", gp = gpar(cex = 1, fontface = "bold", col = "black"))
  ),
  t = 3, l = 10, b = 3, r = 13, name = c("a", "b")
)

ggM10 <- gtable_add_grob(ggM9,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Post-Infection in Latency and Recovered (PSI-LR)", gp = gpar(cex = 1, fontface = "bold", col = "black"))
  ),
  t = 9, l = 10, b = 9, r = 13, name = c("a", "b")
)

ggM11 <- gtable_add_grob(ggM10,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Post-Infection in Latency (PSI-L)", gp = gpar(cex = 1, fontface = "bold", col = "black"))
  ),
  t = 9, l = 4, b = 9, r = 7, name = c("a", "b")
)

ggM12 <- gtable_add_grob(ggM11,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 4, b = 5, r = 4, name = c("a", "b")
)

ggM13 <- gtable_add_grob(ggM12,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 7, b = 5, r = 7, name = c("a", "b")
)

ggM14 <- gtable_add_grob(ggM13,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 10, b = 5, r = 10, name = c("a", "b")
)

ggM15 <- gtable_add_grob(ggM14,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 5, l = 13, b = 5, r = 13, name = c("a", "b")
)

ggM16 <- gtable_add_grob(ggM15,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 4, b = 11, r = 4, name = c("a", "b")
)

ggM17 <- gtable_add_grob(ggM16,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 7, b = 11, r = 7, name = c("a", "b")
)

ggM18 <- gtable_add_grob(ggM17,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Adolescent", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 10, b = 11, r = 10, name = c("a", "b")
)

ggM19 <- gtable_add_grob(ggM18,
  list(
    rectGrob(gp = gpar(col = NA, fill = "gray85", size = .5)),
    textGrob("Older Adult", gp = gpar(cex = 1, fontface = "plain", col = "black", vjust = 0.1))
  ),
  t = 11, l = 13, b = 11, r = 13, name = c("a", "b")
)

grid.newpage()
grid.draw(ggM19)

ggsave("ggM19.pdf", plot = ggM19, width = 16, height = 9, dpi = 120)

ggRI <- ggplot(med_RI2, aes(x = vacc, y = V1)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = med_RI2$V10, ymax = med_RI2$V19, width = 0.25)) +
  theme_classic() +
  annotate("text", x = seq(0.7, 95, length.out = 96), y = -1, label = Dur) +
  annotate("text", x = seq(0.7, 95, length.out = 48), y = -2, label = rep(c("40", "60", "80"), 16)) +
  annotate("text", x = seq(0.7, 95, length.out = 16), y = -3, label = rep(c("30", "70"), 8)) +
  annotate("text", x = seq(0.7, 95, length.out = 8), y = -4, label = med_RI2$Type[c(1, 13, 25, 37, 49, 61, 73, 85)]) +
  annotate("text", x = seq(0.7, 95, length.out = 2), y = -5, label = med_RI2$Age[c(1, 55)]) +
  theme(
    plot.margin = unit(c(1, 1, 4, 1), "lines"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Vaccine", y = "Reduction in TB incidence in 2050 compared to no vaccine scenario (%)") +
  facet_grid(facets, margins = FALSE, scales = "fixed", space = "fixed", shrink = TRUE, labeller = "label_value", as.table = TRUE, drop = TRUE)

g2 <- ggplot_gtable(ggplot_build(ggRI))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.draw(g2)

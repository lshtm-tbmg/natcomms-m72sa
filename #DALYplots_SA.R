vaccout <- here("Vxoutput_M72_refit/1.heatmap")
vacc <- here("Vxoutput_M72_refit")

setwd(vacc)
setwd(vaxfol)
setwd("1.heatmap")
inc1550 <- t(inc1550)
incnum <- t(incnum)
incyrs <- nrow(inc1550)
numvx <- combn * typen
numvx <- 1

med_inc1550 <- matrix(0, incyrs, (3 * (numvx + 1)))
med_incnum <- matrix(0, incyrs, (3 * (numvx + 1)))

for (ii in 1:(numvx + 1)) {
  for (jj in 1:incyrs) {
    seq1 <- seq(ii, (rrun * numvx), (numvx + 1))
    med_inc1550[jj, (3 * ii - 2)] <- median(inc1550[jj, seq1])
    med_inc1550[jj, (3 * ii - 1)] <- min(inc1550[jj, seq1])
    med_inc1550[jj, (3 * ii)] <- max(inc1550[jj, seq1])

    med_incnum[jj, (3 * ii - 2)] <- median(incnum[jj, seq1])
    med_incnum[jj, (3 * ii - 1)] <- min(incnum[jj, seq1])
    med_incnum[jj, (3 * ii)] <- max(incnum[jj, seq1])
  }
}

Year <- seq(2015, 2050, 1)
med_inc1550 <- as.data.frame(cbind(Year, med_inc1550))
med_incnum <- as.data.frame(cbind(Year, med_incnum))

ggincrate <- ggplot(med_inc1550, aes(x = Year, y = med_inc1550[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_inc1550[, 3], ymax = med_inc1550[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_inc1550[, 2]), colour = "black", size = 1) +
  geom_ribbon(aes(ymin = med_inc1550[, 6], ymax = med_inc1550[, 7], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_inc1550[, 5]), colour = "red", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1020), xlim = c(2015, 2052)) +
  labs(x = "Year", y = "All-age incidence rate in South Africa (/100,000 pop)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggincrate

ggsave(paste("ggincrate_SA_", vaxfol, ".png", sep = ""), plot = ggincrate, width = 9, height = 8, dpi = 120)

med_incnum2 <- med_incnum

med_incnum <- cbind(med_incnum[, 1], med_incnum[, 2:7] / 1000)

ggincnum <- ggplot(med_incnum, aes(x = Year, y = med_incnum[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_incnum[, 3], ymax = med_incnum[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_incnum[, 2]), colour = "black", size = 1) +
  geom_ribbon(aes(ymin = med_incnum[, 6], ymax = med_incnum[, 7], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_incnum[, 5]), colour = "red", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 4190), xlim = c(2014, 2052)) +
  labs(x = "Year", y = "Annual incidence number in South Africa (Thousands)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggincnum

ggsave(paste("ggincnum_SA_", vaxfol, ".png", sep = ""), plot = ggincnum, width = 9, height = 8, dpi = 120)

ggincnum <- ggplot(med_incnum, aes(x = Year, y = med_incnum[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_incnum[, 3], ymax = med_incnum[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_incnum[, 2]), colour = "black", size = 1) +
  geom_ribbon(aes(ymin = med_incnum[, 6], ymax = med_incnum[, 7], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_incnum[, 5]), colour = "red", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 510), xlim = c(2014, 2052)) +
  labs(x = "Year", y = "Annual incidence number in South Africa (Thousands)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 22),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggincnum

ggsave(paste("ggincnumsml_SA_", vaxfol, ".png", sep = ""), plot = ggincnum, width = 9, height = 8, dpi = 120)

med_annual2 <- as.data.frame(med_annual, colnames = TRUE)
med_annual_numvxM <- as.data.frame(med_annual2[, 1:4])
med_annual_numvxM2 <- cbind(med_annual_numvxM[, 1], med_annual_numvxM[, 2:4] / 1000000)

ggvaxnum <- ggplot(med_annual_numvxM2, aes(x = Year[11:36], y = med_annual_numvxM2[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxM2[, 3], ymax = med_annual_numvxM2[, 4], width = 3), fill = "grey70") +
  geom_line(aes(y = med_annual_numvxM2[, 2]), colour = "black", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.5), xlim = c(2024, 2052)) +
  labs(x = "Year", y = paste("Number of ", vage, " year olds vaccinated \n per year in South Africa (million)")) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggvaxnum

ggsave(paste("ggvaxnum_SA_", vaxfol, ".png", sep = ""), plot = ggvaxnum, width = 9, height = 8, dpi = 120)

ggdosenum <- ggplot(med_annual_numvxM2, aes(x = Year[11:36], y = (med_annual_numvxM2[, 2] * doses))) +
  theme_classic() +
  geom_ribbon(aes(ymin = (med_annual_numvxM2[, 3] * doses), ymax = (med_annual_numvxM2[, 4] * doses), width = 3), fill = "grey70") +
  geom_line(aes(y = (med_annual_numvxM2[, 2] * doses)), colour = "black", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 3), xlim = c(2024, 2051)) +
  labs(x = "Year", y = "Number of doses delivered \n per year in South Africa (million)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggdosenum

ggsave(paste("ggdosenumsml_SA_", vaxfol, ".png", sep = ""), plot = ggdosenum, width = 9, height = 8, dpi = 120)

ggdosenum <- ggplot(med_annual_numvxM2, aes(x = Year[11:36], y = (med_annual_numvxM2[, 2] * doses))) +
  theme_classic() +
  geom_ribbon(aes(ymin = (med_annual_numvxM2[, 3] * doses), ymax = (med_annual_numvxM2[, 4] * doses), width = 3), fill = "grey70") +
  geom_line(aes(y = (med_annual_numvxM2[, 2] * doses)), colour = "black", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 49), xlim = c(2024, 2051)) +
  labs(x = "Year", y = "Number of doses delivered \n per year in South Africa (million)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggdosenum

ggsave(paste("ggdosenum_SA_", vaxfol, ".png", sep = ""), plot = ggdosenum, width = 9, height = 8, dpi = 120)

med_annual_numvxMdalB <- as.data.frame(med_annual2[, c(1, 77:82)])
med_annual_numvxMdal <- as.data.frame(cbind(med_annual2[, 1], med_annual2[, c(77:82)] / 1000000))

ggDALYnum <- ggplot(med_annual_numvxMdal, aes(x = Year[11:36], y = med_annual_numvxMdal[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMdal[, 3], ymax = med_annual_numvxMdal[, 4], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMdal[, 2]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMdal[, 6], ymax = med_annual_numvxMdal[, 7], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMdal[, 5]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 2.5), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual DALYs averted in South Africa (millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggDALYnum

ggsave(paste("ggDALYnum2_SA_", vaxfol, ".png", sep = ""), plot = ggDALYnum, width = 9, height = 8, dpi = 120)

ggDALYnum <- ggplot(med_annual_numvxMdal, aes(x = Year[11:36], y = med_annual_numvxMdal[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMdal[, 3], ymax = med_annual_numvxMdal[, 4], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMdal[, 2]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMdal[, 6], ymax = med_annual_numvxMdal[, 7], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMdal[, 5]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 0.32), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual DALYs averted in South Africa (millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggDALYnum

ggsave(paste("ggDALYnum2sml_SA_", vaxfol, ".png", sep = ""), plot = ggDALYnum, width = 9, height = 8, dpi = 120)

med_annual_numvxMcost <- as.data.frame(cbind(med_annual2[, 1], (med_annual2[, c(8:10, 14:16, 29:76)] / 1000000)))

ggcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 3], ymax = med_annual_numvxMcost[, 4], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 2]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 6], ymax = med_annual_numvxMcost[, 7], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 5]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 111), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual expenditure on vaccines \n in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggcostnum

ggsave(paste("ggcostnum_SA_", vaxfol, ".png", sep = ""), plot = ggcostnum, width = 9, height = 8, dpi = 120)

ggcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 2])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 3], ymax = med_annual_numvxMcost[, 4], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 2]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 6], ymax = med_annual_numvxMcost[, 7], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 5]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 5.5), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual expenditure on vaccines \n in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggcostnum

ggsave(paste("ggcostnumsml_SA_", vaxfol, ".png", sep = ""), plot = ggcostnum, width = 9, height = 8, dpi = 120)

ggtcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 8])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 9], ymax = med_annual_numvxMcost[, 10], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 8]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 12], ymax = med_annual_numvxMcost[, 13], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 11]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(4, -250), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental TB treatment expenditure \n in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggtcostnum

ggsave(paste("ggtreatcostnum_SA_", vaxfol, ".png", sep = ""), plot = ggtcostnum, width = 9, height = 8, dpi = 120)

ggtcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 8])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 9], ymax = med_annual_numvxMcost[, 10], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 8]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 12], ymax = med_annual_numvxMcost[, 13], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 11]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(4, -15), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental TB treatment expenditure \n in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggtcostnum

ggsave(paste("ggtreatcostnumsml_SA_", vaxfol, ".png", sep = ""), plot = ggtcostnum, width = 9, height = 8, dpi = 120)

ggdcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 14])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 15], ymax = med_annual_numvxMcost[, 16], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 14]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 18], ymax = med_annual_numvxMcost[, 19], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 17]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(4, -250), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental TB diagnostics expenditure \n in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggdcostnum

ggsave(paste("ggdiagcostnum_SA_", vaxfol, ".png", sep = ""), plot = ggdcostnum, width = 9, height = 8, dpi = 120)

ggdcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 14])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 15], ymax = med_annual_numvxMcost[, 16], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 14]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 18], ymax = med_annual_numvxMcost[, 19], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 17]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(4, -15), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental TB diagnostics expenditure \n in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggdcostnum

ggsave(paste("ggdiagcostnumsml_SA_", vaxfol, ".png", sep = ""), plot = ggdcostnum, width = 9, height = 8, dpi = 120)

ggTBcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 20])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 21], ymax = med_annual_numvxMcost[, 22], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 20]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 24], ymax = med_annual_numvxMcost[, 25], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 23]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(4, -250), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental TB programme expenditure \n (treatment & diagnostics) in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBcostnum

ggsave(paste("ggTBprogcostnum_SA_", vaxfol, ".png", sep = ""), plot = ggTBcostnum, width = 9, height = 8, dpi = 120)

ggTBcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 20])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 21], ymax = med_annual_numvxMcost[, 22], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 20]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 24], ymax = med_annual_numvxMcost[, 25], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 23]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(4, -15), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental TB programme expenditure \n (treatment & diagnostics) in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggTBcostnum

ggsave(paste("ggTBprogcostnumsml_SA_", vaxfol, ".png", sep = ""), plot = ggTBcostnum, width = 9, height = 8, dpi = 120)

ggARTcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 26])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 27], ymax = med_annual_numvxMcost[, 28], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 26]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 30], ymax = med_annual_numvxMcost[, 31], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 29]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(111, 0), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental ART \n expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggARTcostnum

ggsave(paste("ggARTcost_SA_", vaxfol, ".png", sep = ""), plot = ggARTcostnum, width = 9, height = 8, dpi = 120)

ggARTcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 26])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 27], ymax = med_annual_numvxMcost[, 28], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 26]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 30], ymax = med_annual_numvxMcost[, 31], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 29]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(5.5, 0), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental ART \n expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggARTcostnum

ggsave(paste("ggARTcostsml_SA_", vaxfol, ".png", sep = ""), plot = ggARTcostnum, width = 9, height = 8, dpi = 120)

ggallcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 32])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 33], ymax = med_annual_numvxMcost[, 34], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 32]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 36], ymax = med_annual_numvxMcost[, 37], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 35]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(100, -160), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (vaccine, TB treatment \n & diagnostics, ART) expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallcostnum

ggsave(paste("ggALLcostnum_SA_", vaxfol, ".png", sep = ""), plot = ggallcostnum, width = 9, height = 8, dpi = 120)

ggallcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 32])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 33], ymax = med_annual_numvxMcost[, 34], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 32]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 36], ymax = med_annual_numvxMcost[, 37], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 35]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(5.5, -15), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (vaccine, TB treatment \n & diagnostics, ART) expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallcostnum

ggsave(paste("ggALLcostnumsml_SA_", vaxfol, ".png", sep = ""), plot = ggallcostnum, width = 9, height = 8, dpi = 120)

ggallNDcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 38])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 39], ymax = med_annual_numvxMcost[, 40], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 38]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 42], ymax = med_annual_numvxMcost[, 43], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 41]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(100, -150), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (vaccine, TB treatment, \n ART) expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallNDcostnum

ggsave(paste("ggALL_nodiag_costnum_SA_", vaxfol, ".png", sep = ""), plot = ggallNDcostnum, width = 9, height = 8, dpi = 120)

ggallNDcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 38])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 39], ymax = med_annual_numvxMcost[, 40], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 38]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 42], ymax = med_annual_numvxMcost[, 43], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 41]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(5.5, -15), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (vaccine,TB treatment, \n ART) expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallNDcostnum

ggsave(paste("ggALL_nodiagsml_costnum_SA_", vaxfol, ".png", sep = ""), plot = ggallNDcostnum, width = 9, height = 8, dpi = 120)

ggallNAcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 44])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 45], ymax = med_annual_numvxMcost[, 46], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 44]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 48], ymax = med_annual_numvxMcost[, 49], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 47]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(100, -160), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (vaccine, TB treatment & diagnostics) \n expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallNAcostnum

ggsave(paste("ggALL_noart_costnum_SA_", vaxfol, ".png", sep = ""), plot = ggallNAcostnum, width = 9, height = 8, dpi = 120)

ggallNAcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 44])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 45], ymax = med_annual_numvxMcost[, 46], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 44]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 48], ymax = med_annual_numvxMcost[, 49], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 47]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(5.5, -15), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (TB treatment & diagnostics) \n expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallNAcostnum

ggsave(paste("ggALL_noartsml_costnum_SA_", vaxfol, ".png", sep = ""), plot = ggallNAcostnum, width = 9, height = 8, dpi = 120)

ggallNADcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 50])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 51], ymax = med_annual_numvxMcost[, 52], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 50]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 54], ymax = med_annual_numvxMcost[, 55], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 53]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(100, -160), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (vaccine&TB treatment) \n expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallNADcostnum

ggsave(paste("ggALL_noartdiag_costnum_SA_", vaxfol, ".png", sep = ""), plot = ggallNADcostnum, width = 9, height = 8, dpi = 120)

ggallNADcostnum <- ggplot(med_annual_numvxMcost, aes(x = Year[11:36], y = med_annual_numvxMcost[, 50])) +
  theme_classic() +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 51], ymax = med_annual_numvxMcost[, 52], width = 3), fill = "red", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 50]), colour = "red", size = 1) +
  geom_ribbon(aes(ymin = med_annual_numvxMcost[, 54], ymax = med_annual_numvxMcost[, 55], width = 3), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = med_annual_numvxMcost[, 53]), colour = "blue", size = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(5.5, -15), xlim = c(2024, 2052)) +
  labs(x = "Year", y = "Annual incremental total (vaccine&TB treatment) \n expenditure in South Africa (USD millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 28), axis.title = element_text(size = 24),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
ggallNADcostnum

ggsave(paste("ggALL_noartdiagsml_costnum_SA_", vaxfol, ".png", sep = ""), plot = ggallNADcostnum, width = 9, height = 8, dpi = 120)

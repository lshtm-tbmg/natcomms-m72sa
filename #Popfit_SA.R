library(plyr)
library(grid)
library(ggplot2)
library(leaflet)
library(gtable)
library(grid)
library(gridExtra)

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
  output <- "/home/lsh355020/SA_Gates/output"
}

if (C == 0) {
  input <- "/Users/lsh355020/R/Gates_sa/Vxoutput"
}
if (C == 1) {
  input <- "/home/lsh355020/SA_Gates/Vxoutput"
}

setwd(home)
setwd("Data")
data_pop <- as.data.frame(read.csv("population_SA.csv", sep = ","))

rrun <- 1

yrpoptot <- matrix(0, 51, rrun)
yrpop014 <- matrix(0, 51, rrun)
yrpop1564 <- matrix(0, 51, rrun)
yrpop65 <- matrix(0, 51, rrun)

for (kkk in 1:rrun) {
  print(kkk)

  Xn <- as.data.frame(Xn)

  yrpoptot[, kkk] <- Xn$PSIZEALL[101:151]
  yrpop014[, kkk] <- Xn$"YearPsize0-14"[101:151] + Xn$"YearPsize0-14H"[101:151]
  yrpop1564[, kkk] <- (Xn$YearPsize15plus[101:151] - Xn[101:151, 70]) + (Xn$"YearPsize15+H"[101:151] - Xn$"YearPsize65+H"[101:151])
  yrpop65[, kkk] <- Xn[101:151, 70] + Xn$"YearPsize65+H"[101:151]
}

med_pop <- matrix(0, 51, 12)

for (jj in 1:51) {
  med_pop[jj, c(1:3)] <- c(median(yrpoptot[jj, ]), min(yrpoptot[jj, ]), max(yrpoptot[jj, ]))
  med_pop[jj, c(4:6)] <- c(median(yrpop014[jj, ]), min(yrpop014[jj, ]), max(yrpop014[jj, ]))
  med_pop[jj, c(7:9)] <- c(median(yrpop1564[jj, ]), min(yrpop1564[jj, ]), max(yrpop1564[jj, ]))
  med_pop[jj, c(10:12)] <- c(median(yrpop65[jj, ]), min(yrpop65[jj, ]), max(yrpop65[jj, ]))
}

setwd(output)

Year <- seq(2000, 2050, 1)
med_pop <- as.data.frame(cbind(Year, med_pop))

write.table(med_pop, "med_pop.csv", sep = ",", row.names = F)

med_pop <- as.data.frame(med_pop)
med_pop2 <- cbind(med_pop[, 1], (med_pop[, 2:13] / 1000))
med_pop <- med_pop2
data_pop <- as.data.frame(data_pop)
data_pop <- cbind(data_pop[, 1], data_pop[, 2:5] / 1000)

setwd(home)
setwd("Demog")

ggpop <- ggplot(med_pop, aes(x = Year, y = med_pop[, 2])) +
  theme_classic() +
  geom_line(aes(y = med_pop[, 2]), colour = "black", size = 1) +
  geom_point(aes(x = data_pop[, 1], y = data_pop[, 5]), colour = "black", size = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Year", y = "All-age population (millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

ggpop014 <- ggplot(med_pop, aes(x = Year, y = med_pop[, 5])) +
  theme_classic() +
  geom_line(aes(y = med_pop[, 5]), colour = "red", size = 1) +
  geom_point(aes(x = data_pop[, 1], y = data_pop[, 2]), colour = "red", size = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Year", y = "Population 0-14 years (millions)") +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

ggpop1564 <- ggplot(med_pop, aes(x = Year, y = med_pop[, 8])) +
  theme_classic() +
  geom_line(aes(y = med_pop[, 8]), colour = "blue", size = 1) +
  geom_point(aes(x = data_pop[, 1], y = data_pop[, 3]), colour = "blue", size = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Year", y = "Population 15-64 years (millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

ggpop65 <- ggplot(med_pop, aes(x = Year, y = med_pop[, 11])) +
  theme_classic() +
  geom_line(aes(y = med_pop[, 11]), colour = "green", size = 1) +
  geom_point(aes(x = data_pop[, 1], y = data_pop[, 4]), colour = "green", size = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Year", y = "Population 65+ years (millions)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

ggpopall <- grid.arrange(ggpop, ggpop014, ggpop1564, ggpop65, ncol = 2)

ggsave("ggpopallTBon.pdf", plot = ggpopall, width = 16, height = 9, dpi = 120)

ggsave("ggpopH.pdf", plot = ggpop, width = 16, height = 9, dpi = 120)
ggsave("ggpop014H.pdf", plot = ggpop014, width = 16, height = 9, dpi = 120)
ggsave("ggpop1564H.pdf", plot = ggpop1564, width = 16, height = 9, dpi = 120)
ggsave("ggpop65H.pdf", plot = ggpop65, width = 16, height = 9, dpi = 120)

setwd(home)
setwd("Data")
hivP <- as.data.frame(read.csv("HIVprevfit.csv", sep = ","))

model_HIVP <- xoutsubH[xoutsubH[, "year"] %in% 1990:2016, ]

med_hivP <- as.data.frame(model_HIVP)
data_hivP <- as.data.frame(hivP)

setwd(home)
setwd("Demog")

gghivp014 <- ggplot(med_hivP, aes(x = med_hivP[, 5], y = med_hivP[, 2])) +
  theme_classic() +
  geom_line(aes(y = med_hivP[, 2]), colour = "red", size = 1) +
  geom_point(aes(x = data_hivP[, 1], y = data_hivP[, 2]), colour = "red", size = 2) +
  coord_cartesian(ylim = c(0, 4)) +
  labs(x = "Year", y = "HIV prevalence (0-14 year olds)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

gghivp1549 <- ggplot(med_hivP, aes(x = med_hivP[, 5], y = med_hivP[, 4])) +
  theme_classic() +
  geom_line(aes(y = med_hivP[, 4]), colour = "blue", size = 1) +
  geom_point(aes(x = data_hivP[, 1], y = data_hivP[, 4]), colour = "blue", size = 2) +
  coord_cartesian(ylim = c(0, 25)) +
  labs(x = "Year", y = "HIV prevalence (15 - 49 year olds)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

gghivp15plus <- ggplot(med_hivP, aes(x = med_hivP[, 5], y = med_hivP[, 3])) +
  theme_classic() +
  geom_line(aes(y = med_hivP[, 3]), colour = "blue", size = 1) +
  geom_point(aes(x = data_hivP[, 1], y = data_hivP[, 3]), colour = "blue", size = 2) +
  coord_cartesian(ylim = c(0, 25)) +
  labs(x = "Year", y = "HIV prevalence (15 plus year olds)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype = "solid")) +
  theme(
    axis.text = element_text(size = 18), axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 0.7, colour = "black", linetype = "solid")
  ) +
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

gghivall <- grid.arrange(gghivp014, gghivp1549, gghivp15plus, ncol = 2)

ggsave("gghivallzerou.pdf", plot = gghivall, width = 16, height = 9, dpi = 120)

ggsave("gghivp014.pdf", plot = gghivp014, width = 16, height = 9, dpi = 120)
ggsave("gghivp1549.pdf", plot = gghivp1549, width = 16, height = 9, dpi = 120)
ggsave("gghivp15plus.pdf", plot = gghivp15plus, width = 16, height = 9, dpi = 120)

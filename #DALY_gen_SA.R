library(purrr)
library(data.table)
data.table::setDTthreads(1)

setwd(home)
setwd("Vxoutput_M72_refit")

dalywt <- 0.333
dalywtH <- 0.408

disc_rate <- 0.03
disc_yr <- 2025

vcost <- 2.5
doses <- 2
vax_cost <- vcost * doses

pcdr <- 0.046
durdr <- 18
tcostds <- 187.4
tcostdr <- 6163.4
tcostdrY1 <- tcostdr / durdr * 12
tcostdrY2 <- tcostdr / durdr * 6

tcostall <- (tcostds * (1 - pcdr)) + (tcostdr * pcdr)
tcostallY1 <- (tcostds * (1 - pcdr)) + (tcostdrY1 * pcdr)
tcostallY2 <- (tcostdrY2 * pcdr)

diag_ratio <- 1.16
diag_cost <- 19.1
dcost <- diag_cost * diag_ratio

acostc <- 284.1
acosta <- 249.2

setwd(home)
setwd("Data")
ledf <- fread("unpd_le_cleaned_SA.csv")
ledf <- ledf[11:36, ]

setwd(home)
setwd("Vxoutput_M72_refit")
setwd(vaxfol)

prev <- fread(paste("TBPneg_age", kkk, ".csv", sep = ""))
prevH <- fread(paste("TBPHIV_age", kkk, ".csv", sep = ""))
prev <- prev[126:151, ]
prevH <- prevH[126:151, ]

prevb <- fread(paste("TBPneg_age_baseline", kkk, ".csv", sep = ""))
prevbH <- fread(paste("TBPHIV_age_baseline", kkk, ".csv", sep = ""))
prevb <- prevb[126:151, ]
prevbH <- prevbH[126:151, ]

deaths <- fread(paste("TBdeaths_age", kkk, ".csv", sep = ""))
deaths <- deaths[126:151, ]

deathsb <- fread(paste("TBdeaths_age_baseline", kkk, ".csv", sep = ""))
deathsb <- deathsb[126:151, ]

YLD <- prev[, .(YLD = rowSums(prev[, 1:101] * dalywt))][, Year := ledf[, 1]]
YLDH <- prevH[, .(YLDH = rowSums(prevH[, 1:101] * dalywtH))][, Year := ledf[, 1]]

YLDb <- prevb[, .(YLDb = rowSums(prevb[, 1:101] * dalywt))][, Year := ledf[, 1]]
YLDbH <- prevbH[, .(YLDbH = rowSums(prevbH[, 1:101] * dalywtH))][, Year := ledf[, 1]]

YLL <- deaths[, map2(deaths[, 1:101], ledf[, -1], `*`)]
YLL <- YLL[, .(YLL = rowSums(YLL[, 1:101]))][, Year := ledf[, 1]]

YLLb <- deathsb[, map2(deathsb[, 1:101], ledf[, -1], `*`)]
YLLb <- YLLb[, .(YLLb = rowSums(YLLb[, 1:101]))][, Year := ledf[, 1]]

DALY <- merge(YLL, YLD)
DALY <- merge(DALY, YLDH)
DALYb <- merge(YLLb, YLDb)
DALYb <- merge(DALYb, YLDbH)

DALY[, DALY := YLL + YLDH + YLD]
DALYb[, DALYb := YLLb + YLDbH + YLDb]

DALY[, d_DALY := DALY / (1 + disc_rate)^(.I - 1)]
DALYb[, d_DALYb := DALYb / (1 + disc_rate)^(.I - 1)]

NumV2 <- as.data.table(NumV[126:151, 1]) * 1000
vax_spend <- NumV2[, .(num_vax = (V1))][, Year := ledf[, 1]]
vax_spend[, vax_cost := vax_cost]
vax_spend[, vax_spend := (num_vax * vax_cost)]
vax_spend[, vax_doses := (num_vax * doses)]

vax_spend[, d_vax_spend := vax_spend / (1 + disc_rate)^(.I - 1)]

treat <- fread(paste("TBtreat", kkk, ".csv", sep = ""))
treatb <- fread(paste("TBtreat_baseline", kkk, ".csv", sep = ""))
t2024 <- treat[1, 1:3]
t2024b <- treatb[1, 1:3]

treat <- as.data.table(treat)
treatb <- as.data.table(treatb)

treat <- treat[2:27, 1]
treatb <- treatb[2:27, 1]

tb_spend <- treat[, treatcostY1 := V1 * tcostallY1]
tb_spend[, treatcost_toY2 := V1 * tcostallY2]
tb_spend[, treatcostY2 := c(0, treatcost_toY2[1:25])]
tb_spend[, treatcost := (treatcostY1 + treatcostY2)]
tb_spend[, d_treat_cost := treatcost / (1 + disc_rate)^(.I - 1)]

tb_spend[, diag_cost := V1 * dcost]
tb_spend[, d_diag_cost := diag_cost / (1 + disc_rate)^(.I - 1)]
tb_spend[, tb_spend := treatcost + diag_cost]
tb_spend[, d_tb_spend := d_treat_cost + d_diag_cost]

tb_spendb <- treatb[, treatcostY1b := V1 * tcostallY1]
tb_spendb[, treatcost_toY2b := V1 * tcostallY2]
tb_spendb[, treatcostY2b := c(0, treatcost_toY2b[1:25])]
tb_spendb[, treatcostb := (treatcostY1b + treatcostY2b)]
tb_spendb[, d_treat_costb := treatcostb / (1 + disc_rate)^(.I - 1)]

tb_spendb[, diag_costb := V1 * dcost]
tb_spendb[, d_diag_costb := diag_costb / (1 + disc_rate)^(.I - 1)]

tb_spendb[, tb_spendb := treatcostb + diag_costb]
tb_spendb[, d_tb_spendb := d_treat_costb + d_diag_costb]

tb_spend_diff <- tb_spend - tb_spendb
tb_spend_diff[, Year := ledf[, 1]]

artpop <- fread(paste("ART", kkk, ".csv", sep = ""))

artpopb <- fread(paste("ART_baseline", kkk, ".csv", sep = ""))

artcost <- c()
artcost <- as.data.table(artcost)
artcostb <- c()
artcostb <- as.data.table(artcostb)

artcostc <- artpop[, 1] * acostc
artcosta <- artpop[, 2] * acosta
artcost <- artcost[, artcost := artcostc + artcosta]
artcost[, d_artcost := artcost / (1 + disc_rate)^(.I - 1)]

artcostcb <- artpopb[, 1] * acostc
artcostab <- artpopb[, 2] * acosta
artcostb <- artcostb[, artcostb := artcostcb + artcostab]
artcostb[, d_artcostb := artcostb / (1 + disc_rate)^(.I - 1)]

art_spend_diff <- artcost - artcostb
art_spend_diff[, Year := ledf[, 1]]

costALL <- merge(vax_spend, tb_spend_diff)
costALL <- merge(costALL, art_spend_diff)
costALL[, costall := vax_spend + tb_spend + artcost]
costALL[, d_costall := d_vax_spend + d_tb_spend + d_artcost]
costALL[, costall_nd := vax_spend + treatcost + artcost]
costALL[, d_costall_nd := d_vax_spend + d_treat_cost + d_artcost]
costALL[, costall_na := vax_spend + tb_spend]
costALL[, d_costall_na := d_vax_spend + d_tb_spend]
costALL[, costall_nad := vax_spend + treatcost]
costALL[, d_costall_nad := d_vax_spend + d_treat_cost]

write.table(costALL, paste("All_costs", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

vcost_DALY <- ledf[, 1]
vcost_DALY[, DALYb := DALYb[, .(DALYb)]]
vcost_DALY[, DALY := DALY[, .(DALY)]]
vcost_DALY[, DALYdiff := DALYb - DALY]
vcost_DALY[, vax_spend := vax_spend[, vax_spend]]
vcost_DALY[, vcost_DALY := (vax_spend / DALYdiff)]

vcost_DALYd <- ledf[, 1]
vcost_DALYd[, d_DALYb := DALYb[, .(d_DALYb)]]
vcost_DALYd[, d_DALY := DALY[, .(d_DALY)]]
vcost_DALYd[, d_DALYdiff := d_DALYb - d_DALY]
vcost_DALYd[, d_vax_spend := vax_spend[, d_vax_spend]]
vcost_DALYd[, d_vcost_DALY := (d_vax_spend / d_DALYdiff)]

vcostDALY_st <- merge(vax_spend[, -c(4, 5)], vcost_DALY)
vcostDALY_all <- merge(vcostDALY_st, vcost_DALYd[, -5])

write.table(vcostDALY_all, paste("vaccinecost_DALY", kkk, ".csv", sep = ""), sep = ",", row.names = FALSE)

vcumDALY_list[[kkk]] <- vcostDALY_all[, lapply(.SD, sum), .SDcol = c("vax_spend", "d_vax_spend", "DALYdiff", "d_DALYdiff")][, cumcostDALY := vax_spend / DALYdiff][, cumcostDALYdis := d_vax_spend / d_DALYdiff][, run := kkk]

cost_DALY <- costALL[, .(Year, costall, d_costall, costall_nd, d_costall_nd, costall_na, d_costall_na, costall_nad, d_costall_nad)]
cost_DALY <- merge(cost_DALY, vcost_DALY[, .(DALYdiff, Year)])
cost_DALY <- merge(cost_DALY, vcost_DALYd[, .(d_DALYdiff, Year)])
cost_DALY[, cost_DALY := (costall / DALYdiff)]
cost_DALY[, d_cost_DALY := (d_costall / d_DALYdiff)]
cost_DALY[, cost_nd_DALY := (costall_nd / DALYdiff)]
cost_DALY[, d_cost_nd_DALY := (d_costall_nd / d_DALYdiff)]
cost_DALY[, cost_na_DALY := (costall_na / DALYdiff)]
cost_DALY[, d_cost_na_DALY := (d_costall_na / d_DALYdiff)]
cost_DALY[, cost_nad_DALY := (costall_nad / DALYdiff)]
cost_DALY[, d_cost_nad_DALY := (d_costall_nad / d_DALYdiff)]

cumDALY_list[[kkk]] <- cost_DALY[, lapply(.SD, sum), .SDcol = c("costall", "d_costall", "costall_nd", "d_costall_nd", "costall_na", "d_costall_na", "costall_nad", "d_costall_nad", "DALYdiff", "d_DALYdiff")][, cumcostDALY := costall / DALYdiff][, cumcostDALYdis := d_costall / d_DALYdiff][, cumcostDALY_nd := costall_nd / DALYdiff][, cumcostDALYdis_nd := d_costall_nd / d_DALYdiff][, cumcostDALY_na := costall_na / DALYdiff][, cumcostDALYdis_na := d_costall_na / d_DALYdiff][, cumcostDALY_nad := costall_nad / DALYdiff][, cumcostDALYdis_nad := d_costall_nad / d_DALYdiff][, run := kkk]

anDAL <- cbind(costALL[, -c(19:26)], cost_DALY[, -1])
annualDALY <- rbind(annualDALY, anDAL)

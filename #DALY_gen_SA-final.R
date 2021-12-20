
setwd(home);setwd("Vxoutput_M72_refit")

#'#############################################################################################################################################################################################
#'
#'Setting up stuff to save files down
#'
#'#############################################################################################################################################################################################

setwd(home)
#source('#vxSA_M72.R')

#pointing this code to Vaccine scenarios 
#vaxfol<-paste(typen,"_",effDis,"_",durs,"_",vage,"_",cover)
setwd(home);setwd("Vxoutput_M72_refit");setwd("duration_scens");setwd(vaxfol)


#'#############################################################################################################################################################################################
#'
#'Econ PSA parameters
#'
#'#############################################################################################################################################################################################

#PSA functions

#list of psa params - not used yet but preparing for loop
psa_params<- c("vcost","tcostds","tcostdr","diag_cost","acostc","acosta","patcost")

#Unit cost parameters 

#'#############################################################################################################################################################################################
#'
#'## Vaccine cost
#'
#'#############################################################################################################################################################################################
doses<-1

#Vaccine costs 
################################################################
#One dose $5 vaccine
#vcost_mean <- 5
#vcost_lci <- 1
#vcost_uci <- 9

#Reviewers requested sensitivity analysis of higher price
#One dose $5 vaccine
vcost_mean <- 7.53
vcost_lci <- 1.53
vcost_uci <- 13.53
##########################
cost_mean<-vcost_mean
cost_lci <-vcost_lci 
cost_uci <-vcost_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

vcost=rgamma(n = 1, shape = shape, scale = scale)
vax_cost<-vcost*doses
#'#############################################################################################################################################################################################
#'
#'## cost of treatment
#'
#'#############################################################################################################################################################################################

#DS treatment costs
################################################################

tcostds_mean<-165.2
tcostds_lci<-92.9
tcostds_uci<-286.7

##########################
cost_mean<-tcostds_mean
cost_lci <-tcostds_lci 
cost_uci <-tcostds_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

tcostds=rgamma(n = 1, shape = shape, scale = scale)

#DR treatment costs
################################################################

tcostdr_mean<-6163.4
tcostdr_lci<-6163.4-(1.96*1847)
tcostdr_uci<-6163.4+(1.96*1847)

##########################
cost_mean<-tcostdr_mean
cost_lci <-tcostdr_lci 
cost_uci <-tcostdr_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

tcostdr=rgamma(n = 1, shape = shape, scale = scale)

#%MDR 
pcdr<-0.046
durdr<-18

tcostdrY1<-tcostdr/durdr*12
tcostdrY2<-tcostdr/durdr*6

# no private in SA

tcostall<-(tcostds*(1-pcdr))+(tcostdr*pcdr)
tcostallY1<-(tcostds*(1-pcdr))+(tcostdrY1*pcdr)
tcostallY2<-(tcostdrY2*pcdr)


#'#############################################################################################################################################################################################
#'
#'## Diagnostic costs
#'
#'##############################################################################################################################################################################################Diagnostic costs

test_diag_ratio<-1.16
diag_cost_mean<-24.42*test_diag_ratio
diag_cost_lci<-24.3*test_diag_ratio
diag_cost_uci<-24.54*test_diag_ratio

##########################
cost_mean<-diag_cost_mean
cost_lci <-diag_cost_lci 
cost_uci <-diag_cost_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

dcost=rgamma(n = 1, shape = shape, scale = scale)

#'#############################################################################################################################################################################################
#'
#'## ART costs
#'
#'#############################################################################################################################################################################################

#ART costs - children 
################################################################

acostc_mean<-284.1
acostc_lci<-acostc_mean*1.2
acostc_uci<-acostc_mean*0.8

##########################
cost_mean<-acostc_mean
cost_lci <-acostc_lci 
cost_uci <-acostc_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

acostc=rgamma(n = 1, shape = shape, scale = scale)

#ART costs - adult
################################################################

acosta_mean<-249.8
acosta_lci<-acostc_mean*1.2
acosta_uci<-acostc_mean*0.8

##########################
cost_mean<-acosta_mean
cost_lci <-acosta_lci 
cost_uci <-acosta_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

acosta=rgamma(n = 1, shape = shape, scale = scale)

#'#############################################################################################################################################################################################
#'
#'## Patient costs
#'
#'##############################################################################################################################################################################################Diagnostic costs


patcost_mean<-112
patcost_lci<-36
patcost_uci<-168

##########################
cost_mean<-patcost_mean
cost_lci <-patcost_lci 
cost_uci <-patcost_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

patcost=rgamma(n = 1, shape = shape, scale = scale)


#'#############################################################################################################################################################################################
#'
#'setting DALY weights
#'
#'#############################################################################################################################################################################################

## Specify DALY weight
dalywt <- 0.333
dalywtH <- 0.408

## Specify discount rate
disc_rate <- 0.03
disc_yr<-2025

## First column of all dataframes is assumed to be "Year"
## The Calendar Years of all dataframes are assumed to be aligned and all data frames have same number of age (in years) groups

## Life Expectancy Data Frame - LEDF - Age-specific Life Expectancy by age and Calendar Year
setwd(home);setwd("Data")
ledf <- fread("unpd_le_cleaned_SA.csv")
ledf <- ledf[11:36,]


#'#############################################################################################################################################################################################
#'
#'Making summary files to calc costs 
#'
#'#############################################################################################################################################################################################


setwd(home);setwd("Vxoutput_M72_refit");setwd("duration_scens");setwd(vaxfol)

#'#############################################################################################################################################################################################
#'
##TB and HIV Prevalence - By Age and by calendar year
#'
#'#############################################################################################################################################################################################


### NEED HIV disagregated data from cfunct

prev <- fread(paste("TBPneg_age",kkk,".csv", sep=""))
prevH <- fread(paste("TBPHIV_age",kkk,".csv", sep=""))
prev <- prev[126:151,]
prevH <- prevH[126:151,]

prevb <- fread(paste("TBPneg_age_baseline",kkk,".csv", sep=""))
prevbH <- fread(paste("TBPHIV_age_baseline",kkk,".csv", sep=""))
prevb <- prevb[126:151,]
prevbH <- prevbH[126:151,]

#'#############################################################################################################################################################################################
#'
##TB deaths - By Age and by calendar year
#'
#'#############################################################################################################################################################################################

#not separate for HIV,as assume LE estimtes account for HIV prevalence in the population
deaths <- fread(paste("TBdeaths_age",kkk,".csv", sep=""))
deaths <- deaths[126:151,]

deathsb <- fread(paste("TBdeaths_age_baseline",kkk,".csv", sep=""))
deathsb <- deathsb[126:151,]

#'#############################################################################################################################################################################################
#'
##DALY calcs
#'
#'#############################################################################################################################################################################################

library(data.table)

#Intervention/vaccine scenario - ------------------------------------------------------------------------------------------------------------------------------------------- 

# Assume that model output is 2025-50
# Discount rate
dr <- 0.03
# Time of analysis for discounting
dyr <- 2018
# Set daly weight
dalywt <- 0.333
# Vaccine interval
vi_start <- 2025
vi_end <- 2050
# Endpoint years
epy <- c(2030, 2035, 2050)

# Clean up prev and deaths
# This code is brittle and was written temporarily during cleanup for my purposes - CW 05-01-2021 - MQ, I have changed 100=>99 to correct rows
prev[, Year := vi_start:vi_end]
prev[, V1 := NULL]
setnames(prev, new = c(0:99, "Year"))

deaths[, Year := vi_start:vi_end]
deaths[, V1 := NULL]
setnames(deaths, new = c(0:99, "Year"))

# Convert to long form
prev_l <- melt(prev, id.vars = "Year", variable.name = "Age", value.name = "Prevalence")
deaths_l <- melt(deaths, id.vars = "Year", variable.name = "Age", value.name = "Mortality")
ledf_l <- melt(ledf, id.vars = "Year", variable.name = "Age", value.name = "CLE") # conditional life expectancy

# Merge YLL working data frame
YLL_working <- ledf_l[deaths_l, on = c("Year", "Age")]
# Raw YLLs
YLL_working[, YLL_raw := Mortality  * CLE]
# Discount YLLs to time of death
YLL_working[, YLL_dtod := Mortality  * (1 - exp(-dr * CLE)) / dr]
# Discount YLLs to time of analysis
YLL_working[, YLL_d := YLL_dtod / (1 + dr)^(Year - dyr)]
# Sum YLLs across ages
YLL <- YLL_working[, .(YLL_nd = sum(YLL_raw), YLL_d = sum(YLL_d)), by = Year]

# YLD calculation
YLD_working <- prev_l[, YLD_nd := Prevalence ]
YLD_working <- prev_l[, YLD_d := Prevalence  / (1 + dr)^(Year - dyr)]
YLD <- YLD_working[, .(YLD_nd = sum(YLD_nd), YLD_d = sum(YLD_d)), by = Year]

# DALY raw
DALY <- merge(YLL, YLD)
DALY[, DALY_nd := YLL_nd + YLD_nd ]
DALY[, DALY_d := YLL_d + YLD_d ]



# Compute summed DALYs to endpoint years
sum_DALY <- rbindlist(
  lapply(epy, function(x) {
    DALY[Year <= x, lapply(X = .SD, sum), .SDcol = !"Year"][, .SD, .SDcol = !patterns("^YL")][, EPY := x]
  }
  )
)

#Renaming to fit into rest of RH code
DALY[, d_DALY := DALY_d]
DALY[, DALY := YLL_d + YLD_d ]

# Summary output:
#sum_DALY


#baseline scenario - -------------------------------------------------------------------------------------------------------------------------------------------

# Clean up prev and deaths
prevb[, Year := vi_start:vi_end]
prevb[, V1 := NULL]
setnames(prevb, new = c(0:99, "Year"))

deathsb[, Year := vi_start:vi_end]
deathsb[, V1 := NULL]
setnames(deathsb, new = c(0:99, "Year"))

# Convert to long form
prev_lb <- melt(prevb, id.vars = "Year", variable.name = "Age", value.name = "Prevalence")
deaths_lb <- melt(deathsb, id.vars = "Year", variable.name = "Age", value.name = "Mortality")
ledf_l <- melt(ledf, id.vars = "Year", variable.name = "Age", value.name = "CLE") # conditional life expectancy

# Merge YLL working data frame - MQ the square brackets relate to working data frame, hence not calling with X$Y
YLL_workingb <- ledf_l[deaths_lb, on = c("Year", "Age")]
# Raw YLLs
YLL_workingb[, YLL_raw := Mortality  * CLE]
# Discount YLLs to time of death
YLL_workingb[, YLL_dtod := Mortality  * (1 - exp(-dr * CLE)) / dr]
# Discount YLLs to time of analysis
YLL_workingb[, YLL_d := YLL_dtod / (1 + dr)^(Year - dyr)]
# Sum YLLs across ages
YLLb <- YLL_workingb[, .(YLL_nd = sum(YLL_raw), YLL_d = sum(YLL_d)), by = Year]

# YLD calculation
YLD_workingb <- prev_lb[, YLD_nd := Prevalence ]
YLD_workingb <- prev_lb[, YLD_d := Prevalence  / (1 + dr)^(Year - dyr)]
YLDb <- YLD_workingb[, .(YLD_nd = sum(YLD_nd), YLD_d = sum(YLD_d)), by = Year]

# DALY raw
DALYb <- merge(YLLb, YLDb)
DALYb[, DALY_nd := YLL_nd + YLD_nd ]
DALYb[, DALY_d := YLL_d + YLD_d ]


# Compute summed DALYs to endpoint years
sum_DALYb <- rbindlist(
  lapply(epy, function(x) {
    DALYb[Year <= x, lapply(X = .SD, sum), .SDcol = !"Year"][, .SD, .SDcol = !patterns("^YL")][, EPY := x]
  }
  )
)

#Renaming to fit into rest of RH code
DALYb[, d_DALY := DALY_d]
DALYb[, DALY := YLL_d + YLD_d ]

# Summary output:
#sum_DALY

#Old RH code below

#YLD <- prev[, .(YLD = rowSums(prev[,1:101] * dalywt))][, Year := ledf[, 1]]
#YLDH <- prevH[, .(YLDH = rowSums(prevH[,1:101] * dalywtH))][, Year := ledf[, 1]]
#
#YLDb <- prevb[, .(YLDb = rowSums(prevb[,1:101] * dalywt))][, Year := ledf[, 1]]
#YLDbH <- prevbH[, .(YLDbH = rowSums(prevbH[,1:101] * dalywtH))][, Year := ledf[, 1]]
#
#YLL <- deaths[, map2(deaths[,1:101], ledf[,-1], `*`)]
#YLL <- YLL[, .(YLL = rowSums(YLL[,1:101]))][, Year := ledf[, 1]]
#
#YLLb <- deathsb[, map2(deathsb[,1:101], ledf[,-1], `*`)]
#YLLb <- YLLb[, .(YLLb = rowSums(YLLb[,1:101]))][, Year := ledf[, 1]]
#
#DALY <- merge(YLL, YLD)
#DALY <- merge(DALY, YLDH)
#DALYb <- merge(YLLb, YLDb)
#DALYb <- merge(DALYb, YLDbH)
#
#DALY[, DALY := YLL + YLDH + YLD]
#DALYb[, DALYb := YLLb + YLDbH + YLDb]
#
## discounting
#DALY[, d_DALY := DALY/(1 + disc_rate)^(.I - 1)]
#DALYb[, d_DALYb := DALYb/(1 + disc_rate)^(.I - 1)]

#'#############################################################################################################################################################################################
#'
##Vaccine cost
#'
#'#############################################################################################################################################################################################

NumV<- fread(paste("number_vaccinated_",kkk,".csv", sep=""))

#NumV<- fread(paste("number_vaccinated_40.csv", sep=""))

NumV2<-as.data.table(NumV[126:151,1])*1000 ## multiply by 1000 as everything is /1000 in model
vax_spend<-NumV2[,.(num_vax = (V1))][, Year := ledf[, 1]]
vax_spend[, vax_cost := vax_cost]
vax_spend[, vax_spend := (num_vax * vax_cost)]
vax_spend[, vax_doses := (num_vax * doses)]
#vax_spend<-NumV2[,.(vax_spend = (V1 * vax_cost))][, Year := ledf[, 1]]

vax_spend[, d_vax_spend := vax_spend/(1 + disc_rate)^(.I - 1)]


#'#############################################################################################################################################################################################
#'
##Treatment and diagnosis cost
#'
#'#############################################################################################################################################################################################


#number treatment initiations 
treat <- fread(paste("TBtreat",kkk,".csv", sep="")) 
treatb <- fread(paste("TBtreat_baseline",kkk,".csv", sep="")) 
t2024<-treat[1,1:3]
t2024b<-treatb[1,1:3]

treat<-as.data.table(treat)
treatb<-as.data.table(treatb)

treat<-treat[2:27,1]
treatb<-treatb[2:27,1]

#treatment
tb_spend<-treat[,treatcostY1 := V1*tcostallY1]
tb_spend[,treatcost_toY2 := V1*tcostallY2]
#tb_spend[,treatcostY2 := c(as.numeric(t2024[,1]*tcostallY2), treatcost_toY2[1:25])]
tb_spend[,treatcostY2 := c(0, treatcost_toY2[1:25])]
tb_spend[,treatcost := (treatcostY1 + treatcostY2)]
tb_spend[, d_treat_cost := treatcost/(1 + disc_rate)^(.I - 1)]

#diagnostics
tb_spend[,diag_cost := V1*dcost] #diagnotic costs 
tb_spend[, d_diag_cost := diag_cost/(1 + disc_rate)^(.I - 1)]
#allTBprogramme
tb_spend[, tb_spend := treatcost + diag_cost]
tb_spend[, d_tb_spend := d_treat_cost + d_diag_cost]
#patient costs
tb_spend[, patient_cost := V1*patcost]
tb_spend[, d_patient_cost := patient_cost/(1 + disc_rate)^(.I - 1)]


###treatment baseline
tb_spendb<-treatb[,treatcostY1b := V1*tcostallY1]
tb_spendb[,treatcost_toY2b := V1*tcostallY2]
#tb_spendb[,treatcostY2b := c(as.numeric(t2024b[,1]*tcostallY2), treatcost_toY2b[1:25])]
tb_spendb[,treatcostY2b := c(0, treatcost_toY2b[1:25])]
tb_spendb[,treatcostb := (treatcostY1b + treatcostY2b)]
tb_spendb[, d_treat_costb := treatcostb/(1 + disc_rate)^(.I - 1)]

#diagnostics baseline
tb_spendb[,diag_costb := V1*dcost] #diagnotic costs - public only
tb_spendb[, d_diag_costb := diag_costb/(1 + disc_rate)^(.I - 1)]
#allTBprogramme
tb_spendb[, tb_spendb := treatcostb + diag_costb]
tb_spendb[, d_tb_spendb := d_treat_costb + d_diag_costb]
#patient costs
tb_spendb[, patient_costb := V1*patcost]
tb_spendb[, d_patient_costb := patient_costb/(1 + disc_rate)^(.I - 1)]


tb_spend_diff<-tb_spend-tb_spendb
tb_spend_diff[, Year := ledf[, 1]]

#'#############################################################################################################################################################################################
#'
##ART cost
#'
#'#############################################################################################################################################################################################

artpop<-fread(paste("ART",kkk,".csv", sep="")) ## 

artpopb<-fread(paste("ART_baseline",kkk,".csv", sep="")) ## 

artcost<-c()
artcost<-as.data.table(artcost)
artcostb<-c()
artcostb<-as.data.table(artcostb)

artcostc<-artpop[,1]*acostc
artcosta<-artpop[,2]*acosta
artcost<-artcost[,artcost:= artcostc+artcosta]
artcost[,d_artcost := artcost/(1 + disc_rate)^(.I - 1)]

artcostcb<-artpopb[,1]*acostc
artcostab<-artpopb[,2]*acosta
artcostb<-artcostb[,artcostb:= artcostcb+artcostab]
artcostb[,d_artcostb := artcostb/(1 + disc_rate)^(.I - 1)]

art_spend_diff<-artcost-artcostb
#art_spend_diff<-art_spend_diff[,-1]
art_spend_diff[, Year := ledf[, 1]]

#'#############################################################################################################################################################################################
#'
#'**
##Total cost i.e. vax costs plus programmatic costs
#'**
#'
#'#############################################################################################################################################################################################

costALL<-merge(vax_spend,tb_spend_diff)
costALL<-merge(costALL,art_spend_diff)
costALL[,costall := vax_spend + tb_spend + artcost]
costALL[,d_costall := d_vax_spend + d_tb_spend + d_artcost]
costALL[,costall_nd := vax_spend + treatcost + artcost] #without diagnostics
costALL[,d_costall_nd := d_vax_spend + d_treat_cost + d_artcost] #without diagnostics
costALL[,costall_na := vax_spend + tb_spend] #without art
costALL[,d_costall_na := d_vax_spend + d_tb_spend] #without art
costALL[,costall_nad := vax_spend + treatcost] #without art or diag
costALL[,d_costall_nad := d_vax_spend + d_treat_cost] #without art or diag
costALL[,costall_pat := vax_spend + tb_spend + artcost + patient_cost] #total incremental costs inc patient costs
costALL[,d_costall_pat := d_vax_spend + d_tb_spend + d_artcost + d_patient_cost] #total incremental costs inc patient costs - discounted

write.table(costALL,paste('All_costs',kkk,'.csv',sep=""),sep=",",row.names=FALSE)

#Delete after artspenddiff

All_costs_list[[kkk]] <- costALL [, run := kkk][,d_vax_spend := d_vax_spend][,d_tb_spend := d_tb_spend][,d_artcost := d_artcost]


#'#############################################################################################################################################################################################
#'
##Cost per DALY
#'
#'#############################################################################################################################################################################################


### Vaccine cost per DALY ### This isn't that useful as doesn't look at costs averted by vax - use overall costs below

#vax_spend<-NumV2[,.(num_vax = (V1))][, Year := ledf[, 1]]

#undiscounted
vcost_DALY<-ledf[,1]
vcost_DALY[,DALYb := DALYb[,.(DALY)]]
vcost_DALY[,DALY := DALY[,.(DALY)]]
vcost_DALY[, DALYdiff := DALYb-DALY]
vcost_DALY[, vax_spend := vax_spend[,vax_spend]]
vcost_DALY[, vcost_DALY := (vax_spend/DALYdiff)]
vcost_DALY[, costall := costALL[,costall]]
vcost_DALY[, costall_pat := costALL[,costall_pat]]


#discounted
vcost_DALYd<-ledf[,1]
vcost_DALYd[, d_DALYb := DALYb[,.(d_DALY)]]
vcost_DALYd[, d_DALY := DALY[,.(d_DALY)]]
vcost_DALYd[, d_DALYdiff := d_DALYb-d_DALY]
vcost_DALYd[, d_vax_spend := vax_spend[,d_vax_spend]]
vcost_DALYd[, d_vcost_DALY := (d_vax_spend/d_DALYdiff)]
vcost_DALY[, d_costall := costALL[,d_costall]]
vcost_DALY[, d_costall_pat := costALL[,d_costall_pat]]

vcostDALY_st<-merge(vax_spend[,-c(4,5)],vcost_DALY)
vcostDALY_all<-merge(vcostDALY_st,vcost_DALYd[,-5])

write.table(vcostDALY_all,paste('vaccinecost_DALY',kkk,'.csv',sep=""),sep=",",row.names=FALSE)

# # cumulative
# 
# cumDALY[kkk,1:5]<-c(kkk,(sum(vcostDALY_all[,vax_spend])),(sum(vcostDALY_all[,d_vax_spend])),(sum(vcostDALY_all[,DALYdiff])),(sum(vcostDALY_all[,d_DALYdiff])))
# cumDALY[kkk,6]<-(cumDALY[kkk,2]/cumDALY[kkk,4])
# cumDALY[kkk,7]<-(cumDALY[kkk,3]/cumDALY[kkk,5])


vcumDALY_list[[kkk]] <- vcostDALY_all[, lapply(.SD, sum), .SDcol = c("vax_spend", "d_vax_spend", "DALYdiff", "d_DALYdiff")][, cumcostDALY := vax_spend / DALYdiff][, cumcostDALYdis := d_vax_spend / d_DALYdiff][, run := kkk]

#vcumDALY_list_test <- vcostDALY_all[, lapply(.SD, sum), .SDcol = c("vax_spend", "d_vax_spend", "DALYdiff", "d_DALYdiff")][, cumcostDALY := vax_spend / DALYdiff][, cumcostDALYdis := d_vax_spend / d_DALYdiff][, run := kkk]



#'#############################################################################################################################################################################################
#'
##Net cost per DALY
#'
#'#############################################################################################################################################################################################


#undiscounted and discounted
cost_DALY<-costALL[,.(Year,costall,d_costall,costall_nd, d_costall_nd,costall_na, d_costall_na,costall_nad, d_costall_nad,costall_pat,d_costall_pat)]
cost_DALY<-merge(cost_DALY,vcost_DALY[,.(DALYdiff, Year)])
cost_DALY<-merge(cost_DALY,vcost_DALYd[,.(d_DALYdiff, Year)])
cost_DALY[, cost_DALY := (costall/DALYdiff)]
cost_DALY[, d_cost_DALY := (d_costall/d_DALYdiff)]
cost_DALY[, cost_nd_DALY := (costall_nd/DALYdiff)]
cost_DALY[, d_cost_nd_DALY := (d_costall_nd/d_DALYdiff)]
cost_DALY[, cost_na_DALY := (costall_na/DALYdiff)]
cost_DALY[, d_cost_na_DALY := (d_costall_na/d_DALYdiff)]
cost_DALY[, cost_nad_DALY := (costall_nad/DALYdiff)]
cost_DALY[, d_cost_nad_DALY := (d_costall_nad/d_DALYdiff)]
#adding patient costs
cost_DALY[, cost_DALY_pat := (costall_pat/DALYdiff)]
cost_DALY[, d_cost_DALY_pat := (d_costall_pat/d_DALYdiff)]


cumDALY_list[[kkk]] <- cost_DALY[, lapply(.SD, sum), .SDcol = c("costall", "d_costall", "costall_pat", "d_costall_pat", "costall_nd", "d_costall_nd","costall_na", "d_costall_na","costall_nad", "d_costall_nad", "DALYdiff", "d_DALYdiff")][, cumcostDALY := costall / DALYdiff][, cumcostDALYdis := d_costall / d_DALYdiff][, cumcostDALY_nd := costall_nd / DALYdiff][, cumcostDALYdis_nd := d_costall_nd / d_DALYdiff][, cumcostDALY_na := costall_na / DALYdiff][, cumcostDALYdis_na := d_costall_na / d_DALYdiff][, cumcostDALY_nad := costall_nad / DALYdiff][, cumcostDALYdis_nad := d_costall_nad / d_DALYdiff][, cumcostDALY_pat := costall_pat / DALYdiff][, cumcostDALYdis_pat := d_costall_pat / d_DALYdiff][, run := kkk]



###  other annual data ####
#anDAL<-cbind(costALL[,-c(19:26)],cost_DALY[,-1])
#annualDALY<-rbind(annualDALY,anDAL)

setwd(home)

#Below lines are from Jack checking gamma dist samples, but will do that after combining kkk files
# Simulations
#sims <- 
# Check the dist looks ok 
#hist(sims)
#Check we can recreate the confidence intervals ok and that the median is reasonable
#quantile(sims, c(0.025, 0.5, 0.975)) 


#taking onemost recent cost kkk outputs and summarising 

cost_summary<-c()
cost_summary$d_vax_spend<-vax_spend$d_vax_spend
cost_summary$d_treat_cost<-   tb_spend$d_treat_cost
cost_summary$d_diag_cost<-   tb_spend$d_diag_cost
cost_summary$d_patient_cost<-   tb_spend$d_patient_cost
cost_summary$d_artcost<-   artcost$d_artcost
cost_summary$d_tb_spend<-   tb_spend$d_tb_spend
cost_summary<-data.frame(cost_summary)

cost_totals<-c()
cost_totals$d_vax_spend<-    sum(cost_summary$d_vax_spend)
cost_totals$d_treat_cost<-   sum(cost_summary$d_treat_cost)
cost_totals$d_diag_cost<-    sum(cost_summary$d_diag_cost)
cost_totals$d_patient_cost<- sum(cost_summary$d_patient_cost)
cost_totals$d_artcost<-      sum(cost_summary$d_artcost)
cost_totals$d_tb_spend<-     sum(cost_summary$d_tb_spend)
cost_totals<-data.frame(cost_totals)        

#add columns for total discounted costs - vax, tb programme (treat,diag,patient), art
cost_totals<-cost_totals%>%
  mutate(
    total_d_cost = d_vax_spend+d_tb_spend+d_artcost,
    total_d_cost_pat = d_vax_spend+d_tb_spend+d_artcost+d_patient_cost,
    total_d_cost_noart = d_vax_spend+d_tb_spend+d_artcost+d_patient_cost,
    prop_tbprog=d_tb_spend/total_d_cost,
    prop_art=d_artcost/total_d_cost,
    prop_vax=d_vax_spend/total_d_cost,
    prop_tbprog_pat=d_tb_spend/total_d_cost_pat,
    prop_art_pat=d_artcost/total_d_cost_pat,
    prop_vax_pat=d_vax_spend/total_d_cost_pat,
    prop_pat=d_patient_cost/total_d_cost_pat,
    prop_tbprog_noart=d_tb_spend/total_d_cost_noart,
    prop_vax_noart=d_artcost/total_d_cost_noart
  )

#check=1
cost_totals<-cost_totals%>%
  mutate(
    check=prop_tbprog+prop_art+prop_vax
  )

#cost_tot[[kkk]] <-cost_totals





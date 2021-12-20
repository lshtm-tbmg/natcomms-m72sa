
#setwd(home);setwd("Vxoutput_M72_rerun_Ind")

#'#############################################################################################################################################################################################
#'
#'Setting up stuff to save files down
#'
#'#############################################################################################################################################################################################


#pointing this code to Vaccine scenarios 
#vaxfol<-paste(typen,"_",effDis,"_",durs,"_",vage,"_",cover)
setwd(home);setwd("Vxoutput_M72_rerun_Ind");setwd(vaxfol)

#test
#kkk<-2

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
nutritional_support=7.49
priv_sector_incentive=3.73
#DS treatment costs - public
################################################################

tcostds_mean<-317.16+nutritional_support
tcostds_lci<- 253.8+nutritional_support
tcostds_uci<- 374.04+nutritional_support

##########################
cost_mean<-tcostds_mean
cost_lci <-tcostds_lci 
cost_uci <-tcostds_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

tcostds=rgamma(n = 1, shape = shape, scale = scale)

#DR treatment costs - public
################################################################

tcostdr_mean<-3891.42+nutritional_support
tcostdr_lci<- 3381.84+nutritional_support
tcostdr_uci<- 4401+nutritional_support

##########################
cost_mean<-tcostdr_mean
cost_lci <-tcostdr_lci 
cost_uci <-tcostdr_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

tcostdr=rgamma(n = 1, shape = shape, scale = scale)


#DS treatment costs - private
################################################################


tcostdsp_mean<-317.16+nutritional_support+priv_sector_incentive
tcostdsp_lci<- 253.8+nutritional_support+priv_sector_incentive
tcostdsp_uci<- 374.04+nutritional_support+priv_sector_incentive

##########################
cost_mean<-tcostds_mean
cost_lci <-tcostds_lci 
cost_uci <-tcostds_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

tcostdsp=rgamma(n = 1, shape = shape, scale = scale)

#DR treatment costs - private #need to update
################################################################
#old costs
#tcostdsp<-31.8
#tcostdrp<-1635.8

tcostdrp_mean<-3891.42+nutritional_support+priv_sector_incentive
tcostdrp_lci<- 3381.84+nutritional_support+priv_sector_incentive
tcostdrp_uci<- 4401+nutritional_support+priv_sector_incentive

##########################
cost_mean<-tcostdr_mean
cost_lci <-tcostdr_lci 
cost_uci <-tcostdr_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

tcostdrp=rgamma(n = 1, shape = shape, scale = scale)

#'#############################################################################################################################################################################################
#'
#'## Diagnostic costs
#'
#'##############################################################################################################################################################################################Diagnostic costs

#for SA test_diag_ratio<-1.16
diag_cost_mean<-55.76
diag_cost_lci<-55.76*.8
diag_cost_uci<-55.76*1.2

##########################
cost_mean<-diag_cost_mean
cost_lci <-diag_cost_lci 
cost_uci <-diag_cost_uci 

cost_se <-     (cost_uci - cost_lci) / (2 * 1.96) 
shape <- (cost_mean ^ 2) / (cost_se ^ 2)
scale <- (cost_se ^ 2) / cost_mean

dcost=rgamma(n = 1, shape = shape, scale = scale)

#from here from RH files------------------------------------------------------------------------------------------------------

## cost of treatment

#%MDR 
pcdr<-0.0619
durdr<-18
#public
#tcostds<-150.7
#tcostdr<-4183.6
tcostdrY1<-tcostdr/durdr*12
tcostdrY2<-tcostdr/durdr*6
#private
#tcostdsp<-31.8
#tcostdrp<-1635.8
tcostdrpY1<-tcostdrp/durdr*12
tcostdrpY2<-tcostdrp/durdr*6

tcost<-(tcostds*(1-pcdr))+(tcostdr*pcdr)
tcostY1<-(tcostds*(1-pcdr))+(tcostdrY1*pcdr)
tcostY2<-(tcostdrY2*pcdr)

tcostp<-(tcostdsp*(1-pcdr))+(tcostdrp*pcdr)
tcostpY1<-(tcostdsp*(1-pcdr))+(tcostdrpY1*pcdr)
tcostpY2<-(tcostdrpY2*pcdr)

tcostall<-0.6*((tcostds*(1-pcdr))+(tcostdr*pcdr)) + 0.4*((tcostdsp*(1-pcdr))+(tcostdrp*pcdr))
tcostallY1<-0.6*((tcostds*(1-pcdr))+(tcostdrY1*pcdr)) + 0.4*((tcostdsp*(1-pcdr))+(tcostdrpY1*pcdr))
tcostallY2<-0.6*((tcostdrY2*pcdr)) + 0.4*(tcostdrpY2*pcdr)

## cost of diagnosis ###
## COSTS and ratio NEED INPUTTING ###
# diag_ratio<-6.48
# diag_cost<-15.5
#Fiamma's estiamte from XTEND is $53.65

dcost<-55.76

#------------------------------------------------------------------------------------------------------

#'#############################################################################################################################################################################################
#'
#'## Patient costs
#'
#'##############################################################################################################################################################################################Diagnostic costs


patcost_mean<-124
patcost_lci<-40 
patcost_uci<-186

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
#'End of cost section------------
#'
#'
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
setwd(home);setwd("Vxoutput_M72_rerun_Ind")
ledf <- fread("unpd_le_cleaned_IN.csv")
ledf <- ledf[11:36,]


#'#############################################################################################################################################################################################
#'
##TB (no HIV) Prevalence - By Age and by calendar year
#'
#'#############################################################################################################################################################################################

setwd(home);setwd("Vxoutput_M72_rerun_Ind");setwd(vaxfol)

## Prevalence - By Age and by calendar year
#will need to cycle over the different runs - integrate in to runs so that is faster? Also will need to run for baseline
prev <- fread(paste("TBP_age",kkk,".csv", sep=""))
prev <- prev[126:151,]

prevb <- fread(paste("TBP_age_baseline",kkk,".csv", sep=""))
prevb <- prevb[126:151,]

#'#############################################################################################################################################################################################
#'
##TB deaths - By Age and by calendar year
#'
#'#############################################################################################################################################################################################
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
# This code is brittle and was written temporarily during cleanup for my purposes - CW 05-01-2021 - MQ, I have changed 100=>99 to correct rows
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


#Old code below here

#YLD <- prev[, .(YLD = rowSums(prev[,1:101] * dalywt))][, Year := ledf[, 1]]
#YLDb <- prevb[, .(YLDb = rowSums(prevb[,1:101] * dalywt))][, Year := ledf[, 1]]
#
#YLL <- deaths[, map2(deaths[,1:101], ledf[,-1], `*`)]
#YLL <- YLL[, .(YLL = rowSums(YLL[,1:101]))][, Year := ledf[, 1]]
#
#YLLb <- deathsb[, map2(deathsb[,1:101], ledf[,-1], `*`)]
#YLLb <- YLLb[, .(YLLb = rowSums(YLLb[,1:101]))][, Year := ledf[, 1]]
#
#
#DALY <- merge(YLL, YLD)
#DALYb <- merge(YLLb, YLDb)
#
#DALY[, DALY := YLL + YLD]
#DALYb[, DALYb := YLLb + YLDb]
#
## discounting
#DALY[, d_DALY := DALY/(1 + disc_rate)^(.I - 1)]
#DALYb[, d_DALYb := DALYb/(1 + disc_rate)^(.I - 1)]

#'#############################################################################################################################################################################################
#'
##Vaccine cost
#'
#'#############################################################################################################################################################################################

#if it falls over here, there is an import from csv in the SA file which is not here - but looks like NumV2 is used
#Yep it fell over

NumV<- fread(paste("number_vaccinated_",kkk,".csv", sep=""))

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

treat<-treat[2:27,]
treatb<-treatb[2:27,]


#treatment
tb_spend<-treat[,treatcostY1 := V1*tcostallY1]
tb_spend[,treatcost_toY2 := V1*tcostallY2]
#tb_spend[,treatcostY2 := c(as.numeric(t2024[,1]*tcostallY2), treatcost_toY2[1:25])]
tb_spend[,treatcostY2 := c(0, treatcost_toY2[1:25])]
tb_spend[,treatcost := (treatcostY1 + treatcostY2)]
tb_spend[, d_treat_cost := treatcost/(1 + disc_rate)^(.I - 1)]
tb_spend[,tcostpubY1 := V2*tcostY1]
tb_spend[,tcostpub_toY2 := V2*tcostY2]
tb_spend[,tcostpubY2 := c(0, tcostpub_toY2[1:25])]
#tb_spend[,tcostpubY2 := c(as.numeric(t2024[,2]*tcostY2), tcostpub_toY2[1:25])]
tb_spend[,tcostpublic := (tcostpubY1 + tcostpubY2)]
tb_spend[, d_tcost_public := tcostpublic/(1 + disc_rate)^(.I - 1)]
tb_spend[,tcostpriY1 := V3*tcostpY1]
tb_spend[,tcostpri_toY2 := V3*tcostpY2]
#tb_spend[,tcostpriY2 := c(as.numeric(t2024[,3]*tcostpY2), tcostpri_toY2[1:25])]
tb_spend[,tcostpriY2 := c(0, tcostpri_toY2[1:25])]
tb_spend[,tcostprivate := (tcostpriY1 + tcostpriY2)]
tb_spend[, d_tcost_private := tcostprivate/(1 + disc_rate)^(.I - 1)]

#diagnostics
tb_spend[,diag_cost := V2*dcost] #diagnotic costs - public only
tb_spend[, d_diag_cost := diag_cost/(1 + disc_rate)^(.I - 1)]
#allTBprogramme
tb_spend[, tb_spend := treatcost + diag_cost]
tb_spend[, d_tb_spend := d_treat_cost + d_diag_cost]
#patient costs
tb_spend[, patient_cost := V2*patcost]
tb_spend[, d_patient_cost := patient_cost/(1 + disc_rate)^(.I - 1)]


#treatment baseline
tb_spendb<-treatb[,treatcostY1b := V1*tcostallY1]
tb_spendb[,treatcost_toY2b := V1*tcostallY2]
#tb_spendb[,treatcostY2b := c(as.numeric(t2024b[,1]*tcostallY2), treatcost_toY2b[1:25])]
tb_spendb[,treatcostY2b := c(0, treatcost_toY2b[1:25])]
tb_spendb[,treatcostb := (treatcostY1b + treatcostY2b)]
tb_spendb[, d_treat_costb := treatcostb/(1 + disc_rate)^(.I - 1)]
tb_spendb[,tcostpubY1b := V2*tcostY1]
tb_spendb[,tcostpub_toY2b := V2*tcostY2]
# tb_spendb[,tcostpubY2b := c(as.numeric(t2024b[,2]*tcostY2), tcostpub_toY2b[1:25])]
tb_spendb[,tcostpubY2b := c(0, tcostpub_toY2b[1:25])]
tb_spendb[,tcostpublicb := (tcostpubY1b + tcostpubY2b)]
tb_spendb[, d_tcost_publicb := tcostpublicb/(1 + disc_rate)^(.I - 1)]
tb_spendb[,tcostpriY1b := V3*tcostpY1]
tb_spendb[,tcostpri_toY2b := V3*tcostpY2]
#tb_spendb[,tcostpriY2b := c(as.numeric(t2024b[,3]*tcostpY2), tcostpri_toY2b[1:25])]
tb_spendb[,tcostpriY2b := c(0, tcostpri_toY2b[1:25])]
tb_spendb[,tcostprivateb := (tcostpriY1b + tcostpriY2b)]
tb_spendb[, d_tcost_privateb := tcostprivateb/(1 + disc_rate)^(.I - 1)]

#diagnostics
tb_spendb[,diag_costb := V2*dcost] #diagnotic costs - public only
tb_spendb[, d_diag_costb := diag_costb/(1 + disc_rate)^(.I - 1)]
#allTBprogramme
tb_spendb[, tb_spendb := treatcostb + diag_costb]
tb_spendb[, d_tb_spendb := d_treat_costb + d_diag_costb]
#patient costs
tb_spendb[, patient_costb := V2*patcost]
tb_spendb[, d_patient_costb := patient_costb/(1 + disc_rate)^(.I - 1)]


tb_spend_diff<-tb_spend-tb_spendb
tb_spend_diff[, Year := ledf[, 1]]

#Loads of stuff here was commented out, so have deleted now

#'#############################################################################################################################################################################################
#'
#'**
##Total cost i.e. vax costs plus programmatic costs
#'**
#'
#'#############################################################################################################################################################################################
### TOTAL incremental costs ### ART commented out, and removed in patient cost sum
costALL<-merge(vax_spend,tb_spend_diff)
#costALL<-merge(costALL,art_spend)
costALL[,costall := vax_spend + tb_spend]
#costALL[,costall := vax_spend + tb_spend + artcost]
costALL[,d_costall := d_vax_spend + d_tb_spend]
#costALL[,d_costall := d_vax_spend + d_tb_spend + d_artcost]
costALL[,costall_nd := vax_spend + treatcost] #without diagnostics
#costALL[,costall_nd := vax_spend + treatcost + d_artcost]
costALL[,d_costall_nd := d_vax_spend + d_treat_cost] #without diagnostics
costALL[,costall_pat := vax_spend + tb_spend  + patient_cost] #total incremental costs inc patient costs
costALL[,d_costall_pat := d_vax_spend + d_tb_spend  + d_patient_cost] #total incremental costs inc patient costs - discounted

write.table(costALL,paste('All_costs',kkk,'.csv',sep=""),sep=",",row.names=FALSE)

All_costs_list[[kkk]] <- costALL [, run := kkk][,vax_spend := vax_spend][,tb_spend := tb_spend]

#'#############################################################################################################################################################################################
#'
##Cost per DALY
#'
#'#############################################################################################################################################################################################

### Vaccine cost per DALY ### This isn't that useful as doesn't look at costs averted by vax - use overall costs below

#undiscounted
vcost_DALY<-ledf[,1]
vcost_DALY[,DALYb := DALYb[,.(DALY)]]
vcost_DALY[,DALY := DALY[,.(DALY)]]
vcost_DALY[, DALYdiff := DALYb-DALY]
vcost_DALY[, vax_spend := vax_spend[,vax_spend]]
vcost_DALY[, vcost_DALY := (vax_spend/DALYdiff)]

#discounted
vcost_DALYd<-ledf[,1]
vcost_DALYd[, d_DALYb := DALYb[,.(d_DALY)]]
vcost_DALYd[, d_DALY := DALY[,.(d_DALY)]]
vcost_DALYd[, d_DALYdiff := d_DALYb-d_DALY]
vcost_DALYd[, d_vax_spend := vax_spend[,d_vax_spend]]
vcost_DALYd[, d_vcost_DALY := (d_vax_spend/d_DALYdiff)]

vcostDALY_st<-merge(vax_spend[,-c(4,5)],vcost_DALY)
vcostDALY_all<-merge(vcostDALY_st,vcost_DALYd[,-5])

write.table(vcostDALY_all,paste('vaccinecost_DALY',kkk,'.csv',sep=""),sep=",",row.names=FALSE)

# cumulative vaccine cost per daly #

# vcumDALY[kkk,1:5]<-c(kkk,(sum(vcostDALY_all[,vax_spend])),(sum(vcostDALY_all[,d_vax_spend])),(sum(vcostDALY_all[,DALYdiff])),(sum(vcostDALY_all[,d_DALYdiff])))
# vcumDALY[kkk,6]<-(vcumDALY[kkk,2]/vcumDALY[kkk,4])
# vcumDALY[kkk,7]<-(vcumDALY[kkk,3]/vcumDALY[kkk,5])

vcumDALY_list[[kkk]] <- vcostDALY_all[, lapply(.SD, sum), .SDcol = c("vax_spend", "d_vax_spend", "DALYdiff", "d_DALYdiff")][, cumcostDALY := vax_spend / DALYdiff][, cumcostDALYdis := d_vax_spend / d_DALYdiff][, run := kkk]

#'#############################################################################################################################################################################################
#'
##Net cost per DALY
#'
#'#############################################################################################################################################################################################
#undiscounted and discounted
cost_DALY<-costALL[,.(Year,costall,d_costall,costall_nd, d_costall_nd,costall_pat,d_costall_pat)]
cost_DALY<-merge(cost_DALY,vcost_DALY[,.(DALYdiff, Year)])
cost_DALY<-merge(cost_DALY,vcost_DALYd[,.(d_DALYdiff, Year)])
cost_DALY[, cost_DALY := (costall/DALYdiff)]
cost_DALY[, d_cost_DALY := (d_costall/d_DALYdiff)]
cost_DALY[, cost_nd_DALY := (costall_nd/DALYdiff)]
cost_DALY[, d_cost_nd_DALY := (d_costall_nd/d_DALYdiff)]
#adding patient costs
cost_DALY[, cost_DALY_pat := (costall_pat/DALYdiff)]
cost_DALY[, d_cost_DALY_pat := (d_costall_pat/d_DALYdiff)]


# cumulative  cost per daly #


# cumDALY[kkk,1:7]<-c(kkk,(sum(cost_DALY[,costall])),(sum(cost_DALY[,d_costall])),(sum(cost_DALY[,costall_nd])),(sum(cost_DALY[,d_costall_nd])),(sum(cost_DALY[,DALYdiff])),(sum(cost_DALY[,d_DALYdiff])))
# cumDALY[kkk,8]<-(cumDALY[kkk,2]/cumDALY[kkk,6])
# cumDALY[kkk,9]<-(cumDALY[kkk,3]/cumDALY[kkk,7])
# cumDALY[kkk,10]<-(cumDALY[kkk,4]/cumDALY[kkk,6])
# cumDALY[kkk,11]<-(cumDALY[kkk,5]/cumDALY[kkk,7])

cumDALY_list[[kkk]] <- cost_DALY[, lapply(.SD, sum), .SDcol = c("costall", "d_costall",  "costall_pat", "d_costall_pat","costall_nd", "d_costall_nd", "DALYdiff", "d_DALYdiff")][, cumcostDALY := costall / DALYdiff][, cumcostDALYdis := d_costall / d_DALYdiff][, cumcostDALY_nd := costall_nd / DALYdiff][, cumcostDALYdis_nd := d_costall_nd / d_DALYdiff][, cumcostDALY_pat := costall_pat / DALYdiff][, cumcostDALYdis_pat := d_costall_pat / d_DALYdiff][, run := kkk]
 


#other annual data
anDAL<-cbind(costALL[,-c(29:32)],cost_DALY[,-1])

#annualDALY<-rbind(annualDALY,anDAL)

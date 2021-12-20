#MQ file to run and synthesise PSA - India

library(purrr)
library(data.table)

home<-"~/GitHub/M72Ind_v2"
setwd(home)

#using the SA file for the below as nothing in this is used in anger, just making a nice big matrix
para<-as.matrix(droplevels(read.csv('~/GitHub/M72SA/Data/fit1000finalSA_191005_1615seeds.csv',header=TRUE,check.names=F)))
para<-as.data.frame(para)
nm<-1:1000

#making lists for things used in this file
rrun<-1000
vcumDALY<-matrix(0,rrun,7);
cumDALY<-matrix(0,rrun,11);
annualDALY<-c()
cumDALY_list <- list()
vcumDALY_list <- list()
inc1550<-c()
incnum<-c()
treatnum<-c()
numart<-c()
CAV2535all<-c()
CAV2550all<-c()
DAV2550all<-c()
NVax2550all<-c()

#Manual debugging


#This is the bit that takes ages as cycles through each of 1000 runs 
setwd(home)

#list of non missing files
list<-1:1000 #exists_list
#ptm <- proc.time()
for (kkk in 1:1000){
  source("C:/Users/lsh390512/Documents/GitHub/M72Ind_v2/#DALY_gen_In_MQ-sa.R")
  print(kkk)
}
#proc.time() - ptm


#'#############################################################################################################################################################################################
#'
#'Now taking code from modelfit_active_GSA_vHratio_M72.R to pull runs together
#'
#'#############################################################################################################################################################################################
vaxfol<-scenario
setwd(home);setwd("Vxoutput_M72_rerun_Ind");setwd("duration_scens");setwd(vaxfol)

# colnames(cumDALY)<-c("run", "cumvaxcost","cumvaxcostdis", "cumDALY","cumDALYdis", "cumcostDALY", "cumcostDALYdis")
# write.table(cumDALY,'cumulativecostDALY.csv',row.names=FALSE)

cumDALY <- rbindlist(cumDALY_list,fill=TRUE)
vcumDALY <- rbindlist(vcumDALY_list)

write.table(vcumDALY,'1_cumulativeVaccinecostDALY.csv',sep=",",row.names=FALSE)

write.table(cumDALY,'1_cumulativecostDALY.csv',sep=",",row.names=FALSE)

#write.table(All_costs_list,'All_costs.csv',sep=",",row.names=FALSE)



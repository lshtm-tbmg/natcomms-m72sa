#To run and synthesise PSA - SA

library(purrr)
library(data.table)

home<-"~/GitHub/M72SA"
setwd(home)

#para<-as.matrix(droplevels(read.csv('fit1000finalSA_191005_1615seeds.csv',header=TRUE,check.names=F)))
#para<-as.data.frame(para)
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
All_costs_list<-c()

#Manual debugging
#kkk<-3
#source('#DALY_gen_SA_MQ.R')

#This is the bit that takes ages as cycles through each of 1000 runs 
setwd(home)

#list of non missing files
list<-1:1000 #exists_list
ptm <- proc.time()
for (kkk in list){
  source('#DALY_gen_SA_MQ.R')
  print(kkk)
}
proc.time() - ptm


#'#############################################################################################################################################################################################
#'
#'Now taking code from modelfit_active_GSA_vHratio_M72.R to pull runs together
#'
#'#############################################################################################################################################################################################
vaxfol<-scenario
setwd(home);setwd("Vxoutput_M72_refit");setwd("duration_scens");setwd(vaxfol)

# colnames(cumDALY)<-c("run", "cumvaxcost","cumvaxcostdis", "cumDALY","cumDALYdis", "cumcostDALY", "cumcostDALYdis")
# write.table(cumDALY,'cumulativecostDALY.csv',row.names=FALSE)

cumDALY <- rbindlist(cumDALY_list,fill=TRUE)
vcumDALY <- rbindlist(vcumDALY_list)
All_costs_list <- rbindlist(All_costs_list)
cost_total_scen<-rbindlist(cost_tot)


write.table(vcumDALY,'1_cumulativeVaccinecostDALY.csv',sep=",",row.names=FALSE)

write.table(cumDALY,'1_cumulativecostDALY.csv',sep=",",row.names=FALSE)

write.table(All_costs_list,'All_costs.csv',sep=",",row.names=FALSE)

write.table(cost_total_scen,'cost_total_scen.csv',sep=",",row.names=FALSE)




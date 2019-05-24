remove(list=ls())
library(sqldf)
set.seed(20181116)

rawM01 <- read.csv("E:/projects/TPD/data/TPDT_SongJian_20181115.csv",header=T,stringsAsFactors=F)
sampleId <-rawM01[,1]
t0 <- rawM01[,-1]

protIds <- colnames(t0)
t0$label <- as.character(sapply(sampleId,function(v){substr(v,1,1)}))
t0$patientId <- as.numeric(sapply(sampleId,function(v){substr(v,2,nchar(v)-1)}))
t0$repc <- sapply(sampleId,function(v){substr(v,nchar(v),nchar(v))})
t0 <- t0[,c("label","patientId","repc",protIds)]

tmp <- t0[,c("label","patientId",protIds[1])]

t1 <- sqldf(sprintf("SELECT label,patientId,avg(%s) %s FROM tmp group by label,patientId",protIds[1],protIds[1]))

indx <-2
while(indx < (length(protIds)-50)){
  avgs <- paste0(as.character(sapply(protIds[indx:(indx+49)],function(v){sprintf("avg(%s) %s",v,v)})),collapse=",")
  tmp <- t0[,c("label","patientId",protIds[indx:(indx+49)])]
  t2 <- sqldf(sprintf("SELECT %s FROM tmp group by label,patientId",avgs))
  t1 <- cbind(t1,t2)
  indx <- indx+50
}
avgs <- paste0(as.character(sapply(protIds[indx:length(protIds)],function(v){sprintf("avg(%s) %s",v,v)})),collapse=",")
tmp <- t0[,c("label","patientId",protIds[indx:length(protIds)])]
t2 <- sqldf(sprintf("SELECT %s FROM tmp group by label,patientId",avgs))
trainM <- cbind(t1,t2)
trainM <- trainM[,-2]
##################################### validation dataset
rawM02 <- read.csv("E:/projects/TPD/data/TPDV_SongJian_20181115.csv",header=T,stringsAsFactors=F)
sampleId <- rawM02[,1]
t0 <- rawM02[,-1]

protIds <- colnames(t0)
t0$patientId <- as.character(sapply(sampleId,function(v){strsplit(v,"_")[[1]][1]}))
t0 <- t0[,c("patientId",protIds)]

tmp <- t0[,c("patientId",protIds[1])]
t1 <- sqldf(sprintf("SELECT patientId,avg(%s) %s FROM tmp group by patientId",protIds[1],protIds[1]))

indx <-2
while(indx < (length(protIds)-50)){
  avgs <- paste0(as.character(sapply(protIds[indx:(indx+49)],function(v){sprintf("avg(%s) %s",v,v)})),collapse=",")
  tmp <- t0[,c("patientId",protIds[indx:(indx+49)])]
  t2 <- sqldf(sprintf("SELECT %s FROM tmp group by patientId",avgs))
  t1 <- cbind(t1,t2)
  indx <- indx+50
}
avgs <- paste0(as.character(sapply(protIds[indx:length(protIds)],function(v){sprintf("avg(%s) %s",v,v)})),collapse=",")
tmp <- t0[,c("patientId",protIds[indx:length(protIds)])]
t2 <- sqldf(sprintf("SELECT %s FROM tmp group by patientId",avgs))
testM <- cbind(t1,t2)
################################################data filter ##########################################################################



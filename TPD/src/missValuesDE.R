remove(list=ls())
library(data.table)

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

tradeOff <- 100
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
SWT <- SWT[,c("patientId",names(missingRate)[missingRate<tradeOff])]
SWT$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
SWT <- SWT[,c("label",names(missingRate)[missingRate<tradeOff])]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
SWV <- SWV[,names(missingRate)[missingRate<tradeOff]]

protIds <- intersect(colnames(SWT),colnames(SWV))
SWT <- SWT[,c("label",protIds)]
SWV <- SWV[,protIds]
#colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
#colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
protIds <- intersect(colnames(SWT),colnames(SWV))





missingRateMarks <- function(label,K,M0){
  M01 <- M0[M0$label==label[1],]
  M02 <- M0[M0$label==label[2],]
  miss01 <- apply(M01,2,function(v){100*sum(is.na(v))/length(v)})
  miss02 <- apply(M02,2,function(v){100*sum(is.na(v))/length(v)})
  miss <- abs(miss01-miss02)
  miss <- miss[order(miss,decreasing=T)]
  data.frame("protId"=names(miss[1:K]),"value"=miss[1:K],stringsAsFactors=F)
}

K=10
t1 <- missingRateMarks(c('N','M'),K,SWT)
t1 <- cbind(t1,missingRateMarks(c('N','A'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('N','C'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('N','P'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('M','A'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('M','C'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('M','P'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('A','C'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('A','P'),K,SWT))
t1 <- cbind(t1,missingRateMarks(c('C','P'),K,SWT))




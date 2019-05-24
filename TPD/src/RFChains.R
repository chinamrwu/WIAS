remove(list=ls())
library("randomForest")
set.seed(20181203)
source("E:/projects/TPD/src/growRF.R")

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
SWT <- SWT[,c("patientId",names(missingRate)[missingRate<tradeOff])]
SWT$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
SWT <- SWT[,c("label",names(missingRate)[missingRate<tradeOff])]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
SWV <- SWV[,names(missingRate)[missingRate<tradeOff]]

protIds <- intersect(colnames(SWT),colnames(SWV))
SWT <- SWT[,c("label",protIds)]
SWT[,-1] <- t(apply(SWT[,-1],1,function(v){v/mean(v,na.rm=T)}))
SWV <- SWV[,protIds]
SWV <- data.frame(t(apply(SWV,1,function(v){v/mean(v,na.rm=T)})),check.names=F)
SWT[is.na(SWT)] <- 0
SWV[is.na(SWV)] <- 0
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
protIds <- intersect(colnames(SWT),colnames(SWV))

LbLs <- c('NM','MA','AC','CP')

L=dim(SWT)[1]
predicted=c()
for(i in 1:L){
   trainM <- SWT[-i,];
   Models <- list()

   for(j in 1:4){
      label <- c(substr(LbLs[j],1,1),substr(LbLs[j],2,2))
      Models[[length(Models)+1]] <- growRF(label,label,trainM)
   }
   names(Models) <- LbLs
   predictions <- sapply(Models,function(model){as.matrix(predict(model$model,SWT[i,],type="prob"))})
   predicted <- rbind(predicted, matrix(predictions,nrow=1,ncol=8))
   print(sprintf('%d/%d',i,L))
 }
 cnames <- c()
 for(model in Models){
  cnames <- c(cnames,colnames(model$model$confusion)[1:2])
 }




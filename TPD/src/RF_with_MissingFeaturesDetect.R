remove(list=ls())
library("randomForest")
library("pROC")
library("data.table")
set.seed(20181203)
source("E:/projects/TPD/src/growRF.R")


getDataMatrix <- function(trainFile,testFile,missingCutOff=100,normalizedBy='row',na.replaceBy=0){
	#SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
	#SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
	SWT <- read.table(trainFile,sep="\t",header=T,stringsAsFactors=F,check.names=F)
	SWV <- read.table(testFile,sep="\t",header=T,stringsAsFactors=F,check.names=F)

	tradeOff <- missingCutOff
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
	SWT[is.na(SWT)] <- na.replaceBy
	SWV[is.na(SWV)] <- na.replaceBy
}
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
protIds <- intersect(colnames(SWT),colnames(SWV))

LbLs <- c('NM','NA','NC','NP','MA','MC','MP','AC','AP','CP')

L=dim(SWT)[1]
predicted=c()
trainM <- SWT
Models <- list()

for(j in 1:length(LbLs)){
      label <- c(substr(LbLs[j],1,1),substr(LbLs[j],2,2))
      Models[[length(Models)+1]] <- growRF(label,label,trainM)
}
names(Models) <- LbLs
 
mRates <- sapply(Models,function(m){
   label <- rownames(m$model$confusion)
   R1 <- t(sapply(label,function(ch){ M <- SWT[SWT$label==ch,m$protIds];apply(M,2,function(v){sum(v==0)/length(v)*100})}))
   R2 <- apply(SWT[SWT$label %in% label,m$protIds],2,function(v){sum(v==0)/length(v)*100})
   R3 <- apply(SWT[,m$protIds],2,function(v){sum(v==0)/length(v)*100})
   rbind(rbind(R1,R2),R3)
})

predictions <- c()
for(lb in c("NM","MA","AC","CP")){
  predictions <- cbind(predictions,as.matrix(predict(Models[[lb]]$model,SWV,type="prob")))
}
predictions <- predictions[,c(2,1,4,3,5,6,7,8)]

# write.table(predictions,file="E:/projects/TPD/results/RF_prediction_with_missingValues.txt",sep="\t",col.names=T,row.names=T,quote=F)

votes <- apply(predictions,1,function(v){
   r1 <- rep(0,5);names(r1) <- c('N','M','A','C','P');
   for(i in seq(1,7,by=2)){ ifelse(v[i] < v[i+1],r1[names(v)[i:(i+1)]] <- r1[names(v)[i:(i+1)]]+c(-1,1),r1[names(v)[i:(i+1)]] <- r1[names(v)[i:(i+1)]]+c(1,-1)) }
   r1
})


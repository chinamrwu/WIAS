remove(list=ls())
library("randomForest")
library("pROC")
library(data.table)
set.seed(20181119)
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
#SWT[,-1] <- t(apply(SWT[,-1],1,function(v){v/mean(v,na.rm=T)}))
SWV <- SWV[,protIds]
#SWV <- data.frame(t(apply(SWV,1,function(v){v/mean(v,na.rm=T)})))
SWT[is.na(SWT)] <- 0
SWV[is.na(SWV)] <- 0
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
protIds <- intersect(colnames(SWT),colnames(SWV))

#######################################################################################################

## featurePower computes the ability to classifying an instance two classes using Random Forest algorithmns 
## this function return a vector, c(AUC,TPR1,TPR2) where TPR equals 'true positive rate'
featurePower <- function(feature,tM){
	tmpRF <- randomForest(as.factor(label) ~ . ,data=tM[,c("label",feature)],importance=T,ntree=500,nodesize=5)
	predictions=as.data.frame(tmpRF$votes)
	clss <- colnames(predictions)
	predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	predictions$observed <- tM$label
	tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
        v1=c(tmpROC$auc,1-tmpRF$confusion[,3])
        names(v1) <- c("AUC",rownames(tmpRF$confusion))
	v1
}




Labels <- c('M','A','C','P')
labelPairs <- c()

powers <- list()
for(i in 1:3){
  for(j in (i+1):4){
     tmp <- SWT[SWT$label %in% c(Labels[i],Labels[j]),]
     powers[[length(powers)+1]] <- t(sapply(protIds,featurePower,tM=tmp))
     labelPairs <- c(labelPairs,paste0(Labels[i],Labels[j]))
  }
}
names(powers) <- labelPairs

comps <- list()
for(i in 1:3){
  for(j in (i+1):4){
   comps[[length(comps)+1]] <- as.matrix(dist(powers[,-1]))
  }
}


   
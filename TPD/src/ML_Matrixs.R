remove(list=ls())
library("data.table")
set.seed(20181203)
getDataMatrix <- function(trainFile,testFile,missingCutOff=100,normalizedBy='none'){
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
	SWV <- SWV[,protIds]

        if(normalizedBy=='row'){
	  SWT[,names(SWT)!='label'] <- t(apply(SWT[,names(SWT)!='label'],1,function(v){v/mean(v,na.rm=T)}));
	  SWV      <- data.frame(t(apply(SWV,1,function(v){v/mean(v,na.rm=T)})),check.names=F);
	}else if(normalizedBy=='column') {
          SWT[,-1] <- apply(   SWT[,-1], 2,function(v){m1 <- mean(v,na.rm=T);dv <- sd(v,na.rm=T); v-m1/dv})
	  SWV <- data.frame(apply(SWV,   2,function(v){m1 <- mean(v,na.rm=T);dv <- sd(v,na.rm=T); v-m1/dv}))
	}
	obj <- list("training"=SWT,"testing"=SWV)
	obj
}

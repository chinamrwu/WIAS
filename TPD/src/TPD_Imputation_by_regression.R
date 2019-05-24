remove(list=ls())
library("randomForest")
library("caret")
library("ggplot2")
set.seed(20181119)
source("E:/projects/TPD/src/growRF.R")
options(width=350)

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
rownames(SWT) <- SWT$patientId
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求

tradeOff <- 100
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
trainM <- SWT[,c('patientId',clnames)]
trainM$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
trainM <- trainM[,c("label",clnames)]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
testM <- SWV[,clnames]
protIds <- intersect(colnames(trainM),colnames(testM))

trainM <- trainM[,c(labelprotIds]
testM <- testM[,protIds]
#######################################################################################################################
lbls <- c('N','M','A','C','P')
innerRate=0.5
for( lb in lbls){
	 tmp <- trainM[trainM$label==lb,-1]
	 L <- dim(tmp)[1]
	 rate <- apply(tmp,2,function(v){sum(is.na(v))});
	 M0 <- tmp[,rate/L<=innerRate]
	 M1 <- tmp[,rate==0]
	 apply(M0,2,function(v){
            indx <- which(!is.na(v))
	    Y=v[indx]
	    cr=cor(Y,M1[indx,])
	    names(cr) <- colnames(M1)
	    cr <- cr[order(cr,decreasing=T)]

	 })



}
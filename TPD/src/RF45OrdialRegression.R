remove(list=ls())
library("sda")
library("pROC")
library("randomForest")
library("data.table")
set.seed(20181203)
source("E:/projects/TPD/src/growRF.R")
##http://strimmerlab.org/software/sda/download/sda-khan-data.pdf
matTrain <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
matV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

protList <- list()
for(tradeOff in rep(5,10)){
	SWT <- matTrain
	SWV <- matV
	missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
	SWT <- SWT[,c("patientId",names(missingRate)[missingRate<tradeOff])]
	SWT$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
	SWT <- SWT[,c("label",names(missingRate)[missingRate<tradeOff])]

	missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
	SWV <- SWV[,names(missingRate)[missingRate<tradeOff]]

	protIds <- intersect(colnames(SWT),colnames(SWV))
	SWT <- SWT[,c("label",protIds)]
	SWT[,-1] <- t(apply(SWT[,-1],2,function(v){m1 <- mean(v,na.rm=T);v/m1}))
	SWV <- SWV[,protIds]
	SWV <- t(apply(SWV,2,function(v){m1 <- mean(v,na.rm=T);v/m1}))
	SWT[is.na(SWT)] <- 0
	SWV[is.na(SWV)] <- 0
	colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
	colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
	protIds <- intersect(colnames(SWT),colnames(SWV))
	############################################################
	rfNM <- growRF(c('N','M'),c('N','M'),SWT)
	rfMA <- growRF(c('M','A'),c('M','A'),SWT)
	rfAC <- growRF(c('A','C'),c('A','C'),SWT)
	rfCP <- growRF(c('C','P'),c('C','P'),SWT)

	t1 <- unique(c(rfNM$protIds,rfMA$protIds,rfAC$protIds,rfCP$protIds))
	protList[[length(protList)+1]] <- t1
	print( sprintf("tradeOff:%d proteinNumber:%d,%f,%f,%f,%f",tradeOff,length(t1),rfNM$AUC,rfMA$AUC,rfAC$AUC,rfCP$AUC))
}

parseGene <- function(protId){
   a <- read.table(sprintf("https://www.uniprot.org/uniprot/%s.fasta",protId),sep="\t")[1,1]
   a <- strsplit(as.character(a)," ")[[1]]
   a <- a[grepl("GN=",a)]
   a <- strsplit(a,"=")[[1]][2]
   a
}
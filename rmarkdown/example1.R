
rm(list=ls())
library(data.table)
library(pROC)
library(randomForest)
library(caret)
library(sqldf)

SDCProbes <- c('cg13096260;cg18719750;cg24732574;cg08979737;cg25070637','cg08979737;cg25070637;cg14538332;cg16935295')
SDCProbes <- c(SDCProbes,'cg14538332;cg16935295')
TFPI2Probes <- c('cg24531255;cg17338208;cg26739865;cg22441533;cg14377593')
TFPI2Probes <- c(TFPI2Probes,'cg12973591;cg22799321')
probes <- unique(as.character(unlist(sapply(c(SDCProbes,TFPI2Probes),function(v){strsplit(v,";")[[1]]}))))

COAD <- fread('F:/projects/allData/TCGA/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/allData/TCGA/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)
mat450 <- mat450[,probes]

mat450 <- mat450[grepl('-01A-|-11A-',rownames(mat450)),] #438 samples:393 cancer vs 45 normal
Label <- as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'Cancer','Normal')}))
mat450$label <- Label
mat450 <- mat450[,c('label',probes)]

clinic <- read.table('F:/projects/allData/TCGA/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
submitters <- as.character(sapply(rownames(mat450),function(v){paste0(strsplit(v,"-")[[1]][1:3],collapse="-")}))
clinic <- clinic[clinic$submitter_id %in% submitters,c('submitter_id','site_of_resection_or_biopsy','tissue_or_organ_of_origin')] 
rownames(clinic) <- clinic$submitter_id
sampleSite <- data.frame('sampleId'=rownames(mat450),stringsAsFactors=F)
sampleSite$site =as.character(sapply(sampleSite$sampleId,function(v){a <- paste0(strsplit(v,"-")[[1]][1:3],collapse="-");
clinic[a,'site_of_resection_or_biopsy']}))
rownames(sampleSite) <- sampleSite[,1]
sampleSite$label <- "NO"
sampleSite$label[grepl('-01A-',rownames(sampleSite))] <- 'cancer'
sampleSite$label[grepl('-11A-',rownames(sampleSite))] <- 'normal'

T1 <- sqldf("SELECT site,count(*) CT FROM sampleSite group by site order by count(*) desc")
T1$ord <- 1:dim(T1)[1]
T2 <- sqldf("SELECT site,label,count(*) CT FROM sampleSite group by site,label")
T3 <- sqldf("SELECT T2.* FROM T1,T2 where T1.site=T2.site order by T1.ord")

#########################################################
getROC <- function(M0,iter=5){
	#M0 <- mat[,c(1:8)]
	set.seed(1)
	iter=5
	predictions <- matrix(0,nrow=dim(M0)[1],ncol=2)
	rownames(predictions) <- rownames(M0)
	for(i in 1:iter){
		folds <- createFolds(M0$label,5)
		for(fold in folds){
			valids <- M0[fold,]
			trains <- M0[setdiff(1:dim(M0)[1],fold),]
			trains$label <- as.factor(trains$label)
			tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
			predicted <- predict(tmpRF,valids,type='prob')
			predictions[rownames(predicted),] <- predictions[rownames(predicted),]+predicted
		}
	}
	colnames(predictions) <- colnames(predicted)
	predicts <- t(apply(predictions,1,function(v){v/sum(v)}))
	colnames(predicts) <- colnames(predicted)
	predicts <- data.frame(predicts,check.names=F)
	predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
	predicts$observed <- M0$label
	ROC <- roc(ifelse(predicts$observed=="Cancer", "Cancer", "Normal"), as.numeric(predicts$Cancer),direction='>')
   return(ROC)
}
normalSamples <- rownames(mat450)[grepl('-11A-',rownames(mat450))]

siteROC <- function(site,DMRs){
	sampleIds <- sampleSite$sampleId[sampleSite$site==site]
	sampleIds <- unique(c(sampleIds,normalSamples))
   
	k0 <- data.frame(sapply(DMRs,function(v){
       pbs <- strsplit(v,";")[[1]]
		 as.numeric(apply(mat450[sampleIds,pbs],1,function(v){mean(v,na.rm=T)}))
	}))
	k0$label <- as.factor(mat450[sampleIds,'label'])

   return(getROC(k0))

}


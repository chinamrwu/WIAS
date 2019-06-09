rm(list=ls())
library(data.table)
library(sqldf)

setwd("F:/WIAS/TPD")
source("src/common.R")
source("src/MachineLearning.R")

matDL <- read.csv('data/TPDDL_prot_190609_to_wu.csv',header=T,stringsAsFactors=F,check.names=F)
matDL <- matDL[grepl('HUMAN',matDL$prot),]
rownames(matDL)<- as.character(sapply(matDL$prot,function(v){strsplit(v,"\\|")[[1]][2]}))
protIds <- rownames(matDL)
matDL <- matDL[,-1]

df0 <- matDL[,!grepl('pool',colnames(matDL))]
fnames <- colnames(df0)
patients <- as.character(sapply(fnames,function(v){a <- strsplit(v,'_DIA_')[[1]][2];b <- strsplit(a,"_")[[1]][1:2];paste0(b,collapse="_")}))
patients <- unique(patients)

colnames(df0) <- as.character(sapply(colnames(df0),function(v){strsplit(v,'_DIA_')[[1]][2]}))
df0 <- 2^df0
df0 <- t(df0)

k0 <- sapply(patients,function(v){
    rname <- c(v,paste0(v,"_repB"))
    v1 <- apply(df0[rname,],2,function(v){log2(mean(v,na.rm=T))})
    v1
})
DL94 <- data.frame(t(k0))
DL94[is.na(DL94)] <- NA

tmp           <- sapply(rownames(DL94),function(v){ a <- strsplit(v,"_")[[1]];c(a[1],a[2])})
tmp           <- data.frame(t(tmp),stringsAsFactors=F)
tmp$X2        <- as.integer(tmp$X2)
tmp$patientId <- rownames(tmp)
tmp <- sqldf('SELECT * FROM tmp order by X1,X2') 

DL94 <- DL94[tmp$patientId,]
#write.table(DL94,file="data/TPD_avg_matDL94.txt",sep="\t",col.names=T,row.names=T,quote=F)

###################################################
B0  <- read.table("data/B0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
B0$label[B0$label %in% c('N','M','A')]  <- 'B'
B0$label[B0$label %in% c('C','P','W')]  <- 'M'

allModels <- read.table("data/AllModels_sta.txt",sep="\t",header=T,stringsAsFactors=F)
allModels <- allModels[order(allModels$avgAcc,decreasing=T),]
allModels <- allModels[allModels$src %in% c('A','B','C'),]
allProts  <- read.table("data/Models_Prots.txt",sep="\t",header=T,stringsAsFactors=F)

k0 <- sqldf("SELECT A.*,B.* FROM allModels A,allProts B where A.ID=B.ID and A.src='B'  and A.avgAcc > 85 and protNumber between 7 and 10 order by accB desc") 
seed <- k0$seed[1]
protId <- strsplit(k0$prots[1],",")[[1]];#"P01011" "P02751" "P09382" "P15090" "P19338" "P39059" "P46783" "Q13642" "Q99584" "Q9Y646"
# ID =11578

##########################

scaleRow <- function(M0){
 tmp <- M0
 tmp[,colnames(tmp)!='label'] <- t(scale(t(tmp[,colnames(tmp)!='label'])))
 tmp[is.na(tmp)] <- 0
 tmp
}


testModel <- function(model,testMat){
	feature    <- rownames(coef(model,s="lambda.1se"))[-1]
	matK       <-  scaleRow(testMat[,feature])
	score      <- predict(model,newx = as.matrix(matK),s = model$lambda.min,type='response')[,1]
	prediction <- data.frame("M"=score,"B"=1-score)
	prediction$predicted <- apply(prediction,1,function(v){names(v)[max(v)==v]})
	prediction
}



#ÑµÁ·Ä£ÐÍ
set.seed(seed)
trainDS <- scaleRow(PseudoSamples(B0[,c('label',protId)])) # normalization dataset
model <-  cv.glmnet(as.matrix(trainDS[,-1]),as.factor(trainDS$label),family='binomial',alpha=0,type.measure='class') ## training logistic regression model

########################################################
DL94Prediction <- testModel(model,DL94)
write.table(DL94Prediction,file='results/DL94_10prot_prediction.txt',sep="\t",col.names=T,row.names=T,quote=F)

########################################################
model6prots <- allModels[allModels$src=='B' & allModels$protNumber < 7,]
prot6 <- c('P01011','P02647','P02751','P39059','Q14103','Q9Y646')
seed <- 677
rate <- 35

set.seed(seed)
trainDS6prots <- scaleRow(PseudoSamples(B0[,c('label',prot6)]))
#trainDS <- scaleRow(kk)
model6prots <-  cv.glmnet(as.matrix(trainDS6prots[,-1]),as.factor(trainDS$label),family='binomial',alpha=0,type.measure='class')
DL94Prot6Prediction <- testModel(model6prots,DL94)
write.table(DL94Prot6Prediction,file='results/DL94_6prot_prediction.txt',sep="\t",col.names=T,row.names=T,quote=F)

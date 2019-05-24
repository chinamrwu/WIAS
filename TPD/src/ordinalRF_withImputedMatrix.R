remove(list=ls())
library("party")
library("randomForest")
library("pROC")
source("E:/projects/TPD/src/growRF.R")

rawMat <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_avgReplicas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(rawMat) <- rawMat[,1]
testData <- read.table("E:/projects/TPD/data/RF_TPDV_openSWATH_avgRepcas_20181108.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")

missingValueCutOff=80

missRate1 <- apply(rawMat[,-c(1:2)],2,function(v){100*sum(v==0)/(dim(rawMat)[1])})
missRate2 <- apply(testData,2,function(v){100*sum(is.na(v))/(dim(testData)[1])})
protIds <- intersect(names(missRate1)[missRate1<missingValueCutOff],names(missRate2)[missRate2<missingValueCutOff])



trainM <- rawMat[,c("label",protIds)]
validM <- testData[,protIds]

missRate1 <- apply(trainM[,-1],1,function(v){100*sum(v==0)/(dim(trainM)[2]-1)})
trainM <- trainM[names(missRate1)[missRate1<missingValueCutOff],]


#### in the condition without dimension reduction
cRF1 <- cforest(label ~ ., data = trainM,scores = list(label = c(12, 37,42,56)), control = cforest_unbiased(ntree = 1000))
validate1<- predict(cRF1,newdata =testData, OOB=TRUE, type = "response")


####################################################

tmp <- trainM[,c("label",protIds)]
cRF2 <- cforest(label ~ ., data = tmp,scores = list(label = c(12, 37,42,56)), control = cforest_unbiased(ntree = 1000))
validate2 <- predict(cRF2,newdata =testData, OOB=TRUE, type = "response")
predictions <- data.frame("ordRF588prot"=validate1,"ordRF45prot"=validate2,stringsAsFactors=F)
write.table(predictions,file="E:/projects/TPD/results/ordRF_validations.txt",sep="\t",col.names=T,row.names=F,quote=F)

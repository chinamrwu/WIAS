remove(list=ls())
library("party")
library("randomForest")
library("pROC")
source("E:/projects/TPD/src/growRF.R")
#rawMat <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_avgReplicas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rawMat <- read.table("E:/projects/TPD/data/TPDT_TRIC_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

rownames(rawMat) <- rawMat[,1]
testData <- read.table("E:/projects/TPD/data/TPDV_TRIC_20181115.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")

tradeOff=40

missCol <- apply(rawMat[,-1],2,function(v){100*sum(is.na(v))/length(v)})
protIds <- names(missCol)[missCol < tradeOff]
trainM <- rawMat[,protIds]
trainM$label <- as.character(sapply(rawMat$patientId,function(v){substr(v,1,1)}))
trainM <- trainM[,c("label",protIds)]
missRow <- apply(trainM[,-1],1,function(v){100*sum(is.na(v))/length(v)})
trainM <- trainM[missRow<tradeOff,]

missCol <- apply(testData,2,function(v){100*sum(is.na(v))/length(v)})
testM <- testData[,names(missCol)[missCol < tradeOff]]

protIds <- intersect(colnames(trainM),colnames(testM))

trainM <- trainM[,c("label",protIds)]
testM <- testM[,protIds]

distM <- as.matrix(dist(trainM[,-1]))

#### feature selection
lbls <- strsplit("NMACP","")[[1]]
clssDist <- sapply(lbls[1:5],function(c1){
    sapply(lbls[1:5],function(c2){
        rn <- rownames(trainM)[trainM$label==c1];
	cn <- rownames(trainM)[trainM$label==c2]
	sum(distM[rn,cn])/(length(rn)*length(cn))
})
})

print(clssDist)

set.seed(20181103)

topModel <- growRF(c("MA","CP"),c("B","M"))
benignModel <- growRF(c('M','A'),c('M','A'))
benignModel <- growRF(c('C','P'),c('C','P'))





////////////////////////////
trainM <- trainM[trainM$label!="N",]
trainM$label <- factor(trainM$label,ordered = TRUE,levels = c("M", "A","C","P"))

#### in the condition without dimension reduction
cRF1 <- cforest(label ~ ., data = trainM,scores = list(label = c(12, 37,42,56)), control = cforest_unbiased(ntree = 1000))
validate1<- predict(cRF1,newdata =testData, OOB=TRUE, type = "response")


####################################################

tmp <- trainM[,c("label",protIds)]
cRF2 <- cforest(label ~ ., data = tmp,scores = list(label = c(12, 37,42,56)), control = cforest_unbiased(ntree = 1000))
validate2 <- predict(cRF2,newdata =testData, OOB=TRUE, type = "response")
predictions <- data.frame("ordRF588prot"=validate1,"ordRF45prot"=validate2,stringsAsFactors=F)
write.table(predictions,file="E:/projects/TPD/results/ordRF_validations.txt",sep="\t",col.names=T,row.names=F,quote=F)

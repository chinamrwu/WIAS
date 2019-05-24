remove(list=ls())
library("party")
library("randomForest")
library("pROC")

rawMat <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
#rawMat <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_avgReplicas_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

rownames(rawMat) <- rawMat[,1]
testData <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")

tradeOff=5

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
trainM[is.na(trainM)] <- 0
testM[is.na(testM)] <- 0


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

protIds <- c()
set.seed(20181103)
for(i in 1:4){
   for(j in (i+1):5){
       
        a <- lbls[i]
        b <- lbls[j]
	lbl=paste0(a,b)

	tmp <- trainM[trainM$label==a | trainM$label==b,]
	fullRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000,nodesize=12,type="regression")

	imps  <- data.frame(importance(fullRF));
	impScore <- imps$MeanDecreaseAccuracy

	imps <- imps[order(impScore,decreasing=T),]
	orderedFeatures <- rownames(imps)
	
	indx <- 1
	featureIndex <- c()
	bestScore <- 0
	while ( indx <= 100 ){
		tmpRF <- NULL
		score <- 0
		currentFeatures <- c("label",orderedFeatures[1:indx])
		tmp <- trainM[trainM$label==a | trainM$label==b,currentFeatures]

		tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
		predictions=as.data.frame(tmpRF$votes)
		clss <- colnames(predictions)
		predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
		predictions$observed <- as.character(sapply(rownames(predictions),function(rn){substr(rn,1,1)}))

		tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
		score <- tmpROC$auc

		if(score > bestScore ){
		        featureIndex <- c(featureIndex,indx)
			bestScore <- score
		        print(sprintf("find better forest for %s with score: %f and %d features",lbl,bestScore,indx))
		}
      
         indx <- indx+1
       }##while()
       protIds <- c(protIds,orderedFeatures[featureIndex])
   }
}
protIds <- unique(protIds)
print(sprintf(".... features selection finished, %d features are found......",length(protIds)))

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

remove(list=ls())
library("party")
library("randomForest")
library("pROC")

rawMat <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_20181022.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(rawMat) <- rawMat[,1]
testData <- read.table("E:/projects/TPD/data/RF_TPDV_openSWATH_20181025.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
protIds <- intersect(colnames(rawMat)[3:dim(rawMat)[2]],colnames(testData))
trainM <- rawMat[,protIds]
missCol <- apply(trainM,2,function(v){100*sum(v==0)/(dim(trainM)[1])})
trainM <- trainM[,missCol<25]
missRow <- apply(trainM,1,function(v){100*sum(v==0)/(dim(trainM)[2]-1)})
trainM <- trainM[missRow<25,]
usedProtIds <- colnames(trainM)
trainM[is.na(trainM)] <- 0
patients <- rownames(trainM)
trainM$label <- rawMat[patients,"label"]
trainM <- trainM[,c("label",usedProtIds)]
#

#trainM$label <- factor(trainM$label,ordered = TRUE,levels = c("N", "M", "A","C","P"))

lbls=strsplit("NMACP","")[[1]]
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
	while ( indx <= 200 ){
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

testData <- testData[,protIds]
testData[is.na(testData)] <- 0

tmp <- trainM[trainM$label!="N",c("label",protIds)]
#tmp <- trainM[trainM$label!="N",c("label",protIds)]
forest <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
predictions <- as.data.frame(predict(forest, testData, type = "prob"))
predictions$predicted <- as.character(apply(predictions,1,function(v){names(v)[v==max(v)]}))

write.table(predictions,file="E:/projects/TPD/results/RF_multi_classes_validate_with_normal.txt",sep="\t",col.names=T,row.names=T,quote=F)
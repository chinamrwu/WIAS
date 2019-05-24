
###### Analysis of the output of RF
remove(list=ls())
library("randomForest")
require(pROC)

trainM <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_20181022.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(trainM) <- trainM[,1]
trainM <- trainM[,-1]

set.seed(1978)
uniqLabels =unique(trainM$label)
k <- length(uniqLabels)
nms <- c()
for(i in 1:(k-1)){
   for(j in (i+1):k){
     nms <- c(nms,paste0(uniqLabels[i],uniqLabels[j]))
   }
}

k=length(uniqLabels)

swaRFs <- list()
bestRFs <- NULL
bestScores <- 0
score <- 0

for(i in 1:(k-1)){
   for(j in (i+1):k){
   tmp <- trainM[trainM$label==uniqLabels[i] | trainM$label==uniqLabels[j],]
   tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)

   predictions=as.data.frame(tmpRF$votes)
   clss <- colnames(predictions)
   predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
   predictions$observed <- as.character(sapply(rownames(predictions),function(rn){substr(rn,1,1)}))
   
   tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
   tmpRF.roc <- tmpROC
   structure(tmpRF)
   score <- score+tmpROC$auc
   swaRFs[[length(swaRFs)+1]] <- tmpRF
}}
names(swaRFs) <- nms
if(score > bestScores){
   bestRFs <- swaRFs
   bestScores <- score
   print(bestScores)
}


orderedProtIds=list()
for(obj in swaRFs){ ## find important proteinId
             
	     imps  <- data.frame(importance(obj));
	     impScore <- imps$MeanDecreaseAccuracy * imps$MeanDecreaseGini

	     imps <- imps[order(impScore,decreasing=T),]
	     orderedProtIds[[length(orderedProtIds)+1]] <- rownames(imps)
}
names(orderedProtIds) <- nms

dlt <- 1

bestFeature <-list();
bestIndx <- 3
bestForest <- list()

for(lbl in nms){
   indx <- 3
   while ( indx < 1000 ){
        
	Forest <- NULL
	score <- 0
	bestScore <- 0
	   a <- substr(lbl,1,1);
	   b <- substr(lbl,2,2);
	   currentFeatures <- c("label",orderedProtIds[[lbl]][1:indx])
	   tmp <- trainM[trainM$label==a | trainM$label==b,currentFeatures]
	   
	   tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)

	   predictions=as.data.frame(tmpRF$votes)
	   clss <- colnames(predictions)
	   predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	   predictions$observed <- as.character(sapply(rownames(predictions),function(rn){substr(rn,1,1)}))
	   
	   RF <- list()
	   tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
	   RF$ROC <- tmpROC
	   RF$label <- lbl
	   RF$tree <- tmpRF
	   score <- tmpROC$auc
	   print(sprintf("%d  %f",indx,score))
	   if(score > bestScore){
		Forest <- RF
		bestScore <- score
		bestFeature <- currentFeatures
		bestIndx <- indx
		print(sprintf("find better tree for %s with score: %f",lbl,bestScores))
	   }
      
     indx <- indx+dlt
    }
     bestForest[[length(bestForst)+1]] <- Forest 
   }
 
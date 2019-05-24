
###### Analysis of the output of RF
remove(list=ls())
library("randomForest")
require(pROC)

trainM <- read.table("E:/projects/TPD/data/RF_TPDT_expert_trainingM_20181024.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
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


dlt <- 1
L <- dim(trainM)[2]-1
classNumber <- length(uniqLabels)

indx <- 50
currentFeatures <- c()
bestFeature <- c();

while ( indx < 100 ){
        swaRFs <- list()
	for( protIds in orderedProtIds){
	   currentFeatures <- c(currentFeatures,protIds[1:indx])
	}
	currentFeatures <- unique(currentFeatures)
	tmpM <- trainM[,c("label",currentFeatures)]

        score <- 0
        for(i in 1:(classNumber-1)){
           for(j in (i+1):classNumber){
	   tmp <- tmpM[tmpM$label==uniqLabels[i] | tmpM$label==uniqLabels[j],]
	   tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)

	   predictions=as.data.frame(tmpRF$votes)
	   clss <- colnames(predictions)
	   predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	   predictions$observed <- as.character(sapply(rownames(predictions),function(rn){substr(rn,1,1)}))
	   
	   RF <- list()
	   tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
	   RF$ROC <- tmpROC
	   RF$label <- paste0(clss,collapse="")
	   RF$tree <- tmpRF
	   score <- score+tmpROC$auc
	   swaRFs[[length(swaRFs)+1]] <- RF
        }
       }
       
       names(swaRFs) <- nms
       print(score)
       if(score > bestScores){
		bestRFs <- swaRFs
		bestScores <- score
		bestFeature <- currentFeatures
		print(paste0("find better: ",bestScores))
        }
       indx <- indx+dlt
 }

############################## deal with best Random Forest


pdf("E:/projects/TPD/results/RF_TPD_expert_Best_ROC_All_in_ONE.pdf")
mt <- matrix(c(1:10,0,0),nrow=4,ncol=3)
layout(mt)
for(i in 1:10){
    obj <- bestRFs[[i]]
    plot.roc(obj$ROC,print.auc=T,col = "blue3",print.thres="best",main=obj$label,legacy.axes = TRUE,print.auc.cex=1.2)
}
dev.off()


##### get best importance features





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
dlt <- 1
bestRFs <- list()

for(lbl in nms){
   
	a <- substr(lbl,1,1);
	b <- substr(lbl,2,2);
	tmp <- trainM[trainM$label==a | trainM$label==b,]
	tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)

	imps  <- data.frame(importance(tmpRF));
	impScore <- imps$MeanDecreaseAccuracy * imps$MeanDecreaseGini

	imps <- imps[order(impScore,decreasing=T),]
	orderedFeatures <- rownames(imps)

	indx <- 3
	bestScore <- 0
	bestForest <- NULL
	while ( indx < 1000 ){
        
		Forest <- NULL
		score <- 0
		
		currentFeatures <- c("label",orderedFeatures[1:indx])
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

		if(score > bestScore){
			bestScore <- score
			bestForest <- RF
			bestFeature <- currentFeatures
		        print(sprintf("find better tree for %s with score: %f and %d features",lbl,bestScore,length(bestFeature)))
		}
      
     indx <- indx+dlt
    }
     bestRFs[[length(bestForest)+1]] <- bestForest 
   }
 
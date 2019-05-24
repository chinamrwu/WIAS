remove(list=ls())
library("randomForest")
library("pROC")

rawMat <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_avgReplicas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(rawMat) <- rawMat[,1]
testData <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
protIds <- intersect(colnames(rawMat)[3:dim(rawMat)[2]],colnames(testData))
trainM <- rawMat[,protIds]
missCol <- apply(trainM,2,function(v){100*sum(v==0)/(dim(trainM)[1])})
trainM <- trainM[,missCol<25]
missRow <- apply(trainM,1,function(v){100*sum(v==0)/(dim(trainM)[2]-1)})
trainM <- trainM[missRow<25,]
usedProtIds <- colnames(trainM)

patients <- rownames(trainM)
trainM$label <- rawMat[patients,"label"]
trainM <- trainM[,c("label",usedProtIds)]


#### feature selection
lbls <- strsplit("NMACP","")[[1]]
protIds <- c()
set.seed(20181103)
for(i in 1:4){
   for(j in (i+1):5){
       
        a <- lbls[i]
        b <- lbls[j]
	lbl=paste0(a,b)

	tmp <- trainM[trainM$label==a | trainM$label==b,]
	fullRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)

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

		tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)
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
M1 <- trainM[,c("label",protIds)]
distM1 <- as.matrix(dist(M1[,-1]))
for(i in 1:dim(distM1)[1]){
   distM1[i,i] <- .Machine$double.xmax
}

pLabels <- rawMat[,c("patientId","label")]

for(k in 1:15){
t1 <- apply(distM1,2,function(v){
  kNNs <- names(v)[order(v)][1:k]
  iLbls <- pLabels$label[match(kNNs,pLabels$patientId)]
  r1 <- sapply(lbls,function(s1){sum(s1==iLbls)})
  mxNumber <- sum(r1==max(r1))
  if(mxNumber > 1 ){ ## in cases where different labels have the max neighbors number,eg [1,1,3,3,1,1],then nearest neighbor used
     nn   <- names(v)[order(v)][1]
     iLbl <- pLabels$label[match(nn,pLabels$patientId)]
     r1 <- rep(0,5)
     names(r1) <- lbls;
     r1[iLbl] <- 1
  }
  r1
})

NNs <- data.frame(t(t1))
NNs$predicted <- unlist(apply(NNs,1,function(v){names(v)[v==max(v)]}))
NNs$observed <- pLabels$label[match(rownames(NNs),pLabels$patientId)]
  print(sprintf("k=%d, accuracy=%f",k,sum(NNs$predicted==NNs$observed)/399.0))
}


##### predict the validation dataset
testData[is.na(testData)] <- 0
testData <- testData[,protIds]
distM2 <- apply(testData,1,function(v){
     apply(trainM[,protIds],1,function(v1){
        sqrt(sum((v-v1)^2))
})
})

k=7

t2 <- apply(distM2,2,function(v){
  kNNs <- names(v)[order(v)][1:k]
  iLbls <- pLabels$label[match(kNNs,pLabels$patientId)]
  r1 <- sapply(lbls,function(s1){sum(s1==iLbls)})
  mxNumber <- sum(r1==max(r1))
  if(mxNumber > 1 ){ ## in cases where different labels have the max neighbors number,eg [1,1,3,3,1,1],then nearest neighbor used
     nn   <- names(v)[order(v)][1]
     iLbl <- pLabels$label[match(nn,pLabels$patientId)]
     r1 <- rep(0,5)
     names(r1) <- lbls;
     r1[iLbl] <- 1
  }
  r1
})

vNNs <- data.frame(t(t2))
vNNs <- data.frame(t(apply(vNNs,1,function(v){v/sum(v)})))
vNNs$predicted <- unlist(apply(vNNs,1,function(v){names(v)[v==max(v)]}))
#write.table(vNNs,file="E:/projects/TPD/results/validation_TPD_SWATH_kNN.txt",sep="\t",col.names = T,row.names = T,quote = F)

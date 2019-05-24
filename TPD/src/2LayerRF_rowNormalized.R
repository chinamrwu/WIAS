remove(list=ls())
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
trainM[is.na(trainM)] <- 0;
usedProtIds <- colnames(trainM)



patients <- rownames(trainM)
trainM$label <- rawMat[patients,"label"]
trainM <- trainM[,c("label",usedProtIds)]
testData <- testData[,protIds]
testData[is.na(testData)] <- 0

set.seed(20181106)
#############################################
 trainModel <- function(a,b,it=100){
	lbls01 <- strsplit(a[1],"")[[1]]
	lbls02 <- strsplit(a[2],"")[[1]]
	M1 <- trainM[trainM$label %in% c(lbls01,lbls02),]
	M1$label[M1$label %in% lbls01] <- b[1]
	M1$label[M1$label %in% lbls02] <- b[2]
	observed <- data.frame("sampleId"=rownames(M1),"label"=M1$label)

	result <- list()
	result$Matrix <- M1

	fullRF <- randomForest(formula=as.factor(M1$label) ~ . ,data=M1,importance=T,ntree=1000,nodesize=12,type="regression")

	imps  <- data.frame(importance(fullRF));
	impScore <- imps$MeanDecreaseAccuracy
	imps <- imps[order(impScore,decreasing=T),]
	orderedFeatures <- rownames(imps)

	indx <- 1
	featureIndex <- c()
	bestScore <- 0
	while ( indx <= it ){
		tmpRF <- NULL
		score <- 0
		currentFeatures <- c("label",orderedFeatures[1:indx])
		tmp <- M1[,currentFeatures]
		tmp$label <- as.factor(tmp$label);

		tmpRF <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
		predictions=as.data.frame(tmpRF$votes)
		clss <- colnames(predictions)
		predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
		predictions$observed <- observed$label

		tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
		score <- tmpROC$auc

		if(score > bestScore ){
			featureIndex <- c(featureIndex,indx)
			bestScore <- score
			#print(sprintf("find better forest with score: %f and %d features",bestScore,indx))
		}

                #print(indx)
		indx <- indx+1
	}##while()
	selected <- orderedFeatures[featureIndex]
	result$protIds <- selected
        print(sprintf("finished features selection and %d proteins were selected",length(selected)))
	tmp <- M1[,c('label',selected)]
	tmp$label <- as.factor(tmp$label);
        result$model <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
	result
}

 topLayer <- trainModel(c("MA","CP"),c("B","M"))
 benign   <- trainModel(c('M','A'),c('M','A'))
 mal      <- trainModel(c('C','P'),c('C','P'))

 topPrediction    <- as.matrix(predict(topLayer$model,testData, type = "prob"))
 benignPrediction <- as.matrix(predict(benign$model, testData, type = "prob"))
 malPrediction    <- as.matrix(predict(mal$model, testData, type = "prob"))

 #predictions<- t(sapply(1:dim(topPrediction)[1],function(i){c(topPrediction[i,1]*benignPrediction[i,],topPrediction[i,2]*malPrediction[i,])}))
 #t1 <- apply(predictions,1,function(v){colnames(predictions)[max(v)==v]})
 #predictions <- data.frame(predictions)
 #predictions$predicted <- t1
 #predictions$patient <- rownames(topPrediction)
 #predictions <- predictions[,c(6,1,2,3,4,5)]
 #write.table(predictions,file="E:/projects/TPD/results/two_layers_RF_validation.txt",sep="\t",col.names=T,row.names=F,quote=F)
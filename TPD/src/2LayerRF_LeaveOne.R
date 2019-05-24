remove(list=ls())
library("randomForest")
library("pROC")

set.seed(20181106)

 trainModel <- function(a,b,M0,it=100){
	lbls01 <- strsplit(a[1],"")[[1]]
	lbls02 <- strsplit(a[2],"")[[1]]
	M1 <- M0[M0$label %in% c(lbls01,lbls02),]
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
        #print(sprintf("finished features selection and %d proteins were selected",length(selected)))
	tmp <- M1[,c('label',selected)]
	tmp$label <- as.factor(tmp$label);
        result$model <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
	result
}

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

trainM <- data.frame(t(apply(trainM,1,function(v){v/(mean(v))}))) #Normalizition

usedProtIds <- colnames(trainM)
patients <- rownames(trainM)
trainM$label <- rawMat[patients,"label"]
trainM <- trainM[,c("label",usedProtIds)]

tmp <- trainM[trainM$label!='N',]
t1 <- sapply(1:dim(tmp)[1],function(i){
    top      <- trainModel(c('MA','CP'),c('B','M'),tmp[-i,])
    benign   <- trainModel(c('M','A'),c('M','A'),  tmp[-i,])
    mal      <- trainModel(c('C','P'),c('C','P'),  tmp[-i,])

    topPredict    <- as.matrix(predict(top$model,tmp[i,], type = 'prob'))
    benignPredict <- as.matrix(predict(benign$model,tmp[i,], type = "prob"))
    malPredict    <- as.matrix(predict(mal$model,tmp[i,], type = "prob"))
    v <- c(topPredict[1,'B']*benignPredict[1,],topPredict[1,'M']*malPredict[1,])
    pred <- names(v)[max(v)==v]
    v <- c(v,substr(rownames(tmp)[i],1,1),pred)
    print(paste0(v,collapse=" "))
    v
})

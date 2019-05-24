 # training a random forest model of binary classification
 #  a: a character vector with 2 elements,each represents a class;for example c("AP","MC") indicates all samples from "A" or "P" 
 #   
 library("randomForest")
 library("data.table")
 library("pROC")
 growRF <- function(a,b,M0,it=100){
	lbls01 <- strsplit(a[1],"")[[1]]
	lbls02 <- strsplit(a[2],"")[[1]]
	M1 <- M0[M0$label %in% c(lbls01,lbls02),]
	M1$label[M1$label %in% lbls01] <- b[1]
	M1$label[M1$label %in% lbls02] <- b[2]
	M1$label <- factor(M1$label)
	observed <- data.frame("sampleId"=rownames(M1),"label"=M1$label)

	result <- list()
	#result$Matrix <- M1

	fullRF <- randomForest(label ~ . ,data=M1,importance=T,ntree=1000,nodesize=12,na.action=na.exclude)

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
        model <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
	result$model <- model
        
	predictions=as.data.frame(model$votes)
	clss <- colnames(predictions)
	predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	predictions$observed <- M1$label
	tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
	result$AUC <- tmpROC$auc
	result
}

naiveRF <-  function(a,b,features,M0,it=100){
	
	lbls01 <- strsplit(a[1],"")[[1]]
	lbls02 <- strsplit(a[2],"")[[1]]
	M1 <- M0[M0$label %in% c(lbls01,lbls02),c("label",features)]
	M1$label[M1$label %in% lbls01] <- b[1]
	M1$label[M1$label %in% lbls02] <- b[2]
	M1$label <- factor(M1$label)


	tmpRF <- randomForest(label ~ . ,data=M1,importance=T,ntree=500,nodesize=5)
	predictions=as.data.frame(tmpRF$votes)
	clss <- colnames(predictions)
	predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	predictions$observed <- M1$label
	tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
	list("model"=tmpRF,"AUC"=tmpROC$auc)
}

#get AUC value of  specific features
featureAUC <- function(features,label,M0){
	tmp <- M0[M0$label %in% label,c("label",features)]
	tmp$label <- as.factor(tmp$label);
	tmp[is.na(tmp)] <- 0
	tmpRF <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
	predictions=as.data.frame(tmpRF$votes)
	clss <- colnames(predictions)
	predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	predictions$observed <- tmp$label

	tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
	tmpROC$auc
}

singleFeatureKNN <- function(feature,label,M0,K=5){
 tmp <- M0[M0$label %in% label,c("label",feature)]
 dm <- as.matrix(dist(tmp[,feature]))
 L=dim(dm)[1]
 for(i in 1:L){
   dm[i,i] <- .Machine$double.xmax
 }

 predictions <- apply(dm,2,function(v){
    kNNLabels <- tmp$label[order(v,decreasing = F)[1:K]]
    votes <-sapply(label,function(lbl){sum(kNNLabels==lbl)})
    names(votes)[votes==max(votes)]
 })
   
   if(any(colnames(tmp)=="label")) {
     indx01 <- which(tmp$label==label[1])
     indx02 <- which(tmp$label==label[2])
     indxT  <- which(tmp$label==predictions)
     performance <- c(length(intersect(indx01,indxT))/length(indx01),length(intersect(indx02,indxT))/length(indx02))
     names(performance) <- label
  }
 list("yes"=indxT,"no"=which(tmp$label!=predictions),"performance"=performance,"prediction"=predictions)
}


### computing the complementarity between two features
## features : string vector containing the names of two features
## label: class labels of the two classes  
complementarity <- function(feature1,feature2,label,M0,K=5){
   result <- 0
   if(feature1!=feature2){
       p1 <- singleFeatureKNN(feature1,label,M0,K)
       p2 <- singleFeatureKNN(feature2,label,M0,K)
       result <- (length(intersect(p1$yes,p2$no))/(1+length(p2$no))+length(intersect(p2$yes,p1$no))/(1+length(p1$no)))
       #unionYes <- c(p1$yes,p2$yes)
       #interYes <- intersect(p1$yes,p2$yes)
       #result <- length(setdiff(unionYes,interYes))/(length(p1$yes)+length(p2$no))
   }
   result
}

multiFeatureKNN <- function(features,label,M0){
   tmp <- M0[M0$label %in% label,]
   result <- list()
   votes <- matrix(0,nrow=dim(tmp)[1],ncol=length(label))
   colnames(votes) <- label

   for(f0 in features){
      sg <- singleFeatureKNN(f0,label,tmp);
      pm <- matrix(0,nrow=dim(tmp)[1],ncol=length(label))
      colnames(pm) <- label
     
      for(ch in label){
          pm[sg$prediction==ch,ch] <- sg$performance[ch]
      }
      votes <- votes+pm
  }
  votes <- data.frame(t(apply(votes,1,function(v){v/sum(v)})))
  names(votes) <- label
  votes$prediction <- apply(votes,1,function(v){names(v)[max(v)==v]})
  result$prediction <- votes$prediction

  if(any(colnames(tmp)=="label")) {
     indx01 <- which(tmp$label==label[1])
     indx02 <- which(tmp$label==label[2])
     indxT  <- which(tmp$label==votes$prediction)
     result$performance <- c(length(intersect(indx01,indxT))/length(indx01),length(intersect(indx02,indxT))/length(indx02))
  }
  result
}

 trainRF <- function(M0,it=100){
        M1 <- M0
	observed <- data.frame("sampleId"=rownames(M1),"label"=M1$label)

	result <- list()
	#result$Matrix <- M1
        M1$label <- factor(M1$label)
	fullRF <- randomForest(label ~ . ,data=M1,importance=T,ntree=1000,nodesize=12,na.action=na.exclude)

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
		#currentFeatures <- c("label",orderedFeatures[1:indx])
		currentFeatures <- c("label",orderedFeatures[c(featureIndex,indx)])
		tmp <- M0[,currentFeatures]
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
			print(sprintf("find better forest with score: %f and %d features",bestScore,indx))
		}

                #print(indx)
		indx <- indx+1
	}##while()
	selected <- orderedFeatures[featureIndex]
	result$protIds <- selected
        #print(sprintf("finished features selection and %d proteins were selected",length(selected)))
	tmp <- M0[,c('label',selected)]
	tmp$label <- factor(tmp$label);
        model <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
	result$model <- model
	predictions=as.data.frame(model$votes)
	clss <- colnames(predictions)
	predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	predictions$observed <- M0$label
	tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
	result$AUC <- tmpROC$auc
	result
}

growRF01 <- function(a,M0,it=100){
	lbls01 <- strsplit(a,"")[[1]]
	print(lbls01)
	M1 <- M0[M0$label %in% lbls01,]
	
	M1$label <- factor(M1$label)
	observed <- data.frame("sampleId"=rownames(M1),"label"=M1$label)

	result <- list()
	#result$Matrix <- M1
	fullRF <- randomForest(label ~ . ,data=M1,importance=T,ntree=1000,nodesize=12,na.action=na.exclude)

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
	print(selected)
        #print(sprintf("finished features selection and %d proteins were selected",length(selected)))
	tmp <- M1[,c('label',selected)]
	tmp$label <- as.factor(tmp$label);
        model <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
	result$model <- model
        
	predictions=as.data.frame(model$votes)
	clss <- colnames(predictions)
	predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
	predictions$observed <- M1$label
	tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
	result$AUC <- tmpROC$auc
	result
}



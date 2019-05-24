remove(list=ls())
library("randomForest")
library("pROC")
set.seed(20181119)
source("E:/projects/TPD/src/growRF.R")
options(width=350)

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
rownames(SWT) <- SWT$patientId
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
trainM <- SWT[,c('patientId',clnames)]
trainM$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
trainM <- trainM[,c("label",clnames)]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
testM <- SWV[,clnames]
protIds <- intersect(colnames(trainM),colnames(testM))

trainM <- trainM[,c("label",protIds)]
trainM[,-1] <- t(apply(trainM[,-1],1,function(v){v/mean(v,na.rm=T)}))

testM <- testM[,protIds]
testM <- data.frame(t(apply(testM,1,function(v){v/mean(v,na.rm=T)})),check.names=F)
trainM[is.na(trainM)] <- 0
testM[is.na(testM)] <- 0
#############################################
M0 <- trainM[trainM$label=='M',];
A0 <- trainM[trainM$label=='A',];
C0 <- trainM[trainM$label=='C',];
P0 <- trainM[trainM$label=='P',];

distMA <- data.frame(as.matrix(dist(rbind(M0[,-1],A0[,-1])))[1:dim(M0)[1],(dim(M0)[1]+1):(dim(M0)[1]+dim(A0)[1])])
distAC <- data.frame(as.matrix(dist(rbind(A0[,-1],C0[,-1])))[1:dim(A0)[1],(dim(A0)[1]+1):(dim(A0)[1]+dim(C0)[1])])
distCP <- data.frame(as.matrix(dist(rbind(C0[,-1],P0[,-1])))[1:dim(C0)[1],(dim(C0)[1]+1):(dim(C0)[1]+dim(P0)[1])])

splitGraph <- function(dMat,pointNumbers){
    L1 <- dim(dMat)[1]
    L2 <- dim(dMat)[2]
    
    t1 <- t(apply(dMat,1,function(v){c(names(v[order(v)[1]]),v[order(v)[1]])}))
    t2 <- data.frame("partA"=rownames(t1),"partB"=t1[,1],"distance"=as.numeric(t1[,2]),stringsAsFactors = F)
    t2 <- t2[order(t2$distance),]

    partA <- c()
    partB <- c()
    index <- 1
    while( length(partA)< pointNumbers[1] & index<L1) {
          if(length(partA)<pointNumbers[1]){partA <- unique(c(partA,t2[index,"partA"]))}
      	  index <- index+1
       
    }
    
    t1 <- t(apply(t(dMat),1,function(v){c(names(v[order(v)[1]]),v[order(v)[1]])}))
    t2 <- data.frame("partB"=rownames(t1),"partA"=t1[,1],"distance"=as.numeric(t1[,2]),stringsAsFactors = F)
    t2 <- t2[order(t2$distance),]
    index <- 1
    while( length(partB)< pointNumbers[2] & index<L2){
          if(length(partB)<pointNumbers[2]){partB <- unique(c(partB,t2[index,"partB"]))}
	  index <- index+1
       
    }

    list("partA"=partA,"partB"=partB)
}

splitMA <- splitGraph(distMA,c(61,43))
splitAC <- splitGraph(distAC,c(0,29))
splitCP <- splitGraph(distCP,c(0,53))

#############################################################################################################

M0 <- trainM[setdiff(rownames(trainM[trainM$label=='M',]),splitMA$partA),];
M1 <- trainM[splitMA$partA,];

A0 <- trainM[splitMA$partB,];
A1 <- trainM[setdiff(rownames(trainM[trainM$label=='A',]),splitMA$partB),];

t1 <- read.table("E:/projects/TPD/data/C_patient_subtypes.txt",sep="\t",header=T,stringsAsFactors=F)
t1 <- as.character(sapply(t1$patientId[t1$subtype=='mC'],function(v){paste0('C',v)}))
rnames <- intersect(t1,splitAC$partB)
## "C13" "C19" "C2"  "C22" "C24" "C25" "C26" "C29" "C33" "C34" 
## "C36" "C37" "C39" "C4"  "C41" "C45" "C47" "C7"  "C8" 
C0 <- trainM[rnames,];
C1 <- trainM[setdiff(rownames(trainM[trainM$label=='C',]),rnames),];

P0 <- trainM[splitCP$partB,];
P1 <- trainM[setdiff(rownames(trainM[trainM$label=='P',]),splitCP$partB),];

M0$label <- rep("M0",dim(M0)[1]);  M1$label <- rep("M1",dim(M1)[1]);
A0$label <- rep("A0",dim(A0)[1]);  A1$label <- rep("A1",dim(A1)[1]);
C0$label <- rep("C0",dim(C0)[1]);  C1$label <- rep("C1",dim(C1)[1]);
P0$label <- rep("P0",dim(P0)[1]);  P1$label <- rep("P1",dim(P1)[1])


trainRF <- function(M0,it=100){
        M1 <- M0
	observed <- data.frame("sampleId"=rownames(M1),"label"=M1$label)

	result <- list()
	#result$Matrix <- M1
        M1$label <- factor(M1$label)
	lbls <- paste0(unique(M1$label),collapse="")
	fullRF <- randomForest(label ~ . ,data=M1,importance=T,ntree=1000,nodesize=12,na.action=na.exclude)

	imps  <- data.frame(importance(fullRF));
	impScore <- imps$MeanDecreaseAccuracy
	imps <- imps[order(impScore,decreasing=T),]
	orderedFeatures <- rownames(imps)

	indx <- 1
	featureIndex <- c()
	bestScore <- 0
	bestModel <- NULL
	bestROC <- NULL
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
			bestModel <- tmpRF
			bestROC <- tmpROC
			print(sprintf("%s find better forest with score: %f when adding %s",lbls,bestScore,orderedFeatures[indx]))
		}

                #print(indx)
		indx <- indx+1
	}##while()
	selected <- orderedFeatures[featureIndex]
	cat("","\n")
	print(sprintf("---------------------------------%d features found  for %s --------------------------------------------------------------",length(selected),lbls))
	print(selected)
	cat("","\n")
	result$protIds <- selected
	result$model <- bestModel
	result$AUC <- bestROC
	result
}

###bestRF: training RF n times and get the best

bestRF <- function(M0,it=100,ntimes=100){
  auc=0;
  bestModel <- NULL
  nFeatures <- dim(M0)[2]

  for(i in 1: ntimes){
    model <- trainRF(M0,it)
    if(model$AUC$auc>auc){
       auc <- model$AUC$auc
       bestModel <- model
    }
  }
  bestModel
}


rfM01    <- bestRF(rbind(M0,M1));
rfM0A0   <- bestRF(rbind(M0,A0));
rfM0A1   <- bestRF(rbind(M0,A1));
rfM1A0   <- bestRF(rbind(M1,A0));
rfM1A1   <- bestRF(rbind(M1,A1));
rfA01    <- bestRF(rbind(A0,A1));
rfA0C0   <- bestRF(rbind(A0,C0));
rfA0C1   <- bestRF(rbind(A0,C1));
rfA1C0   <- bestRF(rbind(A1,C0));
rfA1C1   <- bestRF(rbind(A1,C1));
rfC01    <- bestRF(rbind(C0,C1));
rfC0P0   <- bestRF(rbind(C0,P0));
rfC0P1   <- bestRF(rbind(C0,P1)); R3'
·

rfC1P0   <- bestRF(rbind(C1,P0));
rfC1P1	 <- bestRF(rbind(C1,P1));

rfNM     <- bestRF(trainM[trainM$label %in% c('N','M'),])
rfMA     <- bestRF(trainM[trainM$label %in% c('M','A'),])
rfAC     <- bestRF(trainM[trainM$label %in% c('A','C'),])
rfCP     <- bestRF(trainM[trainM$label %in% c('C','P'),])

models <- list()
modelNames <- ls()[grepl(pattern = "rf",x=ls())]
for(nm in modelNames){
  eval(parse(text=paste0("models[[length(models)+1]] <-",nm)))
  eval(parse(text=sprintf('remove(%s)',nm)))
}
names(models) <- modelNames

###
TM <- rbind(rbind(rbind(M0,M1),rbind(A0,A1)),rbind(rbind(C0,C1),rbind(P0,P1)))
TM$label <- factor(TM$label,ordered = TRUE,levels = c('M0', 'M1','A0','A1','C0','C1','P0','P1'))
TM <- TM[,c("label",features)]

pids <- unique(c(models$rfM0A1$protIds,models$rfM1A1$protIds))
MA0 <- TM[TM$label %in% c('M0','M1','A1'),c("label",pids)]
MA0$label <- factor(MA0$label,ordered=TRUE,levels = c('M0','M1','A1'))
OF1 <- cforest(label ~ ., data <- MA0,scores=list(label=c(1,2,13)),control = cforest_unbiased(ntree = 1000))
vd1 <- predict(OF1,newdata=TM[TM$label=='A0',],OOB=TRUE,type="response")

pids <- unique(c(models$rfM0A1$protIds,models$rfM1A1$protIds))




ordinalRF <- cforest(label ~ ., data <-TM ,scores = list(label = c(1,2,5,6,8,9,12,13)), control = cforest_unbiased(ntree = 1000))
validate1<- predict(ordinalRF,newdata =SWV, OOB=TRUE, type = "response")

remove(list=ls())
library("randomForest")
library("pROC")
set.seed(20181119)
source("E:/projects/TPD/src/growRF.R")

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
SWT <- SWT[,c("patientId",names(missingRate)[missingRate<tradeOff])]
SWT$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
SWT <- SWT[,c("label",names(missingRate)[missingRate<tradeOff])]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
SWV <- SWV[,names(missingRate)[missingRate<tradeOff]]

protIds <- intersect(colnames(SWT),colnames(SWV))
SWT <- SWT[,c("label",protIds)]
SWT[,-1] <- t(apply(SWT[,-1],1,function(v){v/mean(v,na.rm=T)}))
SWV <- SWV[,protIds]
SWV <- data.frame(t(apply(SWV,1,function(v){v/mean(v,na.rm=T)})),check.names=F)
SWT[is.na(SWT)] <- 0
SWV[is.na(SWV)] <- 0
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
protIds <- intersect(colnames(SWT),colnames(SWV))


###########################################################

Labels <- c('NM','MA','AC','CP')
nFolds=10
Models <- list()

for(i in 1:4){
 	label <- c(substr(Labels[i],1,1),substr(Labels[i],2,2))
	bestSample=c()
	bestTest=c()
	bestModel=list()
	bestFeatures=c()
	bestPerformance <- 0;


        for(k in 1:100){
	     indexA <- which(SWT$label==label[1])
	     indx <- cut(1:length(index),breaks=nFolds,labels=F)
	     foldsA <- sapply(1:nFolds,function(i){list(indexA[which(indx==i)])})
		
	     indexB <- which(SWT$label==label[2])
	     indx <- cut(1:length(index),breaks=nFolds,labels=F)
	     foldsB <- sapply(1:nFolds,function(i){list(indexB[which(indx==i)])})

             for(m in 1:nFolds){
	       for(n in 1:nfolds){
                  testIndex <- c(foldsA[[m]],foldsB[[n]])
		  testM <-  SWT[testIndex,]
		  trainM <- SWT[setdiff(c(indexA,indexB),testIndex),]


		  RF1 <- growRF(label,label,trainM)
		  RF2 <- growRF(label,label,testM)

		  prediction1 <-  predict(RF1$model,testM)
		  prediction2 <-  predict(RF2$model,trainM)
		  percent <- (sum(prediction1==testM$label)+sum(prediction2==trainM$label))/(dim(testM)[1]+dim(trainM)[1])
		
		#percents <- c(sum(prediction1==testM$label)/55,sum(prediction2==testM$label)/55)
		print(sprintf("%s%s %d: %f",label[1],label[2],k,percent))
		if(percent > bestPerformance){
		   bestModel$testM  <- testM
		   bestModel$trainM <- trainM
		   bestModel$RF1    <- RF1
		   bestModel$RF2    <- RF2
                   bestModel$accuracy <- percent
		   bestPerformance <- percent
		}
        }##for
        Models[[length(Models)+1]] <- bestModel
    }
  }

names(Models) <- c("MA","MC","MP","AC","AP","CP")

########################### predict with validation dataset

predictions <- sapply(c("MA","AC","CP"),function(lbl){
	RFs <- Models[[lbl]];
	prediction1 <- as.matrix(predict(RFs$RF1$model,SWV,type="prob"))
	prediction1 <- prediction1 *(1.0 - RFs$RF1$model$confusion[,3])
	prediction2 <- as.matrix(predict(RFs$RF2$model,SWV,type="prob"))
	prediction2 <- prediction2 *(1.0 - RFs$RF2$model$confusion[,3])
        t0 <- t(apply(prediction1 + prediction2,1,function(v){v/sum(v)}))
});
pd <- t(sapply(1:181,function(i){
    p1 <- predictions[i,]
    p2 <- predictions[i+181,]
    p0 <- t(rbind(p1,p2))
    p0
}))
pd <- pd[,c(4,1,2,5,3,6)]

colnames(pd) <- c("M","A","A","C","C","P")

decision <- data.frame(t(apply(pd,1,function(v){
  result="Unknown"
  if(v[1]>v[2]){
    result="M"
  }else{
    if(v[3]>v[4]){
      result='A'
    }else{
       if(v[5]>v[6]){
         result='C'
       }else{
         result='P'
       }
    }
  }
  if(result=='M' | result=='A'){  result <- c(result,'benign')}
  else {result <- c(result,'malignant')}
  result
})))
decision$predicted <- apply(decision,1,function(v){names(v)[max(v)==v]})






       
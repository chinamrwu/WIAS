remove(list=ls())
library("randomForest")
library("pROC")
set.seed(20181119)
source("/work/src/growRF.R")
options(width=350)

SWT <- read.table("/work/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("/work/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
missingMarks <- c()
missingMarks <- c('P24821','P09668','Q9NR28','O95861','Q8N6C5')

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
clnames <- unique(c(names(missingRate)[missingRate<tradeOff],missingMarks))
SWT <- SWT[,c('patientId',clnames)]
SWT$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
SWT <- SWT[,c("label",clnames)]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
clnames <- unique(c(names(missingRate)[missingRate<tradeOff],missingMarks))
SWV <- SWV[,clnames]
protIds <- intersect(colnames(SWT),colnames(SWV))

SWT <- SWT[,c("label",protIds)]
SWT[,-1] <- t(apply(SWT[,-1],1,function(v){v/mean(v,na.rm=T)}))
SWV <- SWV[,protIds]
SWV <- data.frame(t(apply(SWV,1,function(v){v/mean(v,na.rm=T)})),check.names=F)
SWT[is.na(SWT)] <- 0
SWV[is.na(SWV)] <- 0



###########################################################

Labels <- "NMACP"
testNumber <- c(61,43,22,53)
names(testNumber) <- c('M','A','C','P')

LbLs <- c('MA','AC','CP')

Models <- list()

for(i in 1:3){
	label <- c(substr(LbLs[i],1,1),substr(LbLs[i],2,2))
	bestSample=c()
	bestTest=c()
	bestModel=list()
	bestFeatures=c()
	bestPerformance <- 0;


        for(k in 1:100){
		index <- which(SWT$label %in% label)
		indexA <- sample(which(SWT$label==label[1]),testNumber[label[1]],replace=F)
		indexC <- sample(which(SWT$label==label[2]),testNumber[label[2]],replace=F)

		testM <-  SWT[c(indexA,indexC),]
		trainM <- SWT[setdiff(index,c(indexA,indexC)),]


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

names(Models) <- LbLs

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
    p0 <- rbind(p1,p2)
    p0
}))
pd <- pd[,c(4,1,2,5,3,6)]

colnames(pd) <- c("N","M","M","A","A","C","C","P")

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






       
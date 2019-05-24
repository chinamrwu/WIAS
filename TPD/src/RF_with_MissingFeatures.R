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
###########################  without deaf features    #########################################
Models <- list()

Models$MA <- growRF(c('M','A'),c('M','A'),trainM)
Models$AC <- growRF(c('A','C'),c('A','C'),trainM)
Models$CP <- growRF(c('C','P'),c('C','P'),trainM)
#Models$CP <- growRF(c('A','P'),c('A','P'),trainM)
#Models$CP <- growRF(c('M','P'),c('M','P'),trainM)


predictions01 <- sapply(c("MA","AC","CP"),function(lbl){
	RF <- Models[[lbl]];
	predict01 <- as.matrix(predict(RF$model,testM,type="prob"))
	predict01
});
pd01 <- t(sapply(1:181,function(i){
    p1 <- predictions01[i,]
    p2 <- predictions01[i+181,]
    p0 <- rbind(p1,p2)
    p0
}))
colnames(pd01) <- c('A','M','A','C','C','P')
pd01 <- pd01[,c(2,1,3,4,5,6)]
############################################# adding deaf features ####################################
effecientProtIds <- unique(c(Models[['MA']]$protIds,Models[['AC']]$protIds,Models[['CP']]$protIds))
deafMarks <- c('P24821','P09668','Q9NR28','O95861','Q8N6C5')

clnames <- c(effecientProtIds,deafMarks)

trainM01 <- SWT[,clnames]
trainM01$label <- trainM$label

trainM01 <- trainM01[,c("label",clnames)]
trainM01[,-1] <- t(apply(trainM01[,-1],1,function(v){v/mean(v,na.rm=T)}))
trainM01[is.na(trainM01)] <- 0

testM01 <- SWV[,clnames]
testM01 <- data.frame(t(apply(testM01,1,function(v){v/mean(v,na.rm=T)})),check.names=F)
testM01[is.na(testM01)] <- 0

Models01 <- list()

M0 <- trainM01[trainM01$label %in% c('M','A'),]
M0$label <- factor(M0$label)
Models01$MA <- randomForest(label ~ . ,data=M0,importance=T,ntree=1000,nodesize=5)

M0 <- trainM01[trainM01$label %in% c('A','C'),]
M0$label <- factor(M0$label)
Models01$AC <- randomForest(label ~ . ,data=M0,importance=T,ntree=1000,nodesize=5)

M0 <- trainM01[trainM01$label %in% c('C','P'),]
M0$label <- factor(M0$label)
Models01$CP <- randomForest(label ~ . ,data=M0,importance=T,ntree=1000,nodesize=5)

predictions02 <- sapply(c("MA","AC","CP"),function(lbl){
	model <- Models01[[lbl]];
	predict02 <- as.matrix(predict(model,testM01,type="prob"))
	predict02
});
pd02 <- t(sapply(1:181,function(i){
    p1 <- predictions02[i,]
    p2 <- predictions02[i+181,]
    p0 <- rbind(p1,p2)
    p0
}))
colnames(pd02) <- c('A','M','A','C','C','P')
pd02 <- pd02[,c(2,1,3,4,5,6)]

###########################################################
decide <- function(v){
  result='Unknown'
  p1 <- v[1]>=v[2]
  p2 <- v[3]>=v[4]
  p3 <- v[5]>=v[6]
  if(!p1 & !p2 & !p3) result <- 'P'#1
  if(!p1 & !p2 & p3)  result <- 'C'#2
  if(!p1 &  p2 & !p3) result <- 'A'#3
  if(p1 &  !p2 & !p3) result <- 'P'#4
  if(!p1 &  p2 & p3) result <- 'A'#5
  if(p1 &  !p2 & p3) result <- 'C'#6
  if(p1 &  p2 & !p3) result <- 'M'#7
  if(p1 &  p2 &  p3) result <- 'M'#7
  result
}
votes <- function(v){
   names(v) <- c('M','A','A','C','C','P')
   vote <- rep(0,6)
   names(vote) <- c('M','A','A','C','C','P')
   ifelse(v[1]>v[2],vote[1:2] <- c(1,-1),vote[1:2] <- c(-1,1))
   ifelse(v[3]>v[4],vote[3:4] <- c(1,-1),vote[3:4] <- c(-1,1))
   ifelse(v[5]>v[6],vote[5:6] <- c(1,-1),vote[5:6] <- c(-1,1))
   vote <- c(vote[1],sum(vote[2:3]),sum(vote[4:5]),vote[6])
   names(vote) <- c('M','A','C','P')
   result <- names(vote)[max(vote)==vote][1]
   if(length(result)>1){ ifelse(v[1]>v[6],result<-'M',result <- 'P')}
   result
}

######################################
predicted <- data.frame("patientId"=1:dim(pd01)[1],"RF"=apply(pd01,1,decide),"RFWithMissed"=apply(pd02,1,decide))
write.table(predicted,file="/work/predicted_RF_with_MissingFeatures.txt",sep="\t",col.names=T,row.names=F,quote=F)







       
remove(list=ls())
library("sda")
library("pROC")
library("data.table")
set.seed(20181203)
source("E:/projects/TPD/src/growRF.R")
##http://strimmerlab.org/software/sda/download/sda-khan-data.pdf
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
SWT[,-1] <- t(apply(SWT[,-1],2,function(v){m1 <- mean(v,na.rm=T);dv <- sd(v,na.rm=T); v-m1/dv}))
SWV <- SWV[,protIds]
SWV <- t(apply(SWV,2,function(v){m1 <- mean(v,na.rm=T);dv <- sd(v,na.rm=T); v-m1/dv}))
SWT[is.na(SWT)] <- 0
SWV[is.na(SWV)] <- 0
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
protIds <- intersect(colnames(SWT),colnames(SWV))
#############################################################################################################################
predfun01 <-  function(Xtrain, Ytrain, Xtest, Ytest, numVars, diagonal=FALSE){
	# estimate ranking and determine the best numVars variables
	ra = sda.ranking(Xtrain, Ytrain, verbose=FALSE, diagonal=diagonal, fdr=FALSE)
	selVars = ra[,"idx"][1:numVars]
	# fit and predict
	sda.out = sda(Xtrain[, selVars, drop=FALSE], Ytrain, diagonal=diagonal, verbose=FALSE)
	ynew = predict(sda.out, Xtest[, selVars, drop=FALSE], verbose=FALSE)$class
	# count false and true positives/negatives
	negative = levels(Ytrain)[2] # "healthy"
	cm = confusionMatrix(Ytest, ynew, negative=negative)
	return(cm)
}

predfun02 <- function(Xtrain, Ytrain, Xtest, Ytest, numVars, diagonal=FALSE,ranking.score="entropy"){
	# estimate ranking and determine the best numVars variables
	ra = sda.ranking(Xtrain, Ytrain, verbose=FALSE, diagonal=diagonal,fdr=FALSE, ranking.score=ranking.score)
	selVars = ra[,"idx"][1:numVars]
	# fit and predict
	sda.out = sda(Xtrain[, selVars, drop=FALSE], Ytrain, diagonal=diagonal,verbose=FALSE)
	ynew = predict(sda.out, Xtest[, selVars, drop=FALSE], verbose=FALSE)$class
	# compute accuracy
	acc = mean(Ytest == ynew)

	negative = levels(Ytrain)[2] # "healthy"
	cm = confusionMatrix(Ytest, ynew, negative=negative)
	return(c(acc,cm))
}


library("crossval")
nFolds = 10 # number of folds
nRepeat = 20 # number of repetitions
set.seed(12345)
#trainX <- as.matrix(SWT[SWT$label %in% c('N','M'),-1])
#trainY <- factor(SWT$label[SWT$label %in% c('N','M')])


#for(varNumbers in 121:200){
#  cv.LDA10= crossval(predfun, trainX, trainY, K=nFolds, B=nRepeat, numVars=varNumbers, diagonal=FALSE, verbose=FALSE)
#  print(c(varNumbers,diagnosticErrors(cv.LDA10$stat),cv.LDA10$stat))
#}

trainX <- as.matrix(SWT[,-1])
trainY <- factor(SWT$label)

rankAvg <- sda.ranking(trainX, trainY, fdr=TRUE, plot.fdr=F, ranking.score="avg")
rankMax <- sda.ranking(trainX, trainY, fdr=TRUE, plot.fdr=F, ranking.score="max")
rankEntropy <- sda.ranking(trainX, trainY, fdr=TRUE, plot.fdr=F, ranking.score="entropy")
sum(rankAvg[,'lfdr']< 0.8)
sum(rankMax[,'lfdr']< 0.8)
sum(rankEntropy[,'lfdr']< 0.8)

cv.Entropy <- crossval(predfun02, trainX, trainY, K=nFolds, B=nRepeat, numVars=100, diagonal=FALSE, verbose=FALSE)
cv.Avg     <- crossval(predfun02, trainX, trainY, K=nFolds, B=nRepeat, numVars=9, diagonal=FALSE, verbose=FALSE,ranking.score='avg')
cv.Max     <- crossval(predfun02, trainX, trainY, K=nFolds, B=nRepeat, numVars=5, diagonal=FALSE, verbose=FALSE,ranking.score='max')

diagnosticErrors(cv.Entropy$stat)
diagnosticErrors(cv.Avg$stat)
diagnosticErrors(cv.Max$stat)


#########################


Xtrain
rankEntropy <- sda.ranking(trainX, trainY, fdr=TRUE, plot.fdr=F, ranking.score="entropy")
selVars = rankEntropy[,"idx"][1:19]
model = sda(Xtrain[, selVars, drop=FALSE], Ytrain, diagonal=diagonal,verbose=FALSE)

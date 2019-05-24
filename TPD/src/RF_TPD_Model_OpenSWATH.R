
###### Analysis of the output of RF
remove(list=ls())
library("randomForest")
require(pROC)

trainM <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_20181022.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(trainM) <- trainM[,1]
trainM <- trainM[,-1]

set.seed(1978)

uniqLabels =unique(trainM$label)
k=length(uniqLabels)
for(i in 1:(k-1)){
   for(j in (i+1):k){
   tmp <- trainM[trainM$label==uniqLabels[i] | trainM$label==uniqLabels[j],]
   tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)
   assign(paste0("model.",uniqLabels[i],uniqLabels[j]),tmpRF)
 }}

#plot(rt.NC.ROC,legend.title="NC")
#auc(rt.NC.ROC)

#RF.output <- randomForest(formula=as.factor(trainM$label) ~ . ,data=trainM,importance=T)
for(i in 1:(k-1)){
   for(j in (i+1):k){
    txt=paste0("print(model.",uniqLabels[i],uniqLabels[j],")")
    eval(parse(text=txt))
   }
}

uniqLabels =unique(trainM$label)

RFList <- list()
for(i in 1:(k-1)){
   for(j in (i+1):k){
     lb=paste0(uniqLabels[i],uniqLabels[j])
     txt=paste0("RFList[[length(RFList)+1]]"," <- model.",lb)
     eval(parse(text=txt))
}}

lbs <-c()
for(i in 1:(k-1)){
   for(j in (i+1):k){
     lbs <- c(lbs,paste0(uniqLabels[i],uniqLabels[j]))
   }
}
names(RFList) <- lbs
############################################################ Draw ROC Plot
ROCList <- list()
for(obj in RFList){
   predictions=as.data.frame(obj$votes)
   clss <- colnames(predictions)
   predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
   predictions$observed <- as.character(sapply(rownames(predictions),function(rn){substr(rn,1,1)}))
   ROCList[[length(ROCList)+1]] <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
}



### one plot per classifier
#mt <- matrix(c(1:10,0,0),nrow=3,ncol=4)
#layout(mt)
for(i in 1:10){
	nam<-names(RFList)[i]
	pdf(paste0("E:/projects/TPD/results/RF_TPD_openSWATH_ROC_",nam,".pdf"))
	plot.roc(ROCList[[i]],print.auc=T,col = "blue3",print.thres="best",main=names(RFList)[i],legacy.axes = TRUE,print.auc.cex=1.2)
	dev.off()
}

pdf("E:/projects/TPD/results/RF_TPD_openSWATH_ROC_All_in_ONE.pdf")
mt <- matrix(c(1:10,0,0),nrow=4,ncol=3)
layout(mt)
for(i in 1:10){
    plot.roc(ROCList[[i]],print.auc=T,col = "blue3",print.thres="best",main=names(RFList)[i],legacy.axes = TRUE,print.auc.cex=1.2)
}
dev.off()

###### find important proteins
importantProts <- list()
importantMatrix <- c();
indx=1
for( obj in RFList){
     imps  <- data.frame(importance(obj));
     score <- imps$MeanDecreaseAccuracy * imps$MeanDecreaseGini

     imps <- imps[order(score,decreasing=T),]
     imps$protId <- rownames(imps)
     imps$score <- score[order(score,decreasing=T)]
     imps <-  imps[1:50,c("protId","score")];
     imps$label <- rep(names(RFList)[indx],50)
     indx <- indx+1
     importantProts[[length(importantProts)+1]] <- imps
     importantMatrix <- rbind(importantMatrix,imps)
}
rownames(importantMatrix) <- 1:dim(importantMatrix)[1]
write.table(importantMatrix,file="E:/projects/TPD/results/RF_TPD_openSWATH_top50_important_proteins.txt",sep="\t",col.names=T,row.names=F,quote=F)

#names(impVars) <- names(RFList)

#find common 
cmProts <-c()
for(v in impVars){
  cmProts <- c(cmProts,v$protId)
}

coreProts <-unique(cmProts)

for(obj in impVars){
  coreProts <- intersect(obj$protId,coreProts)
  print(coreProts)
}

overlaps <- matrix(nrow=10,ncol=10,data=NA)
k=length(impVars)
for(i in 1:(k-1)){
  for(j in (i+1):k){
   overlaps[i,j] <- length(intersect(impVars[[i]]$protId,impVars[[j]]$protId))
  }
}
colnames(overlaps) <- names(RFList)
rownames(overlaps) <- names(RFList)
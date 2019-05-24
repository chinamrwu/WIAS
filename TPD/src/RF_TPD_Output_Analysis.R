
###### Analysis of the output of RF

trainM <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_20181022.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(trainM) <- trainM[,1]
trainM <- trainM[,-1]

library("randomForest")
require(pROC)
set.seed(1978)

uniqLabels =unique(trainM$label)
k=length(uniqLabels)
for(i in 1:(k-1)){
   for(j in (i+1):k){
   tmp <- trainM[trainM$label==uniqLabels[i] | trainM$label==uniqLabels[j],]
   tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)
   assign(paste0("RF.",uniqLabels[i],uniqLabels[j]),tmpRF)
   roc <- roc(tmp$label,tmpRF$votes[,2])
   assign(paste0("ROC.",uniqLabels[i],uniqLabels[j]),roc)
 }}

#plot(rt.NC.ROC,legend.title="NC")
#auc(rt.NC.ROC)

#RF.output <- randomForest(formula=as.factor(trainM$label) ~ . ,data=trainM,importance=T)
for(i in 1:(k-1)){
   for(j in (i+1):k){
    txt=paste0("print(RF.",uniqLabels[i],uniqLabels[j],")")
    eval(parse(text=txt))
   }
}

uniqLabels =unique(trainM$label)

RFList <- list()
for(i in 1:(k-1)){
   for(j in (i+1):k){
     lb=paste0(uniqLabels[i],uniqLabels[j])
     txt=paste0("RFList[[length(RFList)+1]]"," <- RF.",lb)
     eval(parse(text=txt))
}}

lbs <-c()
for(i in 1:(k-1)){
   for(j in (i+1):k){
     lbs <- c(lbs,paste0(uniqLabels[i],uniqLabels[j]))
   }
}
names(RFList) <- lbs

impVars <- list()
for( obj in RFList){
     imps  <- data.frame(importance(obj));
     score <- imps$MeanDecreaseAccuracy * imps$MeanDecreaseGini

     imps <- imps[order(score,decreasing=T),]
     imps$protId <- rownames(imps)
     imps$score <- score[order(score,decreasing=T)]
     impVars[[length(impVars)+1]] <- imps[1:10,c("protId","score")]
}

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
remove(list=ls())
library("randomForest")
#From local directory,load the models,which are trainned on trainning steps

#testData <- read.table("E:/projects/TPD/data/RF_TPDV_openSWATH_20181025.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
#testData[is.na(testData)] <- 0;
testData <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_20181022.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
load(file = "E:/projects/TPD/data/RF_TPD_Models_SWATH_v3.Rdata")

models <- bestRFs
remove(bestRFs)


modelNames <- c()
for(model in models){
  modeNames <- c(modeNames,model$label)
}

names(models) <- modelNames

 



lbls <- strsplit("NMACP","")[[1]]
predictions <- data.frame(matrix(0,nrow=dim(testData)[1],ncol=5))
colnames(predictions) <- lbls
rownames(predictions) <- rownames(testData)

details=matrix(0,nrow=dim(testData)[1],ncol=1)

for(model in models){
   forest <- model$tree
   predicted <- as.data.frame(predict(forest, testData, type = "prob"))
   details <- cbind(details,predicted)
   predictions[,colnames(predicted)] <- predictions[,c(colnames(predicted))] + predicted
}

predictions <- data.frame(t(apply(predictions,1,function(v){v/sum(v)})),stringsAsFactors=F)
predictions <- cbind(details,predictions)
predictions <- predictions[,-1]
predictions$label <- testData$label
write.table(predictions,"E:/projects/TPD/results/test_SWATH_no_weight_details.txt",col.names=T,row.names=T,sep="\t",quote=F)

weightPredictions <- data.frame(matrix(0,nrow=dim(testData)[1],ncol=5))
colnames(weightPredictions) <- lbls
rownames(weightPredictions) <- rownames(testData)

for(model in models){
   forest <- model$tree
   weight <- 1- forest$confusion[,3]
   predicted <- as.data.frame(predict(forest, testData, type = "prob"))
   predicted[,1] <- predicted[,1]*weight[1]
   predicted[,2] <- predicted[,2]*weight[2]
   weightPredictions[,colnames(predicted)] <- weightPredictions[,c(colnames(predicted))] + predicted
}
weightPredictions <- data.frame(t(apply(weightPredictions,1,function(v){v/sum(v)})),stringsAsFactors=F)
weightPredictions$label <- testData$label
write.table(weightPredictions,"E:/projects/TPD/results/test_SWATH_weight.txt",col.names=T,row.names=T,sep="\t",quote=F)

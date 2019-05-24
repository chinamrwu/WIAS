remove(list=ls())
library("randomForest")
#From local directory,load the models,which are trainned on trainning steps

#testData <- read.table("E:/projects/TPD/data/RF_TPDV_expert_20181024.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
#testData[is.na(testData)] <- 0;
testData <- read.table("E:/projects/TPD/data/RF_TPDT_expert_trainingM_20181024.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")

load(file = "E:/projects/TPD/data/RF_TPD_Models_expert_v3.Rdata")
models <- bestRFs
remove(bestRFs)

lbls <- strsplit("NMACP","")[[1]]
predictions <- data.frame(matrix(0,nrow=dim(testData)[1],ncol=5))
colnames(predictions) <- lbls
rownames(predictions) <- rownames(testData)

for(model in models){
   forest <- model$tree
   predicted <- as.data.frame(predict(forest, testData, type = "prob"))
   predictions[,colnames(predicted)] <- predictions[,c(colnames(predicted))] + predicted
}

predictions <- data.frame(t(apply(predictions,1,function(v){v/sum(v)})),stringsAsFactors=F)
predictions$label <- testData$label
write.table(predictions,"E:/projects/TPD/results/test_expert_no_weight.txt",col.names=T,row.names=T,sep="\t",quote=F)





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
write.table(weightPredictions,"E:/projects/TPD/results/test_expert_weight.txt",col.names=T,row.names=T,sep="\t",quote=F)

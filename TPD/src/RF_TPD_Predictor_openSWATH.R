remove(list=ls())
library("randomForest")
#From local directory,load the models,which are trainned on trainning steps

testData <- read.table("E:/projects/TPD/data/RF_TPDV_openSWATH_20181025.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
testData[is.na(testData)] <- 0;

load(file = "E:/projects/TPD/data/RF_TPD_Models_SWATH_v3.Rdata")

models <- bestRFs
remove(bestRFs)

modelNames <- c()
for(model in models){
  modelNames <- c(modelNames,model$label)
}

names(models) <- modelNames


lbls <- strsplit("NMACP","")[[1]]
predictions <- data.frame(matrix(0,nrow=dim(testData)[1],ncol=5))
colnames(predictions) <- lbls
rownames(predictions) <- rownames(testData)

details=matrix(0,nrow=dim(testData)[1],ncol=1)
for(model in models){
   forest <- model$tree
   predicted <- as.matrix(predict(forest, testData, type = "prob"))
   details <- cbind(details,predicted)
   predictions[,colnames(predicted)] <- predictions[,c(colnames(predicted))] + predicted
}

predictions <- t(apply(predictions,1,function(v){v/sum(v)}))
predictions <- data.frame(predictions,stringsAsFactors=F)
predictions$predicted <- as.character(apply(predictions,1,function(v){names(v)[which(v==max(v))[1]]}))
predictions <- cbind(details,predictions)
predictions <- data.frame(predictions,stringsAsFactors=F)
write.table(predictions,"E:/projects/TPD/results/validation_SWATH_no_weight_details.txt",col.names=T,row.names=T,sep="\t",quote=F)

########################################################## weighted prediction
weightPredictions <- data.frame(matrix(0,nrow=dim(testData)[1],ncol=5))
colnames(weightPredictions) <- lbls
rownames(weightPredictions) <- rownames(testData)

details=matrix(0,nrow=dim(testData)[1],ncol=1)
for(model in models){
   forest <- model$tree
   weight <- data.frame(forest$confusion)

   predicted <- as.matrix(predict(forest, testData, type = "prob"))
   predicted[,1] <- predicted[,1]*(1-weight[colnames(predicted)[1],3])
   predicted[,2] <- predicted[,2]*(1-weight[colnames(predicted)[2],3])
   details <- cbind(details,predicted)
   weightPredictions[,colnames(predicted)] <- weightPredictions[,c(colnames(predicted))] + predicted
}
weightPredictions <- t(apply(weightPredictions,1,function(v){v/sum(v)}))
weightPredictions <- data.frame(weightPredictions,stringsAsFactors=F)
weightPredictions$predicted <- as.character(apply(weightPredictions,1,function(v){names(v)[which(v==max(v))[1]]}))
weightPredictions <- cbind(details,weightPredictions)
weightPredictions <- data.frame(weightPredictions,stringsAsFactors=F)
write.table(weightPredictions,"E:/projects/TPD/results/validation_SWATH_weight_details.txt",col.names=T,row.names=T,sep="\t",quote=F)

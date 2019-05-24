remove(list=ls())
library("randomForest")
#From local directory,load the models,which are trainned on trainning steps

testData <- read.table("E:/projects/TPD/data/RF_TPDV_expert_20181024.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
testData[is.na(testData)] <- 0;

load(file = "E:/projects/TPD/data/RF_TPD_Models_expert_v5.Rdata")

models <- bestRFs
remove(bestRFs)

lbls <- strsplit("NMACP","")[[1]]
predictions <- data.frame(matrix(0,nrow=dim(testData)[1],ncol=5))
colnames(predictions) <- lbls
rownames(predictions) <- rownames(testData)

for(model in models){
   forest <- model$forest
   predicted <- as.data.frame(predict(forest, testData, type = "prob"))
   predictions[,colnames(predicted)] <- predictions[,c(colnames(predicted))] + predicted
}

predictions <- t(apply(predictions,1,function(v){v/sum(v)}))
predictions <- data.frame(predictions,stringsAsFactors=F)
predictions$predicted <- as.character(apply(predictions,1,function(v){names(v)[which(v==max(v))[1]]}))

write.table(predictions,"E:/projects/TPD/results/validation_expert_no_weight_v5.txt",col.names=T,row.names=T,sep="\t",quote=F)

print("prediction without weight done")
########################################################## weighted prediction
weightPredictions <- data.frame(matrix(0,nrow=dim(testData)[1],ncol=5))
colnames(weightPredictions) <- lbls
rownames(weightPredictions) <- rownames(testData)

for(model in models){
   forest <- model$forest
   weight <- data.frame(forest$confusion)

   predicted <- as.data.frame(predict(forest, testData, type = "prob"))
   predicted[,1] <- predicted[,1]*(1-weight[colnames(predicted)[1],3])
   predicted[,2] <- predicted[,2]*(1-weight[colnames(predicted)[2],3])
   weightPredictions[,colnames(predicted)] <- weightPredictions[,c(colnames(predicted))] + predicted
}
weightPredictions <- t(apply(weightPredictions,1,function(v){v/sum(v)}))
weightPredictions <- data.frame(weightPredictions,stringsAsFactors=F)
weightPredictions$predicted <- as.character(apply(weightPredictions,1,function(v){names(v)[which(v==max(v))[1]]}))
write.table(weightPredictions,"E:/projects/TPD/results/validation_expert_weight_v5.txt",col.names=T,row.names=T,sep="\t",quote=F)

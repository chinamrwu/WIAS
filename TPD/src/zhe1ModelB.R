resps1 <- c()
for(model in models){
  feature    <- rownames(coef(model,s="lambda.1se"))[-1]
  zM <- scaleRow(zhe1[,feature])
  resps1 <-cbind(resps1, predict(model,newx = as.matrix(zM),s = model$lambda.min,type='response')[,1])
}

zhe1PredictB <- data.frame('patientId'=rownames(resps1),'M'=resps1[,2],'B'=1-resps1[,2],'observed'=as.character(sapply(rownames(resps1),function(v){ifelse(substr(v,1,1)=='A','B','M')})),
'predictd'=as.character(sapply(resps1[,2],function(v){ifelse(v>=0.5,'M','B')})),stringsAsFactors=F)

################

obj <- zhe1PredictB
ROC <- roc(ifelse(obj$observed=="M", "M", "B"), as.numeric(obj$M))
pdf("ROC_zy54_ModelB_20190424.pdf")
plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",	
  main="model B",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()

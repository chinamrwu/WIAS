vd0 <- read.csv("data/TPDZE_predict_190421_to_wu.csv",check.names = F,header = T,stringsAsFactors = F)
protIds <- as.character(sapply(colnames(vd0)[-1],function(v){a=strsplit(v,"\\|")[[1]][2]}))
colnames(vd0)[-1] <- protIds

k0 <- vd0[!grepl('pool',vd0$Sample),]
rownames(k0) <- k0$Sample
k1 <- k0$Sample[grepl('repB',k0$Sample)]

k2 <- sapply(k1,function(v){
    a <- strsplit(v,"_")[[1]][1]
    tmp <- apply(2^k0[c(a,v),-1],2,function(v1){log2(mean(v1,na.rm=T));})
    tmp
    })
k2 <- t(k2)
k2[is.na(k2)] <- NA
rownames(k2) <- sapply(rownames(k2),function(v){strsplit(v,"_")[[1]][1]})
k2 <- data.frame(k2)
R0 <- apply(k2,2,function(v){sum(is.na(v))})

zhe2 <- k2[,features]
#zhe2[is.na(zhe2)] <- 0

resps <- c()
for(model in models){
  feature    <- rownames(coef(model,s="lambda.1se"))[-1]
  zM <- scaleRow(zhe2[,feature])
  resps <-cbind(resps, predict(model,newx = as.matrix(zM),s = model$lambda.min,type='response')[,1])
}
scores <- apply(resps,1,mean)
zhe2Predict <- data.frame("B"=as.numeric(1-scores),"M"=as.numeric(scores))
rownames(zhe2Predict) <- names(scores)
zhe2Predict$predicted <- as.character(apply(zhe2Predict,1,function(v){names(v)[max(v)==v]}))
zhe2Predict$sampleId <- rownames(zhe2Predict)
zhe2Predict <- zhe2Predict[,c('sampleId','B','M','predicted')]

################################
resps <- data.frame(resps)
t1 <- sapply(1:3,function(i){
 
 data.frame('M'=resps[,i],'B'=1-resps[,i],'predicted'=as.character(sapply(resps[,i],function(v){ifelse(v>=0.5,'M','B')})),stringsAsFactors=F);
}
)
predictA <- data.frame(t1[,1]);
predictA$patientId <- rownames(resps)
predictA <- predictA[,c('patientId','M','B','predicted')]

predictB <- data.frame(t1[,2]);
predictB$patientId <- rownames(resps)
predictB <- predictB[,c('patientId','M','B','predicted')]

predictC <- data.frame(t1[,3]);
predictC$patientId <- rownames(resps)
predictC <- predictC[,c('patientId','M','B','predicted')]

if(F){
write.table(zhe2Predict,file="zhe2_150Samples_prediction_0422.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(predictA,file="predictA_150samples.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(predictB,file="predictB_150samples.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(predictC,file="predictC_150samples.txt",sep="\t",col.names=T,row.names=F,quote=F)


####################### ROC 
obj <- predictB
ROC1 <- roc(ifelse(obj$observed=="M", "M", "B"), as.numeric(obj$M))
pdf("ROC1_zhe150_ModelB_20190424.pdf")
plot.roc(ROC1,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",	
  main="model B",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()



}


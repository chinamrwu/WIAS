
setwd("F:/projects/TPD")
matDL <- read.csv('data/TPDDL_prot190503_for_predict.csv',header=T,stringsAsFactors=F,check.names=F)
matDL <- matDL[,c(1,which(grepl('HUMAN',colnames(matDL))))]
colnames(matDL)[-1] <- as.character(sapply(colnames(matDL)[-1],function(v){strsplit(v,"\\|")[[1]][2]}))
sampleIds <- matDL$sample[!grepl('rep',matDL$sample)]
rownames(matDL) <- matDL$sample
matDL <- matDL[,-1]
matDL <- 2^matDL

k0 <- sapply(sampleIds,function(v){
    rname <- c(v,paste0(v,"repB"))
    v1 <- apply(matDL[rname,],2,function(v){log2(mean(v,na.rm=T))})
    v1
})
DL <- data.frame(t(k0))
DL[is.na(DL)] <- NA


#############################################
matZE <- read.csv('data/TPDZE_prot190503_for_predict.csv',header=T,stringsAsFactors=F,check.names=F)
matZE <- matZE[,c(1,which(grepl('HUMAN',colnames(matZE))))]
colnames(matZE)[-1] <- as.character(sapply(colnames(matZE)[-1],function(v){strsplit(v,"\\|")[[1]][2]}))
sampleIds <- matZE$sample[!grepl('rep',matZE$sample)]
rownames(matZE) <- matZE$sample
matZE <- matZE[,-1]
matZE <- 2^matZE

k0 <- sapply(sampleIds,function(v){
    rname <- c(v,paste0(v,"repB"))
    v1 <- apply(matZE[rname,],2,function(v){log2(mean(v,na.rm=T))})
    v1
})

ZE <- data.frame(t(k0))
ZE[is.na(ZE)] <- NA


###################################################################################
model <- models[[2]]

matTest <- DL
feature  <- rownames(coef(model,s="lambda.1se"))[-1]
matK     <-  scaleRow(matTest[,feature])
score    <- predict(model,newx = as.matrix(matK),s = model$lambda.min,type='response')[,1]

predDL <- data.frame("M"=score,"B"=1-score)
predDL$predicted <- apply(predDL,1,function(v){names(v)[max(v)==v]})


matTest <- ZE
feature  <- rownames(coef(model,s="lambda.1se"))[-1]
matK     <-  scaleRow(matTest[,feature])
score    <- predict(model,newx = as.matrix(matK),s = model$lambda.min,type='response')[,1]

predZE <- data.frame("M"=score,"B"=1-score)
predZE$predicted <- apply(predZE,1,function(v){names(v)[max(v)==v]})
write.table(predZE,file='results/predict_Zhe2_149_0503.txt',sep="\t",col.names=T,row.names=T,quote=F)
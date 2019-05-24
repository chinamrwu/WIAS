rm(list=ls())
library(sqldf)
setwd("F:/projects/TPD")
source("src/common.R")
df0 <- read.table("data/TPD_RPM_SG_20190430.txt",sep="\t",header = T,stringsAsFactors = F)
df0 <- df0[grepl("_HUMAN",df0$ProteinName),-2]
df0$ProteinName <- as.character(sapply(df0$ProteinName,function(v){strsplit(v,"\\|")[[1]][2]}))
rownames(df0) <- as.character(sapply(1:dim(df0)[1],function(i){paste0('pepitde',i)}))
df0 <- df0[,-1]
df0 <- data.frame(t(df0))
batchId <- as.character(sapply(rownames(df0),function(v){strsplit(v,"_rep")[[1]][1]}))
dfA <- df0[grepl("repA",rownames(df0)),]
dfB <- df0[grepl("repB",rownames(df0)),]
rownames(dfA) <- as.character(sapply(rownames(dfA),function(v){strsplit(v,"_rep")[[1]][1]}))
rownames(dfB) <- as.character(sapply(rownames(dfB),function(v){strsplit(v,"_rep")[[1]][1]}))
clinic <- read.table("data/PRM_label_SG_20190430.txt",sep = "\t",header = T,stringsAsFactors = F)
#index <- match(rownames(dfA),clinic$number)nomatch=-1) # b62_6

rnames <- intersect(rownames(dfA),clinic$number)
dfA <- log2(dfA[rnames,])
index <- match(rownames(dfA),clinic$number,nomatch=-1);
clnames <- colnames(dfA)
dfA$label <- clinic$label[index]

########################################################### classification between A and C
library(glmnet)
library(umap)

AC <- dfA[dfA$label %in% c('A','C'),]
ump0 <- drawUMAP(AC,c('A'='blue','C'='red'), rowNormalization=F,strTitle="UMAP:A v C with orignal proteins ")

#tmp0 <- as.matrix(AC[,colnames(AC)!='label'])
#tmp0 <- apply(AC[,colnames(AC)!='label'],2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})

Y0   <- as.factor(AC$label)
selects <- list()
plots <- list()

tmp0 <- t(apply(AC[,colnames(AC)!='label'],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0

cvfit  <- cv.glmnet(tmp0,Y0,family='binomial',alpha=1,type.measure='class')
cf <- coef(cvfit,s='lambda.min')
selection <- rownames(cf)[cf[,1]!=0][-1]
length(selection)
selects[[length(selects)+1]] <- selection
plots[[length(plots)+1]] <- drawUMAP(AC[,c('label',selection)],c('A'='blue','C'='red'), rowNormalization=T,strTitle=sprintf("UMAP:%d proteins",length(selection)))


tmp0 <- t(apply(AC[,selection],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0
cvfit  <- cv.glmnet(tmp0,Y0,family='binomial',alpha=1,type.measure='class')
cf <- coef(cvfit,s='lambda.min')
selection <- rownames(cf)[cf[,1]!=0][-1]
length(selection)
selects[[length(selects)+1]] <- selection

plots[[length(plots)+1]] <-  drawUMAP(AC[,c('label',selection)],c('A'='blue','C'='red'), rowNormalization=T,strTitle=sprintf("UMAP:%d proteins",length(selection)))

tmp0 <- t(apply(AC[,selection],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0
cvfit  <- cv.glmnet(tmp0,Y0,family='binomial',alpha=1,type.measure='class')
cf <- coef(cvfit,s='lambda.min')
selection <- rownames(cf)[cf[,1]!=0][-1]
length(selection)
selects[[length(selects)+1]] <- selection
plots[[length(plots)+1]] <- drawUMAP(AC[,c('label',selection)],c('A'='blue','C'='red'), rowNormalization=T,strTitle=sprintf("UMAP:%d proteins",length(selection)))



######################################################### B M
BM <- dfA
BM$label[BM$label %in% c('N','M','A')] <- 'B'
BM$label[BM$label %in% c('C','P','W')] <- 'M'
ubmp <- drawUMAP(BM,c('B'='blue','M'='red'));

Y0 <- as.factor(BM$label)
selection <- colnames(BM)[colnames(BM) !='label']

ubmp0 <- drawUMAP(BM[,c('label',selection)],c('B'='blue','M'='red'),strTitle=sprintf("UMAP:%d features",length(selection)));

tmp0 <- t(apply(BM[,selection],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0
cvfit  <- cv.glmnet(tmp0,Y0,family='binomial',alpha=1,type.measure='class')
cf <- coef(cvfit,s='lambda.min')
selection <- rownames(cf)[cf[,1]!=0][-1]
length(selection)
ubmp1 <- drawUMAP(BM[,c('label',selection)],c('B'='blue','M'='red'),strTitle=sprintf("UMAP:%d features",length(selection)));

tmp0 <- t(apply(BM[,selection],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0
cvfit  <- cv.glmnet(tmp0,Y0,family='binomial',alpha=1,type.measure='class')
cf <- coef(cvfit,s='lambda.min')
selection <- rownames(cf)[cf[,1]!=0][-1]
length(selection)
ubmp2 <- drawUMAP(BM[,c('label',selection)],c('B'='blue','M'='red'),strTitle=sprintf("UMAP:%d features",length(selection)));

tmp0 <- t(apply(BM[,selection],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0
cvfit  <- cv.glmnet(tmp0,Y0,family='binomial',alpha=1,type.measure='class')
cf <- coef(cvfit,s='lambda.min')
selection <- rownames(cf)[cf[,1]!=0][-1]
length(selection)
ubmp3 <- drawUMAP(BM[,c('label',selection)],c('B'='blue','M'='red'),strTitle=sprintf("UMAP:%d features",length(selection)));

tmp0 <- t(apply(BM[,selection],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0
cvfit  <- cv.glmnet(tmp0,Y0,family='binomial',alpha=1,type.measure='class')
cf <- coef(cvfit,s='lambda.1se')
selection <- rownames(cf)[cf[,1]!=0][-1]
length(selection)
ubmp4 <- drawUMAP(BM[,c('label',selection)],c('B'='blue','M'='red'),strTitle=sprintf("UMAP:%d features",length(selection)));

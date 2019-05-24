rm(list=ls())
library(ggplot2)
library(caret)
library(umap)
library(tsne)
library(glmnet)


##############################
setwd("E:/projects/TPD")
source("src/common.R")
set.seed(190309)
pool <- read.csv('data/TPD_SG579_116poolProt_matrix_190304.csv',header=T,stringsAsFactors = F)
sampleInf <- read.table('data/AllSampleInformation_20190304.txt',sep="\t",header=T,stringsAsFactors = F)

machines <- as.character(sapply(colnames(pool)[-1],function(v){substr(v,1,1)}))
prots <- sapply(pool$prot,function(v){strsplit(v,"\\|")[[1]][2]})
pool <- data.frame(t(pool[,-1]))
colnames(pool) <- prots

rownames(pool) <- as.character(sapply(rownames(pool),function(v){a=strsplit(v,"sunyt_TPD_DIA_")[[1]];b=strsplit(a[2],"_pool")[[1]][1];paste0(a[1],"_",b)}))
pool$MS <- machines
pool$label <- machines
pool <- pool[,c('label',prots)]
R1 <- apply(pool,2,function(v){sum(is.na(v))/length(v)*100})
tmp0 <- pool[,R1 <= 50]
table(pool$label)
color3=c(A="red",B="black",C="blue")
color2=c(M='red',B='blue')
poolPCA50 <- drawPCA(tmp0,ptColors = color3,rowNormalization = T,colNormalization = T,strTitle="PCA:pool samples with 50%% missingRate")
print(poolPCA50)
tmp0 <- pool[,R1<=5]
poolPCA5 <- drawPCA(tmp0,ptColors = color3,rowNormalization = T,colNormalization = T,strTitle="PCA for pool samples with 5%% missingRate")
print(poolPCA5)

poolPCA5 <- poolPCA5 +scale_x_continuous(breaks = round(seq(min(poolPCA5$data$PC1), max(poolPCA5$data$PC1), by = 1),1)) +
     scale_y_continuous(breaks = round(seq(min(poolPCA5$data$PC2), max(poolPCA5$data$PC2), by = 1),1)) 
rowsA <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -10.5 & poolPCA5$data$PC2 < 3.9 & poolPCA5$data$PC1 < -3.5 & poolPCA5$data$label=='A'] # 25 batches
rowsB <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -3 & poolPCA5$data$PC1> 10 & poolPCA5$data$PC2 < 6 & poolPCA5$data$label=='A'] # 35batches

batchA <- data.frame(t(sapply(rowsA,function(v){a <- strsplit(v,"_")[[1]];b1=substr(a[1],1,1);b2=substr(a[1],2,nchar(a[1]));c(b1,b2,a[2])})),stringsAsFactors=F)
batchB <- data.frame(t(sapply(rowsB,function(v){a <- strsplit(v,"_")[[1]];b1=substr(a[1],1,1);b2=substr(a[1],2,nchar(a[1]));c(b1,b2,a[2])})),stringsAsFactors=F)
batchA$batchId <- as.character(sapply(rownames(batchA),function(v){a=strsplit(v,"_")[[1]];paste0(a[1:2],collapse="_")}))
batchB$batchId <- as.character(sapply(rownames(batchB),function(v){a=strsplit(v,"_")[[1]];paste0(a[1:2],collapse="_")}))

reps0 <- read.table('data/TPD_prot_matrix_rep_20190304.txt',sep="\t",header=T,stringsAsFactors = F)
R0 <- apply(reps0,2,function(v){sum(is.na(v))/length(v)*100})
for(missingRate in seq(5,100,by=5)){ 
set.seed(190309)

reps <- reps0[,R0 <= missingRate]
repInf <- reps[,1:3]
repInf$batchId <- as.character(sapply(repInf$filename,function(v){
    a=strsplit(v,"sunyt_TPD_DIA_")[[1]];
    b=strsplit(a[2],"_")[[1]][1];
    paste0(c(a[1],b),collapse="_")}))

#########
tmp <- sqldf("SELECT filename FROM repInf R,batchA A where A.batchId=R.batchID")
matA <- reps[match(tmp[,1],reps$filename),]
specimens <- unique(matA$SpecimenID)
tmp <- c()
for(sp in specimens){
  k0 <- apply(matA[matA$SpecimenID==sp,-c(1:3)],2,function(v){log2(mean(2^v,na.rm=T))})
  tmp <- rbind(tmp,k0)
}
A <- tmp
A[is.na(A)] <- NA
rownames(A) <- specimens
protA <- colnames(A)


########
tmp <- sqldf("SELECT filename FROM repInf R,batchB B where B.batchId=R.batchID")
matB <- reps[match(tmp[,1],reps$filename),]
specimens <- unique(matB$SpecimenID)
tmp <- c()
for(sp in specimens){
  k0 <- apply(matB[matB$SpecimenID==sp,-c(1:3)],2,function(v){log2(mean(2^v,na.rm=T))})
  tmp <- rbind(tmp,k0)
}
B <- tmp 
B[is.na(B)] <- NA
rownames(B) <- specimens
protB <- colnames(B)

tmp <- reps[,c("SpecimenID","label")]
tmp <- sqldf("SELECT distinct * from tmp")
rownames(tmp) <- tmp$SpecimenID

A <- data.frame(A)##225
B <- data.frame(B)##140
print(dim(A))
print(dim(B))
print(table(A$label))
print(table(B$label))

A$label <- tmp[rownames(A),"label"]
A <- A[,c("label",protA)]
B$label <- tmp[rownames(B),"label"]
B <- B[,c("label",protB)]

AB <- rbind(A,B)
print(dim(AB))
print(table(AB$label))
###############################################################
AB$label[AB$label %in% c('N','M','A')] <- 'B'
AB$label[AB$label %in% c('C','P','W')] <- 'M'

labelAB <- as.factor(AB$label)
print(table(AB$label))

color2=c(M='red',B='blue')
############################################################
AM <- as.matrix(AB[,-1]);AM[is.na(AM)] <- 0;
cvfitA0 <- cv.glmnet(AM,labelAB,family='binomial',alpha=0.5,type.measure='class')
cf <- coef(cvfitA0,s="lambda.min")
selectedA0 <- rownames(cf)[cf[,1] !=0][-1]
pAB1 <-  drawPCA(AB[,c('label',selectedA0)],ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset from machineA without normalization")
print(pAB1)

#############################################
AM <- t(apply(AB[,-1],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}));AM[is.na(AM)] <- 0;
cvfitA1 <- cv.glmnet(AM,labelAB,family='binomial',alpha=0.5,type.measure='class')
cf <- coef(cvfitA1,s="lambda.min")
selectedA1 <- rownames(cf)[cf[,1] !=0][-1]
pAB2 <-  drawPCA(AB[,c('label',selectedA1)],ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:machineA z-scored")
print(pAB2)

selection0 <- colnames(AB)[-1]
hitFitAB <- NULL
mse=100
flg=T
fitsAB <- list()
while(flg){
	tmp <- AB[,selection0]
	tmp    <- t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	tmp[is.na(tmp)] <- 0
	
	cvfit <- cv.glmnet(tmp,labelAB,family='binomial',alpha=0.5,type.measure='class')
	fitsAB[[length(fitsAB)+1]] <- cvfit
        cf <- coef(cvfit, s = "lambda.1se")
        selection1 <- rownames(cf)[cf[,1]!=0][-1]
        if(min(cvfit$cvm) <= mse){
           mse    <-  min(cvfit$cvm);
	   hitFitAB <-  cvfit
	}
	if(length(selection1)>= 4){ ## 4 proteins at least 
		    flg <- length(selection1) != length(selection0)
                    selection0 <- selection1
	}else{  flg=F	 }

}#flg

stsA <- sapply(fitsAB,function(obj){
     cf <- coef(obj,s = "lambda.1se");
     a <- sum(cf[,1]!= 0)-1;
      b=min(obj$cvm);
      c(a,b)
})
print("lasso result")
print(stsA)

obj <- hitFitAB
cf <- coef(obj, s = "lambda.1se")
selectionA0 <- rownames(cf)[cf[,1]!=0][-1]
print(sprintf("%d selected proteins from dataset A are:%s",length(selectionA0),paste0(selectionA0,collapse=",")))
strTit <- sprintf("PCA for machineA % d prots by %d iterations of ElasticNet",length(selectionA0),dim(stsA)[2])

pA3 <-  drawPCA(AB[,c('label',selectionA0)],ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=strTit)
print(pA3)



###################################################
avgMat <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(avgMat) <- avgMat$SpecimenID
avgMat <- avgMat[,-2]

testAB  <-  avgMat[setdiff(rownames(avgMat),rownames(AB)),selectionA0]
testAB    <- t(apply(testAB,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
testAB[is.na(testAB)] <- 0
lblA <- avgMat[rownames(testAB),"label"]
lblA[lblA %in% c('N','M','A')] <- 'B'
lblA[lblA %in% c('C','P','W')]  <- 'M'

tmpA <- AB[,selectionA0]
tmpA    <- t(apply(tmpA,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmpA[is.na(tmpA)] <- 0

cvfitAB <- glmnet(tmpA,labelAB,family='binomial',alpha=0,lambda=hitFitAB$lambda.min)
predA <- data.frame(predict(cvfitAB,newx = testAB,s = cvfitAB$lambda.min,type="class"))
predA$observed <- lblA
colnames(predA) <- c("predicted","observed")
print(sprintf("predicted £ºacc=%4.3f,missing cutoff=%d ",sum(predA$predicted==predA$observed)/dim(testAB)[1]*100,missingRate ))
}
#library(randomForest)
#RF <- 

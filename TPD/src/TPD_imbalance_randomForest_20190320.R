rm(list=ls())
library(ggplot2)
library(caret)
library(umap)
library(tsne)
library(glmnet)
library(sqldf)
library(caret)
library(pROC)
library(randomForest)

##############################
setwd("E:/projects/TPD")
source("src/common.R")
set.seed(190311)
pool <- read.csv('data/TPD_SG579_116poolProt_matrix_190304.csv',header=T,stringsAsFactors = F)
sampleInf <- read.table('data/AllSampleInformation_20190304.txt',sep="\t",header=T,stringsAsFactors = F)
reps <- read.table('data/TPD_prot_matrix_rep_20190304.txt',sep="\t",header=T,stringsAsFactors = F)
repInf <- reps[,1:3]
repInf$batchId <- as.character(sapply(repInf$filename,function(v){a=strsplit(v,"sunyt_TPD_DIA_")[[1]];b=strsplit(a[2],"_")[[1]][1];paste0(c(a[1],b),collapse="_")}))

rowNorm <- T
scaleRow <- function(M0){
  tmp <- M0;
  if(rowNorm){ tmp[,colnames(tmp)!='label'] <- t(scale(t(tmp[,colnames(tmp)!='label'])))}
  tmp[is.na(tmp)] <- 0
  tmp
}
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
color3=c(A="red",B="black",C="blue")
color2=c(M='red',B='blue')
poolPCA50 <- drawPCA(tmp0,ptColors = color3,rowNormalization = T,colNormalization = T,strTitle="PCA:pool samples with 50%% missingRate")
print(poolPCA50)
tmp0 <- pool[,R1<=5]
poolPCA5 <- drawPCA(tmp0,ptColors = color3,rowNormalization = T,colNormalization = T,strTitle="PCA for pool samples with 5%% missingRate")
print(poolPCA5)

poolPCA5 <- poolPCA5 +scale_x_continuous(breaks = round(seq(min(poolPCA5$data$PC1), max(poolPCA5$data$PC1), by = 1),1))+ 
     scale_y_continuous(breaks = round(seq(min(poolPCA5$data$PC2), max(poolPCA5$data$PC2), by = 1),1)) 
print(poolPCA5)
################################################
rowsA <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -10.5 & poolPCA5$data$PC2  < 0 & poolPCA5$data$PC1 < -3.5 & poolPCA5$data$label=='A'] # 25 batches
rowsB <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -3    &  poolPCA5$data$PC2 < 6 & poolPCA5$data$PC1> 11 ] # 35batches
rowsC <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -1.5  & poolPCA5$data$PC2  < 3.5]#51batch

getDataFrame <- function(rows){
   batchA <- data.frame(t(sapply(rows,function(v){a <- strsplit(v,"_")[[1]];b1=substr(a[1],1,1);b2=substr(a[1],1,nchar(a[1]));c(b1,paste0(b2,"_",a[2]))})),stringsAsFactors=F)
   tmp <- merge(repInf,batchA,by.x ='batchId',by.y="X2")["filename"]

	#tmp <- sqldf("SELECT filename FROM repInf,batchA  where batchA.X2=repInf.batchID")
   matA <- reps[match(tmp[,1],reps$filename),]
   specimens <- unique(matA$SpecimenID)
   tmp <- c()
   print("processing replicas .....")
   for(sp in specimens){
      k0 <- apply(matA[matA$SpecimenID==sp,-c(1:3)],2,function(v){log2(mean(2^v,na.rm=T))})
      tmp <- rbind(tmp,k0)
   }
   matA <- tmp
   matA[is.na(matA)] <- NA
   rownames(matA) <- specimens
   protA <- colnames(matA)

	tmp <- reps[,c("SpecimenID","label")]
   tmp <- sqldf("SELECT distinct * from tmp")
   rownames(tmp) <- tmp$SpecimenID

	matA <- data.frame(matA)

   matA$label <- tmp[rownames(matA),"label"]
   matA <- matA[,c("label",protA)]
   return(matA)
}
A  <- getDataFrame(rowsA)
B  <- getDataFrame(rowsB)
C1 <- getDataFrame(rowsC)
########
print("samples distribution in A:")
print(table(A$label))
print("samples distribution in B:")
print(table(B$label))

###############################################################
A$label[A$label %in% c('N','M','A')] <- 'B'
A$label[A$label %in% c('C','P','W')] <- 'M'
print("samples distribution in A:")
print(table(A$label))

B$label[B$label %in% c('N','M','A')] <- 'B'
B$label[B$label %in% c('C','P','W')] <- 'M'
print("samples distribution in B:")
print(table(B$label))

C1$label[C1$label %in% c('N','M','A')] <- 'B'
C1$label[C1$label %in% c('C','P','W')] <- 'M'
print("samples distribution in C1:")
print(table(C1$label))


color2=c(M='red',B='blue')
pA <- drawPCA(A,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset A with full features")
pA1 <- drawPCA(A,ptColors = color2,rowNormalization = F,colNormalization = F,strTitle="PCA:dataset A with full features")

pB <- drawPCA(B,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset B with full features")
pC <- drawPCA(C1,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset C with full features")

############################################################


###################################################
avgMat <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(avgMat) <- avgMat$SpecimenID
avgMat <- avgMat[,-2]

######################################################
p20s <- read.table("data/better6combins.txt",sep=" ",stringsAsFactors=F)


plotsA <- apply(p20s,1,function(v){
  tmp <- A[,c('label',as.character(v[1:6]))]
  p0  <- drawPCA(tmp,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=sprintf("A:%s",paste0(v[1:6],collapse=",")))
  p0
})

plotsB <- apply(p20s,1,function(v){
  tmp <- B[,c('label',as.character(v[1:6]))]
  p0  <- drawPCA(tmp,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=sprintf("A:%s",paste0(v[1:6],collapse=",")))
  p0
})
plotsC <- apply(p20s,1,function(v){
  tmp <- C1[,c('label',as.character(v[1:6]))]
  p0  <- drawPCA(tmp,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=sprintf("A:%s",paste0(v[1:6],collapse=",")))
  p0
})
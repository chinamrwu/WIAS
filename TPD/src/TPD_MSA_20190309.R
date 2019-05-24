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
reps <- read.table('data/TPD_prot_matrix_rep_20190304.txt',sep="\t",header=T,stringsAsFactors = F)
repInf <- reps[,1:3]
repInf$batchId <- as.character(sapply(repInf$filename,function(v){a=strsplit(v,"sunyt_TPD_DIA_")[[1]];b=strsplit(a[2],"_")[[1]][1];paste0(c(a[1],b),collapse="_")}))


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

rowsA <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -10.5 & poolPCA5$data$PC2 < 0 & poolPCA5$data$PC1 < -3.5 & poolPCA5$data$label=='A'] # 25 batches
rowsB <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -3   &  poolPCA5$data$PC2 < 6 & poolPCA5$data$PC1> 11 ] # 35batches

getDataFrame <- function(rows){
   batchA <- data.frame(t(sapply(rows,function(v){a <- strsplit(v,"_")[[1]];b1=substr(a[1],1,1);b2=substr(a[1],2,nchar(a[1]));c(b1,b2,a[2])})),stringsAsFactors=F)
   tmp <- sqldf("SELECT filename FROM repInf R,batchA A where A.batchId=R.batchID")
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

A <- getDataFrame(rowsA)
B <- getDataFrame(rowsB)

########
print("samples distribution in A:")
print(table(A$label))
print("samples distribution in B:")
print(table(B$label))

###############################################################
A$label[A$label %in% c('N','M','A')] <- 'B'
A$label[A$label %in% c('C','P','W')] <- 'M'

B$label[B$label %in% c('N','M','A')] <- 'B'
B$label[B$label %in% c('C','P','W')] <- 'M'

print("samples distribution in A:")
print(table(A$label))
print("samples distribution in B:")
print(table(B$label))


color2=c(M='red',B='blue')
pA <- drawPCA(A,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset A with full features")
pB <- drawPCA(B,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset B with full features")
############################################################
if(F){
	AM <- as.matrix(A[,-1]);AM[is.na(AM)] <- 0;
	BM <- as.matrix(B[,-1]);BM[is.na(BM)] <- 0;

	cvfitA0 <- cv.glmnet(AM,as.factor(A$label),family='binomial',alpha=0.5,type.measure='class')
	cvfitB0 <- cv.glmnet(BM,as.factor(B$label),family='binomial',alpha=0.5,type.measure='class')

	cf <- coef(cvfitA0,s="lambda.min")
	selectedA0 <- rownames(cf)[cf[,1] !=0][-1]
	cf <- coef(cvfitB0,s="lambda.min")
	selectedB0 <- rownames(cf)[cf[,1] !=0][-1]
	str1 <- "PCA: %d prots from %s unormlized by 1 iteration of elasticNet"
	pA1 <-  drawPCA(A[,c('label',selectedA0)],ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=sprintf(str1,length(selectedA0),"A"))
	pB1 <-  drawPCA(B[,c('label',selectedB0)],ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=sprintf(str1,length(selectedB0),"B"))

	#############################################
	AM <- t(apply(A[,-1],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}));AM[is.na(AM)] <- 0;
	BM <- t(apply(B[,-1],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}));BM[is.na(BM)] <- 0;

	cvfitA1 <- cv.glmnet(AM,as.factor(A$label),family='binomial',alpha=0.5,type.measure='class')
	cvfitB1 <- cv.glmnet(BM,as.factor(B$label),family='binomial',alpha=0.5,type.measure='class')

	cf <- coef(cvfitA1,s="lambda.min")
	selectedA1 <- rownames(cf)[cf[,1] !=0][-1]
	cf <- coef(cvfitB1,s="lambda.min")
	selectedB1<- rownames(cf)[cf[,1] !=0][-1]
	pA2 <-  drawPCA(A[,c('label',selectedA1)],ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset A samples z-scored")
	pB2 <-  drawPCA(B[,c('label',selectedB1)],ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="PCA:dataset B samples z-scored")
}
#############################################################################################################################################################
iterate <- function(M1){
     selection0 <- colnames(A)
     hitFit <- NULL
     mse=100
     flg=T

	  result <- list()
     fits <- list()
     while(flg){
		  tmp <- M1[,selection0]
		  tmp    <- t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
		  tmp[is.na(tmp)] <- 0
	
	     cvfit <- cv.glmnet(tmp,as.factor(M1$label),family='binomial',alpha=0.5,type.measure='class')
	     fits[[length(fits)+1]] <- cvfit
        cf <- coef(cvfit, s = "lambda.1se")
        selection1 <- rownames(cf)[cf[,1]!=0][-1]
        if(min(cvfit$cvm) <= mse){
           mse    <-  min(cvfit$cvm);
	        hitFit <-  cvfit
	     }
	     if(length(selection1)>= 4){ ## 4 proteins at least 
		    flg <- length(selection1) != length(selection0)
          selection0 <- selection1
	      }else{  flg=F	 }
     }#flg
     result$fits <- fits
     sts <- sapply(fitsA,function(obj){
			  cf <- coef(obj,s = "lambda.1se");
			  a <- sum(cf[,1]!= 0)-1;
			  b=min(obj$cvm);
			  c(a,b)
     })
	  result$sts <- sts
	  result$hitFit <- hitFit
     print("search result:")
     print(stsA)
     obj <- hitFit
     cf <- coef(obj, s = "lambda.1se")
     selectionA0 <- rownames(cf)[cf[,1]!=0][-1]
	  result$selected <- selected
     print(sprintf("%d selected proteins from dataset A are:%s",length(selectionA0),paste0(selectionA0,collapse=",")))
	  result
  }##function iterate

itA <- iterate(A)
itB <- iterate(B)


#################################### predict each other between dataset A and B ###############

tmpA <- A[,selectionA0]
tmpA    <- t(apply(tmpA,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmpA[is.na(tmpA)] <- 0

testA <- B[,selectionA0]
testA    <- t(apply(testA,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
testA[is.na(testA)] <- 0

cvfitA <- glmnet(tmpA,labelA,family='binomial',alpha=0,lambda=hitFitA$lambda.min)

predA <- data.frame(predict(cvfitA,newx = testA,s = cvfitA$lambda.min,type="class"))
predA$observed <- as.character(labelB)
colnames(predA) <- c("predicted","observed")
print(sprintf("predicted by A£º%4.3f",sum(predA$predicted==predA$observed)/dim(testA)[1]*100))


tmpB <- B[,selectionB0]
tmpB    <- t(apply(tmpB,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmpB[is.na(tmpB)] <- 0

testB <- A[,selectionB0]
testB    <- t(apply(testB,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
testB[is.na(testB)] <- 0

cvfitB <- glmnet(tmpB,labelB,family='binomial',alpha=0,lambda=hitFitB$lambda.min)

predB <- data.frame(predict(cvfitB,newx = testB,s = cvfitB$lambda.min,type="class"))
predB$observed <- as.character(labelA)
colnames(predB) <- c("predicted","observed")
print(sprintf("predicted by B£º%4.3f",sum(predA$predicted==predA$observed)/dim(testB)[1]*100))

###################################################
avgMat <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(avgMat) <- avgMat$SpecimenID
avgMat <- avgMat[,-2]

testA  <-  avgMat[setdiff(rownames(avgMat),rownames(A)),selectionA0]
testA    <- t(apply(testA,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
testA[is.na(testA)] <- 0
lblA <- avgMat[rownames(testA),"label"]
lblA[lblA %in% c('N','M','A')] <- 'B'
lblA[lblA %in% c('C','P','W')]  <- 'M'

predA <- data.frame(predict(cvfitA,newx = testA,s = cvfitA$lambda.min,type="class"))
predA$observed <- lblA
colnames(predA) <- c("predicted","observed")
print(sprintf("predicted by A£º%4.3f",sum(predA$predicted==predA$observed)/dim(testA)[1]*100))

testB  <-  avgMat[setdiff(rownames(avgMat),rownames(B)),selectionB0]
testB    <- t(apply(testB,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
testB[is.na(testB)] <- 0
lblB <- avgMat[rownames(testB),"label"]
lblB[lblB %in% c('N','M','A')] <- 'B'
lblB[lblB %in% c('C','P','W')]  <- 'M'

predB  <- data.frame(predict(cvfitB,newx = testB,s = cvfitB$lambda.min,type="class"))
predB$observed <- lblB
colnames(predB) <- c("predicted","observed")
print(sprintf("predicted by B£º%4.3f",sum(predA$predicted==predA$observed)/dim(testB)[1]*100))


################## 
tmp <- avgMat[,c('label',selectionA0)]
tmp$label[tmp$label %in% c('N','M','A')] <- 'B'
tmp$label[tmp$label %in% c('C','P','W')]  <- 'M'
pAB0 <- drawPCA(tmp,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="AB" )

tmp <- avgMat[,c('label',selectionB0)]
tmp$label[tmp$label %in% c('N','M','A')] <- 'B'
tmp$label[tmp$label %in% c('C','P','W')]  <- 'M'
pAB1 <- drawPCA(tmp,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle="AB" )

##########################################################

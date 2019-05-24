options("width"=500)
if(T){
	#rm(list=ls())
	library(ggplot2)
	library(caret)
	library(umap)
	library(tsne)
	library(glmnet)
	library(sqldf)
	library(caret)
	setwd("F:/projects/TPD")
	#setwd("/work")
	source("src/common.R")
	source("src/MachineLearning.R")
	set.seed(190311)
	#############################################################################################################
	color3=c(A="red",B="black",C="blue")
	color2=c(M='red',B='blue')

	#####################################  Loading all needed dataset  ############################################
	print("Loading datasets ......")

   ################################
	A0  <- read.table("data/A0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
	B0  <- read.table("data/B0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
	C0  <- read.table("data/C0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
	df0 <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
	rownames(df0) <- df0$SpecimenID
	df0 <- df0[,-2]
	print("distribution A:");print(table(A0$label))
	print("distribution B:");print(table(B0$label))
	print("distribution C:");print(table(C0$label))

	A0$label[A0$label %in% c('N','M','A')]  <- 'B'
	A0$label[A0$label %in% c('C','P','W')]  <- 'M'
	B0$label[B0$label %in% c('N','M','A')]  <- 'B'
	B0$label[B0$label %in% c('C','P','W')]  <- 'M'
	C0$label[C0$label %in% c('N','M','A')]  <- 'B'
	C0$label[C0$label %in% c('C','P','W')]  <- 'M'

	df0$label[df0$label %in% c('N','M','A')]  <- 'B'
	df0$label[df0$label %in% c('C','P','W')]  <- 'M'

	print("distribution A:");print(table(A0$label))
	print("distribution B:");print(table(B0$label))
	print("distribution C:");print(table(B0$label))


	zhe1 <- read.table("data/zhe1_54samples_protMat.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
	zhe1$label <- sapply(rownames(zhe1),function(v){ch=substr(v,1,1);ifelse(ch=='A','B','M')})
	zhe1 <- zhe1[,c(dim(zhe1)[2],1:(dim(zhe1)[2]-1))]

	TestA <- df0[setdiff(rownames(df0),rownames(A0)),]
	TestB <- df0[setdiff(rownames(df0),rownames(B0)),]
	TestC <- df0[setdiff(rownames(df0),rownames(C0)),]
	rowNorm <- T

	missingA <- apply(A0,2,function(v){sum(is.na(v))/length(v)*100})
	missingB <- apply(B0,2,function(v){sum(is.na(v))/length(v)*100})
	missingC <- apply(C0,2,function(v){sum(is.na(v))/length(v)*100})

	##############################################################################################################
	allModels <- read.table("data/AllModels_sta.txt",sep="\t",header=T,stringsAsFactors=F)
	allModels <- allModels[order(allModels$avgAcc,decreasing=T),]
	allModels <- allModels[allModels$src %in% c('A','B','C'),]
	allProts  <- read.table("data/Models_Prots.txt",sep="\t",header=T,stringsAsFactors=F)
}
##########################################################################################################################

scaleRow <- function(M0){
 tmp <- M0
 tmp[,colnames(tmp)!='label'] <- t(scale(t(tmp[,colnames(tmp)!='label'])))
 tmp[is.na(tmp)] <- 0
 tmp
}

################# 根据 6 proteins 模型 在B0数据集上性能较好，但对benign样本错误率较高，做以下过滤
k0 <- sqldf("SELECT A.*,B.* FROM allModels A,allProts B where A.ID=B.ID and A.src='B'  and A.avgAcc > 85 and protNumber between 7 and 10 order by accB desc") 

seed <- k0$seed[1]
protId <- strsplit(k0$prots[1],",")[[1]];

#训练模型
set.seed(seed)
trainDS <- scaleRow(PseudoSamples(B0[,c('label',protId)])) # normalization dataset
model <-  cv.glmnet(as.matrix(trainDS[,-1]),as.factor(trainDS$label),family='binomial',alpha=0,type.measure='class') ## training logistic regression model

###### the following function test the model 
testModel <- function(model,testMat){
	feature    <- rownames(coef(model,s="lambda.1se"))[-1]
	matK       <-  scaleRow(testMat[,feature])
	score      <- predict(model,newx = as.matrix(matK),s = model$lambda.min,type='response')[,1]
	prediction <- data.frame("M"=score,"B"=1-score)
	prediction$predicted <- apply(prediction,1,function(v){names(v)[max(v)==v]})
	prediction
}

### Test in SG307
sg307Prediction <- testModel(model,TestB)
sg307Prediction$observed  <- TestB$label
n <- sum(sg307Prediction$observed==sg307Prediction$predicted)
print(sprintf("%d samples are predicted, %d correct ,accuracy is %4.3f",dim(TestB)[1],n,n/dim(TestB)[1]))

### Test in Zhe1 54 samples

zhe54Prediction <- testModel(model, zhe1)
zhe54Prediction$observed <- zhe1$label
n <- sum(zhe54Prediction$observed==zhe54Prediction$predicted)
print(sprintf("%d samples are predicted, %d correct ,accuracy is %4.3f",dim(zhe1)[1],n,n/dim(zhe1)[1]))
##########################  DL53 samples

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
DL53 <- data.frame(t(k0))
DL53[is.na(DL53)] <- NA

DL53Prediction <- testModel(model,DL53)
##### Test in Zhe2 150 samples
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

ZE150 <- data.frame(t(k0))
ZE150[is.na(ZE150)] <- NA

ZE150Prediction <- testModel(model,ZE150)

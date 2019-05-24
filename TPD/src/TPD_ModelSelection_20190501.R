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

	#################################################################################
	print("Loading datasets ......")


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

allModels <- allModels[allModels$src %in% c('A','B','C'),]

getModel <- function(modelId){
	  #modelId=5744
	  if(sum(allModels$ID==modelId) == 0){print(sprintf('No model with ID %d found!',modelId));return(NULL)}
	  protId <- c('label',strsplit(allProts$prots[allProts$ID==modelId],",") [[1]])
	  seed <- allModels$seed[allModels$ID==modelId]
	  src  <- allModels$src[allModels$ID==modelId]
	  
	  kk=NULL
	  ifelse(src=='A',kk <- A0[,protId],ifelse(src=='B',kk <- B0[,protId],kk <- C0[,protId]))

	  set.seed(seed)
	  trainDS <- scaleRow(PseudoSamples(kk))
	  model <-  cv.glmnet(as.matrix(trainDS[,-1]),as.factor(trainDS$label),family='binomial',alpha=0,type.measure='class')
	  return(list('mod'=model,'ID'=modelId,'src'=src))
}

TestIt <- function(model,matV=NULL){
        if(is.null(matV)){ return(NULL)}
	feature    <- rownames(coef(model,s="lambda.1se"))[-1]
	cfeatures  <- colnames(matV)[colnames(matV) !='label']
	if(! all(feature %in% cfeatures))
	 {
	     print(sprintf('need %d, found: %d',length(feature),sum(feature %in% cfeatures)))
	      return(NULL)
	 }


	Y1         <- matV$label
	#print(sprintf('%d features extracted:%s',length(feature),paste0(feature,collapse=",")))
	tmp1 <- scaleRow(matV[,feature])
	response <- predict(model,newx = as.matrix(tmp1),s = model$lambda.1se,type='response')
	clss     <- predict(model,newx = as.matrix(tmp1),s = model$lambda.1se,type='class')
	pr0 <- data.frame("B"=1-response[,1],"M"=response[,1],"predicted"=clss[,1],"observed"=Y1,stringsAsFactors=F)
	acc <- sum(pr0$predicted==pr0$observed)/dim(pr0)[1]*100
	#predicts1[[length(predicts1)+1]] <- pr0
	if(acc >= 85 ){
	   #print(sprintf('%d samples  acc=%4.3f;%d features:%s',dim(matV)[1],acc,length(feature),paste0(feature,collapse=",")))
	}
	acc/100
}

####################################################################################

zhe1 <- read.table("data/zhe1_54samples_protMat.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
zhe1$label <- sapply(rownames(zhe1),function(v){ch=substr(v,1,1);ifelse(ch=='A','B','M')})
zhe1 <- zhe1[,c(dim(zhe1)[2],1:(dim(zhe1)[2]-1))]

zhe2 <- read.table('data/zhe2_150_protMat_0501.txt',sep="\t",header=T,stringsAsFactors=F)
zhe2$label <- as.character(sapply(zhe2$label,function(v){ifelse(v %in% c('N','M','A'),'B','M')}))

selectedModels <- list()
for(ID in allModels$ID[allModels$protNumber <= 18 & allModels$avgAcc > 85]){
  selectedModels[[length(selectedModels)+1]] <- getModel(ID)
}

index=0
accs <- c()
for(model in selectedModels){
   acc1=TestIt(model$mod,zhe1)
   acc2=TestIt(model$mod,zhe2)
   acc <- (acc1*dim(zhe1)[1]+acc2*dim(zhe2)[1])/(dim(zhe1)[1]+dim(zhe2)[1])*100
   accs <- c(accs,acc)

   index <- index +1
   print(sprintf('%d ========FROM %s :%d ================ overall %4.3f ',index,model$src,model$ID,acc))
}

model <- selectedModels[[which(accs==max(accs))]]
model <- model$mod

resps1 <- c() ### zhe2 150
feature    <- rownames(coef(model,s="lambda.1se"))[-1]
zM <- scaleRow(zhe2[,feature])
resps1 <-cbind(resps1, predict(model,newx = as.matrix(zM),s = model$lambda.min,type='response')[,1])

scores <- apply(resps1,1,mean)
zhe2Predict <- data.frame("B"=as.numeric(1-scores),"M"=as.numeric(scores))
rownames(zhe2Predict) <- names(scores)
zhe2Predict$predicted <- as.character(apply(zhe2Predict,1,function(v){names(v)[max(v)==v]}))
zhe2Predict$sampleId <- zhe2$sampleId
zhe2Predict <- zhe2Predict[,c('sampleId','B','M','predicted')]
rownames(zhe2Predict) <- rownames(zhe2)
##### zhe1 54
resps2 <- c() ### zhe1 54
feature    <- rownames(coef(model,s="lambda.1se"))[-1]
zM <- scaleRow(zhe1[,feature])
resps2 <-cbind(resps2, predict(model,newx = as.matrix(zM),s = model$lambda.min,type='response')[,1])

scores <- apply(resps2,1,mean)
zhe1Predict <- data.frame("B"=as.numeric(1-scores),"M"=as.numeric(scores))
rownames(zhe1Predict) <- rownames(resps2)
zhe1Predict$predicted <- as.character(apply(zhe1Predict,1,function(v){names(v)[max(v)==v]}))
zhe1Predict$sampleId <- rownames(zhe1Predict)
zhe1Predict <- zhe1Predict[,c('sampleId','B','M','predicted')]
zhe1Predict$observed <- zhe1$label
write.table(zhe1Predict,file="zhe2_54samplesBy10Proteins_0502.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################## Zhe Er 149
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

matTest <- ZE
feature  <- rownames(coef(model,s="lambda.1se"))[-1]
matK     <-  scaleRow(matTest[,feature])
score    <- predict(model,newx = as.matrix(matK),s = model$lambda.min,type='response')[,1]

predZE <- data.frame("M"=score,"B"=1-score)
predZE$predicted <- apply(predZE,1,function(v){names(v)[max(v)==v]})
write.table(predZE,file='results/predictBy10proteins_Zhe2_149_0503.txt',sep="\t",col.names=T,row.names=T,quote=F)

############################ DL 53

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

matTest <- DL
feature  <- rownames(coef(model,s="lambda.1se"))[-1]
matK     <-  scaleRow(matTest[,feature])
score    <- predict(model,newx = as.matrix(matK),s = model$lambda.min,type='response')[,1]

predDL <- data.frame("M"=score,"B"=1-score)
predDL$predicted <- apply(predDL,1,function(v){names(v)[max(v)==v]})
#write.table(predDL,file='results/predictBy10proteins_DL53_0503.txt',sep="\t",col.names=T,row.names=T,quote=F)
################################################# SG-307

matTest <- TestB
feature  <- rownames(coef(model,s="lambda.1se"))[-1]
matK     <-  scaleRow(matTest[,feature])
score    <- predict(model,newx = as.matrix(matK),s = model$lambda.min,type='response')[,1]
predSG307 <- data.frame("M"=score,"B"=1-score)
predSG307$predicted <- apply(predSG307,1,function(v){names(v)[max(v)==v]})
predSG307$observed  <- TestB$label
write.table(predSG307,file='results/predictBy10proteins_SG307_0505.txt',sep="\t",col.names=T,row.names=T,quote=F)
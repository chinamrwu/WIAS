options("width"=500)
if(T){
rm(list=ls())
library(ggplot2)
library(caret)
library(umap)
library(tsne)
library(glmnet)
library(sqldf)
library(caret)
setwd("D:/westlake/projects/TPD")
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
allProts  <- read.table("data/Models_Prots.txt",sep="\t",header=T,stringsAsFactors=F)
}
######################################
getProtIds <- function(IDs){
    tmp <- allProts[allProts$ID %in% IDs,"prots"]
	 pids <- c()
	 for(obj in tmp){pids <- unique(c(pids,strsplit(obj,",")[[1]]))}
	 print(sprintf("%d proteins total",length(pids)))
	 pids
}

################################ try to use best model ###############
print("------------------------------------------------------")
modelsA <- allModels[allModels$src=='A',]
modelsB <- allModels[allModels$src=='B',]
modelsC <- allModels[allModels$src=='C',]

modelsA <- modelsA[modelsA$protNumber < 7,]
modelsB <- modelsB[modelsB$protNumber < 7,]
modelsC <- modelsC[modelsC$protNumber < 7,]

selectedModels <- rbind(modelsA[1,],modelsB[1,])
selectedModels <- rbind(selectedModels,modelsC[1,])
selectedModels <- merge(x=selectedModels,y=allProts,by='ID')
selectedProts  <- unique(unlist(sapply(selectedModels$prots,function(v){strsplit(v,",")[[1]]})))
models <- list()
validSet <- list();
for(i in 1:dim(selectedModels)[1]){
    protId <- c('label',strsplit(selectedModels$prots[i],",")[[1]])
	 #protId  <- c('label',selectedProts)
	 rate   <- selectedModels$rate[i]
	 seed   <- selectedModels$seed[i]
	 ch     <- selectedModels$src[i]
		kk=NULL
		ifelse(ch=='A',kk <- A0[,protId],ifelse(ch=='B',kk <- B0[,protId],kk <- C0[,protId]))
		ifelse(ch=='A',validSet[[length(validSet)+1]] <- TestA,ifelse(ch=='B',validSet[[length(validSet)+1]] <- TestB,validSet[[length(validSet)+1]] <- TestC))
		set.seed(seed)
		trainDS <- scaleRow(PseudoSamples(kk))
      #trainDS <- scaleRow(kk)
		model <-  cv.glmnet(as.matrix(trainDS[,-1]),as.factor(trainDS$label),family='binomial',alpha=0,type.measure='class')
		models[[length(models)+1]] <- model
}

features <- c()

rows  <- unique(c(rownames(TestA),rownames(TestB),rownames(TestC)))
TestD <- df0[setdiff(rownames(df0),rows),]

predictions <- data.frame(matrix(0,nrow=dim(df0)[1],ncol=2))
rownames(predictions) <- rownames(df0)
predicts1 <- list()

for(i in 1:length(models)){
    Y1    <- validSet[[i]]$label
	 model <- models[[i]]
	 feature    <- rownames(coef(model,s="lambda.1se"))[-1]
	 features   <- unique(c(features,feature))
    tmp1 <- scaleRow(validSet[[i]][,feature])
    response <- predict(models[[i]],newx = as.matrix(tmp1),s = model$lambda.1se,type='response')
	 clss <- predict(models[[i]],newx = as.matrix(tmp1),s = model$lambda.1se,type='class')
	 pr0 <- data.frame("B"=1-response[,1],"M"=response[,1],"predicted"=clss[,1],"observed"=Y1,stringsAsFactors=F)
	 
	 TPV <- sum(pr0$M > pr0$B & pr0$observed=='M')/sum(pr0$predicted=='M') 
	 NPV <- sum(pr0$M < pr0$B & pr0$observed=='B')/sum(pr0$predicted=='B')
	 acc <- sum(pr0$predicted==pr0$observed)/dim(pr0)[1]*100
	 predicts1[[length(predicts1)+1]] <- pr0
	 print(sprintf("acc=%4.3f",acc))
	 predictions[rownames(pr0),] <-  predictions[rownames(pr0),] + pr0[,1:2]*c(NPV,TPV)
    
	 tmp1 <- scaleRow(TestD[,feature])
	 Y2   <- TestD$label
	 response <- predict(models[[i]],newx = as.matrix(tmp1),s = model$lambda.1se,type='response')
	 clss <- predict(models[[i]],newx = as.matrix(tmp1),s = model$lambda.1se,type='class')
	 pr1 <- data.frame("B"=1-response[,1],"M"=response[,1],"predicted"=clss[,1],"observed"=Y2,stringsAsFactors=F)
    predictions[rownames(pr1),] <- pr1[,1:2]*c(NPV,TPV) + predictions[rownames(pr1),]

}

predictions <- data.frame(t(apply(predictions,1,function(v){v/sum(v)})))
colnames(predictions) <- c("B","M")
predictions$predicted <- apply(predictions,1,function(v){names(v)[max(v)==v]})
predictions$observed=df0$label

predicts1[[length(predicts1)+1]] <- predictions


names(predicts1) <- c("modelA","modelB","modelC","Ensemble") ### computing sensitivity and specificity
performances <- sapply(predicts1,function(obj){
   TP <- sum(obj$observed=='M' &  obj$predicted=='M' )
	FP <- sum(obj$observed=='B' &  obj$predicted=='M' )
	TN <- sum(obj$observed=='B' &  obj$predicted=='B' )
	FN <- sum(obj$observed=='M' &  obj$predicted=='B' )
	c(TP/(TP+FN),TN/(TN+FP),(TP+TN)/dim(obj)[1],TP/(TP+FP),TN/(TN+FN))
})
colnames(performances) <- c("modelA","modelB","modelC","Ensemble")
rownames(performances) <- rownames(performances) <- c("Sensitivity","Specificity","Accuracy","PPV","NPV")
prevalences <- c(sum(TestA$label=='M')/dim(TestA)[1],sum(TestB$label=='M')/dim(TestB)[1])
prevalences <- c(prevalences,sum(TestC$label=='M')/dim(TestC)[1],sum(df0$label=='M')/dim(df0)[1])
correctedPPVs <- sapply(1:length(prevalences),function(i){
  v <- performances[,i]
  prev <- prevalences[i]
  c(v[1]*prev/(v[1]*prev+(1-v[2])*(1-prev)),(v[2]*(1-prev))/((1-v[1])*prev+v[2]*(1-prev)))
})
rownames(correctedPPVs) <- c("PPV","NPV")
colnames(correctedPPVs) <- c("modelA","modelB","modelC","Ensemble")

resps <- c()
for(model in models){
  feature    <- rownames(coef(model,s="lambda.1se"))[-1]
  zM <- scaleRow(zhe1[,feature])
  resps <-cbind(resps, predict(model,newx = as.matrix(zM),s = model$lambda.min,type='response')[,1])
}
scores <- apply(resps,1,mean)
zPred <- data.frame("B"=as.numeric(1-scores),"M"=as.numeric(scores))
rownames(zPred) <- names(scores)
zPred$observed  <- as.character(sapply(rownames(zPred),function(v){a=substr(v,1,1);ifelse(a=='A','B','M')}))
zPred$predicted <- apply(zPred,1,function(v){v1=v[1:2];names(v1)[max(v1)==v1]})
print(sprintf("accuracy for 54 samples:%4.3f",sum(zPred$observed==zPred$predicted)/dim(zPred)[1]*100))
if(F){
#################################################### ROC
obj <- predicts1[[1]]
ROC <- roc(ifelse(obj$observed=="M", "M", "B"), as.numeric(obj$M))
pdf("ROC_ModelA_9proteins.pdf")
plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",	
  main="model A",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
##################
obj <- predicts1[[2]]
ROC <- roc(ifelse(obj$observed=="M", "M", "B"), as.numeric(obj$M))
pdf("ROC_ModelB_9proteins.pdf")
plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",	
  main="model B",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
##################
obj <- predicts1[[3]]
ROC <- roc(ifelse(obj$observed=="M", "M", "B"), as.numeric(obj$M))
pdf("ROC_ModelC_9proteins.pdf")
plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",	
  main="model C",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
#################
obj <- predicts1[[4]]
ROC <- roc(ifelse(obj$observed=="M", "M", "B"), as.numeric(obj$M))
pdf("ROC_Model_Ensemble_9proteins.pdf")
plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",	
  main="Ensemble 3 models",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
}

############ the codes below extract misclassified samples 
ump1 <- drawUMAP(df0[,c('label',features)],color2,rowNormalization=T,colNormalization=T,strTitle="UMAP:579 sampels with 9 proteins")
ump1 <- ump1+scale_x_continuous(breaks = round(seq(min(ump1$data$X), max(ump1$data$X), by = 0.5),1))+ 
     scale_y_continuous(breaks = round(seq(min(ump1$data$Y), max(ump1$data$Y), by = 0.5),1)) 
dat <- ump1$data
wrongA <- rownames(dat)[dat$Y > 0 & dat$X < 1 & dat$label=='M'] #
wrongB <- rownames(dat)[dat$X > 1.3  & dat$label=='M']
wrongC <- rownames(dat)[dat$Y < -1.6 & dat$label=='B']

#########################
if(F){
df0 <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(df0) <- df0$SpecimenID
df0 <- df0[,-2]

sup1 <- data.frame("specimenId"=wrongA,src='A',"subType"=df0[wrongA,'label'],stringsAsFactors=F)
sup2 <- data.frame("specimenId"=wrongB,src='B',"subType"=df0[wrongB,'label'],stringsAsFactors=F)
sup3 <- data.frame("specimenId"=wrongC,src='C',"subType"=df0[wrongC,'label'],stringsAsFactors=F)
allSup <- rbind(rbind(sup1,sup2),sup3)
write.csv(allSup,file="results/TPD_suspicious_samples.csv",col.names=T,row.names=F,quote=F)
}
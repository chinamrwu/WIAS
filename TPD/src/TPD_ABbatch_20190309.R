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
rowsA <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -10.5 & poolPCA5$data$PC2 < 0 & poolPCA5$data$PC1 < -3.5 & poolPCA5$data$label=='A'] # 25 batches
rowsB <-  rownames(poolPCA5$data)[poolPCA5$data$PC2 > -3   &  poolPCA5$data$PC2 < 6 & poolPCA5$data$PC1> 11 ] # 35batches

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
print("samples distribution in A:")
print(table(A$label))

B$label[B$label %in% c('N','M','A')] <- 'B'
B$label[B$label %in% c('C','P','W')] <- 'M'
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
iterate <- function(train,validation){
     selection0 <- colnames(train)[-1]
     hitFit <- NULL
     mse=100
     flg=T
     fits <- list()
     while(flg){
		  tmp <- scaleRow(train[,selection0])
	
	     cvfit <- cv.glmnet(as.matrix(tmp),as.factor(train$label),family='binomial',alpha=0.5,type.measure='class')
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
     sts <- sapply(fits,function(obj){
			  cf <- coef(obj,s = "lambda.1se");
			  a <- sum(cf[,1]!= 0)-1;
			  b=min(obj$cvm);
			  c(a,b)
     })
	  #print("search result:")
     #print(sts)
     obj <- hitFit
     cf <- coef(obj, s = "lambda.1se")
     selectionA0 <- rownames(cf)[cf[,1]!=0][-1]
  
	  tmp1 <- scaleRow(train[,selectionA0])
	  cvfit <- glmnet(as.matrix(tmp1),as.factor(train$label),family='binomial',alpha=0,lambda=hitFit$lambda.min)

	  tmp2 <- scaleRow(validation[,selectionA0])
     
	  pred <- data.frame(predict(cvfit,newx = as.matrix(tmp2),s = cvfit$lambda.1se,type="class"),stringsAsFactors=F)
	  pred$observed <- as.character(validation$label)
	  colnames(pred) <- c("predicted","observed")
	  acc <- sum(pred$predicted==pred$observed)/dim(tmp2)[1]*100
	  result <- sprintf("%s %4.3f",paste0(sort(selectionA0),collapse=","),acc)
     result
  }##function iterate

#itA <- iterate(A)
#itB <- iterate(B)

itorTimes <- 20
print("iteration for searching proteins for A ......")
R0 <- apply(A,2,function(v){sum(is.na(v))/length(v)*100})
resultA <- c()
for( miss in seq(5,100,by=5)){
	 for(i in 1:itorTimes){
       resultA <- rbind(resultA,iterate(A[,R0 <= miss],B)) 
	 }
}
resultA <- t(apply(resultA,1,function(v){strsplit(v," ")[[1]]}))
resultA <- data.frame(resultA,stringsAsFactors=F)
resultA <- resultA[order(resultA$X2,decreasing=T),]
rownames(resultA) <- 1:dim(resultA)[1]
resultA$protNumber <- as.integer(apply(resultA,1,function(v){length(strsplit(v,",")[[1]])}))

indxA <- c()
for(i in 1:dim(resultA)[1]){indxA <- c(indxA,which(resultA$protNumber > resultA[i,3] & resultA$X2 <= resultA[i,2]))}
optionA <- unique(resultA[- indxA,])
colnames(optionA) <- c('prots','accuracy','protNumber')
optionA <- sqldf("SELECT prots,max(accuracy) accuracy,max(protNumber) protNumber FROM optionA group by prots order by max(accuracy) desc")

protA <- c()
for(i in 1:dim(optionA)[1]){protA <- unique(c(protA,strsplit(optionA[i,1],",")[[1]]))}


print("iteration for searching proteins for B ......")
R0 <- apply(B,2,function(v){sum(is.na(v))/length(v)*100})
resultB <- c()
for( miss in seq(5,100,by=5)){
	 for(i in 1:itorTimes){
       resultB <- rbind(resultB,iterate(B[,R0 <= miss],A)) 
	 }
}

resultB <- t(apply(resultB,1,function(v){strsplit(v," ")[[1]]}))
resultB <- data.frame(resultB,stringsAsFactors=F)
resultB <- resultB[order(resultB$X2,decreasing=T),]
rownames(resultB) <- 1:dim(resultB)[1]
resultB$protNumber <- as.integer(apply(resultB,1,function(v){length(strsplit(v,",")[[1]])}))

indxB <- c()
for(i in 1:dim(resultB)[1]){indxB <- c(indxB,which(resultB$protNumber > resultB[i,3] & resultB$X2 <= resultB[i,2]))}
optionB <- unique(resultB[- indxB,])
colnames(optionB) <- c('prots','accuracy','protNumber')
optionB <- sqldf("SELECT prots,max(accuracy) accuracy,max(protNumber) protNumber FROM optionB group by prots order by max(accuracy) desc")

protB <- c()
for(i in 1:dim(optionB)[1]){protB <- unique(c(protB,strsplit(optionB[i,1],",")[[1]]))}
selectedA <- strsplit(optionA[1,1],",")[[1]]
selectedB <- strsplit(optionA[3,1],",")[[1]]

#################################### predict each other between dataset A and B ###############


Model <- function(M1,measureType='class'){
   mat <- scaleRow(M1[,colnames(M1)!='label'])
   cvfit <- cv.glmnet(as.matrix(mat),as.factor(M1$label),family='binomial',alpha=0,type.measure=measureType,standardize = F)
   cvfit
}
Test <- function(model,M1,predType='class'){
   mat <- scaleRow(M1[,colnames(M1)!='label'])
	pred <- NULL

	if(predType=='class'){
		pred <- data.frame(predict(model,newx = as.matrix(mat),s = model$lambda.min,type=predType))
		pred$observed <- as.character(M1$label)
		colnames(pred) <- c("predicted","observed")
	}else if(predType=='response'){
      lbls <- levels(as.factor(M1$label))[c(2,1)]
	   pred  <-  data.frame(t(apply(predict(model,newx = as.matrix(mat),s = model$lambda.min,type="response"),1,function(v){c(v,1-v)})))
	   colnames(pred) <- lbls
		pred$predicted <- apply(pred,1,function(v){names(v)[max(v)==v]})
      pred$observed <- M1$label
	}
	L <- dim(M1)[1]
   print(sprintf(" %d samples predicted £ºaccruacy=%4.3f",L, sum(pred$predicted==pred$observed)/L*100))
	pred
}
modelA <- Model(A[,c('label',selectedA)])
modelB <- Model(B[,c('label',selectedB)])

testA <- Test(modelA,B[,c('label',selectedA)])
testB <- Test(modelB,A[,c('label',selectedB)])

###################################################
avgMat <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(avgMat) <- avgMat$SpecimenID
avgMat <- avgMat[,-2]

avgMatA <-  avgMat[setdiff(rownames(avgMat),rownames(A)),c('label',selectedA)]
avgMatA$label[avgMatA$label %in% c('N','M','A')]  <- 'B'
avgMatA$label[avgMatA$label %in% c('C','P','W')]  <- 'M'

avgMatB <-  avgMat[setdiff(rownames(avgMat),rownames(B)),c('label',selectedB)]
avgMatB$label[avgMatB$label %in% c('N','M','A')]  <- 'B'
avgMatB$label[avgMatB$label %in% c('C','P','W')]  <- 'M'

avgTestA <- Test(modelA,avgMatA)
avgTestB <- Test(modelB,avgMatB)
###########################################################################
crossAssess <- function(features,mat1=NULL,mat2=NULL){
   if(!('label'  %in% colnames(mat1) & 'label' %in% colnames(mat2))){
      print('factor column named label not exists')
		return(NULL)
	}
	if(! all(features %in% colnames(mat1))){
	  printf(sprintf("features %s not exist in the first matrix",paste0(features[ features %in% colnames(mat1)],collapse=","))) 
	  return(NULL)    
	}
   if(! all(features %in% colnames(mat2))){
	  printf(sprintf("features %s not exist in the second matrix",paste0(features[ features %in% colnames(mat2)],collapse=","))) 
	  return(NULL)    
	}
	lblA <- as.factor(mat1$label)
	lblB <- as.factor(mat2$label)
   
	M1 <- scaleRow(mat1[,features])
	M2 <- scaleRow(mat2[,features])

   cvfit1 <- cv.glmnet(as.matrix(M1),lblA,family='binomial',alpha=0,type.measure="class")
   cvfit2 <- cv.glmnet(as.matrix(M2),lblB,family='binomial',alpha=0,type.measure="class")

   pred1 <- data.frame(predict(cvfit1,newx = as.matrix(M2),s = cvfit1$lambda.min,type="class"))
   pred1$observed <- as.character(lblB)
   colnames(pred1) <- c("predicted","observed")
	pred2 <- data.frame(predict(cvfit2,newx = as.matrix(M1),s = cvfit2$lambda.min,type="class"))
   pred2$observed <- as.character(lblA)
   colnames(pred2) <- c("predicted","observed")


   L1 <- dim(M1)[1]
	L2 <- dim(M2)[1]

	acc1 <- sum(pred1$predicted==pred1$observed)/L2
	acc2 <- sum(pred2$predicted==pred2$observed)/L1
   #print(sprintf("%d %d %4.4f  %4.4f",L1,L2,acc1,acc2))
	return(c(acc1,acc2,(acc1*L2/(L1+L2)+acc2*L1/(L1+L2)))*100)
}
###########################################
featuresA <- sapply(optionA$prots,function(v){strsplit(v,",")[[1]]})
featuresB <- sapply(optionB$prots,function(v){strsplit(v,",")[[1]]})
assess1 <- t(sapply(featuresA,crossAssess,mat1=A,mat2=B))
assess2 <- t(sapply(featuresB,crossAssess,mat1=B,mat2=A))

assessment <- rbind(assess1,assess2)
assessment <- assessment[order(assessment[,3],decreasing=T),]
optionAll <- data.frame(assessment,stringsAsFactors=F)
colnames(optionAll) <- c('accuracy1','accuracy2','avgAcc')
optionAll$prots <- rownames(optionAll)
rownames(optionAll) <- 1:dim(optionAll)[1]
optionAll$protNumber <- as.integer(sapply(optionAll$prots,function(v){length(strsplit(v,",")[[1]])}))
optionAll <- optionAll[,c('avgAcc','accuracy1','accuracy2','protNumber','prots')]
optionAll <- optionAll[optionAll$avgAcc > 70,]
###################################
###################### ROC and AUC for eache protein group
print("Now do 100 times 10-folds cross validation on all samples for ROC")
matAll <- avgMat
matAll$label[matAll$label %in% c('N','M','A')] <- 'B'
matAll$label[matAll$label %in% c('C','P','W')] <- 'M'


cvRigePredict <- function(features,mat=NULL,itors=10){
   lbls <- levels(as.factor(mat$label))[c(2,1)]
   predictions <- data.frame(matrix(0,nrow=dim(mat)[1],ncol=2));
   rownames(predictions) <- rownames(mat)
   colnames(predictions) <- lbls
   for(i in 1:itors){
     folds <- createFolds(mat$label,k=10)
     for(j in 1:10){
			train <- mat[unlist(folds[c(1:10)[-j]]),c('label',features)]
			valid <- scaleRow(mat[folds[[j]],features])
			model <- Model(train)
			pred  <-  data.frame(t(apply(predict(model,newx = as.matrix(valid),s = model$lambda.min,type="response"),1,function(v){c(v,1-v)})))
			colnames(pred) <- lbls
			predictions[rownames(pred),] <- predictions[rownames(pred),] + pred
     }
  }
  predictions <- predictions/itors
  predictions$predicted <- apply(predictions,1,function(v){names(v)[max(v)==v]})
  predictions$observed <- mat$label
  acc <- sum(predictions$predicted==predictions$observed)/dim(mat)[1]
  print(sprintf("%d predicted,acc=%4.3f %s",dim(mat)[1],acc,paste0(features,collapse=",")))
  predictions
}

cvRandomForest <- function(features,mat=NULL,itors=10){
   lbls <- levels(as.factor(mat$label))
   predictions <- data.frame(matrix(0,nrow=dim(mat)[1],ncol=2));
   rownames(predictions) <- rownames(mat)
   colnames(predictions) <- lbls

   for(i in 1:itors){
     folds <- createFolds(mat$label,k=10)
     for(j in 1:10){
	      indx   <-  unlist(folds[c(1:10)[-j]])
			labelA <- mat$label[indx]
         
			train  <- scaleRow(mat[indx,features])
			train$label <- as.factor(labelA)
			
			valid <- scaleRow(mat[folds[[j]],features])
         model <- randomForest(label ~ . ,data=train,importance=T,ntree=1000,nodesize=5)
	      pred <- predict(model,valid,type='prob')
	      predictions[rownames(pred),] <-  predictions[rownames(pred),] + pred
     }
  }
  predictions <- predictions/itors
  predictions$predicted <- apply(predictions,1,function(v){names(v)[max(v)==v]})
  predictions$observed <- mat$label
  acc <- sum(predictions$predicted==predictions$observed)/dim(mat)[1]
  print(sprintf("%d predicted,acc=%4.3f %s",dim(mat)[1],acc,paste0(features,collapse=",")))
  predictions
}
set.seed(1)
allSelected <- sapply(optionAll$prots,function(v){strsplit(v,",")[[1]]})
#rgPredicts <- list()
#for(i in 1:length(allSelected)){ rgPredicts[[length(rgPredicts)+1]] <- cvRigePredict(allSelected[[i]],mat=matAll)}
#t1 <- sapply(features,cvPredict,mat=mat)

#set.seed(1)
#rfPredicts <- list()
#for(i in 1:length(allSelected)){ rfPredicts[[length(rfPredicts)+1]] <- cvRandomForest(allSelected[[i]],mat=matAll)}
#source("src/zhe1Validation_glm_20190313.R")

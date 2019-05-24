rm(list=ls())
library(ggplot2)
library(caret)
library(umap)
library(tsne)
library(glmnet)


##############################
setwd("E:/projects/TPD")
source("src/common.R")
set.seed(190310)
allRep     <- read.table('data/TPD_prot_matrix_rep_20190304.txt',sep="\t",header=T,stringsAsFactors = F)
sampleInf <- read.table('data/AllSampleInformation_20190304.txt',sep="\t",header=T,stringsAsFactors = F)

pattern <- 'C201'
repTrain <- allRep[grepl(pattern,allRep$filename),]
specimens <- unique(repTrain$SpecimenID)
print(sprintf("processing replicas from machine %s.....",substr(pattern,1,1)))
tmp <- c()
for(sp in specimens){
  k0 <- apply(repTrain[repTrain$SpecimenID==sp,-c(1:3)],2,function(v){log2(mean(2^v,na.rm=T))})
  tmp <- rbind(tmp,k0)
}

tmp[is.na(tmp)] <- NA
rownames(tmp) <- specimens
df0 <- data.frame(tmp)
prots <- colnames(df0)

tmpLabel <- unique(repTrain[,c('SpecimenID','label')])
rownames(tmpLabel) <- tmpLabel$SpecimenID
label01 <- tmpLabel[rownames(df0),"label"]
print(table(label01))
label01[label01 %in% c('N','M','A')] <- 'B'
label01[label01 %in% c('C','P','W')] <- 'M'
print(table(label01))
label01 <- as.factor(label01)

repValidation <- allRep[!grepl(pattern,allRep$filename),]
specimens <- setdiff(unique(repValidation$SpecimenID),rownames(df0))

print("processing replicas from other machines .....")
tmp <- c()
for(sp in specimens){
  k0 <- apply(repValidation[repValidation$SpecimenID==sp,-c(1:3)],2,function(v){log2(mean(2^v,na.rm=T))})
  tmp <- rbind(tmp,k0)
}
rownames(tmp) <- specimens
tmp[is.na(tmp)] <- 0
df1 <- data.frame(tmp);

tmpLabel <- unique(repValidation[,c('SpecimenID','label')])
rownames(tmpLabel) <- tmpLabel$SpecimenID
label02 <- tmpLabel[rownames(df1),"label"]
print(table(label02))
label02[label02 %in% c('N','M','A')] <- 'B'
label02[label02 %in% c('C','P','W')] <- 'M'
print(table(label02))
label02 <- as.factor(label02)
######################################################################
color3=c(A="red",B="black",C="blue")
color2=c(M='red',B='blue')

#machine3 <- rbind(machineA,df0)

R0 <- apply(df0,2,function(v){sum(is.na(v))/length(v)*100})
print("feature selection iteration....")

iterations <- 10
selecteds <- list()

for(missingRate in seq(5,100,by=5)){ 
    #set.seed(190310)
    training <- t(apply(df0[,R0 <= missingRate],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}));
    training[is.na(training)] <- 0;

	 itorAcc <- 0
	 itorSelected <- c()
    for(itor in 1:iterations){

       selection0 <- colnames(training)
       hitFit <- NULL
       mse=1000
       flg=T
       fits <- list()
       while(flg){
				tmp <- df0[,selection0]
				tmp    <- t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
				tmp[is.na(tmp)] <- 0

				cvfit <- cv.glmnet(tmp,label01,family='binomial',alpha=0.5,type.measure='class')
				fits[[length(fits)+1]] <- cvfit
				cf <- coef(cvfit, s = "lambda.1se")
				selection1 <- rownames(cf)[cf[,1]!=0][-1]
				if(min(cvfit$cvm) <= mse){
				    mse     <-  min(cvfit$cvm);
				    hitFit <-  cvfit
				}
				if(length(selection1)>= 4){ ## 4 proteins at least 
				flg <- length(selection1) != length(selection0)
				selection0 <- selection1
				}else{  flg <- F	 }

		 }#flg
       sts <- sapply(fits,function(obj){
				cf <- coef(obj,s = "lambda.1se");
				a <- sum(cf[,1]!= 0)-1;
				b=min(obj$cvm);
				c(a,b)
       }) 

       #print(stsA)
       obj <- hitFit
       cf <- coef(obj, s = "lambda.1se")
       selected <- rownames(cf)[cf[,1]!=0][-1]
		 
	    #print(sprintf("%d proteins selected with missing rate=%d",length(selected),missingRate))

       tmpA <- df0[,selected]
	    tmpA    <- t(apply(tmpA,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	    tmpA[is.na(tmpA)] <- 0
	    cvfit <- glmnet(tmpA,label01,family='binomial',alpha=0,lambda=hitFit$lambda.min)

			
		 test  <-  df1[,selected]
		 test    <- t(apply(test,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
		 test[is.na(test)] <- 0
	
		 pred <- data.frame(predict(cvfit,newx = test,s = cvfit$lambda.1se,type="class"),stringsAsFactors=F)
		 pred$observed <- as.character(label02)
		 colnames(pred) <- c("predicted","observed")
		 acc <- sum(pred$predicted==pred$observed)/dim(test)[1]*100
		 itorSelected <- c(itorSelected,paste0(sort(selected),collapse=","))
		 if(length(selected) <=25){ print(sprintf("%s %4.3f %d",paste0(sort(selected),collapse=","),acc,missingRate))}
		 #print(sprintf("predicted £ºacc=%4.3f,missing cutoff=%d ",acc,missingRate ))
       itorAcc <- itorAcc + acc
	}
	   print(sprintf("************************************************* average acc %4.3f with missing rate %d",itorAcc/iterations,missingRate))
	   #strMain <- sprintf('PCA:  %4.3f missing %d;%d proteins',acc,missingRate,length(selected))
	   #tmp <- df0[,selected]
	   #tmp$label <- labelA
      #pcaA <-  drawPCA(tmp,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=strMain)
      #print(pcaA)

}
#library(randomForest)
#RF <- 

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
reps0 <- read.table('data/TPD_prot_matrix_rep_20190304.txt',sep="\t",header=T,stringsAsFactors = F)
sampleInf <- read.table('data/AllSampleInformation_20190304.txt',sep="\t",header=T,stringsAsFactors = F)

repsA <- reps0[grepl('A2018',reps0$filename),]
specimens <- unique(repsA$SpecimenID)
print("processing replicas from machine A.....")
tmp <- c()
for(sp in specimens){
  k0 <- apply(repsA[repsA$SpecimenID==sp,-c(1:3)],2,function(v){log2(mean(2^v,na.rm=T))})
  tmp <- rbind(tmp,k0)
}
tmp[is.na(tmp)] <- NA
rownames(tmp) <- specimens
machineA <- data.frame(tmp)
prots <- colnames(machineA)

specimenLabelsA <- unique(repsA[,c('SpecimenID','label')])
rownames(specimenLabelsA) <- specimenLabelsA$SpecimenID
labelA <- specimenLabelsA[rownames(machineA),"label"]
labelA[labelA %in% c('N','M','A')] <- 'B'
labelA[labelA %in% c('C','P','W')]     <- 'M'
labelA <- as.factor(labelA)

repsBC <- reps0[!grepl('A201',reps0$filename),]
specimensBC <- setdiff(unique(repsBC$SpecimenID),rownames(machineA))

print("processing replicas from machine B C .....")
tmp <- c()
for(sp in specimensBC){
  k0 <- apply(repsBC[repsBC$SpecimenID==sp,-c(1:3)],2,function(v){log2(mean(2^v,na.rm=T))})
  tmp <- rbind(tmp,k0)
}
rownames(tmp) <- specimensBC
tmp[is.na(tmp)] <- 0
machineBC <- data.frame(tmp);

specimenLabelsBC <- unique(repsBC[,c('SpecimenID','label')])
rownames(specimenLabelsBC) <- specimenLabelsBC$SpecimenID
labelBC <- specimenLabelsBC[rownames(machineBC),"label"]
labelBC[labelBC %in% c('N','M','A')] <- 'B'
labelBC[labelBC %in% c('C','P','W')] <- 'M'
labelBC <- as.factor(labelBC)
######################################################################
color3=c(A="red",B="black",C="blue")
color2=c(M='red',B='blue')

#machine3 <- rbind(machineA,machineB)


R0 <- apply(machineA,2,function(v){sum(is.na(v))/length(v)*100})
print("feature selection iteration....")

iterations <- 10
selecteds <- list()

for(missingRate in seq(5,100,by=5)){ 
    #set.seed(190310)
    training <- t(apply(machineA[,R0 <= missingRate],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}));
    training[is.na(training)] <- 0;

	 itorAcc <- 0
	 itorSelected <- c()
    for(itor in 1:iterations){

       selection0 <- colnames(training)
       hitFitA <- NULL
       mse=1000
       flg=T
       fitsA <- list()
       while(flg){
				tmp <- machineA[,selection0]
				tmp    <- t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
				tmp[is.na(tmp)] <- 0

				cvfit <- cv.glmnet(tmp,labelA,family='binomial',alpha=0.5,type.measure='class')
				fitsA[[length(fitsA)+1]] <- cvfit
				cf <- coef(cvfit, s = "lambda.1se")
				selection1 <- rownames(cf)[cf[,1]!=0][-1]
				if(min(cvfit$cvm) <= mse){
				mse    <-  min(cvfit$cvm);
				hitFitA <-  cvfit
				}
				if(length(selection1)>= 4){ ## 4 proteins at least 
				flg <- length(selection1) != length(selection0)
				selection0 <- selection1
				}else{  flg=F	 }

		 }#flg
       stsA <- sapply(fitsA,function(obj){
				cf <- coef(obj,s = "lambda.1se");
				a <- sum(cf[,1]!= 0)-1;
				b=min(obj$cvm);
				c(a,b)
       })
       #print(stsA)
       obj <- hitFitA
       cf <- coef(obj, s = "lambda.1se")
       selectedA <- rownames(cf)[cf[,1]!=0][-1]
		 
	    #print(sprintf("%d proteins selected with missing rate=%d",length(selectedA),missingRate))

       tmpA <- machineA[,selectedA]
	    tmpA    <- t(apply(tmpA,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	    tmpA[is.na(tmpA)] <- 0
	    cvfitA <- glmnet(tmpA,labelA,family='binomial',alpha=0,lambda=hitFitA$lambda.min)

			
		 test  <-  machineBC[,selectedA]
		 test    <- t(apply(test,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
		 test[is.na(test)] <- 0
	
		 pred <- data.frame(predict(cvfitA,newx = test,s = cvfitA$lambda.1se,type="class"),stringsAsFactors=F)
		 pred$observed <- as.character(labelBC)
		 colnames(pred) <- c("predicted","observed")
		 acc <- sum(pred$predicted==pred$observed)/dim(test)[1]*100
		 itorSelected <- c(itorSelected,paste0(sort(selectedA),collapse=","))
		 if(length(selectedA) <=25){ print(sprintf("%s %4.3f %d",paste0(sort(selectedA),collapse=","),acc,missingRate))}
		 #print(sprintf("predicted £ºacc=%4.3f,missing cutoff=%d ",acc,missingRate ))
       itorAcc <- itorAcc + acc
	}
	   print(sprintf("************************************************* average acc %4.3f with missing rate %d",itorAcc/iterations,missingRate))
	   #strMain <- sprintf('PCA:  %4.3f missing %d;%d proteins',acc,missingRate,length(selectedA))
	   #tmp <- machineA[,selectedA]
	   #tmp$label <- labelA
      #pcaA <-  drawPCA(tmp,ptColors = color2,rowNormalization = T,colNormalization = T,strTitle=strMain)
      #print(pcaA)

}
#library(randomForest)
#RF <- 

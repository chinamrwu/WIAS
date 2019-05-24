remove(list=ls())
library(randomForest)
library(caret)
library(data.table)
options(width=350)
set.seed(20190122)


df0 <- read.table("E:/projects/TPD/data/TPD_prot_matrix_avg_20190115.txt",sep="\t",header=T,stringsAsFactors=F)

t0 <- df0[df0$label %in% c('M','A'),-1]

R1 <- apply(t0,2,function(v){100*sum(is.na(v))/length(v)}) ## remove  columns with more than 5% missing values
t0 <- t0[,R1<=5]
R0 <- apply(t0[,-1],1,function(v){100*sum(is.na(v))/length(v)}) ## remove rows with more than 5% missing values
t0 <- t0[R0<=5,] ### 200X331
#########################  imputation by kNN conditioned given class   #########################################


K=3
t1 <- sapply(1:dim(t0)[1],function(i) {
     
     lbl <- t0[i, 1]
     v0  <- t0[i,-1]
     feature0 <- names(v0)[is.na(v0)]
     if(length(feature0)>0){
	     k0  <- t0[-i,]
	     k0  <- k0[k0$label==lbl,-1]
	     k0  <- rbind(v0,k0)
	     k0 <- apply(k0,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
             k0[is.na(k0)] <- 0

	     D0  <- as.matrix(dist(k0))
	     d0  <- order(D0[1,])[2:(K+1)]
	     if(length(feature0)==1){
	        v0[feature0] <- mean(t0[d0,feature0],na.rm=T)
	     }else{
                v0[feature0] <- apply(t0[d0,feature0],2,function(v){mean(v,na.rm=T)})
	     }
      }
     v0
})

t1 <- apply(t(t1),2,unlist)    
rownames(t1) <- rownames(t0)
t1 <- data.frame(t1)
protIds <- colnames(t1)
t1$label <- t0$label
t1 <- t1[,c('label',protIds)]
################  RF Feature Selection #######################

RFScore <- function(feature,M1=t1,nFolds=5,itNumber=5){
   tmp <- data.frame(t(apply(M1[,feature],1,function(v){(v-mean(v))/sd(v)})))
   tmp$label <- M1$label
   tmp <- tmp[,c('label',feature)]
   
   avgAcc <- 0;
   for(i in 1:itNumber){
      acc <- 0;
      folds <- createFolds(tmp$label,nFolds)
      for(fold in folds){
        valids <- tmp[fold,]
        trains <- tmp[setdiff(1:dim(tmp)[1],fold),]
	trains$label <- as.factor(trains$label)
        tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=7)
	predicted <- predict(tmpRF,valids,type='response')
        acc <- acc+sum(tmp$label[fold]==predicted)
      }
      acc <- acc/dim(tmp)[1]*100
      avgAcc <- avgAcc+acc
   #print(sprintf("%s %4.3f",paste0(feature,collapse=","),avgAcc/10))
   }
   avgAcc <- avgAcc /itNumber
   rtv <- sprintf("%s %4.3f",paste0(feature,collapse=" "),avgAcc)
   print(rtv)
   rtv
}
#############################################################################
tmp <- t1;
tmp$label <- as.factor(tmp$label)
RF <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=5)
importance <- RF$importanceSD
importance <- importance[order(importance[,3],decreasing=T),]
orderProtIds <- rownames(importance)
triplets <- combn(orderProtIds,3)
result <- c()
L <- dim(triplets)[2]
#print("*************************************************  M vs A**************************************")
#for(i in 1:L){result <- rbind(result,RFScore(triplets[,1:]))}
for(i in seq(L,1,by=-1)){result <- rbind(result,RFScore(triplets[,i]))}

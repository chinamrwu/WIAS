remove(list=ls())
library(randomForest)
library(caret)
library(data.table)

df0 <- read.table("E:/projects/TPD/data/TPD_prot_matrix_avg_20190115.txt",sep="\t",header=T,stringsAsFactors=F)

t0 <- df0[df0$label %in% c('N','M'),-1]

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

RFScore <- function(feature,M1=t1,nFolds=10){
   tmp <- data.frame(t(apply(M1[,feature],1,function(v){(v-mean(v))/sd(v)})))
   tmp$label <- M1$label
   tmp <- tmp[,c('label',feature)]
   clsses <- unique(tmp$label)
   rownames(tmp) <- sapply(1:dim(tmp)[1],function(v){paste0('R',v)})
   sampleNumber <- sapply(clsses,function(v){sum(tmp$label==v)})
   maxLabel <- names(sampleNumber)[max(sampleNumber)==sampleNumber]
   minLabel <- names(sampleNumber)[min(sampleNumber)==sampleNumber]
   size <- min(sampleNumber)
   
   avgAcc <- 0;
   for(i in 1:100){
      indx1 <- sample(rownames(tmp)[tmp$label==maxLabel],size,replace=F)
      indx2 <- rownames(tmp)[tmp$label==minLabel]

      tmp1 <- tmp[c(indx1,indx2),]
      
      acc <- 0;

      folds <- createFolds(tmp1$label,nFolds)
      for(fold in folds){
        valids <- tmp1[fold,]
        trains <- tmp1[setdiff(1:dim(tmp1)[1],fold),]
	trains$label <- as.factor(trains$label)
        tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=500,nodesize=5)
	predicted <- predict(tmpRF,valids,type='response')
        acc <- acc+sum(tmp1$label[fold]==predicted)
      }
      acc <- acc/dim(tmp1)[1]*100
      avgAcc <- avgAcc+acc
   #print(sprintf("%s %4.3f",paste0(feature,collapse=","),avgAcc/10))
   }
   avgAcc <- avgAcc /100
   rtv <- sprintf("%s %4.3f",paste0(feature,collapse=" "),avgAcc)
   print(rtv)
   rtv
   
}


features <- colnames(t1)[-1]
M0 <- t1[t1$label=='N',-1]
#M0 <- data.frame(t(apply(M0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
L0 <- dim(M0)[1]
M1 <- t1[t1$label=='M',-1]
#M1 <- data.frame(t(apply(M1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
L <- dim(t1)[1]

tmp <- t1
tmp$label <- as.factor(t1$label)

fullRF <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=7)
importance <- fullRF$importanceSD
importance <- importance[order(importance[,3],decreasing=T),]

f3s <- combn(rownames(importance),3)



 f1scores <-  FisherScore(features,1)
 f2scores <-  FisherScore(features,2)
 f3scores <-  FisherScore(f1scores$F1[1:50],3)
 f4scores <-  FisherScore(f1scores$F1[1:50],4)
 f5scores <-  FisherScore(f1scores$F1[1:50],5)
 f6scores <-  FisherScore(f1scores$F1[1:30],6)

##############################################################


############################ random record cmds
plots=list()
for(i in 1:30){
    plots[[length(plots)+1]] <-  drawPCA(t1[, c('label',unique(as.vector(unlist(f6scores[i,1:6]))))],'Fisher Scores',T,T)
    plots[[length(plots)+1]] <- plotTSNE(t1[,c('label',unique(as.vector(unlist(f6scores[i,1:6]))))],'Fisher Scores',T,T)
}

pdf("E:/projects/TPD/results/feature_selection_Fisher_6_proteins_PCA_tSNE.pdf",width=12,height=5)
for(i in seq(1,60,by=2)){
    grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 3)));
    print(plots[[i]],vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(plots[[i+1]],vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
}
dev.off()
#################################  cross validation for Random Forest
indxA <- which(t1$label=='A')
indxC <- which(t1$label=='C')


avgAcc=c()
for(p in 1:500){
   predictions <- matrix(0,nrow=dim(t1)[1],ncol=2)
   rownames(predictions) <- 1:dim(t1)[1]
   colnames(predictions) <- c('A','C')

   predictions <- data.frame(predictions)
   acc <- matrix(0,nrow=100,ncol=2)
   for(i in 1:100){
     tA <- sample(indxA,48,replace=F)
     vA <- setdiff(indxA,tA)
     tC <- sample(indxC,48,replace=F)
     vC <- setdiff(indxC,tC)
     dsT <- t1[c(tA,tC),c('label',as.character(f6scores[p,1:6]))]
     dsV <- t1[c(vA,vC),as.character(f6scores[p,1:6])]
     rownames(dsV) <- c(vA,vC)
     dsT$label <- as.factor(dsT$label)
     RF <- randomForest(label ~ . ,data=dsT,importance=T,ntree=1000,nodesize=5)
     predicted <- predict(RF,dsV,type='response')
     acc[i,1] <- sum(t1[c(vA,vC),'label']==predicted)/length(predicted)*100
     for(j in 1:length(predicted)){predictions[names(predicted[j]),as.character(predicted[j])] <- predictions[names(predicted[j]),as.character(predicted[j])]+1}
     
     K=5
     mn <-  apply(dsT[,-1],2,function(v){mean(v,na.rm=T)})
     sdn <- apply(dsT[,-1],2,function(v){sd(v,na.rm=T)})
     tmp <- (dsV-mn)/sdn

     dsT[,-1] <- apply(dsT[,-1],2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})

     L=dim(dsV)[1]
     D0 <- as.matrix(dist(rbind(tmp,dsT[,-1])))
     D0 <- D0[1:L,(L+1):dim(D0)[2]]
     nn <- apply(D0,1,function(v){
        lbls <- dsT$label[order(v)[1:K]];
	ifelse(sum(lbls=='A')>sum(lbls=='C'),'A','C')
     })
     acc[i,2] <- sum(t1[c(vA,vC),'label']==nn)/length(nn)*100

    
 }
  sts=apply(acc,2,mean)
  avgAcc <- rbind(avgAcc,sts)
  print(sprintf("%s %4.3f  %4.3f",paste0(f6scores[p,1:6],collapse=","),sts[1],sts[2]))
}
#################################################  sampling for RF
sampleRF <- function(M1,feats,foldNumber=5,itNumber=100,rowNormalizatin=F){
   if(!is.element('label',colnames(M1))){
      print(" A column named label is required to indicate the classes of samples")
      return(NULL)
   }
   if(! all(feats %in% colnames(M1))) {
     unfound <- setdiff(feats,colnames(M1))
     print(sprintf('found undefined features:%s',paste0(unfound,collapse=",")))
     return(NULL)
   }
   clsses <- unique(M1$label)
   
   predictions <- matrix(0,nrow=dim(M1)[1],ncol=2)
   rownames(predictions) <- rownames(M1)
   colnames(predictions) <- clsses
   predictions <- data.frame(predictions)


   sampleIndex   <- sapply(clsses,function(v){which(v==M1$label)})
   sampleNumber  <- sapply(sampleIndex,length)
   maxLabel <- names(sampleNumber)[sampleNumber==max(sampleNumber)]
   minLabel <- names(sampleNumber)[sampleNumber==min(sampleNumber)]
   sampleSize    <- min(sampleNumber)
   
   
   for(i in 1:itNumber){
     spIndx  <- sample(sampleIndex[[maxLabel]],sampleSize,replace=F,prob=sampleProb)
     spIndx0 <- setdiff(sampleIndex[[maxLabel]],spIndx)
     sampleProb[match(spIndx,sampleIndex[[maxLabel]])] <- sampleProb[match(spIndx,sampleIndex[[maxLabel]])]-1

     rowIndx <- c(sampleIndex[[minLabel]],spIndx)
     tmp <- M1[rowIndx,]

     folds <- createFolds(tmp$label,foldNumber)
     for(fold in folds){
        DSTest  <- rbind(M1[spIndx0,feats],tmp[fold,feats])
	observed <- c(M1$label[spIndx0],tmp$label[fold])
	DSTrain <- tmp[setdiff(1:dim(tmp)[1],fold),c('label',feats)]
	if(rowNormalization){
          DSTest   <- t(apply(DSTest,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	  DSTrain[,-1]  <- t(apply(DSTrain[,-1],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	}
	DSTrain$label <- as.factor(DSTrain$label)
	RF <- randomForest(label ~ . ,data=DSTrain,importance=T,ntree=1000,nodesize=5)
        predicted <- predict(RF,DSTest,type='response')
        for(j in 1:length(predicted)){pred[names(predicted[j]),as.character(predicted[j])] <- pred[names(predicted[j]),as.character(predicted[j])]+1}
     }
     predictions <- predictions + pred
     pred1 <- apply(pred,1,function(v){ifelse(v[1]>0 & v[1]>v[2],names(v)[1],names(v)[2])})
     #print(sprintf("%3d  acc=%4.4f",i,sum(pred1==M1$label)/dim(M1)[1]*100))
  }
   acc <- apply(predictions,1,function(v){ifelse(v[1]>v[2],names(v)[1],names(v)[2])})
   acc <- 100*sum(acc==M1$label)/dim(M1)[1]
   print(sprintf('%s: acc=%4.4f',paste0(feats,collapse=" "),acc))
   predictions
  }


##########  for 
## 
prot11 <- c('Q14103','O00468','Q12906','Q6PCB0','P13489','P08133','P22314','P02765','P09960','P62750','P07108')
cmb6 <-  combn(prot11,6)
cmb7 <-  combn(prot11,7)
cmb8 <-  combn(prot11,8)
cmb9 <-  combn(prot11,9)
cmb10 <- combn(prot11,10)

tmp <- sampleRF(t1,prot11,itNumber=10)
for(i in 1:dim(cmb10)[2]){ tmp <- sampleRF(t1,cmb10[,i],itNumber=10)}
for(i in 1:dim(cmb9)[2]){ tmp <- sampleRF(t1,cmb9[,i],itNumber=10)}
for(i in 1:dim(cmb8)[2]){ tmp <- sampleRF(t1,cmb8[,i],itNumber=10)}
for(i in 1:dim(cmb7)[2]){ tmp <- sampleRF(t1,cmb7[,i],itNumber=10)}
 

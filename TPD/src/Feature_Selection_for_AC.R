remove(list=ls())
library(randomForest)
library(caret)
library(data.table)

df0 <- read.table("E:/projects/TPD/data/TPD_prot_matrix_avg_20190115.txt",sep="\t",header=T,stringsAsFactors=F)
df0$label[df0$label=='W'] <- 'C'

R1 <- apply(df0,2,function(v){100*sum(is.na(v))/length(v)})
t0 <- df0[,R1<=5]
##################################################################
kNN <- function(M,features,K=7){
    tmp <- data.frame(M[,features])
    colnames(tmp) <- features
    tmp <- apply(tmp,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
    minv=min(tmp)
    tmp[is.na(tmp)] <- minv
    DM  <- as.matrix(dist(tmp))
    mxv=max(DM)
    for(i in 1:dim(DM)[1]){DM[i,i] <- mxv}
    clss <- unique(M$label)
    p1 <- apply(DM,2,function(v){  
      neighborLablels <- M$label[order(v)[1:K]]
      lblCounts <- sapply(clss,function(ch){ sum(neighborLablels==ch)})
      rtv       <- lblCounts;
      a         <- sum(lblCounts==max(lblCounts));
      ifelse(a==1,rtv <- c(rtv,names(lblCounts)[lblCounts==max(lblCounts)]),rtv <- c(rtv,neighborLablels[1]))
      rtv
    })
   p1 <- data.frame(t(p1))
   colnames(p1) <- c(clss,'predicted')
   p1$observed <- M$label
   p2 <- sapply(clss,function(ch){indx <- which(p1$observed==ch);sum(p1$observed[indx]==p1$predicted[indx])})
   p2 <- rbind(p2,sapply(clss,function(ch){sum(M$label==ch)}))
   p2 <- data.frame(p2)
   p2$total <- apply(p2,1,sum)
   p2 <- rbind(p2,apply(p2,2,function(v){100*v[1]/v[2]}))
   result=list()
   result$detail <- p1;
   result$stat <- p2
   result
}
#################################################################
growRF <- function(M,it=100){
   M1 <- M
   M1[is.na(M1)] <- 0
   M1$label <- as.factor(M1$label)

   fullRF <- randomForest(label ~ . ,data=M1,importance=T,ntree=1000,nodesize=7)
   imps  <- data.frame(importance(fullRF));
   impScore <- imps$MeanDecreaseAccuracy
   imps <- imps[order(impScore,decreasing=T),]
   orderedFeatures <- rownames(imps)
   
   result <- list()

   indx <- 1
   featureIndex <- c()
   bestScore <- 0
   while ( indx <= it ){
      score <- 0
      currentFeatures <- c("label",orderedFeatures[1:indx])
      kk <- kNN(M,orderedFeatures[1:indx],K=5)
      score <- kk$stat[3,3]
      if(score>bestScore){
          bestScore <- score
          featureIndex <- c(featureIndex,indx)
      }
      indx <- indx+1
   }
    
   result$protIds <- orderedFeatures[featureIndex]
   result$score <- bestScore
   result
}

#########################  imputation by kNN conditioned given class   #########################################
indxA <- which(t0$label=='A')
indxC <- which(t0$label=='C')

R0 <- apply(t0[,-1],1,function(v){100*sum(is.na(v))/length(v)})
t0 <- t0[R0<=5,] ### 200X331

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
################  Fisher Feature Selection #######################

features <- colnames(t1)[-1]
M0 <- t1[t1$label=='A',-1]
#M0 <- data.frame(t(apply(M0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
L0 <- dim(M0)[1]
M1 <- t1[t1$label=='C',-1]
#M1 <- data.frame(t(apply(M1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
L <- dim(t1)[1]

FisherScore <- function(mat,featureNames,K){
  lbls <- unique(mat$label)
  f <- colnames
  M0 <- mat[mat$label==lbls[1],featureNames]
  M1 <- mat[mat$label==lbls[2],featureNames] 
  L0 <- dim(M0)[1]
  L  <- dim(mat)[1]

  fcmb <- combn(featureNames,K);
  if(dim(fcmb)[1]==1){ colnames(fcmb) <- featureNames }
  
  scores <- apply(fcmb,2,function(f1){
    w=length(f1)
    tmp <- c()
    if(w==1){ tmp  <- c(M0[,f1],M1[,f1])    }
    else{     tmp <- rbind(M0[,f1],M1[,f1]) }
    
    D0 <- as.matrix(dist(tmp))/sqrt(w)
    mx <- max(D0)+1
    for(i in 1:dim(D0)[1]){D0[i,i] <- mx}
    D1 <- mean(apply(D0[1:L0,(L0+1):L],1,min))
    D2 <- mean(apply(D0[1:L0,(L0+1):L],2,min))
    score <- (D1+D2)/(mean(apply(D0[1:L0,1:L0],1,min))+mean(apply(D0[(L0+1):L,(L0+1):L],2,min)))
    score
  })
  fcmb        <- data.frame(t(fcmb),stringsAsFactors=F);
  fcmb$scores <- scores
  fcmb <- fcmb[order(scores,decreasing=T),]
  colnames(fcmb)[1:(dim(fcmb)[2]-1)] <- as.character(sapply(1:(dim(fcmb)[2]-1),function(i){paste0("F",i)}))
  fcmb
 }

 f1scores <-  FisherScore(features,1)
 f2scores <-  FisherScore(features,2)
 f3scores <-  FisherScore(f1scores$F1[1:30],3)
 f4scores <-  FisherScore(f1scores$F1[1:30],4)
 f5scores <-  FisherScore(f1scores$F1[1:30],5)
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
 

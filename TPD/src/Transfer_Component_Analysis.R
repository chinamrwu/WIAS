remove(list=ls())
library("data.table")
library("R.matlab")
library("mice")
options(width=350)
set.seed(20181224)

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
rownames(SWT) <- SWT$patientId
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
trainM <- SWT[,c('patientId',clnames)]
trainM$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
trainM <- trainM[,c("label",clnames)]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
testM <- SWV[,clnames]
protIds <- intersect(colnames(trainM),colnames(testM))

trainM <- trainM[,protIds]
testM  <-  testM[,protIds]


########################### imputation using mice
M0 =trainM

imp01 <- mice(data = M0, m = 5, method = "pmm", maxit = 50, seed = 500)
comp01 <- complete(imp01,1)
mean01 <- apply(comp01,2,mean)
sd01   <- apply(comp01,2,sd)
comp01 <- apply(comp01,2,function(v){(v-mean(v))/sd(v)})


M1 =testM
imp02 <- mice(data = M1, m = 5, method = "pmm", maxit = 50, seed = 500)
comp02 <- complete(imp02,1)
comp02 <- (comp02-mean01)/sd01

writeMat("E:/projects/TPD/data/trainM_5percentage_mice_imputate.mat",train=as.matrix(comp01))
writeMat("E:/projects/TPD/data/testM_5percentage_mice_imputate.mat",test=as.matrix(comp02))

total <- rbind(comp01,comp02)
writeMat("E:/projects/TPD/data/TPD_5percentage_mice_imputate.mat",TPD=as.matrix(total))


label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))


labeled <- matrix(0,nrow=580,ncol=2)
labeled[1:399,1] <- 1

labeled[which(label=='N'),2] <- 1
labeled[which(label=='M'),2] <- 2
labeled[which(label=='A'),2] <- 3
labeled[which(label=='C'),2] <- 4
labeled[which(label=='P'),2] <- 5

writeMat("E:/projects/TPD/data/TPD_5percentage_mice_imputate_labels.mat",Labels=labeled)


################################################################# the following codes 
df0 <- read.table("E:/projects/TPD/results/TPD_TCA_transformed.txt",sep="\t",header = F,stringsAsFactors = F)
df0$label=rep('training',dim(df0)[1])
df0$label[400:580] <- rep('validation',181)
drawPCA(df0,'PCA after domain adaption')

total$label <- df0$label
drawPCA(df0,'PCA after domain adaption')

###############################################################################################
df0 <- read.table("E:/projects/TPD/results/TPD_TCA_transformed.txt",sep="\t",header = F,stringsAsFactors = F)
df0 <- data.frame(t(apply(df0,1,function(v){(v-mean(v))/sd(v)})))
colnames(df0) <- colnames(trainM)


M0 <- df0[1:399,]
M0 <- apply(M0,2,function(v){(v-mean(v))/sd(v)})
DM0 <- as.matrix(dist(M0))

L=399
mx=max(DM0)+1
for(i in 1:L){ DM0[i,i] <- mx}


lbls <- strsplit('NMACP',"")[[1]]
K=7
kNN0 <- apply(DM0,2,function(v){  
  tmp <- label[order(v)[1:K]]
  lblCounts <- sapply(lbls,function(ch){ sum(tmp==ch)})
  rtv       <- lblCounts;
  a         <- sum(lblCounts==max(lblCounts));
  ifelse(a==1,rtv <- c(rtv,names(lblCounts)[lblCounts==max(lblCounts)]),rtv <- c(rtv,tmp[1]))
  rtv
})
kNN0 <- data.frame(t(kNN0))
kNN0$observed <- label
kNN0$patientId <- SWT$patientId
colnames(kNN0) <- c(lbls,c('predicted','observed','patientId'))
kNN0 <- kNN0[,c(8,1:5,7,6)]
write.table(kNN0,file="E:/projects/TPD/results/kNN_imputed_218proteins_TCAed.txt",sep="\t",col.names=T,row.names=F,quote=F)

M1 <- apply(trainM,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
DM1 <- as.matrix(dist(M1))
L=399
mx=max(DM1)+1
for(i in 1:L){ DM1[i,i] <- mx}


kNN1 <- apply(DM1,2,function(v){  
  tmp <- label[order(v)[1:K]]
  lblCounts <- sapply(lbls,function(ch){ sum(tmp==ch)})
  rtv       <- lblCounts;
  a         <- sum(lblCounts==max(lblCounts));
  ifelse(a==1,rtv <- c(rtv,names(lblCounts)[lblCounts==max(lblCounts)]),rtv <- c(rtv,tmp[1]))
  rtv
})
kNN1 <- data.frame(t(kNN1))
kNN1$observed <- label
kNN1$patientId <- SWT$patientId
colnames(kNN1) <- c(lbls,c('predicted','observed','patientId'))
kNN1 <- kNN1[,c(8,1:5,7,6)]
write.table(kNN1,file="E:/projects/TPD/results/kNN_imputed_218proteins.txt",sep="\t",col.names=T,row.names=F,quote=F)




VM0 <- df0[400:580,]
VD0 <- as.matrix(dist(rbind(M1,M0)))
##############################################  random forest  ########################################################### 
df0 <- read.table("E:/projects/TPD/results/TPD_TCA_transformed.txt",sep="\t",header = F,stringsAsFactors = F)
df0 <- data.frame(t(apply(df0,1,function(v){(v-mean(v))/sd(v)})))
colnames(df0) <- colnames(trainM)


M0 <- df0[1:399,]

M0$label <- label
library(randomForest)

source("E:/projects/TPD/src/growRF.R")
loopSize=100
modelNM=list();
modelMA=list();
modelAC=list();
modelCP=list();
for(i in 1:loopSize){  modelNM[[length(modelNM)+1]] <-  growRF01('NM',M0)}
for(i in 1:loopSize){  modelMA[[length(modelMA)+1]]  <- growRF01('MA',M0)}
for(i in 1:loopSize){  modelAC[[length(modelAC)+1]]  <- growRF01('AC',M0)}
for(i in 1:loopSize){  modelCP[[length(modelCP)+1]]  <- growRF01('CP',M0)}


models <- list()
models$NM <- modelNM
models$MA <- modelMA
models$AC <- modelAC
models$CP <- modelCP

protIds <- colnames(df0)
protWeights <- matrix(0,nrow=length(protIds),ncol=4)
colnames(protWeights) <- c('NM','MA','AC','CP')
rownames(protWeights) <- protIds
protWeights <- data.frame(protWeights)
for(cl in colnames(protWeights)){
      for(md in models[[cl]]){
         protWeights[md$protIds,cl] <- protWeights[md$protIds,cl]+1
  }
}
protIds <- rownames(protWeights[protWeights$NM>=20,])
protIds <- c(protIds,rownames(protWeights[protWeights$MA>=20,]))
protIds <- c(protIds,rownames(protWeights[protWeights$AC>=20,]))
protIds <- c(protIds,rownames(protWeights[protWeights$CP>=20,]))
protIds <- unique(protIds)

sm <- apply(protWeights,1,sum)
sm <- sm[sm>0]
sm <- sort(sm,decreasing=T)
protWeights <- protWeights[names(sm),]





M0   <- M0[,protIds]
M0 <- apply(M0,2,function(v){(v-mean(v))/sd(v)})
DM0 <- as.matrix(dist(M0))

L=399
mx=max(DM0)+1
for(i in 1:L){ DM0[i,i] <- mx}


lbls <- strsplit('NMACP',"")[[1]]
K=7
kNN0 <- apply(DM0,2,function(v){  
  tmp <- label[order(v)[1:K]]
  lblCounts <- sapply(lbls,function(ch){ sum(tmp==ch)})
  rtv       <- lblCounts;
  a         <- sum(lblCounts==max(lblCounts));
  ifelse(a==1,rtv <- c(rtv,names(lblCounts)[lblCounts==max(lblCounts)]),rtv <- c(rtv,tmp[1]))
  rtv
})
kNN0 <- data.frame(t(kNN0))
kNN0$observed <- label
kNN0$patientId <- SWT$patientId
colnames(kNN0) <- c(lbls,c('predicted','observed','patientId'))
kNN0 <- kNN0[,c(8,1:5,7,6)]
write.table(kNN0,file="E:/projects/TPD/results/kNN_imputed_50proteins_TCAed.txt",sep="\t",col.names=T,row.names=F,quote=F)


M1   <- trainM[,protIds]
M1   <- apply(M1,2,function(v){(v-mean(v))/sd(v)})
DM1  <- as.matrix(dist(M1))

L=399
mx=max(DM1)+1
for(i in 1:L){ DM1[i,i] <- mx}


lbls <- strsplit('NMACP',"")[[1]]
K=7
kNN1 <- apply(DM1,2,function(v){  
  tmp <- label[order(v)[1:K]]
  lblCounts <- sapply(lbls,function(ch){ sum(tmp==ch)})
  rtv       <- lblCounts;
  a         <- sum(lblCounts==max(lblCounts));
  ifelse(a==1,rtv <- c(rtv,names(lblCounts)[lblCounts==max(lblCounts)]),rtv <- c(rtv,tmp[1]))
  rtv
})
kNN1 <- data.frame(t(kNN1))
kNN1$observed <- label
kNN1$patientId <- SWT$patientId
colnames(kNN1) <- c(lbls,c('predicted','observed','patientId'))
kNN1 <- kNN1[,c(8,1:5,7,6)]
write.table(kNN1,file="E:/projects/TPD/results/kNN_imputed_50proteins",sep="\t",col.names=T,row.names=F,quote=F)


#########################################
NM0=list();
MA0=list();
AC0=list();
CP0=list();
M1=trainM
M1[is.na(M1)] <- 0;
M1$label=label
for(i in 1:loopSize){  NM0[[length(NM0)+1]] <-  growRF01('NM',M1)}
for(i in 1:loopSize){  MA0[[length(MA0)+1]]  <- growRF01('MA',M1)}
for(i in 1:loopSize){  AC0[[length(AC0)+1]]  <- growRF01('AC',M1)}
for(i in 1:loopSize){  CP0[[length(CP0)+1]]  <- growRF01('CP',M1)}

PW <- matrix(0,nrow=dim(trainM)[2],ncol=4)
rownames(PW) <- colnames(trainM);
colnames(PW) <- c('NM','MA','AC','CP')
models0 <- list()
models0$NM <- NM0
models0$MA <- MA0
models0$AC <- AC0
models0$CP <- CP0
for(lb in c('NM','MA','AC','CP')){
   for(md in models0[[lb]]){ PW[md$protIds,lb] <- PW[md$protIds,lb]+1}
}
PW <- data.frame(PW)
sm <- apply(PW,1,sum)
sm <- sm[sm>0]
sm <- sort(sm,decreasing=T)
PW <- PW[names(sm),]
write.table(PW,file="E:/projects/TPD/results/proteins_used_in_100times_RF.txt",sep="\t",col.names = T,row.names = T,quote=F)

###################################################################################################
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
##############################################################################################################
M2 <- trainM
M2$label <- label
for(i in 1:50){
  t1 <- kNN(M2,rownames(PW)[1:i])
  print(sprintf("%d,%f",i,t1$stat[3,6]))
}
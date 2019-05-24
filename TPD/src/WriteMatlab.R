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
####################################################imputation

rate   <- apply(trainM,2,function(v){sum(is.na(v))})
M1 <- trainM[,which(rate!=0)]
M0 <- trainM[,which(rate ==0)]

c1 <- c()

t1 <- for(cname in colnames(M0)){
      tmp <- cbind(M0[,cname],M1)
      imp1=amelia(tmp, m=1, parallel = "no")
      c1 <- cbind(c1,imp1$imputations$imp1[,1])
 }

writeMat("E:/projects/TPD/data/trainM_perc5Missing.mat",train=as.matrix(trainM))
writeMat("E:/projects/TPD/data/testM_perc5Missing.mat",test=as.matrix(testM))

########################### imputation using mice
M0 =trainM
imp01 <- mice(data = M0, m = 5, method = "pmm", maxit = 50, seed = 500)
comp01 <- complete(imp01,1)

M1 =testM
imp02 <- mice(data = M1, m = 5, method = "pmm", maxit = 50, seed = 500)
comp02 <- complete(imp02,1)

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


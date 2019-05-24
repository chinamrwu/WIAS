remove(list=ls())
library(sqldf)

sampleT <- read.table("E:/projects/TPD/data/RF_TPDT_Sample_Label.txt",sep="\t",header=T,stringsAsFactors=F)
sampleV <- read.table("E:/projects/TPD/data/RF_TPDV_patient_sample.txt",sep="\t",header=T,stringsAsFactors=F)
T0 <- read.table("E:/projects/TPD/data/TPD_Win600_protMatrix_20181115.csv",sep=",",header=T,stringsAsFactors=F)
tmp <- T0[,c(2:4)]
tmp <- sqldf("SELECT * FROM tmp where prot like '1/sp%'")
indx <- match(tmp$prot,T0$prot)
T0 <- T0[indx,-1]
colnames(T0) <- as.character(sapply(colnames(T0),function(v){strsplit(v,"_with_")[[1]][1]}))
T1 <- t(T0[,-1])
colnames(T1) <- as.character(sapply(T0$prot,function(v){strsplit(v,"\\|")[[1]][2]}))
protIds  <- colnames(T1)
sampleId <- rownames(T1)
T1 <- data.frame(T1)
T1$sampleId <-sampleId
T1 <- T1[,c('sampleId',protIds)]
rownames(T1) <- 1:dim(T1)[1]
T1$patientId <- sampleT[match(T1$sampleId,sampleT$sampleId),1]
T1 <- T1[,c('patientId',protIds)]
T1 <- T1[!grepl('e',T1$patientId),]
patientId <- unique(sapply(T1$patientId,function(v){substr(v,1,(nchar(v)-1))}))
rownames(T1) <- T1$patientId
T1 <- T1[order(T1$patientId),]

T2 <- c()
for(pId in patientId){
     sp <- as.character(sapply(c('a','b','c'),function(v){paste0(pId,v)}))
     tmp  <- T1[sp,-1]
     tmp1 <- apply(tmp,2,function(v){v[is.na(v)] <- mean(v[!is.na(v)]);v})
     T2 <- rbind(T2,tmp1)
 }
T2        <- T2[rownames(T1),]
T2        <- t(apply(T2,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
T2        <- data.frame(apply(T2,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))

T2$label  <- sapply(rownames(T2),function(v){substr(v,1,1)})
########## The following codes  imputes missing values by kNN
cutOff <- 60
R1 <- apply(T2,2,function(v){sum(is.na(v))})
R1 <- 100*R1/dim(T2)[1]
T3 <- T2[,R1<=cutOff]
T3[is.na(T3)] <- NA

K=5
t2 <- sapply(rownames(T3),function(rn){
   v1 <- T3[rn,]
   lbl  <- T3[rn,'label']
   feature0 <- names(v1)[is.na(v1)]
   feature1 <- names(v1)[!is.na(v1)]
   feature1 <- feature1[feature1!='label']

   tmp <- T3[T3$label==lbl,colnames(T3)!='label']
   tmp <-  tmp[rownames(tmp)!=rn,]
   tmp1 <- tmp[,feature1]
   tmp1 <- rbind(v1[feature1],tmp1)
   tmp1[is.na(tmp1)] <- min(tmp1)-1000
   D1 <- as.matrix(dist(tmp1))
   D1[1,1] <- max(D1)+1

   d1 <- tmp[rownames(D1)[order(D1[1,])[1:K]],feature0]
   d1[is.na(d1)] <- NA
   v1[1,feature0] <- as.numeric(apply(d1,2,function(v){mean(v,na.rm=T)}))
   as.numeric(v1[1,])
})
t2 <- data.frame(t(t2))

colnames(t2) <- colnames(T3)
R3 <- apply(t2,2,function(v){sum(is.na(v))*100/length(v)})
t2 <- t2[,R3<=5]

patientId <- unique(as.character(sapply(rownames(t2),function(v){substr(v,1,(nchar(v)-1))})))

T2 <- c()
for(pId in patientId){
     sp <- as.character(sapply(c('a','b','c'),function(v){paste0(pId,v)}))
     tmp  <- t2[sp,]
     tmp1 <- apply(tmp,2,function(v){v[is.na(v)] <- mean(v[!is.na(v)]);v})
     T2 <- rbind(T2,tmp1)
 }
T2  <- T2[rownames(t2),]
T2  <- t(apply(T2,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
T2  <- data.frame(apply(T2,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
T2  <- T2[,colnames(t2)]
T2$label <- as.character(sapply(rownames(T2),function(v){substr(v,1,1)}))

 #  source("E:/projects/TPD/src/aPatient3Sample.R")   
 

indx <-      c(379,184,278,38,662,265,255,314,516,237,98,74,406,495,111,103,183,537,588,238)
indx <- c(indx, 95,316,407,513,522,259,109,102,527,656,88,136,557,675,219,451,476,181,23,128)
indxAC <- c(646,274,211,652,350,406,638,128,244,440,44,424,97,674,658,230,217,226,170,242,71,43,256,181,557,141,676,495,413,401)
pcaAC <- list()
for(i in 2:length(indxAC)){pcaAC[[length(pcaAC)+1]] <- drawPCA(T2[T2$label %in% c('A','C'),c(680,indxAC[1:i])],sprintf("PCA for AC with %d proteins",i))}
tsneAC <- list()
for(i in 2:length(indxAC)){tsneAC[[length(tsneAC)+1]] <- plotTSNE(T2[T2$label %in% c('A','C'),c(680,indxAC[1:i])],sprintf("tSNE for AC with %d proteins",i))}


protAC <- colnames(T2)[indxAC]
t1 <- combn(protAC,6)
colnames(t1) <- 1:dim(t1)[2]
acc    <- -seq(9,0,by=-1)
result <- rep(0,10)
AC <- T2[T2$label %in% c('A','C'),]
for(i in 1:dim(t1)[2]){
  kk <- kNN(AC,t1[,i])
  m1 <- min(acc)
  if(m1 < kk$stat[3,3]){
     acc[acc==m1] <- kk$stat[3,3];
     result[acc==m1] <- i;
     print(sprintf("%d -- %4.3f %s",i,kk$stat[3,3],paste0(t1[,i],collapse=",")))
  }
}

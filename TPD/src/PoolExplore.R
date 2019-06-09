rm(list=ls())
library(ggplot2)
library(caret)
library(umap)
library(sqldf)
setwd("F:/WIAS/TPD")
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
machines <- as.character(sapply(colnames(pool)[-1],function(v){substr(v,1,1)}))
prots <- sapply(pool$prot,function(v){strsplit(v,"\\|")[[1]][2]})
pool <- data.frame(t(pool[,-1]))
colnames(pool) <- prots
rownames(pool) <- as.character(sapply(rownames(pool),function(v){a=strsplit(v,"sunyt_TPD_DIA_")[[1]];b=strsplit(a[2],"_pool")[[1]][1];paste0(a[1],"_",b)}))
pool$MS <- machines
pool$label <- machines
pool <- pool[,c('label',prots)]
color3=c(A="red",B="black",C="blue")
color2=c(M='red',B='blue')
############################################################################################################
R0 <- apply(pool,2,function(v){sum(is.na(v))/length(v)*100})
ump <- drawUMAP(pool[,R0 <=10],color3,strTitle=sprintf("UMAP:pool samples with %d missing rate",10),rowNormalization=T,colNormalization=T)
ump <- ump + scale_x_continuous(breaks = round(seq(min(ump$data$X), max(ump$data$X), by = 0.5),1))+ 
     scale_y_continuous(breaks = round(seq(min(ump$data$Y), max(ump$data$Y), by = 0.5),1)) 
dat <- ump$data
rowsA <- rownames(dat)[dat$X > 5]
rowsB <- rownames(dat)[dat$Y < -4.7 & dat$X > -4.7 & dat$X < -2.2]
rowsC <- rownames(dat)[dat$Y < -5 & dat$X < -5.2]

A0 <- getDataFrame(rowsA)
B0 <- getDataFrame(rowsB)
C0 <- getDataFrame(rowsC)
print("Samples distribution in A:");print(table(A0$label));
print("Samples distribution in B:");print(table(B0$label));
print("Samples distribution in C:");print(table(C0$label));
if(F){
	write.table(A0,file="data/A0.txt",sep="\t",col.names=T,row.names=T,quote=F)
	write.table(B0,file="data/B0.txt",sep="\t",col.names=T,row.names=T,quote=F)
	write.table(C0,file="data/C0.txt",sep="\t",col.names=T,row.names=T,quote=F)
}

A0$label[A0$label %in% c('N','M','A')] <- 'B'
A0$label[A0$label %in% c('C','P','W')] <- 'M'
B0$label[B0$label %in% c('N','M','A')] <- 'B'
B0$label[B0$label %in% c('C','P','W')] <- 'M'

C0$label[C0$label %in% c('N','M','A')] <- 'B'
C0$label[C0$label %in% c('C','P','W')] <- 'M'

print("Samples distribution in A:");print(table(A0$label));
print("Samples distribution in B:");print(table(B0$label));
print("Samples distribution in C:");print(table(C0$label));




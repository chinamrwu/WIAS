remove(list=ls())
library("sqldf")


lbl <- read.table("E:/projects/TPD/data/RF_TPDT_Sample_Label.txt",sep="\t",stringsAsFactors = F,header = T,check.names =F)
lbl <- sqldf("SELECT * FROm lbl where ID2 not like '%ouse%' and ID2 not like '%ool%'") # remove mouseliver and pool
sampleInfo <- t(apply(lbl,1,function(v){
	s1=v[1];
	k=nchar(s1);
	c(substr(s1,1,1),substr(s1,k,k),substr(s1,2,k-1),v[2])
	}))
colnames(sampleInfo) <- c("label","repc","patientId","sampleId")
sampleInfo <- data.frame(sampleInfo)
sampleInfo <- sqldf("SELECT * FROM sampleInfo WHERE sampleId !='' and label in ('A','C','M','N','P')")

mat <- read.csv("E:/projects/TPD/data/tpdt_dia_expert_prot_181024.csv",stringsAsFactors = F,header = T,check.names =F)
mat <-mat[grep("^1/sp*",mat$prot,perl=T),]

sampleIds <- as.character(sapply(colnames(mat)[3:dim(mat)[2]],function(s){a=strsplit(s,"_area")[[1]][1]}));

mat1 <- t(mat[,-c(1:2)])
colnames(mat1) <- as.character(sapply(mat$prot,function(s1){strsplit(s1,"\\|")[[1]][2]}))
#colnames(mat1) <- mat$prot
rownames(mat1) <- sampleIds
mat1 <- data.frame(mat1)
mat1$sampleId <- sampleIds

M <- merge(mat1,sampleInfo)
M[is.na(M)]=0;
M$patientId <- as.character(M$patientId)
w=dim(M)[2]

tmp <- M[,c((w-2):w)]
lbls <- as.character(unique(M$label))
trainM=c();
for( lbl in lbls){
   t1 <- M[M$label==lbl,]
   patients <- unique(t1$patientId)
   for(pid in patients){
      t2 <- t1[t1$patientId==pid,c(2:(w-3))]
      v1 <- apply(t2,2,mean)
      trainM <- rbind(trainM,c(paste0(lbl,pid),v1,lbl))
   }

}
colnames(trainM)[1] <- "patientId"
colnames(trainM)[dim(trainM)[2]] <- "label"
w <- dim(trainM)[2]
trainM <- trainM[,c(1,w,2:(w-1))]
write.table(trainM,file="E:/projects/TPD/data/RF_TPDT_expert_trainingM_20181024.txt",sep="\t",col.names=T,row.names=F,quote=F)

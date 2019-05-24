remove(list=ls())
library("sqldf")

sampleInfo <- read.table("E:/projects/TPD/data/RF_TPDV_patient_sample.txt",sep="\t",stringsAsFactors = F,header = T,check.names =F)
sampleInfo <- sqldf("SELECT * FROM sampleInfo WHERE patientId>0 and sampleId not like '%pool%' and sampleId not like '%ml%' order by patientId ")
#############

mat <- read.csv("E:/projects/TPD/data/tpdv_openswath_prot_181108.csv",stringsAsFactors = F,header = T,check.names =F)
mat <-mat[grep("^1/sp*",mat$prot,perl=T),]
sampleIds <- as.character(sapply(colnames(mat)[3:dim(mat)[2]],function(s){a=strsplit(s,"_with_dscore")[[1]][1]}));
colnames(mat)[3:dim(mat)[2]] <- sampleIds

indx <- match(sampleInfo$sampleId,colnames(mat),nomatch=-1)
sampleInfo <- sampleInfo[-which(indx<0),]
mat1 <- mat[,indx[which(indx>0)]]

rownames(mat1) <- mat$prot



patientIds <- unique(sampleInfo$patientId)
t1 <- sapply(patientIds,function(pid){
   sampleId <- sampleInfo$sampleId[sampleInfo$patientId==pid]
   tmp <- mat1[,sampleId]
   v <- apply(tmp,1,mean,na.rm=T)
   v[is.nan(v)] <- NA
   v
})

colnames(t1) <- as.character(sapply(patientIds,function(id){paste0("P",id)}))
rownames(t1) <- as.character(sapply(rownames(t1),function(s1){a=strsplit(s1,"1/sp\\|")[[1]][2];strsplit(a,"_")[[1]][1]}))
t1 <- data.frame(t(t1),stringsAsFactors=F,check.names=F)
write.table(t1,file="E:/projects/TPD/data/RF_TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",col.names=T,row.names=T,quote=F)

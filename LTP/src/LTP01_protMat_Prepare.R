rm(list=ls())
library(data.table)


setwd('E:/WIAS/LTP')
df0 <- read.csv("data/LTP_protMat_20190413.csv",header = T,stringsAsFactors = F)
tmp <- data.frame(t(df0[,-c(1,2)]))
colnames(tmp) <- df0$prot
tmp <- tmp[,grepl('HUMAN',colnames(tmp))]
colnames(tmp) <- as.character(sapply(colnames(tmp),function(v){a <- strsplit(v,"\\|")[[1]][2]}))

matPool <- tmp[grepl('pool',rownames(tmp)),]
rownames(matPool) <- as.character(sapply(rownames(matPool),function(v){
  a <- strsplit(v,"caix_wangly_LTP_DIA_")[[1]];
  b <- a[1]
  b1 <- strsplit(a[2],"_with_dscore_filtered")[[1]][1]
  s1 <- gsub("_pool","",paste0(b,"_",b1))
}))
protIds <- colnames(matPool)
lbls    <- substr(rownames(matPool),1,1)
matPool$label <- lbls
matPool <- matPool[,c('label',protIds)]

write.table(matPool,file="data/matPool.txt",sep="\t",col.names=T,row.names=T,quote=F)

###############################################################
matProt <- tmp[!grepl('_pool',rownames(tmp)),]

sampleInf <- read.table("data/sample_info_20190409.txt",sep="\t",header=T,stringsAsFactors=F)
patientIds <- unique(sampleInf$PatientId)
k0 <- c()
for(patientId in patientIds){
    v0 <- apply(matProt[sampleInf$rawFile[sampleInf$PatientId==patientId],],2,function(v){log2(mean(2^v,na.rm=T))})
	 k0 <- rbind(k0,v0)
}
rownames(k0) <- patientIds
k0[is.na(k0)] <- NA
R0 <- apply(k0,2,function(v){sum(is.na(v))/length(v)*100})
k0 <- data.frame(k0[,R0!=100])

subtypes <- unique(sampleInf[,c('PatientId','PType')])
rownames(subtypes) <- subtypes$PatientId

protIds  <- colnames(k0)
k0$label <- subtypes[rownames(k0),'PType']
k0 <- k0[,c('label',protIds)]

write.table(k0,file='data/LTP_protMat_avg_20190526.txt',sep="\t",col.names=T,row.names=T,quote=F)
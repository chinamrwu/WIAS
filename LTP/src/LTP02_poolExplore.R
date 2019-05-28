rm(list=ls())
library(data.table)
library(sqldf)


source("src/common.R")

df0 <- read.table('data/matPool.txt',sep="\t",header=T,row.names=1,stringsAsFactors=F)

R0 <- apply(df0,2,function(v){sum(is.na(v))})
color2 <- c('A'='red','B'='blue')
ump1 <- drawUMAP(df0[,R0<=20],color2,strTitle=sprintf("UMAP1:%d proteins of LTP",sum(R0<=20)-1),rowNormalization=T,colNormalization=F)
ump2 <- drawUMAP(df0[,R0<=15],color2,strTitle=sprintf("UMAP2:%d proteins of LTP",sum(R0<=15)-1),rowNormalization=T,colNormalization=F)
ump3 <- drawUMAP(df0[,R0<=10],color2,strTitle=sprintf("UMAP3:%d proteins of LTP",sum(R0<=10)-1),rowNormalization=T,colNormalization=F)
ump4 <- drawUMAP(df0[,R0<=5],color2,strTitle=sprintf("UMAP4:%d proteins of LTP", sum(R0<=5)-1),rowNormalization=T,colNormalization=F)
ump5 <- drawUMAP(df0[,R0<=3],color2,strTitle=sprintf("UMAP5:%d proteins of LTP", sum(R0<=3)-1),rowNormalization=T,colNormalization=F)
ump6 <- drawUMAP(df0[,R0<=2],color2,strTitle=sprintf("UMAP6:%d proteins of LTP", sum(R0<=2)-1),rowNormalization=T,colNormalization=F)
ump7 <- drawUMAP(df0[,R0<=1],color2,strTitle=sprintf("UMAP7:%d proteins of LTP", sum(R0<=1)-1),rowNormalization=T,colNormalization=F)


dat1 <- ump1$data
batchA <- as.character(sapply(rownames(dat1)[dat1$X > -1],function(v){strsplit(v,"_")[[1]][2]}))
batchB <- as.character(sapply(rownames(dat1)[dat1$X < -1],function(v){strsplit(v,"_")[[1]][2]}))

df1 <- read.csv("data/LTP_protMat_20190413.csv",header = T,stringsAsFactors = F)
tmp <- data.frame(t(df1[,-c(1,2)]))
colnames(tmp) <- df1$prot
tmp <- tmp[,grepl('HUMAN',colnames(tmp))]
colnames(tmp) <- as.character(sapply(colnames(tmp),function(v){a <- strsplit(v,"\\|")[[1]][2]}))

rn <- as.character(sapply(rownames(tmp),function(v){a <- strsplit(v,'_DIA_')[[1]][2];b <- strsplit(a,'_with_dscore')[[1]][1];b}))
rn <- as.character(sapply(rn,function(v){strsplit(v,"_")[[1]][1]}))

sampleInf <- read.table("data/sample_info_20190409.txt",sep="\t",header=T,stringsAsFactors=F)
patientIds <- unique(sampleInf$PatientId)

#######################################################
getDataSet <- function(batchIds){
		tA <- tmp[which(rn %in% batchIds),]
		tA <- tA[!grepl('pool',rownames(tA)),] # 198X4238
		sampInfA <- sampleInf[sampleInf$rawFile %in% rownames(tA),]
		patients <- unique(sampInfA$PatientId) # 153 patients belongs to batchA
		k0 <- c()
		for(patientId in patients){
			fnames <- sampInfA$rawFile[sampInfA$PatientId==patientId]
			v0  <- apply(tA[fnames,],2,function(v){log2(mean(2^v,na.rm=T))})
			k0  <- rbind(k0,v0)
		}

		rownames(k0) <- patients
		k0[is.na(k0)] <- NA
		R0 <- apply(k0,2,function(v){sum(is.na(v))/length(v)*100})
		k0 <- data.frame(k0[,R0!=100])

		subtypes <- unique(sampInfA[,c('PatientId','PType')])
		rownames(subtypes) <- subtypes$PatientId
		protIds  <- colnames(k0)
		k0$label <- subtypes[rownames(k0),'PType']
		k0 <- k0[,c('label',protIds)]
		k0
}
matA <- getDataSet(batchA)
matB <- getDataSet(batchB)

write.table(matA,file='data/matA.txt',sep="\t",col.names=T,row.names=T,quote=F)
write.table(matB,file='data/matB.txt',sep="\t",col.names=T,row.names=T,quote=F)

benign <- c('B','Q','P')
malign <- c('C','S','M','G')
normal <- c('N')

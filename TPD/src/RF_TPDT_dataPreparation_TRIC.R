remove(list=ls())
library("sqldf")
library("data.table")


sampleLabels <- read.table("E:/projects/TPD/data/RF_TPDT_Sample_Label.txt",sep="\t",stringsAsFactors = F,header = T,check.names =F)
sampleLabels <- sqldf("SELECT * FROm sampleLabels where ID2 not like '%ouse%' and ID2 not like '%ool%'") # remove mouseliver and pool

sampleInfo <- t(apply(sampleLabels,1,function(v){
	s1=v[1];
	k=nchar(s1);
	c(substr(s1,1,1),substr(s1,k,k),substr(s1,2,k-1),v[2])
	}))
colnames(sampleInfo) <- c("label","repc","patientId","sampleId")
sampleInfo <- data.frame(sampleInfo)
sampleInfo <- sqldf("SELECT * FROM sampleInfo WHERE sampleId !='' and label in ('A','C','M','N','P')")

#mat <- read.csv("E:/projects/TPD/data/tpdt_prot_matrix_181106_byliucan.csv",stringsAsFactors = F,header = T,check.names =F)
mat <- read.csv("E:/projects/TPD/data/tpdt_tric_protmatrix_181104.csv",stringsAsFactors = F,header = T,check.names =F)
sampleIds <- as.character(sapply(colnames(mat)[3:dim(mat)[2]],function(s){a=strsplit(s,"_with_")[[1]][1]}));

mat1 <- t(mat[,-c(1:2)])
colnames(mat1) <- mat$prot
mat1 <- mat1[,colnames(mat1) %like% '1/sp']
colnames(mat1) <- as.character(sapply(colnames(mat1),function(v){strsplit(v,"1/sp\\|")[[1]][2]}))
colnames(mat1) <- as.character(sapply(colnames(mat1),function(v){strsplit(v,"_")[[1]][1]}))
rownames(mat1) <- sampleIds
mat1 <- data.frame(mat1,check.names=F,stringsAsFactors=F)
mat1$sampleId <- sampleIds

M <- merge(mat1,sampleInfo)

protIds <- colnames(M)[2:(dim(M)[2]-3)]
tmp <- M[,c("label","patientId","sampleId",protIds[1])]
sql <- sprintf("SELECT label || patientId patientId,avg([%s]) [%s] FROM tmp GROUP BY label,patientId",protIds[1],protIds[1])
t0 <- sqldf(sql);
t1 <- t0[,2]

indx <- 2
while(indx+49 < length(protIds)){

   cnames <- as.character(sapply(protIds[indx:(indx+49)],function(v){sprintf("avg([%s]) [%s]",v,v)}));
   tmp <- M[,c("label","patientId","sampleId",protIds[indx:(indx+49)])]
   t1  <- cbind(t1,sqldf(sprintf("SELECT %s FROM tmp GROUP BY label,patientId",paste0(cnames,collapse=","))))   
   indx <- indx+50
}

cnames <- as.character(sapply(protIds[indx:length(protIds)],function(v){sprintf("avg([%s]) [%s]",v,v)}));
tmp <- M[,c("label","patientId","sampleId",protIds[indx:length(protIds)])]
t1  <- cbind(t1,sqldf(sprintf("SELECT %s FROM tmp GROUP BY label,patientId",paste0(cnames,collapse=","))))   
colnames(t1)[1] <- protIds[1]
t1$patientId <- t0$patientId
t1 <- t1[,c("patientId",protIds)]

write.table(t1,file="E:/projects/TPD/data/TPDT_TRIC_20181115.txt",sep="\t",col.names=T,row.names=F,quote=F)


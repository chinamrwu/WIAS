remove(list=ls())
library("sqldf")
library("FactoMineR")
library("factoextra")
library("grid")
library("data.table")
library("ggplot2")

sampleInfo <- read.table("E:/projects/TPD/data/RF_TPDV_patient_sample.txt",sep="\t",stringsAsFactors = F,header = T,check.names =F)
sampleInfo <- sqldf("SELECT * FROM sampleInfo WHERE  sampleId  like '%pool%' or sampleId  like '%ml%' order by patientId ")
#############

mat <- read.csv("E:/projects/TPD/data/tpdv_openswath_prot_181108.csv",stringsAsFactors = F,header = T,check.names =F)
mat <-mat[grep("^1/sp*",mat$prot,perl=T),]
sampleIds <- as.character(sapply(colnames(mat)[3:dim(mat)[2]],function(s){a=strsplit(s,"_with_dscore")[[1]][1]}));
colnames(mat)[3:dim(mat)[2]] <- sampleIds

indx <- match(sampleInfo$sampleId,colnames(mat),nomatch=-1)
sampleInfo <- sampleInfo[-which(indx<0),]
mat1 <- mat[,indx[which(indx>0)]]
mat1 <- data.frame(t(mat1))

sampleNames <- as.character(sapply(rownames(mat1),function(v){strsplit(v,"_DIA_")[[1]][2]}))
protIds <- sapply(mat$prot,function(v){strsplit(v,"\\|")[[1]][2]})
colnames(mat1) <- protIds
rownames(mat1) <- 1:dim(mat1)[1]
mat1$sampleId <- sampleNames
mat1 <- mat1[,c("sampleId",protIds)]

poolV <- mat1[grepl('pool',mat1$sampleId),]
liverV   <- mat1[grepl('ml',mat1$sampleId),] 
##############################################################################################################################################
sampleLabels <- read.table("E:/projects/TPD/data/RF_TPDT_Sample_Label.txt",sep="\t",stringsAsFactors = F,header = T,check.names =F)
sampleLabels <- sqldf("SELECT * FROm sampleLabels where sampleId  like '%ml%' or sampleId like '%ool%'") # remove mouseliver and pool

mat <- read.csv("E:/projects/TPD/data/tpdt_prot_matrix_181106_byliucan.csv",stringsAsFactors = F,header = T,check.names =F)
sampleIds <- as.character(sapply(colnames(mat)[3:dim(mat)[2]],function(s){a=strsplit(s,"_with_")[[1]][1]}));
protIds <- mat$prot
colnames(mat)[3:dim(mat)[2]] <- sampleIds
#sampleIds <- as.character(sapply(sampleIds,function(v){a=strsplit(v,"sunyt_TPD_DIA_")[[1]];paste0(a[1],"_",a[2])}))
sampleIds <- intersect(sampleIds,sampleLabels$sampleId)
mat <- data.frame(t(mat[,sampleIds]))
colnames(mat) <- protIds
mat <- mat[,grepl("1/sp",colnames(mat))]
colnames(mat) <- as.character(sapply(colnames(mat),function(v){strsplit(v,"\\|")[[1]][2]}))

poolT  <- mat[grepl('pool',rownames(mat)),]
liverT <- mat[grepl('ml',rownames(mat)),]
################################################################################################################################################


drawPCA <- function(df0,strTitle){ ## M is a matrix or dataframe, rows are samples and columns are features, rownames are sample names
   M <- df0[,colnames(df0)!='label'] 
   m1 <- prcomp(M,F);
   Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
   Y  <- Y[,c(1,2)]
   
   Y <- data.frame(Y,df0$label);
   colnames(Y) <- c("PC1","PC2","label")

   eigs <- m1$sdev^2
   percentages <- eigs[1:2] / sum(eigs)
   

   p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=2.5)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank())
   p <- p +  labs(x = sprintf("PC1(%4.2f%%)",percentages[1]*100),
                  y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
	      title =strTitle)
  
   p
}
##############################
protIds <- intersect(colnames(poolT),colnames(poolV))
poolT   <- poolT[,protIds]
poolV   <- poolV[,protIds]


protIds <- intersect(colnames(liverT),colnames(liverV))
liverT  <- liverT[,protIds]
liverV  <- liverV[,protIds]

rate <- 5
R1 <- apply(poolT,2,function(v){sum(is.na(v))/length(v)*100})
R2 <- apply(poolV,2,function(v){sum(is.na(v))/length(v)*100})
protIds <- intersect(names(R1)[R1<rate],names(R2)[R2<rate])

M0 <- rbind(poolT[,protIds],poolV[,protIds]);
M0[is.na(M0)] <-0
M0$label <- c(rep("Training",dim(poolT)[1]),rep("Testing",dim(poolV)[1]))
drawPCA(M0,"PCA without normalization");

M0 <- rbind(poolT[,protIds],poolV[,protIds]);
M0 <- data.frame(t(apply(M0,1,function(v){v/mean(v,na.rm=T)}))) # by row
M0$label <- c(rep("Training",dim(poolT)[1]),rep("Testing",dim(poolV)[1]))
M0[is.na(M0)] <- 0; 
drawPCA(M0,"PCA for pool samples normalized by row")

M0 <- rbind(poolT[,protIds],poolV[,protIds]);
M0 <- data.frame(apply(M0,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))##by column
M0$label <- c(rep("Training",dim(poolT)[1]),rep("Testing",dim(poolV)[1]))
M0[is.na(M0)] <- 0; 
drawPCA(M0,"PCA for pool samples normalized by column")

M0 <- rbind(poolT[,protIds],poolV[,protIds]);
M0 <- data.frame(t(apply(M0,1,function(v){v/mean(v,na.rm=T)}))) # by row
M0 <- data.frame(t(apply(M0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))##by column
M0$label <- c(rep("Training",dim(poolT)[1]),rep("Testing",dim(poolV)[1]))
M0[is.na(M0)] <- 0; 
drawPCA(M0,"PCA for pool samples normalized by row and column")


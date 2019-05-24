remove(list=ls())
library(data.table)
library(sqldf)
library(ggplot2)
library(grid)


df0  <- read.table('E:/projects/TPD/data/TPD_prot_matrix_20190130.csv',sep=",",header=T,stringsAsFactors=F)
colnames(df0) <- as.character(sapply(colnames(df0),function(v){strsplit(v,"_with_")[[1]][1]}))
df1 <- read.table("E:/projects/TPD/data/AllSampleInformation_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
indx <- match(df1$sampleId,colnames(df0),nomatch=-1)
indx <- indx[indx>0]
df0 <- df0[,c('peptide_group_label','prot',colnames(df0)[indx])]
t0 <- data.frame(t(df0[,-c(1:2)]))
colnames(t0) <- as.character(sapply(df0$prot,function(v){strsplit(v,"\\|")[[1]][2]}))
t0 <- t0[,!grepl('iRT',colnames(t0))]
protIds <- colnames(t0)

indx <- match(rownames(t0),df1$sampleId)
t0$SpecimenID <- df1$SpecimenID[indx]
t0$label <- df1$label[indx]
t0$MS <- as.character(sapply(rownames(t0),function(v){substr(v,1,1)}))
t0$Date <- sapply(rownames(t0),function(v){
    a <- strsplit(v,'suny')[[1]][1];
    b <- substr(a,2,nchar(a));
    b <- strsplit(b,'')[[1]];
    b <-c(paste0(b[1:4],collapse = ""),paste0(b[5:6],collapse = ""),paste0(b[7:8],collapse = ""));
    b <- paste0(b,collapse="-")
    b
   })
t0$Year <- sapply(t0$Date,function(v){as.integer(strsplit(v,"-")[[1]][1])})
t0$Month <- sapply(t0$Date,function(v){as.integer(strsplit(v,"-")[[1]][2])})
t0$Day <- sapply(t0$Date,function(v){as.integer(strsplit(v,"-")[[1]][3])})
t0 <- t0[,c('SpecimenID','label','MS','Year','Month','Day',protIds)]
write.table(t0,file="E:/projects/TPD/data/TPD_prot_matrix_rep_20190304.txt",sep="\t",col.names=T,row.names=T,quote=F)

Specimens <- unique(t0$SpecimenID)
t1 <- sapply(Specimens,function(v){
    tmp <- t0[t0$SpecimenID==v,-c(1:6)]
    tmp1 <- apply(tmp,2,function(v1){v0 <- 2^v1;log2(mean(v0,na.rm=T))})
    tmp1
   })
t1 <- data.frame(t(t1))
R1 <- apply(t1,2,function(v){sum(is.na(v))/length(v)*100})
t1 <- t1[,R1<100]
protIds <- colnames(t1)



tmp <- sqldf("SELECT distinct SpecimenID,label FROM df1")

t1$label <- tmp$label[match(rownames(t1),tmp$SpecimenID)]
t1$SpecimenID <- rownames(t1)
L=dim(t1)[2]
t1 <- t1[,!grepl('iRT',colnames(t1))] ## È¥µôiRTµ°°×
t1 <- t1[,c('label','SpecimenID',protIds)]

write.table(t1,file="E:/projects/TPD/data/TPD_prot_matrix_avg_20190304.txt",sep="\t",col.names=T,row.names=F,quote=F)

#########################################################################################################
tmp <- sqldf("SELECT distinct SpecimenID,groups FROM df1")
t2 <- t1[,-c(1:2)]
R2 <- apply(t2,2,function(v){100*sum(is.na(v))/length(v)})
t2 <- t2[,R2<=5]
protIds <- colnames(t2)
t2[is.na(t2)] <- 0
t2$label <- tmp$groups[match(t1$SpecimenID,tmp$SpecimenID)]
t2 <- t2[,c('label',protIds)]
rownames(t2) <- 1:dim(t2)[1]

p1 <- drawPCA(t2,sprintf('PCA: %d proteins without normalization ',dim(t2)[2]-1),rowNormalization=F,colNormalization=F)
p2 <- drawPCA(t2,sprintf('PCA: %d proteins row normalization only',dim(t2)[2]-1),rowNormalization=T,colNormalization=F)
p3 <- drawPCA(t2,sprintf('PCA: %d proteins column normalization only',dim(t2)[2]-1),rowNormalization=F,colNormalization=T)
p4 <- drawPCA(t2,sprintf('PCA: %d proteins row and column normalization ',dim(t2)[2]-1),rowNormalization=T,colNormalization=T)

pdf("E:/projects/TPD/results/PCA_SWATH_278_protein_20190115.pdf",width=10,height=5)
 grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 2)));
    print(p2,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p3,vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    #print(p3,vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
    #print(p4,vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()
#####################################################################################

tmp <- sqldf("SELECT distinct SpecimenID,label FROM df1")
t2 <- t1[,-c(1:2)]
R2 <- apply(t2,2,function(v){100*sum(is.na(v))/length(v)})
t2 <- t2[,R2<=5]
protIds <- colnames(t2)
t2[is.na(t2)] <- 0
t2$label <- tmp$label[match(t1$SpecimenID,tmp$SpecimenID)]
t2 <- t2[,c('label',protIds)]
rownames(t2) <- 1:dim(t2)[1]

p5 <- drawPCA(t2,sprintf('PCA: %d proteins row normalization only',dim(t2)[2]-1),rowNormalization=T,colNormalization=F)

pdf("E:/projects/TPD/results/PCA_SWATH_278_protein_20190115_6class.pdf",width=10,height=5)
 grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 2)));
    print(p2,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p3,vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    #print(p3,vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
    #print(p4,vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

###################################################
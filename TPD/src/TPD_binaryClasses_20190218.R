rm(list=ls())
library(randomForest)
library(data.table)
library(grid)
setwd("E:/projects/TPD")
source('src/common.R')


df0 <- read.table('data/TPD_prot_matrix_avg_20190131.txt',sep='\t',header=T,stringsAsFactors=F)
df0$label[df0$label=='W'] <- 'C'
df0$label[df0$label %in% c('N','M','A')] <- 'B'
df0$label[df0$label %in% c('C','P')]     <- 'M'

rownames(df0) <- df0$SpecimenID
df0 <- df0[,-1]
df1 <- getSubset(c('B','M'),df0,30,F) ## 30% cutOff for missing value rate

##########################  select proteins by volcano plot #############
lblColors <- c(B='black',M='red')
p0 <-  drawPCA(df1,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA: with all %d proteins",(dim(df1)[2]-1)))
t0 <- drawTSNE(df1,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("tSNE: with all %d proteins",(dim(df1)[2]-1)))

theFC=2.0
thePV=0.05
vc <- Volcano(df1,sprintf('TPD_binaryClasses_volcano_FC%3.1f_Pval%4.3f_20190219.pdf',theFC,thePV),outFile='tpd_binary_volcano_20190219.txt',thresholdFC=theFC,thresholdPValue=thePV,strTitle="volcano plot")
protIds <- as.character(vc$prot)
p1 <- drawPCA(df1[c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA:%d proteins from volcano with FC=%3.1f pval=%4.3f",dim(vc)[1],theFC,thePV))
t1 <- drawTSNE(df1[c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA:%d proteins from volcano with FC=%3.1f pval=%4.3f",dim(vc)[1],theFC,thePV))

theFC=2.5
thePV=0.01
vc <- Volcano(df1,sprintf('TPD_binaryClasses_volcano_FC%3.1f_Pval%4.3f_20190219.pdf',theFC,thePV),outFile='tpd_binary_volcano_20190219.txt',thresholdFC=theFC,thresholdPValue=thePV,strTitle="volcano plot")
protIds <- as.character(vc$prot)
p2 <- drawPCA(df1[c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA:%d proteins from volcano with FC=%3.1f pval=%4.3f",dim(vc)[1],theFC,thePV))
t2 <- drawTSNE(df1[c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA:%d proteins from volcano with FC=%3.1f pval=%4.3f",dim(vc)[1],theFC,thePV))

theFC=3.0
thePV=0.005
vc <- Volcano(df1,sprintf('TPD_binaryClasses_volcano_FC%3.1f_Pval%4.3f_20190219.pdf',theFC,thePV),outFile='tpd_binary_volcano_20190219.txt',thresholdFC=theFC,thresholdPValue=thePV,strTitle="volcano plot")
protIds <- as.character(vc$prot)
p3 <- drawPCA(df1[c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA:%d proteins from volcano with FC=%3.1f pval=%4.3f",dim(vc)[1],theFC,thePV))
t3 <- drawTSNE(df1[c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA:%d proteins from volcano with FC=%3.1f pval=%4.3f",dim(vc)[1],theFC,thePV))

pdf('binary/PCA_tSNE_by_Volcano.pdf',height=24,width=12)
grid.newpage();
pushViewport(viewport(layout = grid.layout(4, 2)));
  print(p0,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(t0,vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p1,vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(t1,vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(p2,vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(t2,vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
  print(p3,vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
  print(t3,vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
dev.off()

########################### feature selection based on FC=2.5 and p-value=0.05
theFC=2.5
thePV=0.01
vc <- Volcano(df1,sprintf('TPD_binaryClasses_volcano_FC%3.1f_Pval%4.3f_20190219.pdf',theFC,thePV),outFile='tpd_binary_volcano_20190219.txt',thresholdFC=theFC,thresholdPValue=thePV,strTitle="volcano plot")
protIds <- as.character(vc$prot)

bp2 <- combn(protIds,2)

library(caret)
L = dim(bp2)[2]

df1[is.na(df1)] <- 0;
result <- c()

for(i in 69:L){result <- rbind(result,RFScore1(df1[c('label',bp2[,i])],itors=3))}

L=dim(bp3)[2]
result3 <- c()
for(i in 1:L){result3 <- rbind(result3,RFScore1(df1[c('label',bp3[,i])],itors=3))}

L=dim(bp4)[2]
result4 <- c()
for(i in 1:L){result4 <- rbind(result4,RFScore1(df1[c('label',bp4[,i])],itors=3))}

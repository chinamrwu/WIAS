rm(list=ls())
library(ggplot2)
library(data.table)
library(sqldf)
library(glmnet)
setwd("E:/projects/TPD")
source("src/common.R")
df0 <- read.table('data/TPD_prot_matrix_avg_20190131.txt',sep='\t',header=T,stringsAsFactors=F)
df0$label[df0$label=='W'] <- 'C'
rownames(df0) <- df0$SpecimenID
df0 <- df0[,-1]
#################################################################################

R0 <- apply(df0,2,function(v){sum(is.na(v))/length(v)*100})
missingCutOff=80


df1 <- df0[,R0<=missingCutOff] 
df1$label[df1$label %in% c('N','M','A')] <- 'B'
df1$label[df1$label %in% c('C','P')]     <- 'M'
label <- as.factor(sapply(df1$label,function(ch){ifelse(ch=='B',0,1)}))
lblColors <- c(B='#56B4E9',M='#CC79A7')


################### 20190303


tmp0   <- df1[,colnames(df1)!='label']
tmp0   <- t(apply(tmp0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0
selection0 <- colnames(tmp0)[colnames(tmp0) != 'label']
p0 <- drawPCA(df1,lblColors,rowNormalization=T,colNormalization=T)


flg=T
while(flg){
  #for(a in seq(0.8,1.0,by=0.01)) {
     set.seed(1)
     cvfit0 <- cv.glmnet(tmp0,label,family='binomial',alpha=0.5)
     cf1se <- coef(cvfit0, s = "lambda.1se")
     selection1 <- rownames(cf1se)[cf1se[,1]!=0][-1]
     #print(sprintf('%d proteins selected with alpha=%4.4f',length(selection1),a))
  #}

   print(sprintf('%d original prots ,%d selected,%d intersected',length(selection0),length(selection1),length(intersect(selection1,selection0))))
   flg <- length(selection1) != length(selection0)
   selection0 <- selection1

   tmp0   <- df1[,selection0]
   tmp0   <- t(apply(tmp0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
   tmp0[is.na(tmp0)] <- 0
}

t1 <- sort(abs(cf1se[-1,]),decreasing=T) #remove interception
selection0 <- names(t1)

pts <- list()
for(i in 5:length(selection0)){
 pts[[length(pts)+1]] <- drawPCA(df1[, c('label',selection0[1:i])],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:%d',(i-1)))
 pts[[length(pts)+1]] <- drawTSNE(df1[,c('label',selection0[1:i])],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('tSNE:%d',(i-1)))
}

library(grid)
pdf("PCA_TSNE_LASSO_BM_20190303_3.pdf",width=12,height=6)
  for(i in seq(1,length(pts)-1,by=2)){
    grid.newpage();    pushViewport(viewport(layout = grid.layout(1, 2)));
    print(pts[[i]],vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(pts[[i+1]],vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
 }
dev.off()


##################### try multiple-response mgaussian with PCA : principal component regression
p1 <- drawPCA(df1[,c('label',selection0)],lblColors,rowNormalization=T,colNormalization=T)
t1 <- drawTSNE(df1[,c('label',selection0)],lblColors,rowNormalization=T,colNormalization=F)

#tmp0 <- df1[,selection0]

tmp0   <- df1[,colnames(df1)!='label']
#tmp0   <- df1[,selection0]
tmp0   <- t(apply(tmp0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp0[is.na(tmp0)] <- 0

mfit = cv.glmnet(x=as.matrix(tmp0), y=as.matrix(p1$data[,1:2]), family = "mgaussian",alpha=0.8)
cf1se <- coef(mfit, s = 'lambda.1se')
cf <- cf1se$PC1
sPC1 <- rownames(cf)[cf[,1]!=0][-1]
cf <- cf1se$PC2
sPC2 <- rownames(cf)[cf[,1]!=0][-1]
################################################################
label <- as.factor(sapply(df1$label,function(ch){ifelse(ch=='B',0,1)}))

missingRates <- seq(85,5,by=-5)

cvfits <- list()
for(missR in missingRates){
     tmp1 <- df0[,R0<=missR] 
     tmp1 <- tmp1[,colnames(tmp1)!='label']
     tmp1   <- t(apply(tmp1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
     tmp1[is.na(tmp1)] <- 0

     selection0 <- colnames(tmp1)

     flg=T
     while(flg){
	      set.seed(1)
	      cvfit0 <- cv.glmnet(tmp1,label,family='binomial',alpha=0.5)
	      cf1se <- coef(cvfit0, s = "lambda.1se")
	      selection1 <- rownames(cf1se)[cf1se[,1]!=0][-1]
	      print(sprintf('%d original prots ,%d selected,%d intersected',length(selection0),length(selection1),length(intersect(selection1,selection0))))
	      flg <- length(selection1) != length(selection0)
	      selection0 <- selection1

	      tmp1   <- df1[,selection0]
	      tmp1   <- t(apply(tmp1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	      tmp1[is.na(tmp1)] <- 0
     }
     cvfits[[length(cvfits)+1]] <- cvfit0

}
names(cvfits) <- missingRates
plots <- list()
for(nm in names(cvfits)){
  mfit <- cvfits[[nm]]
  cf1se <- coef(mfit, s = 'lambda.1se')
  protIds <- rownames(cf1se)[cf1se[,1]!=0][-1]
  
  plots[[length(plots)+1]] <- drawPCA(df1[, c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:%s%% missingRate cutoff,%d protein',nm,length(protIds)))
  plots[[length(plots)+1]] <- drawTSNE(df1[, c('label',protIds)],lblColors,rowNormalization=T,colNormalization=F,strTitle=sprintf('tSNE:%s%% missingRate cutoff,%d protein',nm,length(protIds)))

}

library(grid)
pdf("PCA_TSNE_LASSO_BM_20190304_1.pdf",width=12,height=6)
  for(i in seq(1,length(plots)-1,by=2)){
    grid.newpage();    pushViewport(viewport(layout = grid.layout(1, 2)));
    print(plots[[i]],vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(plots[[i+1]],vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
 }
dev.off()
########################################## in case of 5 classes
label5 <- as.factor(df0$label)
missingRates <- seq(85,5,by=-5)

cvfits5 <- list()
for(missR in missingRates){
     tmp1 <- df0[,R0<=missR] 
     tmp1 <- tmp1[,colnames(tmp1)!='label']
     tmp1   <- t(apply(tmp1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
     tmp1[is.na(tmp1)] <- 0

     selection0 <- colnames(tmp1)

     flg=T
     while(flg){
	      set.seed(1)
	      cvfit0 <- cv.glmnet(tmp1,label5,family='multinomial',type.multinomial = 'grouped',alpha=0.5,parallel = T)
	      cf1se <- coef(cvfit0, s = "lambda.1se")
	      selection1 <- rownames(cf1se[[1]])[cf1se[[1]][,1]!=0][-1]
	      print(sprintf('%d original prots ,%d selected,%d intersected',length(selection0),length(selection1),length(intersect(selection1,selection0))))
	      flg <- length(selection1) != length(selection0)
	      selection0 <- selection1

	      tmp1   <- df1[,selection0]
	      tmp1   <- t(apply(tmp1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	      tmp1[is.na(tmp1)] <- 0
     }
     cvfits5[[length(cvfits5)+1]] <- cvfit0

}

lblColors <- c(N='black',M='#e17832',A='#f5b819',C='#5992c6',P='#537e35')
names(cvfits5) <- missingRates
plots5 <- list()
for(nm in names(cvfits5)){
  mfit <- cvfits[[nm]]
  cf1se <- coef(mfit, s = 'lambda.1se')
  protIds <- rownames(cf1se)[cf1se[,1]!=0][-1]
  
  plots5[[length(plots5)+1]] <- drawPCA(df0[, c('label',protIds)],lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:%s%% missingRate cutoff,%d protein',nm,length(protIds)))
  plots5[[length(plots5)+1]] <- drawTSNE(df0[, c('label',protIds)],lblColors,rowNormalization=T,colNormalization=F,strTitle=sprintf('tSNE:%s%% missingRate cutoff,%d protein',nm,length(protIds)))

}

library(grid)
pdf("PCA_TSNE_LASSO_BM_20190304_3.pdf",width=12,height=6)
  for(i in seq(1,length(plots5)-1,by=2)){
    grid.newpage();    pushViewport(viewport(layout = grid.layout(1, 2)));
    print(plots5[[i]],vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(plots5[[i+1]],vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
 }
dev.off()


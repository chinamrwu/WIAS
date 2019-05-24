rm(list=ls())
library(glmnet)
library(caret)
library(ggplot2)
library(pROC)
library(umap)
source("src/common.R")

set.seed(190305)
options("width"=500)
setwd("E:/projects/TPD")
matAll <- read.table('data/TPD_prot_matrix_avg_20190304.txt',sep='\t',header=T,stringsAsFactors=F)
table(matAll$label)
matAll$label[matAll$label=='W'] <- 'C'
table(matAll$label)
rownames(matAll) <- matAll$SpecimenID
matAll <- matAll[,-2]
benign <- c('N','M','A')
malign <- c('C','P')

matAll$label[matAll$label %in% benign] <- 'B'
matAll$label[matAll$label %in% malign] <- 'M'
df0 <- matAll

label <- as.factor(matAll$label)
matAll <- matAll[,-1]
R0 <- apply(matAll,2,function(v){sum(is.na(v))/length(v)*100})
color2 <- c(B='#56B4E9',M='#CC79A7')

rate=65
selection0 <- colnames(matAll)[R0 <= rate]
hitFit <- NULL
mse=100
flg=T
fits <- list()

while(flg){
	tmp <- matAll[,selection0]
	tmp    <- t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	tmp[is.na(tmp)] <- 0
	
	cvfit <- cv.glmnet(tmp,label,family='binomial',alpha=0.5,type.measure='class')
	fits[[length(fits)+1]] <- cvfit
        cf <- coef(cvfit, s = "lambda.1se")
        selection1 <- rownames(cf)[cf[,1]!=0][-1]
        if(min(cvfit$cvm) <= mse){
           mse    <-  min(cvfit$cvm);
	   hitFit <-  cvfit
	}


	if(length(selection1)>= 4){ ## 4 proteins at least 
		    flg <- length(selection1) != length(selection0)
                    selection0 <- selection1
	}else{  flg=F	 }

}#flg
sts <- sapply(fits,function(obj){
     cf <- coef(obj,s = "lambda.min");
     a <- sum(cf[,1]!= 0)-1;
      b=min(obj$cvm);
      c(a,b)
 })

 plots <- list()
 for(obj in fits){
	cf <- coef(obj, s = "lambda.min")
	prots <- rownames(cf)[cf[,1]!=0][-1]
	tmp <- matAll[,prots]
	tmp$label <- label
	p1 <- drawPCA(tmp,color2,rowNormalization=T,colNormalization=T,strTitle=sprintf("PCA:65%% missingRate;%d prots;%3.4f",length(prots),min(obj$cvm)))
	t1 <- drawTSNE(tmp,color2,rowNormalization=T,colNormalization=T,strTitle=sprintf("tSNE:65%% missingRate;%d prots;%3.4f",length(prots),min(obj$cvm)))
	u1 <- drawUMAP(tmp,color2,rowNormalization=T,colNormalization=T,strTitle=sprintf("UMAP:65%% missingRate;%d prots;%3.4f",length(prots),min(obj$cvm)))
	plots[[length(plots)+1]] <- list("pca"=p1,"tsne"=t1,"umap"=u1)
 }

if(F){
   pdf("PCA_TSNE_UMAP_missing65_LASSO_0308.pdf",width=15,height=6)
   for(pts in plots){
     grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 3)));
    print( pts$pca,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print( pts$tsne,vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(pts$umap,vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
   }
   dev.off()
}

finePlots=list()
for(i in 1:length(plots)){
  obj <- plots[[i]]
  p1 <- obj$pca
  p1 <- p1 +scale_x_continuous(breaks = round(seq(min(p1$data$PC1), max(p1$data$PC1), by = 0.5),1)) +
            scale_y_continuous(breaks = round(seq(min(p1$data$PC2), max(p1$data$PC2), by = 0.5),1)) 
  obj$pca <- p1
  
  p1 <- obj$tsne
  p1 <- p1 +scale_x_continuous(breaks = round(seq(min(p1$data$X), max(p1$data$X), by = 6),1)) +
            scale_y_continuous(breaks = round(seq(min(p1$data$Y), max(p1$data$Y), by = 6),1)) 
  obj$tsne <- p1

   p1 <- obj$umap
   p1 <- p1 +scale_x_continuous(breaks = round(seq(min(p1$data$X), max(p1$data$Y), by = 0.5),1)) +
            scale_y_continuous(breaks = round(seq(min(p1$data$Y), max(p1$data$Y), by = 0.5),1)) 
  obj$umap <- p1
  finePlots[[length(finePlots)+1]] <- obj
}
obj <- finePlots[[2]]
p1 <- obj$pca
t1 <- obj$tsne
u1 <- obj$umap



###########################################################
uM <- rownames(u1$data)[u1$data$Y > 3.5 & u1$data$label=='M']
uB <- rownames(u1$data)[u1$data$Y < -1.0 & u1$data$label=='B']
trs <- c(uM,uB)
other <- setdiff(rownames(matAll),trs)
lbls <- as.factor(u1$data[trs,'label'])

rate=65
selection0 <- colnames(matAll)[R0 <= rate]
hitFit <- NULL
mse=100
flg=T
fits1 <- list()

while(flg){
	tmp <- matAll[trs,selection0]
	tmp    <- t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	tmp[is.na(tmp)] <- 0
	
	cvfit <- cv.glmnet(tmp,lbls,family='binomial',alpha=0.5,type.measure='class')
	fits1[[length(fits1)+1]] <- cvfit
        cf <- coef(cvfit, s = "lambda.1se")
        selection1 <- rownames(cf)[cf[,1]!=0][-1]
        if(min(cvfit$cvm) <= mse){
           mse    <-  min(cvfit$cvm);
	   hitFit <-  cvfit
	}


	if(length(selection1)>= 4){ ## 4 proteins at least 
		    flg <- length(selection1) != length(selection0)
                    selection0 <- selection1
	}else{  flg=F	 }

}#flg

sts1 <- sapply(fits1,function(obj){
     cf <- coef(obj,s = "lambda.min");
     a <- sum(cf[,1]!= 0)-1;
      b=min(obj$cvm);
      c(a,b)
 })

obj   <- fits1[[2]]
cf    <- coef(obj,s = "lambda.min");
prots <- rownames(cf)[cf[,1] !=0][-1]
tmp <- matAll[trs,prots]
tmp    <- t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
tmp[is.na(tmp)] <- 0

cvfit <- glmnet(tmp,lbls,family='binomial',alpha=0,lambda=obj$lambda.min)

test <- matAll[other,prots]
test <- t(apply(test,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
test[is.na(test)] <- 0
pred <- data.frame(predict(cvfit,newx = test,s = cvfit$lambda.min,type="class"))
pred$observed <- u1$data[other,"label"]
colnames(pred) <- c("predicted","observed")

obj <- fits[[2]]
cf    <- coef(obj,s = "lambda.min");
prots <- rownames(cf)[cf[,1] !=0][-1]


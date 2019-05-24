remove(list=ls())
library("randomForest")
library("caret")
library("ggplot2")
set.seed(20181119)
source("E:/projects/TPD/src/growRF.R")
options(width=350)

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
rownames(SWT) <- SWT$patientId
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求

Models=list()
for(tradeOff in c(5,30,50)){
	missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
	clnames <- names(missingRate)[missingRate<tradeOff]
	trainM <- SWT[,c('patientId',clnames)]
	trainM$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
	trainM <- trainM[,c("label",clnames)]

	missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
	clnames <- names(missingRate)[missingRate<tradeOff]
	testM <- SWV[,clnames]
	protIds <- intersect(colnames(trainM),colnames(testM))

	trainM <- trainM[,c('label',protIds)]
	testM <- testM[,protIds]

	######################################################################################
	models=list()
	trainM[is.na(trainM)] <- 0
        for(lbl in c('NM','MA','AC','CP')){
          cnt <- rep(0,length(protIds))
	  names(cnt) <- protIds
          for(k in 1:100){
	    model=growRF01(lbl,trainM)
	    cnt[model$protIds] <- cnt[model$protIds]+1
	  }
	  models[[length(models)+1]] <- cnt
       }
       Models[[length(Models)+1]] <- models
 }

protIds <- colnames(SWT)[-1]
t1 <- c()
for(models in Models){
   for (model in models){
     v <- rep(0,length(protIds))
     names(v) <- protIds
     v[names(model)] <- model
     t1 <- cbind(t1,v)
}
}

f1 <- function(v){ sum(v)!=0}
M0 <- SWT[,-1]
M0$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
f2 <- function(lbl){tmp <- M0[M0$label %in% lbl,]; a <- apply(tmp,2,function(v){sum(is.na(v))/length(v)*100});a}
missN  <- f2(c('N'))
missM  <- f2(c('M'))
missA  <- f2(c('A'))
missC  <- f2(c('C'))
missP  <- f2(c('P'))

missNM <- f2(c('N','M'))
missMA <- f2(c('M','A'))
missAC <- f2(c('A','C'))
missCP <- f2(c('C','P'))

modelNM <- t1[,c(1,5,9)]
colnames(modelNM) <- c(5,30,50)
modelNM <- modelNM[apply(modelNM,1,f1),]
sm <- apply(modelNM,1,sum)
sm <- sm[order(sm,decreasing=T)]
modelNM <- data.frame(modelNM[names(sm),])
modelNM$missNM <- missNM[rownames(modelNM)]

modelMA <- t1[,c(2,6,10)]
colnames(modelMA) <- c(5,30,50)
modelMA <- modelMA[apply(modelMA,1,f1),]
sm <- apply(modelMA,1,sum)
sm <- sm[order(sm,decreasing=T)]
modelMA <- data.frame(modelMA[names(sm),])
modelMA$missMA <- missMA[rownames(modelMA)]


modelAC <- t1[,c(3,7,11)]
colnames(modelAC) <- c(5,30,50)
modelAC <- modelAC[apply(modelAC,1,f1),]
sm <- apply(modelAC,1,sum)
sm <- sm[order(sm,decreasing=T)]
modelAC <- data.frame(modelAC[names(sm),])
modelAC$missAC <- missAC[rownames(modelAC)]

modelCP <- t1[,c(4,8,12)]
colnames(modelCP) <- c(5,30,50)
modelCP <- modelCP[apply(modelCP,1,f1),]
sm <- apply(modelCP,1,sum)
sm <- sm[order(sm,decreasing=T)]
modelCP <- data.frame(modelCP[names(sm),])
modelCP$missCP <- missCP[rownames(modelCP)]


f3 <- function(M2,strTitle){
    tmp <- M2[M2[,3]>0,c(4,3)]
    colnames(tmp) <- c('X','Y')
    p <- ggplot(tmp, aes(x=X, y=Y)) +geom_point()+
    xlab("missing rate (percentage)")+ylab("occurrence times in random forest model")+ggtitle(strTitle)
    p
  }



library(ggplot2)
library(grid)
pdf("E:/projects/TPD/results/missingRate_5_30_50.pdf")
grid.newpage();
pushViewport(viewport(layout = grid.layout(2, 2)));
    print(f3(modelNM,'NM'),vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(f3(modelMA,'MA'),vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(f3(modelAC,'AC'),vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(f3(modelCP,'CP'),vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()
 








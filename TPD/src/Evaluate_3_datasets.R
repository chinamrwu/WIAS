### Goal: evaluate the effects of transferring componen


### The following codes load data 
remove(list=ls())
library("data.table")
library("R.matlab")
library("mice")
options(width=350)
set.seed(20181224)

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
rownames(SWT) <- SWT$patientId
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
trainM <- SWT[,c('patientId',clnames)]
trainM$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
trainM <- trainM[,c("label",clnames)]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
testM <- SWV[,clnames]
protIds <- intersect(colnames(trainM),colnames(testM))

M0 <-  trainM[,protIds]   ## raw training matrix
T0  <-  testM[,protIds]   ## raw validation matrix

tmp  <- read.table("E:/projects/TPD/data/TPD_5percentage_mice_imputate.txt",sep=",",header=F)
colnames(tmp) <- colnames(M0)
M1    <- tmp[1:399,]   ### imputed by MICE
T1    <- tmp[400:580,]

tmp <- read.table("E:/projects/TPD/data/TPD_TCA_transformed.txt",sep="\t",header=F)
colnames(tmp) <- colnames(M0)
M2    <- tmp[1:399,]   ### transformed by TCA
T2    <- tmp[400:580,]
#################################################
prots0 <- read.table("E:/projects/TPD/results/proteins_used_in_100_RF.txt",sep="\t",header=T,stringsAsFactors=F)
prots1 <- read.table("E:/projects/TPD/results/proteins_used_in_100_RF_MICE_imputed.txt",sep=" ",header=T,stringsAsFactors=F)
prots2 <- read.table("E:/projects/TPD/results/proteins_used_in_100_RF_TCAed.txt",sep="\t",header=T,stringsAsFactors=F)

###################### kKN chains ############################
kNN <- function(M,features,K=7){
    tmp <- data.frame(M[,features])
    colnames(tmp) <- features
    tmp <- apply(tmp,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
    minv=min(tmp)
    tmp[is.na(tmp)] <- minv
    DM  <- as.matrix(dist(tmp))
    mxv=max(DM)
    for(i in 1:dim(DM)[1]){DM[i,i] <- mxv}
    clss <- unique(M$label)
    p1 <- apply(DM,2,function(v){  
      neighborLablels <- M$label[order(v)[1:K]]
      lblCounts <- sapply(clss,function(ch){ sum(neighborLablels==ch)})
      rtv       <- lblCounts;
      a         <- sum(lblCounts==max(lblCounts));
      ifelse(a==1,rtv <- c(rtv,names(lblCounts)[lblCounts==max(lblCounts)]),rtv <- c(rtv,neighborLablels[1]))
      rtv
    })
   p1 <- data.frame(t(p1))
   colnames(p1) <- c(clss,'predicted')
   p1$observed <- M$label
   p2 <- sapply(clss,function(ch){indx <- which(p1$observed==ch);sum(p1$observed[indx]==p1$predicted[indx])})
   p2 <- rbind(p2,sapply(clss,function(ch){sum(M$label==ch)}))
   p2 <- data.frame(p2)
   p2$total <- apply(p2,1,sum)
   p2 <- rbind(p2,apply(p2,2,function(v){100*v[1]/v[2]}))
   result=list()
   result$detail <- p1;
   result$stat <- p2
   result
}

library(data.table)

tmp <- M0;tmp$label <- label;

selectedProts=list()
for(lblPair in c('NM','MA','AC','CP')){
   lbl <- strsplit(lblPair,"")[[1]];
   tmpProts <- c();
   acc <- 0;
   indx <- order(prots0[,lblPair],decreasing = T)
   orderProts <- prots0$protId[indx]
   mAcc <- data.frame('protId'=orderProts,'Accuracy'=0)
   rownames(mAcc) <- orderProts;
   for(i in 1:length(orderProts)){
      t1   <- kNN(tmp[tmp$label %in% lbl,],orderProts[1:i]);
      print(sprintf("%s   %d  %f",lblPair,i,t1$stat[3,3]))
      mAcc$Accuracy[i] <- t1$stat[3,3]
   }
   selectedProts[[length(selectedProts)+1]] <- mAcc
}
names(selectedProts) <- c('NM','MA','AC','CP');

L=length(protIds)
topPair=list()
for(lblPair in c('NM','MA','AC','CP')){
   lbl <- strsplit(lblPair,"")[[1]];
   topAcc=0;
   bestPair=c();

   for(i in 1:(L-1)){
     for(j in (i+1):L){
        t1 <- kNN(tmp[tmp$label %in% lbl,],protIds[c(i,j)]);
	acc <- t1$stat[3,3]
	if(acc >=topAcc){
	    if(acc==topAcc){
	       bestPair <- c(bestPair,paste0(protIds[i],'_',protIds[j]));
	    }else{
	       bestPair <- paste0(protIds[i],'_',protIds[j])
	    }
	    topAcc=acc;
	    print(sprintf('%s :%s+%s  Acc=%f ',lblPair,protIds[i],protIds[j],topAcc))
	}
     }
   }
   topPair[[length(topPair)+1]] <- bestPair;
}


bestAcc=c(82.716049,80.769231,80.152672,89.403974)
names(bestAcc) <- c('NM','MA','AC','CP')
names(topPair) <- c('NM','MA','AC','CP');

bPairs <- read.table('E:/projects/TPD/data/bestPairProteins.txt',sep=" ",header=F,stringsAsFactors=F)
top3 <- apply(bPairs,1,function(v){
     pids <- setdiff(protIds,v[2:3])
     t1   <-  sapply(pids,function(v1){
          lbl <- strsplit(v[1],"")[[1]];
	  t2  <- kNN(tmp[tmp$label %in% lbl,],c(v[2:3],v1));
	  acc <- t2$stat[3,3]
	  
     })
     t1 <- data.frame('id'=pids,'acc'=t1)
     t1 <- t1[t1$acc >= v[4],]
     t1[t1$acc==max(t1$acc),]

   })
  #####################################

  plotTSNE <- function(DF,perplexity=10){
    M <- DF
    M[is.na(M)] <- 0
    indx <- match('label',colnames(M))

    tsn = tsne(M[,-indx],perplexity =perplexity)
    cnames <- colnames(M)[-indx]
    tsn <- data.frame(tsn,M[,indx])
    colnames(tsn) <- c("X","Y","label")
   #lblColors=c(M = "forestgreen", N = "gray0", P="orange3", C= "red",A="blue")
   lblColors <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   p <- ggplot(tsn, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
   P <- p+ xlab(paste0(cnames,collapse=","))

   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank())
   p <- p +   scale_colour_manual(values=lblColors)
   p
}

####################
#NM O75396 P02765 P06744 85.80247
#MA P07305 P21397 Q13263 82.69231
#AC O75083 P02765 P10809 81.67939
#CP P04083 P30044 P10809 91.39073
pr1 <- c('O75396','P02765','P06744')
pr2 <- c('P07305','P21397','Q13263')
pr3 <- c('O75083','P02765','P10809')
pr4 <- c('P04083','P30044','P10809')
ten <- unique(c(pr1,pr2,pr3,pr4))

p1 <- plotTSNE(tmp[,c('label',pr1)])
p2 <- plotTSNE(tmp[,c('label',pr1)])
p3 <- plotTSNE(tmp[,c('label',pr1)])
p4 <- plotTSNE(tmp[,c('label',pr1)])
p5 <- plotTSNE(tmp[,c('label',unique(c(pr1,pr2,pr3,pr4)))])
p6 <- plotTSNE(tmp[,c('label',unique(c(pr1,pr2,pr4)))])
p7 <- plotTSNE(tmp[,c('label',unique(c(pr1,pr3,pr4)))])
p8 <- plotTSNE(tmp[,c('label',unique(c(pr2,pr3,pr4)))])
p9 <- plotTSNE(tmp[,c('label',unique(c(pr1,pr2,pr3)))])

##################################
spC=c("C13","C2","C22", "C24","C25","C26","C29","C33","C34","C36","C37","C39","C4","C41","C45","C47","C7","C8")
indx <- match(spC,rownames(tmp))
tmp$label[indx] <- 'C1'
tmp$label[tmp$label=='C'] <- 'W'
tmp$label[tmp$label=='C1'] <- 'C'

pACW <- plotTSNE(tmp[tmp$label %in% c('A','C','W'),c('label',unique(c(pr1,pr2,pr3,pr4)))])

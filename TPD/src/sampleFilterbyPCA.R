remove(list=ls())
library(randomForest)
library(caret)
library(data.table)
source('E:/projects/TPD/src/plotTools.R')
df0 <- read.table("E:/projects/TPD/data/TPD_prot_matrix_avg_20190115.txt",sep="\t",header=T,stringsAsFactors=F)
df1 <- read.table("E:/projects/TPD/data/AllSampleInformation_20190115.txt",sep="\t",header=T,stringsAsFactors=F)

df0$label[df0$label=='W'] <- 'C'
t0 <- df0[,-1]
rownames(t0) <- df0[,1]

protNM <- c('Q6NZI2','P06702','P07355','P02765')
protAC <- c('Q14103','O00468','Q12906','Q6PCB0','P13489','P08133')
protCP <- c('Q14764','P02751','P14543')
#protMA <- c('P21397','P35268','P12110','P08123')
#protMA <- c('P21397','P35268','Q9Y646','P37802')
# protMA <- c('P21397','P35268','Q9Y646','P37802','Q13263')
#P21397 P35268 Q9Y646 Q09666 Q13263
protMA <- c('P21397','P35268','Q9Y646','Q09666','Q13263','P37802','P12110')
sigProts <- unique(c(protNM,protAC,protCP,protMA))
###########################################################################
getSubset <- function(Labels,missRate){
   tmp  <- t0[t0$label %in% Labels,]
   R0   <- apply(tmp,2,function(v){sum(is.na(v))/length(v)*100})
   tmp  <- tmp[,R0<=missRate]

   K=3
   t1 <- sapply(1:dim(tmp)[1],function(i) {
     lbl <- tmp[i, 'label']
     v0  <- tmp[i, colnames(tmp) != 'label']
     feature0 <- names(v0)[is.na(v0)]
     if(length(feature0)>0){
	     k0  <- tmp[-i,]
	     k0  <- k0[k0$label==lbl,-1]
	     k0  <- rbind(v0,k0)
	     k0 <- apply(k0,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
             k0[is.na(k0)] <- 0

	     D0  <- as.matrix(dist(k0))
	     d0  <- order(D0[1,])[2:(K+1)]
	     if(length(feature0)==1){
	        v0[feature0] <- mean(tmp[d0,feature0],na.rm=T)
	     }else{
                v0[feature0] <- apply(tmp[d0,feature0],2,function(v){mean(v,na.rm=T)})
	     }
      }
     v0
   })

   t1 <- apply(t(t1),2,unlist)    
   rownames(t1) <- rownames(tmp)
   t1 <- data.frame(t1)
   clnames <- colnames(t1)
   t1$label <- tmp$label
   t1 <- t1[,c('label',clnames)]
   print(sprintf("get subset for %s with %d rows and %d features",paste0(Labels,collapse=","),dim(t1)[1],dim(t1)[2]))
   t1
}
############################################################
NM <- getSubset(c('N','M'),5)[,c('label',protNM)]
AC <- getSubset(c('A','C'),5)[,c('label',protAC)]
CP <- getSubset(c('C','P'),5)[,c('label',protCP)]
MA <- getSubset(c('M','A'),5)[,c('label',protMA)]
#################################################################
#lblColors <- c(training='#56B4E9',validation='#CC79A7')
p1 <-  drawPCA(CP,c(P='#537e35',C='#f5b819'),colNormalization=T)
p2 <-  drawPCA(AC,c(A='#537e35',C='#f5b819'),colNormalization=T)
p3 <-  drawPCA(MA,c(A='#537e35',M='#f5b819'),colNormalization=T)
p4 <-  drawPCA(NM,c(N='#537e35',M='#f5b819'),colNormalization=T)

filterCP <- rownames(p1$data)[p1$data$PC1>=0.5 & p1$data$label=='P']
filterCP <- c(filterCP,rownames(p1$data)[p1$data$PC1<=0.5 &p1$data$PC2 <= -2 & p1$data$label=='C'])
CP0 <- CP[setdiff(rownames(CP),filterCP),]
pCP0 <- drawPCA(CP0,c(P='#537e35',C='#f5b819'),colNormalization=T)

AC0 <- AC[setdiff(rownames(AC),filterCP),]
pAC0 <- drawPCA(AC0,c(A='#537e35',C='#f5b819'),colNormalization=T)
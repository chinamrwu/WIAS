remove(list=ls())
library(randomForest)
library(caret)
library(data.table)
df0 <- read.table("E:/projects/TPD/data/TPD_prot_matrix_avg_20190115.txt",sep="\t",header=T,stringsAsFactors=F)
df0$label[df0$label=='W'] <- 'C'
t0 <- df0[,-1]
rownames(t0) <- df0[,1]
################# get specified label subset 
getSubset <- function(lbl,missRate){
   tmp  <- t0[t0$label %in% lbl,]
   R0   <- apply(tmp,2,function(v){sum(is.na(v))/length(v)*100})
   tmp  <- tmp[,R0<=missRate]

   K=3
   t1 <- sapply(1:dim(tmp)[1],function(i) {
     lbl <- tmp[i, 1]
     v0  <- tmp[i,-1]
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
   print(sprintf("get subset for %s-%s with %d rows and %d features",lbl[1],lbl[2],dim(t1)[1],dim(t1)[2]))
   t1
}
########################################################################
# NM: Q6NZI2 P06702 P07355 P02765  80%
#MA: 
#AC: Q14103 O00468 Q12906 Q6PCB0 P13489 P08133 70%
#CP: Q14764 P02751 P14543 90%

protNM <- c('Q6NZI2','P06702','P07355','P02765')
protAC <- c('Q14103','O00468','Q12906','Q6PCB0','P13489','P08133')
protCP <- c('Q14764','P02751','P14543')
#protMA <- c('P21397','P35268','P12110','P08123')
#protMA <- c('P21397','P35268','Q9Y646','P37802')
# protMA <- c('P21397','P35268','Q9Y646','P37802','Q13263')
#P21397 P35268 Q9Y646 Q09666 Q13263
protMA <- c('P21397','P35268','Q9Y646','Q09666','Q13263','P37802','P12110')

allProts <- unique(c(protNM,protAC,protCP,protMA))

NM <- getSubset(c('N','M'),5)[,c('label',protNM)]
AC <- getSubset(c('A','C'),5)[,c('label',protAC)]
CP <- getSubset(c('C','P'),5)[,c('label',protCP)]
MA <- getSubset(c('M','A'),5)[,c('label',protMA)]

plotPCA <- list()
plotPCA[[length(plotPCA)+1]] <- drawPCA(NM,c(N='#537e35',M='#e17832'),colNormalization=T)
plotPCA[[length(plotPCA)+1]] <- drawPCA(AC,c(A='#f5b819',C='#5992c6'),colNormalization=T)
plotPCA[[length(plotPCA)+1]] <- drawPCA(CP,c(P='#537e35',C='#f5b819'),colNormalization=T)
plotPCA[[length(plotPCA)+1]] <- drawPCA(MA,c(A='#537e35',M='#e17832'),colNormalization=T)

ptSNE <- list()
ptSNE[[length(ptSNE)+1]] <- drawTSNE(NM,c(N='#537e35',M='#e17832'),colNormalization=T)
ptSNE[[length(ptSNE)+1]] <- drawTSNE(AC,c(A='#f5b819',C='#5992c6'),colNormalization=T)
ptSNE[[length(ptSNE)+1]] <- drawTSNE(CP,c(P='#537e35',C='#f5b819'),colNormalization=T)
ptSNE[[length(ptSNE)+1]] <- drawTSNE(MA,c(A='#537e35',M='#e17832'),colNormalization=T)

library(ggplot2)
library(grid)
pdf("E:/projects/TPD/results/feature_selection_visualization.pdf")


lblColors <- c(N='black',M='#e17832',A='#f5b819',C='#5992c6',P='#537e35')


plots1 <- list()
for(obj in proteins){
  #obj <- proteins[[2]][1:2,]
  m=dim(obj)[1]
  n=dim(obj)[2]
  for(i in 1:m) {
    v=obj[i,]
    acc=sprintf("%4.3f",v[n])
    tmp <- getSubset(c('M','A'),5)[,c('label',as.character(v[1:(n-1)]))]
    plots[[length(plots)+1]] <- drawPCA(tmp,c(A='#537e35',M='#e17832'),rowNormalization=F,colNormalization=T,acc=acc)
    plots[[length(plots)+1]] <- drawTSNE(tmp,c(A='#537e35',M='#e17832'),rowNormalization=F,colNormalization=T,acc=acc)
  }
 }
pdf("E:/projects/TPD/results/feature_selection_4MA.pdf",width=10,height=5)
    for(i in seq(1,length(plots),by=2)) {
    grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 2)));
      print(plots[[i]],vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(plots[[i+1]],vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    }
dev.off()

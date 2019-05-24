rm(list=ls())
library(ggplot2)
library(data.table)
library(sqldf)
##############################################
drawPCA <- function(DF,ptColors,pointSize=4,rowNormalization=F,colNormalization=F,strTitle='PCA'){ 
## M is a matrix or dataframe, rows are samples and columns are features, rownames are sample names
   M <- DF[,colnames(DF)!='label']
   if(rowNormalization==T){
      M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
      #M <- data.frame(t(apply(M,1,function(v){v/mean(v,na.rm=T)})))
      #print('Normalization by row is done!')
   }
   if(colNormalization==T){
      M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
   }
   clnames <- colnames(DF)[colnames(DF)!='label']
   M[is.na(M)] <- 0
   m1 <- prcomp(M,colNormalization);
   Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
   Y  <- Y[,c(1,2)]
   
   Y <- data.frame(Y,DF$label);
   colnames(Y) <- c("PC1","PC2","label")

   eigs <- m1$sdev^2
   percentages <- eigs[1:2] / sum(eigs)
   # lblColors <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   #lblColors <- c(training='#537e35',validation='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=pointSize, shape=16, stroke = 2)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    plot.title   = element_text(size=16),
	    axis.line.x = element_line(color="black", size = 1),
	    axis.line.y = element_line(color="black", size = 1),
	    panel.background = element_blank())
            
   strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
   p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
	      title =strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
  
   p
}
getSubset <- function(Labels,M0,missRate,imputation=F){
   lbls <- Labels
   lbls <- if(length(Labels)==1){lbls <- unlist(strsplit(lbls,""))}else{lbls <- Labels}
   tmp  <- M0[M0$label %in% lbls,]
   R0   <- apply(tmp,2,function(v){sum(is.na(v))/length(v)*100})
   tmp  <- tmp[,R0<=missRate]

   #if(!imputation){tmp[is.na(tmp)] <- 0}
   if(imputation){ #### kNN imputation 
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
      #print(sprintf("get subset for %s with %d rows and %d features",paste0(lbls,collapse=","),dim(t1)[1],dim(t1)[2]))
      tmp <- t1
   }
   return(tmp)
}


#################################################
lblColors <- c(B='#56B4E9',M='#CC79A7')
###################################################################################
df0 <- read.table('data/TPD_prot_matrix_avg_20190131.txt',sep='\t',header=T,stringsAsFactors=F)
rownames(df0) <- df0$SpecimenID
df0 <- df0[,-1]
sampleInf <- read.table("data/AllSampleInformation_20190115.txt",sep="\t",header = T,stringsAsFactors = F)[,c(2,4,6)]
sampleInf <- sqldf("SELECT distinct * from sampleInf")
rownames(sampleInf) <- sampleInf[,1]
df0$label <- sampleInf[rownames(df0),2]
###################################################
ds0 <- getSubset(c('training','validation'),df0,0,F)
ds5 <- getSubset(c('training','validation'),df0,5,F)
ds10 <- getSubset(c('training','validation'),df0,10,F)
ds20 <- getSubset(c('training','validation'),df0,20,F)
ds30 <- getSubset(c('training','validation'),df0,30,F)
ds40 <- getSubset(c('training','validation'),df0,40,F)
ds50 <- getSubset(c('training','validation'),df0,50,F)
ds100 <- getSubset(c('training','validation'),df0,100,F)
sigProts <- c('P04083','P00338','Q99584','P40121','P78417','P31949','P16403','O60437','P04004','Q02880','P02765','Q08380')

p0 <- drawPCA(ds0,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA without missing value %d proteins',(dim(ds0)[2]-1)))
p5 <- drawPCA(ds5,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:5%% missing %d proteins',(dim(ds5)[2]-1)))
p10 <- drawPCA(ds10,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:10%% missing %d proteins',(dim(ds10)[2]-1)))
p20 <- drawPCA(ds20,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:20%% missing %d proteins',(dim(ds20)[2]-1)))
p30 <- drawPCA(ds30,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:30%% missing %d proteins',(dim(ds30)[2]-1)))
p40 <- drawPCA(ds40,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:40%% missing %d proteins',(dim(ds40)[2]-1)))
p50 <- drawPCA(ds40,lblColors,rowNormalization=T,colNormalization=T,strTitle=sprintf('PCA:50%% missing %d proteins',(dim(ds50)[2]-1)))
p100 <- drawPCA(ds100,lblColors,rowNormalization=T,colNormalization=F,strTitle=sprintf('PCA: all %d proteins',(dim(ds100)[2]-1)))

p12 <- drawPCA(df0[,c('label',sigProts)],lblColors,rowNormalization=T,colNormalization=F,strTitle=sprintf("PCA:%s",paste0(sigProts,collapse=",")))

pdf("PCA_training_validation_20190227.pdf",width=18,height=8)
 grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 3)));
    print( p0,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print( p5,vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p10,vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
 grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 3)));
    print( p20,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print( p30,vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p10,vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
 grid.newpage();
    pushViewport(viewport(layout = grid.layout(1, 3)));
    print( p0,vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print( p5,vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p10,vp=viewport(layout.pos.row = 1, layout.pos.col = 3))


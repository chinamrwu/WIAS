library("FactoMineR")
library("factoextra")
library("grid")
library("data.table")
library("ggplot2")
library(pheatmap);
library("RColorBrewer");
library("tsne");


###tSNE plot DE: dataframe with rows as sample and columns as features,a column with name 'label' is required, represents the label of samples
# lblColors: a vector containing colors for each label, like  lblColors=c(M = "forestgreen", N = "gray0", P="firebrick", C= "red",A="blue")
drawTSNE <- function(DF,ptColors,rowNormalization=F,colNormalization=F,perplexity=10,acc=""){
    M <- DF[,colnames(DF)!='label']
    if(rowNormalization){M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))}
    if(colNormalization){M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})}
    M[is.na(M)] <- 0
    indx <- match('label',colnames(M))
    clnames <- colnames(DF)[colnames(DF)!='label']
    strTitle <- sprintf("tSNE:%s  %s",paste0(clnames,collapse=','),acc)
    tsn = tsne(M,perplexity =perplexity)
    cnames <- colnames(M)
    tsn <- data.frame(tsn,DF$label)
    colnames(tsn) <- c("X","Y","label")
    tsn <- tsn[-c(which(tsn$X==min(tsn$X)),which(tsn$Y==min(tsn$Y))),]
     tsn <- tsn[-c(which(tsn$X==max(tsn$X)),which(tsn$Y==max(tsn$Y))),]
   #lblColors <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')

   #lblColors <- c(A='#537e35',M='#e17832',N='#f5b819',B='#5992c6',C='#282f89',W='mediumorchid3')
   p <- ggplot(tsn, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	     plot.title = element_text(size=8),
	    panel.background = element_blank())
   p <-  p +  labs(x =paste0(clnames,collapse=","),y='',title=strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
   p
}

drawPCA <- function(DF,ptColors,rowNormalization=F,colNormalization=F,acc=""){ 
## M is a matrix or dataframe, rows are samples and columns are features, rownames are sample names
   M <- DF[,colnames(DF)!='label']
   if(rowNormalization){
      M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
      #print('Normalization by row is done!')
   }
   clnames <- colnames(DF)[colnames(DF)!='label']
   strTitle <- sprintf("PCA: %s    %s",paste0(clnames,collapse=','),acc)
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
   p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=5)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    plot.title = element_text(size=8),
	    panel.background = element_blank())
            
   strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
   p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
	      title =strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
  
   p
}

# parse gene symbol from uniprot website
# acc: protein access Id like ' "O00264" "O00468" "O14773" "O14979" "O15230" "O43175" "O43707"
# return: the gene symbol for the acc
geneName <- function(acc){
   kk <- read.table(sprintf('https://www.uniprot.org/uniprot/%s.fasta',acc),sep="\n",stringsAsFactors=F,header=F)
   a <- strsplit(kk[1,]," ")[[1]]
   b <- a[grepl("GN=",a)]
   strsplit(b,"=")[[1]][2]
}

# This file includes PCA UMAP TSNE and Vocano plot 
# copyright: wuzhicheng@fudan.edu.cn
#
library("FactoMineR")
library("factoextra")
library("grid")
library("data.table")
library("ggplot2")
library(pheatmap);
library("RColorBrewer");
library("tsne");
library(umap)

# draw PCA plot 

#
#DF is a dataframe, rows are samples and columns are proteins, rownames are sample names and DF must have a column
#named label to indicate which class the corresponding row below to.
#ptColors£º point colors example as : c('B'='red','M'='black'). the number of colors equals with the number of distinct classes 
# strTitle: title text showed on the plot.
# 

drawPCA<- function(DF,ptColors,rowNormalization=T,colNormalization=T,strTitle=NULL){ 
   M <- DF[,colnames(DF)!='label']
   if(rowNormalization){
      M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
      #print('Normalization by row is done!')
   }
   if(colNormalization){
      M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
   }
   clnames <- colnames(DF)[colnames(DF)!='label']
   M[is.na(M)] <- 0
   m1 <- prcomp(M,colNormalization);
   Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
   Y  <- Y[,c(1,2)]
   
   Y <- data.frame(Y,DF$label);
   colnames(Y) <- c("PC1","PC2","label")
   if(is.null(strTitle)){
      strTitle <- sprintf("PCA:%d features",length(clnames))
	}
   eigs <- m1$sdev^2
   percentages <- eigs[1:2] / sum(eigs)
   # lblColors <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   #lblColors <- c(training='#537e35',validation='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    plot.title   = element_text(size=16),
	    panel.background = element_blank())
            
   strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
   p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
	      title =strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
  
   p
}

###tSNE plot DE: dataframe with rows as sample and columns as features,a column with name 'label' is required, represents the label of samples
# lblColors: a vector containing colors for each label, like  lblColors=c(M = "forestgreen", N = "gray0", P="firebrick", C= "red",A="blue")
drawTSNE <- function(DF,ptColors,rowNormalization=F,colNormalization=F,perplexity=10,strTitle='tSNE'){
    M <- DF[,colnames(DF)!='label']
    if(rowNormalization){M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))}
    if(colNormalization){M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})}
    M[is.na(M)] <- 0
    indx <- match('label',colnames(M))
    clnames <- colnames(DF)[colnames(DF)!='label']
    tsn = tsne(M,perplexity =perplexity)
    cnames <- colnames(M)
    tsn <- data.frame(tsn,DF$label)
    colnames(tsn) <- c("X","Y","label")
    rownames(tsn) <- rownames(M)
    #tsn <- tsn[-c(which(tsn$X==min(tsn$X)),which(tsn$Y==min(tsn$Y))),]
    #tsn <- tsn[-c(which(tsn$X==max(tsn$X)),which(tsn$Y==max(tsn$Y))),]
   #lblColors <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')

   #lblColors <- c(A='#537e35',M='#e17832',N='#f5b819',B='#5992c6',C='#282f89',W='mediumorchid3')
   p <- ggplot(tsn, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	     plot.title = element_text(size=15),
	    axis.line.x = element_line(color="black", size = 0.5),
	    axis.line.y = element_line(color="black", size = 0.5),
	    panel.background = element_blank())
   p <-  p +  labs(title=strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
   p
}

# parse gene symbol from uniprot website
# acc: protein access Id like ' "O00264" "O00468" "O14773" "O14979" "O15230" "O43175" "O43707"
# return: the gene symbol for the acc
geneName <- function(acc){kk <- read.table(sprintf('https://www.uniprot.org/uniprot/%s.fasta',acc),sep="\n",stringsAsFactors=F,header=F)
   a <- strsplit(kk[1,]," ")[[1]]
   b <- a[grepl("GN=",a)]
   strsplit(b,"=")[[1]][2]
}
###############################################################################################################
drawVolcano <- function(df1,pdfPath,strTitle="volcano plot",outFile=NULL){ 
   label <- df1$label
   lbls = unique(df1$label)
   table(df1$label) # A143,C75
   fc <- apply(2^df1[,colnames(df1)!='label'],2,function(x){
       log2( mean(na.omit(x[label==lbls[1]])) /mean(na.omit(x[label==lbls[2]])))
    })
   
   df1[is.na(df1)] <- 0
   pValue <- apply(df1[,colnames(df1)!='label'], 2, function(v) {
       p1 <- t.test(v[label == lbls[1]], v[label == lbls[2]], paired = F, var.equal = F)$p.value
       p1 <-  p.adjust(p1,method="BH")
       p1
    }) 

  pdf(pdfPath)
     plot(fc, -log10(pValue), col = '#00000033', pch = 19,xlab = 'log2(FC)', ylab = '-log10(p-value)', main = strTitle)
     abline(h = 1.3, v = c(-log2(1.5),log2(1.5)), lty = 2, lwd = 1)
  
     up  <- fc >= log2(1.5) & pValue <= 0.05
     points(fc[up], -log10(pValue[up]), col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
 
     down <- fc <= -log2(1.5) & pValue <= 0.05
     points(fc[down], -log10(pValue[down]), col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 2)
  dev.off()
  
  name = list(up = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[up],fc = fc[up], p_value = pValue[up],
                            type = rep("Upregulation",sum(up)),stringsAsFactors=F),
            down = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[down],fc = fc[down], p_value = pValue[down],
                              type = rep("Downregulation",sum(down))),stringsAsFactors=F)
   name1 = rbind(name[[1]],name[[2]])
   rownames(name1) <- 1:dim(name1)[1]
   
   if( !is.null(outFile) ){
     write.table(name1,file=outFile,sep="\t",col.names=T,row.names=F,quote=F)
   }
   #write.xlsx(name1, "Table_1_TPD_volcano_prot_AC10_fc1.5_190131.xlsx",showNA = T,row.names = F)
   sum(up)
   sum(down)
   return(name1)
}
# see drawPCA
drawUMAP <- function(M1,ptColors, strTitle="UMAP",rowNormalization=T,colNormalization=T){
     
     if(!'label' %in% colnames(M1)){
        print('A column with named label must existed in data frame')
	return(NULL)
     }

     tmp <- M1[,colnames(M1)!='label']
     if(rowNormalization){
        tmp <- data.frame(t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})),stringsAsFactors=F)
        rownames(tmp) <- rownames(M1)
     }
     if(colNormalization){    tmp <- apply(tmp,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})   }
     	tmp[is.na(tmp)] <- 0

     obj = umap(d=tmp,method='naive')
     clnames <- colnames(tmp)
     df1 <- data.frame(obj$layout)
     df1$label <- M1$label
     colnames(df1) <- c('X','Y','label')
     
     p <- ggplot(df1, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
     p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    plot.title   = element_text(size=16),
	    axis.line.x = element_line(color="black", size = 0.5),
	    axis.line.y = element_line(color="black", size = 0.5),

	    panel.background = element_blank())
            
   #strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
   p <- p +  labs(title =strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
   p
 }
Volcano <- function(df1,pdfPath,outFile=NULL,thresholdFC=1.5,thresholdPValue=0.05,strTitle="volcano plot"){ 
   label <- df1$label
   lbls = unique(df1$label)
   table(df1$label) # A143,C75
   fc <- apply(2^df1[,colnames(df1)!='label'],2,function(x){
       log2( mean(na.omit(x[label==lbls[1]])) /mean(na.omit(x[label==lbls[2]])))
    })
   
   df1[is.na(df1)] <- 0
   pValue <- apply(df1[,colnames(df1)!='label'], 2, function(v) {
       p1 <- t.test(v[label == lbls[1]], v[label == lbls[2]], paired = F, var.equal = F)$p.value
       p1 <-  p.adjust(p1,method="BH")
       p1
    }) 
  if(!is.null(pdfPath)){
	   pdf(pdfPath)
	     plot(fc, -log10(pValue), col = '#00000033', pch = 19,xlab = 'log2(FoldChange)', ylab = '-log10(P-value)', main = strTitle)
	     abline(h = -log10(thresholdPValue), v = c(-log2(thresholdFC),log2(thresholdFC)), lty = 2, lwd = 1)
	  
	     up  <- fc >= log2(thresholdFC) & pValue <= thresholdPValue
	     points(fc[up], -log10(pValue[up]), col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
	 
	     down <- fc <= -log2(thresholdFC) & pValue <= thresholdPValue
	     points(fc[down], -log10(pValue[down]), col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 2)
	  dev.off()
  }

  name = list(up = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[up],fc = fc[up], p_value = pValue[up],
                            type = rep("Upregulation",sum(up)),stringsAsFactors=F),
            down = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[down],fc = fc[down], p_value = pValue[down],
                              type = rep("Downregulation",sum(down))),stringsAsFactors=F)
   name1 = rbind(name[[1]],name[[2]],stringsAsFactors=F)
   rownames(name1) <- 1:dim(name1)[1]
   name1 <- name1[order(abs(name1$fc),decreasing=T),]
   if( !is.null(outFile) ){
     write.table(name1,file=outFile,sep="\t",col.names=T,row.names=F,quote=F)
   }
   #write.xlsx(name1, "Table_1_TPD_volcano_prot_AC10_fc1.5_190131.xlsx",showNA = T,row.names = F)
   sum(up)
   sum(down)
   return(name1)
 }

 #
 #
 #
 #
 #

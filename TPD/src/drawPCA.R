library("FactoMineR")
library("factoextra")
library("grid")
library("data.table")
library("ggplot2")

drawPCA <- function(df0,strTitle){ ## M is a matrix or dataframe, rows are samples and columns are features, rownames are sample names
   M <- df0[,colnames(df0)!='label'] 
   m1 <- prcomp(M,F);
   Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
   Y  <- Y[,c(1,2)]
   
   Y <- data.frame(Y,df0$label);
   colnames(Y) <- c("PC1","PC2","label")

   eigs <- m1$sdev^2
   percentages <- eigs[1:2] / sum(eigs)
   

   p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=2.5)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank())
   p <- p +  labs(x = sprintf("PC1(%4.2f%%)",percentages[1]*100),
                  y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
	      title =strTitle)
  
   p
}
## study  kNN algorithmns 
remove(list=ls())
library(data.table)

kNNPredict <- function(features,DS,K=5){
   tmp <- apply(DS[,features],2,function(v){mv <- mean(v,na.rm=T);(v-mv)/sd(v,na.rm=T)})
   DM  <- as.matrix(dist(tmp))
   L <- dim(DS)[1]
   Labels <- unique(DS$label)

   for(i in 1:L){  DM[i,i] <- .Machine$double.xmax}
   prb <- apply(DM,2,function(v){
	 nns <- DS$label[order(v)][1:K]
	 votes <- sapply(Labels,function(lbl){ sum(nns==lbl)/K})
	 votes
   }                                                                                                                                            
   
   colnames(prb) <- Labels
   prb <- data.frame(prb)
   prb$predicted <- apply(prb,1,function(v){names(v)[max(v)==v]})
   prb   
}

predictionStatistic <- function(predictions,realLabels){
   predLabels <- predictions[,"predicted"]
   Labels <- unique(realLabels)
   indexP <- which(realLabels==Labels[1])
   indexN <- which(realLabels==Labels[2])
   

}

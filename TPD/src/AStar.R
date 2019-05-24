### AStar algorithm + random forest for features selection
library(randomForest)
library(data.table)
df0 <- read.table('E:/projects/TPD/data/TPD_prot_matrix_avg_20190131.txt',sep='\t',header=T,stringsAsFactors=F)[,-1]
missingRate <- 5
imputation  <- F

getSubset <- function(Labels,M0,missRate){
   lbls <- Labels
   lbls <- ifelse(length(Labels)==1,lbls <- strsplit(lbls,"")[[1]],Labels)
   tmp  <- M0[M0$label %in% lbls,]
   R0   <- apply(tmp,2,function(v){sum(is.na(v))/length(v)*100})
   tmp  <- tmp[,R0<=missRate]

   if(!imputation){tmp[is.na(tmp)] <- 0}
   else{ #### kNN imputation 
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

searchFeatures <- function(M0,level=6){
   RF <- randomForest(label ~ . ,data=M0,importance=T,ntree=1000,nodesize=5)
   importance <- RF$importanceSD
   importance <- importance[order(importance[,3],decreasing=T),]
   importance <- rownames(importance)

   L <- length(importance)
   bestScore <- 0
   index     <- 1
   
   while(index<L){
     currentFeatures <- importance[index]
     


   }




}
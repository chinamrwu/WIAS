remove(list=ls())
library("randomForest")
library("pROC")
set.seed(20181119)
source("E:/projects/TPD/src/growRF.R")

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
SWT <- SWT[,c("patientId",names(missingRate)[missingRate<tradeOff])]
SWT$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
SWT <- SWT[,c("label",names(missingRate)[missingRate<tradeOff])]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
SWV <- SWV[,names(missingRate)[missingRate<tradeOff]]

protIds <- intersect(colnames(SWT),colnames(SWV))
SWT <- SWT[,c("label",protIds)]
#SWT[,-1] <- t(apply(SWT[,-1],1,function(v){v/mean(v,na.rm=T)}))
SWV <- SWV[,protIds]
#SWV <- data.frame(t(apply(SWV,1,function(v){v/mean(v,na.rm=T)})))
SWT[is.na(SWT)] <- 0
SWV[is.na(SWV)] <- 0
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
protIds <- intersect(colnames(SWT),colnames(SWV))

#######################################################################################################
Labels <- c('M','A','C','P')
lblPairs <- c()

featureAbility <- list()
for(i in 1:3){
   for(j in (i+1):4){
       label <- c(Labels[i],Labels[j])
       t1 <- t(sapply(protIds,function(id){
         t0 <- singleFeatureKNN(id,label,SWT)
	 c(length(t0$yes)/(length(t0$yes)+length(t0$no)),t0$performance)
       }))
       colnames(t1) <- c("accuracy",label)
       featureAbility[[length(featureAbility)+1]] <- t1
       lblPairs <- c(lblPairs,paste0(label,collapse=""))
      }
}
names(featureAbility) <- lblPairs

if(F){
complementaryAbility <- list()
for(i in 1:3){
   for(j in (i+1):4){
       label <- c(Labels[i],Labels[j])
       
	cmp <- matrix(0,nrow=length(protIds),ncol=length(protIds))
	colnames(cmp) <- protIds
	rownames(cmp) <- protIds
	for(m in 1:(length(protIds)-1)){
		for(n in (m+1):length(protIds)){
			cmp[m,n] <- complementarity(protIds[m],protIds[n],label,SWT)
			cmp[n,m] <- cmp[m,n]
		}
	}
    complementaryAbility[[length(complementaryAbility)+1]] <- cmp
   }
}
names(complementaryAbility) <- lblPairs
}
comp2 <- list()
for(i in 1:3){
   for(j in (i+1):4){
        label <- paste0(c(Labels[i],Labels[j]),collapse="")
        ability <- featureAbility[[label]]
	ds <- as.matrix(dist(ability[,-1]))
        comp2[[length(comp2)+1]] <- ds 
   }
}
names(comp2) <- lblPairs

label <- c('M','A')
fpm <- featureAbility[[paste0(label,collapse="")]]
fpm <- fpm[order(fpm[,1],decreasing=T),]
tM  <- SWT[SWT$label %in% label,]
L <- length(protIds)

K=5;


reduceK <- function(K){
   t1 <- sapply(K,function(i){
	    f1 <- rownames(fpm)[i] ### best feature
	    p1 <- singleFeatureKNN(f1,label,tM)
	    p1T <- sum(p1$prediction==tM$label)
	    p1T
   })
   p2 <- multiFeatureKNN(rownames(fpm)[K],label,tM)
   c(t1,sum(p2$prediction==tM$label))
}

combin <- function(f1,label,M0){
       tM <- M0[M0$label %in% label,]
       p1 <- singleFeatureKNN(f1,label,tM)
       p1T <- sum(p1$prediction==tM$label)
       fpm <- featureAbility[[paste0(label,collapse="")]]

       v <- comp2[[paste0(label,collapse="")]][f1,]
       v <- fpm[names(v),1] * v
       cf <- names(v)[order(v,decreasing = T)][1]
       p2 <- singleFeatureKNN(cf,label,tM)
       p2T <- sum(p2$prediction==tM$label)

       p3<-multiFeatureKNN(c(f1,cf),label,tM)
       p3T <- sum(p3$prediction==tM$label)
       rts <- c(p1T,p2T,p3T) 
       names(rts) <- c(f1,cf,"both")
       rts
}





for(i in 1:10){
   for(j in (i+1):L){
    
    f1 <- rownames(fpm)[i] ### best feature
    p1 <- singleFeatureKNN(f1,label,tM)
    p1T <- sum(p1$prediction==tM$label)

    f2 <- rownames(fpm)[j] ### best feature
    p2 <- singleFeatureKNN(f2,label,tM)
    p2T <- sum(p2$prediction==tM$label)

    p3<-multiFeatureKNN(c(f1,f2),label,tM)
    p3T <- sum(p3$prediction==tM$label)

   # print(c(p1$performance,p2$performance,p3$performance))
    if(p3T < max(c(p1T,p2T,p3T))){
      print(c(p1T,p2T,p3T))
    }
 }
 }


   fcm <- comp2[[paste0(label,collapse="")]];
   v <- fcm[f1,]
   v <- v[order(v,decreasing=T)]
   f2 <- names(v)[1]
   p2 <- multiFeatureKNN(c(f1,f2),label,SWT[SWT$label %in% label,])





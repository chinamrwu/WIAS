remove(list=ls())
library(randomForest)
library(data.table)
library(caret)
remove(list=ls())
library(randomForest)
library(caret)
library(data.table)
#################################################################################################################
protNM <- c('Q6NZI2','P06702','P07355','P02765')
protAC <- c('Q14103','O00468','Q12906','Q6PCB0','P13489','P08133')
protCP <- c('Q14764','P02751','P14543')
#protMA <- c('P21397','P35268','P12110','P08123')
#protMA <- c('P21397','P35268','Q9Y646','P37802')
# protMA <- c('P21397','P35268','Q9Y646','P37802','Q13263')
#P21397 P35268 Q9Y646 Q09666 Q13263
protMA <- c('P21397','P35268','Q9Y646','Q09666','Q13263','P37802','P12110')
sigProts <- unique(c(protNM,protAC,protCP,protMA))


#############################################################################################
zhe1 <- read.table("E:/projects/TPD/data/tpdzy_prot_matrix_20190122.csv",sep=",",header=T,stringsAsFactors=F)[,-1]
zheSamples <- sapply(colnames(zhe1)[-1],function(v){a=strsplit(v,"_with")[[1]][1];strsplit(a,"_DIA_")[[1]][2]})
zheProtId <- as.character(sapply(zhe1$prot,function(v){strsplit(v,"\\|")[[1]][2]}))
zhe1 <- data.frame(t(zhe1[,-1]))
colnames(zhe1) <- zheProtId
rownames(zhe1) <- zheSamples
zhe1 <- zhe1[sort(zheSamples),]
patientIds <- unique(as.character(sapply(rownames(zhe1),function(v){strsplit(v,"_")[[1]][1]})))

zhe0 <- c()
for(pid in patientIds){
  rname <- c(paste0(pid,"_repA"),paste0(pid,"_repB"))
  tmp <- zhe1[rname,]
  zhe0 <- rbind(zhe0,apply(tmp,2,function(v){v0 <- 2^v;log2(mean(v0,na.rm=T))}))
}
zhe0 <- data.frame(zhe0)
zhe0[is.na(zhe0)] <- NA
rownames(zhe0) <- patientIds
zheV <- data.frame(zhe0[,sigProts])
R1 <- apply(zheV,1,function(v){sum(is.na(v))})
validation0 <- zheV[R1==0,]
validation1 <- data.frame(t(apply(zheV[R1==0,],1,function(v){(v-mean(v))/sd(v)})))

################################## training dataset   #################################################
df0  <- read.table("E:/projects/TPD/data/TPD_prot_matrix_avg_20190115.txt",sep="\t",header=T,stringsAsFactors=F)
df0$label[df0$label=='W'] <- 'C'
t0 <- df0[,-1]
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
########################################################################
NM <- getSubset(c('N','M'),5)[,c('label',protNM)]
AC <- getSubset(c('A','C'),5)[,c('label',protAC)]
CP <- getSubset(c('C','P'),5)[,c('label',protCP)]
MA <- getSubset(c('M','A'),5)[,c('label',protMA)]
ACP <- getSubset(c('A','C','P'),10)[,c('label',unique(c(protAC,protCP)))]

M20 <- t0[,c('label',sigProts)]
M20[M20$label %in% c('N','M'),protNM] <- NM[,protNM]
M20[M20$label %in% c('A','C'),protAC] <- AC[,protAC]
M20[M20$label %in% c('C','P'),protCP] <- CP[,protCP]
M20[M20$label %in% c('M','A'),protMA] <- MA[,protMA]


NM$label <- as.factor(NM$label)
AC$label <- as.factor(AC$label)
CP$label <- as.factor(CP$label)
MA$label <- as.factor(MA$label)

if(F){
	M20[,-1] <- data.frame(t(apply(M20[,-1],1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
	M20[is.na(M20)] <- min(M20[,-1],na.rm=T)
	model <- randomForest(label ~ . ,data=M20,importance=T,ntree=1000,nodesize=7)

	modelNM <-  randomForest(label ~ . ,data=NM,importance=T,ntree=1000,nodesize=7)
	modelAC <-  randomForest(label ~ . ,data=AC,importance=T,ntree=1000,nodesize=7)
	modelCP <-  randomForest(label ~ . ,data=CP,importance=T,ntree=1000,nodesize=7)
	modelMA <-  randomForest(label ~ . ,data=MA,importance=T,ntree=1000,nodesize=7)

}

 BinaryScore<- function(V,M1,it=10,rowNormalization=T){
   features <- colnames(M1)
   features <- features[features!='label']
   testData <- V
   if(all(features %in% names(testData))) {
     testData <- testData[features]
     print(sprintf("%d instances to be predicted",dim(testData)[1]))
   }else{
       print("supplied test data missed some features")
       print(sprintf("Please supply the values of these proteins:%s as orderly",paste0(features,collapse=","))) 
       #return(NULL)
   }
   lbls <- unique(M1$label);
   sampleNumbers <- sapply(lbls,function(k){sum(M1$label==k)})
   names(sampleNumbers) <- lbls
   sampleSize <- min(sampleNumbers)
   sampleLabel <- names(sampleNumbers)[max(sampleNumbers)==sampleNumbers]
   L=0;
   ifelse(is.vector(testData),L <- 1,L <- dim(testData)[1])
   pred <- matrix(0,nrow=L,ncol=length(lbls))
   for(i in 1:it){
	    indx <- sample(which(M1$label==sampleLabel),sampleSize,replace=F)
	    indx <- c(indx,which(M1$label!=sampleLabel))
	    tmp  <- M1[indx,features]
	    if(rowNormalization){
	       tmp  <- data.frame(t(apply(tmp,1,function(v){(v-mean(v))/sd(v)})))
	       testData <- data.frame(t(apply(testData,1,function(v){(v-mean(v))/sd(v)})))
	    }
	    tmp$label <- as.factor(M1$label[indx])
	    md <- randomForest(label ~ . ,data=tmp,importance=T,ntree=1000,nodesize=7)
	    predicted <- as.matrix(predict(md,testData,type='prob'))
	    pred <- pred+predicted
           
   }
   pred <- pred/it
   ifelse(L==1,names(pred) <- names(predicted),colnames(pred) <- colnames(predicted))
   pred
}



predicTPD <- function(V,sampleTimes){
    


}

###################################################################
testAC0 <- validation0[!grepl('B',rownames(validation0)),]
testAC1 <- validation1[!grepl('B',rownames(validation1)),]
testCP0 <- validation0[!grepl('A',rownames(validation0)),]
testCP1 <- validation1[!grepl('A',rownames(validation1)),]

scoreAC0 <- BinaryScore(testAC0,AC)
scoreAC1 <- BinaryScore(testAC1,AC,rowNormalization = T)
scoreCP0 <- BinaryScore(testCP0,CP)
scoreCP1 <- BinaryScore(testCP1,CP,rowNormalization = T)

scores0 <- BinaryScore(validation0,AC,it=100,rowNormalization = F)
scores1 <- BinaryScore(validation1,AC,it=100,rowNormalization = T)


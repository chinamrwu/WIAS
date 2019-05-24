#用浙一医的54个样本验证两分类（良性、恶性）
#12个特征蛋白从火山图中挑选出来

remove(list=ls())
library(randomForest)
library(data.table)
library(caret)
remove(list=ls())
library(randomForest)
library(caret)
library(data.table)
library(pROC)
set.seed(2019227)
#################################################################################################################


#############################################################################################
zhe0 <- read.table("E:/projects/TPD/data/tpdzy_prot_20190127.csv",sep=",",header=T,stringsAsFactors=F)[,-1]
zheSamples <- sapply(colnames(zhe0)[-1],function(v){a=strsplit(v,"_with")[[1]][1];strsplit(a,"_DIA_")[[1]][2]})
zheProtId <- as.character(sapply(zhe0$prot,function(v){strsplit(v,"\\|")[[1]][2]}))
zhe0 <- data.frame(t(zhe0[,-1]))
colnames(zhe0) <- zheProtId
rownames(zhe0) <- zheSamples
zhe0 <- zhe0[sort(zheSamples),]
patientIds <- unique(as.character(sapply(rownames(zhe0),function(v){strsplit(v,"_")[[1]][1]})))

zhe1 <- c()
for(pid in patientIds){
  rname <- c(paste0(pid,"_repA"),paste0(pid,"_repB"))
  tmp <- zhe0[rname,]
  zhe1 <- rbind(zhe1,apply(tmp,2,function(v){v0 <- 2^v;log2(mean(v0,na.rm=T))}))
}
zhe1 <- data.frame(zhe1)
zhe1[is.na(zhe1)] <- NA
rownames(zhe1) <- patientIds


################################## training dataset   #################################################
df0 <- read.table('data/TPD_prot_matrix_avg_20190131.txt',sep='\t',header=T,stringsAsFactors=F)
df0$label[df0$label=='W'] <- 'C'
rownames(df0) <- df0$SpecimenID
df0 <- df0[,-1]

df1 <- df0
df1$label[df1$label %in% c('N','M','A')] <- 'B'
df1$label[df1$label %in% c('C','P')]     <- 'M'




###############################################################################
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

#############################################################################
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

###################################################################
sigProts=c('P04083','P60174','P49773','Q13151','P51149','P62861','P31946','P09651')
sigProts=c('P04083','Q14974','P60174','P31946')
sigProts <- c('P04083','P00338','Q99584','P40121','P78417','P31949','P16403','O60437','P04004','Q02880','P02765','Q08380')

sigProts1=c('P04083','Q14974','Q9Y3Y2','P49773','P62826')
sigProts2=c('P04083','Q14974','P60174','P31946','P62266')

sigProts=c('P04083','Q14974','P60174','P42677','P31946')

sigProts=c('P04083','Q14974','P60174','P49773','P05198')
testD <- zhe1[,sigProts]

testD <- data.frame(t(apply(testD,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
testD[is.na(testD)] <- 0

M1 <- df1[,sigProts]
M1 <- data.frame(t(apply(M1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
M1$label <- df1$label
M1[is.na(M1)] <- 0
result <- BinaryScore(testD,M1)

errows=c('A11','A16','A19','A3','A4','A5','A6','A7','C12','C15','C6','C8')
errD <- zhe1[errows,sigProts]
accD <- zhe1[setdiff(rownames(zhe1),errows),sigProts]


############################  ROC for all samples



RFScore1 <- function(M1,nFolds=10,itors=10,index=1){
   if(! 'label' %in% colnames(M1)){
     print('A column with name "label" is required for labeling the class of each sample')
     return(NULL)
   }

   clsses <- unique(M1$label)
   #rownames(tmp) <- sapply(1:dim(tmp)[1],function(v){paste0('R',v)})
   sampleNumber <- sapply(clsses,function(v){sum(M1$label==v)})
   names(sampleNumber) <- clsses
   maxLabel <- names(sampleNumber)[max(sampleNumber)==sampleNumber]
   minLabel <- names(sampleNumber)[min(sampleNumber)==sampleNumber]
   size <- min(sampleNumber)
   
   probs <- rep(0,max(sampleNumber))
   names(probs) <- rownames(M1)[M1$label==maxLabel]
   
   indx2 <- rownames(M1)[M1$label==minLabel]

   prediction <- matrix(0,nrow=dim(M1)[1],ncol=2)
   rownames(prediction) <- rownames(M1)

   while(sum(probs==0)>0){
      indx1 <- names(sort(probs))[1:size]
      probs[indx1] <- probs[indx1]+1
      
      tmp1 <- M1[c(indx1,indx2),]
      for(i in 1:itors){
          folds <- createFolds(tmp1$label,nFolds)
          for(fold in folds){
              valids <- tmp1[fold,]
	      rownames(valids) <- rownames(tmp1)[fold]
              trains <- tmp1[setdiff(1:dim(tmp1)[1],fold),]
	      trains$label <- as.factor(trains$label)
              tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
	      predicted <- predict(tmpRF,valids,type='prob')
	      prediction[rownames(predicted),] <- predicted+prediction[rownames(predicted),]
         }
      }
     }#while
      colnames(prediction) <- colnames(predicted)
      prediction <- data.frame(t(apply(prediction,1,function(v){v/sum(v)})))
      prediction$predicted <- as.character(apply(prediction,1,function(v){names(v)[v==max(v)]}))
      prediction$observed  <- M1$label
      feature <- colnames(M1)[colnames(M1) !='label']
      rtv <-sprintf("%d %s %4.3f",index,paste0(feature,collapse=" "),sum(prediction$observed==prediction$predicted)/dim(M1)[1]*100)
      print(rtv)
      prediction
      #rtv
}


sigProts <- c('P04083','P00338','Q99584','P40121','P78417','P31949','P16403','O60437','P04004','Q02880','P02765','Q08380')

prot01=c('P04083','Q14974','P60174','P31946','P62266')
prot02=c('P04083','Q14974','Q9Y3Y2','P49773','P62826')

M1 <- df1[,unique(c(prot01,prot02))]
M1 <- data.frame(t(apply(M1,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
M1$label <- df1$label
M1[is.na(M1)] <- 0

scores <- RFScore1(M1,itors=10)
predicts <- scores
ROC <- roc(ifelse(predicts$observed=="B", "B", "M"), as.numeric(predicts$B))
pdf("TPD_binary_ROC_12protein_20190228.pdf")
plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",	
  main="ROC with 12 selected proteins",legacy.axes = TRUE,print.auc.cex=1.5)
dev.off()

write.table(scores,file="TPD_binary_ROC_score_with_12proteins_20190228.txt",sep="\t",col.name=T,row.names=T,quote=F)
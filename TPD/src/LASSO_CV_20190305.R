rm(list=ls())
library(glmnet)
library(caret)
library(ggplot2)
set.seed(190305)
options("width"=500)
setwd("E:/projects/TPD")
matAll <- read.table('data/TPD_prot_matrix_avg_20190304.txt',sep='\t',header=T,stringsAsFactors=F)
table(matAll$label)
matAll$label[matAll$label=='W'] <- 'C'
table(matAll$label)
rownames(matAll) <- matAll$SpecimenID
matAll <- matAll[,-2]
benign <- c('N','M','A')
malign <- c('C','P')

allProts <- colnames(matAll)[colnames(matAll) != 'label']
matB <- matAll[matAll$label %in% benign,] ## matrix for benign
matM <- matAll[matAll$label %in% malign,] ## matrix for malign

lblColors <- c(B='#56B4E9',M='#CC79A7')

R0 <- apply(matAll[,-1],2,function(v){sum(is.na(v))/length(v)*100})


        
cmbB <- combn(1:5,2) #
cmbM <- combn(1:5,4)
L1   <- dim(cmbB)[2]
L2   <- dim(cmbM)[2]
predicts <- list()
selected <- list()

#rates <- seq(5,100,by=5)
rates <- c(5,10,70,90)
for(missingCutOff in rates){

   selections  <- data.frame('prot'=allProts,'coefSum'=rep(0,length(allProts)),'number'=rep(0,length(allProts)))
   rownames(selections) <- allProts
   predictions <- data.frame('B'=rep(0,dim(matAll)[1]),'M'=rep(0,dim(matAll)[1]));
   rownames(predictions) <- rownames(matAll)

  
  for(itor in 1:5){
     index =1
     benignFolds <- createFolds(matB$label,k=5)  ## 
     malignFolds <- createFolds(matM$label,k=5) ##
     
    for(i in 1:L1){ ## choose 4 folds from 10 folds of benign dataset to combine a balanced training dataset
	  benignT <- matB[as.integer(unlist(benignFolds[cmbB[,i]])),names(R0)[R0 <= missingCutOff]]
	  benignT$label <- rep('B',dim(benignT)[1])
	  for(j in 1:L2){
	     malT         <-  matM[as.integer(unlist(malignFolds[cmbM[,j]])),names(R0)[R0 <= missingCutOff]]
	     malT$label   <- rep('M',dim(malT)[1])
	     
	     training     <- rbind(benignT,malT)
             label        <- as.factor(training$label)
	     training     <- training[,colnames(training)!='label']

             selection0   <- colnames(training)
             
             
	     flg=T
             while(flg){
	        
		 tmp0 <- training[,selection0]
	         tmp0 <- t(apply(tmp0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
                 tmp0[is.na(tmp0)] <- 0

                 cvfit <- cv.glmnet(tmp0,label,family='binomial',alpha=0.5)
                 cf1se <- coef(cvfit, s = "lambda.1se")
                 selection1 <- rownames(cf1se)[cf1se[,1]!=0][-1]
                 print(sprintf('%d- %d- %d - %d original prots ,%d selected,%d intersected',missingCutOff,itor,index,length(selection0),
		 length(selection1),length(intersect(selection1,selection0))))
                
		 if(length(selection1)>= 4){ ## 4 proteins at least 
		    flg <- length(selection1) != length(selection0)
                    selection0 <- selection1
		 }else{
		  flg=F
		  #selection0 <- selection1
		 }
		 
             }#while
	     selections[selection0,'number']  <- selections[selection0,'number']+1
             selections[selection0,'coefSum'] <- selections[selection0,'coefSum']+abs(cf1se[selection0,1])
             cvfit <- cv.glmnet(tmp0,label,family='binomial',alpha=0.5)

	     validation   <- matAll[setdiff(rownames(matAll),rownames(training)),selection0]
	     validation   <- t(apply(validation,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
	     validation[is.na(validation)] <- 0
             pred  <- predict(cvfit,newx = validation,s = cvfit$lambda.1se,type="class")
	     predictions[rownames(pred)[pred[,1]=='B'],'B'] <- predictions[rownames(pred)[pred[,1]=='B'],'B']+1
             predictions[rownames(pred)[pred[,1]=='M'],'M'] <- predictions[rownames(pred)[pred[,1]=='M'],'M']+1
	     
	     index <- index+1
	  }
   }

   }
   selected[[length(selected)+1]] <- selections
   predicts[[length(predicts)+1]] <- predictions
}
names(selected) <- as.character(rates)
names(predicts) <- as.character(rates)


allLabels <- matAll$label
allLabels[matAll$label %in% benign] <- 'B'
allLabels[matAll$label %in% malign] <- 'M'

pds <- list()
for(obj in predicts){ 
   obj0 <-  data.frame(t(apply(obj,1,function(v){v/sum(v)})))
   obj0$predicted <- as.character(apply(obj,1,function(v){names(v)[v==max(v)]}))
   obj0$observed  <- as.character(allLabels)
   pds[[length(pds)+1]] <- obj0
}
names(pds) <- names(predicts)
accs <- sapply(pds,function(obj){sum(obj$predicted==obj$observed)/dim(obj)[1]*100})## Ô¤²â×¼È·ÂÊ
names(accs) <- names(predicts)

rocs <- list()
for( obj in pds){
  rocs[[length(rocs)+1]] <- roc(ifelse(obj$observed=="B", "B", "M"), as.numeric(obj$B))
}
names(rocs) <- names(predicts)
aucs <- sapply(rocs,function(v){v$auc})
names(aucs) <- names(predicts)

#rates <- c(seq(10,90,by=5),5,95,100)
for(i in 1:20){
   rate <- rates[i]
     pdf(sprintf("ROC_%dNA_20190307.pdf",rate),width=6,height=6);
     plot.roc(rocs[[i]],print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best", main=sprintf("ROC with %d%% missingRate",rate),legacy.axes = TRUE,print.auc.cex=1.3)
   dev.off()
}

##################################################################################################

overlaps <- rownames(selected[[1]])[selected[[1]]$number>0]
for(obj in selected){ overlaps <- intersect(overlaps,rownames(obj)[obj$number>0])}




# generating pseudo samples for sample size balance.
# when required sample size less than the parameter size,this function generates pseudo samples,each of which
# is defined by two points P1 and P2. P1 
# averaging two samples randomly picked from the one class,and adding a noise, which is determined by
# random selected two from the samples of all the adversarial classes,the final pseudo sample is defined
#
scaleRow <- function(M0){
  tmp <- M0;
  tmp[,colnames(tmp)!='label'] <- t(scale(t(tmp[,colnames(tmp)!='label'])))
  tmp[is.na(tmp)] <- 0
  tmp
}
pseudoMixSamples <- function(M0,mixDegree=0.1,size=500){
    #M0 <- matA;mixDegree=0.1;size=200
	 tb <- table(M0$label)
	 features <-  colnames(M0)[colnames(M0)!='label']
	 tmp      <-  M0[,features]
    
	 rtv <- c()
	 for(nm in names(tb)){
        index0 <- which(M0$label==nm)
		  index1 <- which(M0$label!=nm)
		  sampleNumber <- size-tb[nm]
		  if(sampleNumber > 0){ 
				  indexA <- rbind(sample(index0,sampleNumber,replace=T), sample(index0,sampleNumber,replace=T))
				  indexB <- rbind(sample(index1,sampleNumber,replace=T), sample(index1,sampleNumber,replace=T))
				  index  <- rbind(indexA,indexB)
				  t0 <- apply(index,2,function(idx){
					  p1 <- apply(tmp[idx[1:2],],2,function(v){mean(2^v,na.rm=T)})
					  p2 <- apply(tmp[idx[3:4],],2,function(v){mean(2^v,na.rm=T)})
					  log2(p1+mixDegree*(p2-p1))
				  })
				  t0 <- data.frame(t(t0))
				  rownames(t0) <- paste0(nm,'#',1:dim(t0)[1])
				  tmp1 <- rbind(t0,tmp[index0,])
				  tmp1$label <- nm
				  tmp1 <- tmp1[,c('label',features)]
				  rtv <- rbind(rtv,tmp1)
		   }else {
				index <- setdiff(index0,sample(index0,tb[nm]-size,replace=F))
				tmp1   <- tmp[index,]
				tmp1$label <- nm;
				tmp1 <- tmp1[,c('label',features)]
				rtv  <- rbind(rtv,tmp1)
		  }
	 }
	 return(rtv)
}

###########################################


glmFeatures <- function(TM,VM,VM2=NULL){
     	Y0 <- TM$label
		Y1 <- VM$label
		#print(levels(as.factor(Y1)))
      hitFit <- NULL
		fits   <- list()
		acc=0
		flg=T
		selected   <- c()
		selection0 <- colnames(TM)[colnames(TM) != 'label']
      while(flg){
		  if(length(selection0) >=3 ) {
				tmp0   <- scaleRow(TM[,selection0])
				cvfit  <- cv.glmnet(as.matrix(tmp0),Y0,family='multinomial',alpha=1,type.measure='class')
				cf     <- coef(cvfit,s='lambda.1se')
				selection1 <- rownames(cf)[cf[,1]!=0][-1]
            if(length(selection1) >=3){
					tmp1   <- scaleRow(VM[,selection1])
					model  <- cv.glmnet(as.matrix(scaleRow(TM[,selection1])),Y0,family='multinomial',alpha=1,type.measure='class')
					pred   <- data.frame(predict(model,newx = as.matrix(tmp1),s = model$lambda.1se,type='class'))
					acc <- as.numeric(sapply(levels(as.factor(Y1)),function(ch){sum(pred[which(Y1==ch),1]==ch)/sum(Y1==ch)}))
					acc <- round(100*c(sum(pred==Y1)/dim(VM)[1],acc),digits=3)
               
					PPV  <- sum(Y1=='B' & pred=='B')/sum(pred=='B')*100
					NPV  <- sum(Y1=='M' & pred=='M')/sum(pred=='M')*100
               strShow <- sprintf("%d selected:%4.3f %4.3f %4.3f-%4.3f %4.3f",length(selection1),acc[1],acc[2],acc[3],PPV,NPV)

					
					if(length(selection1) <=35){ print(strShow);selected <- rbind(selected,paste0(strShow,"|",paste0(selection1,collapse=",")))}
					flg <- length(selection1)!=length(selection0)
				}
				selection0 <- selection1
			}else{flg <- F}
	   }
		return(data.frame(selected,stringsAsFactors=F))
}



glmIterate <- function(TM){
     	Y0 <- as.factor(TM$label)
    if(F){
		hitFit <- NULL
		fits   <- list()
		flg=T
		selected   <- c()
		selection0 <- colnames(TM)[colnames(TM) != 'label']
      while(flg){
		  if(length(selection0) >=3){
				tmp0   <- scaleRow(TM[,selection0])
				cvfit  <- cv.glmnet(as.matrix(tmp0),Y0,family='binomial',alpha=1,type.measure='class')
				cf <- coef(cvfit,s='lambda.1se')
				selection1 <- rownames(cf)[cf[,1]!=0][-1]
				selected <- rbind(selected,paste0(selection1,collapse=","))
				flg <- length(selection1)!=length(selection0)
				selection0 <- selection1
            print(sprintf("%d proteins selected",length(selection0)))
			}else{flg <- F}
	   }
		
		  #return(data.frame(selected,stringsAsFactors=F))
		}

		###########################
		cvfit  <- cv.glmnet(as.matrix(TM[,colnames(TM) !='label'],Y0,family='binomial',alpha=1,type.measure='class'))
		cf <- coef(cvfit,s='lambda.1se')
		selection1 <- rownames(cf)[cf[,1]!=0][-1]
		return(paste0(paste0(selection1,collapse=",")))
}

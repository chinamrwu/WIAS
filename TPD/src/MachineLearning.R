weightAdd <- function(v1,v2,weight){
	v01 <- weight* v1
	v02 <- (1- weight) * v2
	if(weight>0.5){v02[is.na(v02)] <- 0}
	else {v01[is.na(v01)] <- 0 }
	return(v01+v02)
}

scaleRow <- function(M0){
  tmp <- M0;
  if(rowNorm){ tmp[,colnames(tmp)!='label'] <- t(scale(t(tmp[,colnames(tmp)!='label'])))}
  tmp[is.na(tmp)] <- 0
  tmp
}

PseudoSamples <- function(M1,weight=1,k=2,size=1000){
	if(!'label' %in% colnames(M1)){
		 print('label column is required');
		 return(NULL)
	}
	#print("beginning generating Pseudo samples......")
	sampleTable  <-  table(M1$label)
	label2       <-  names(sampleTable)
	rows01       <-  rownames(M1)[M1$label==label2[1]]
	rows02       <-  rownames(M1)[M1$label==label2[2]]

	f1    <- colnames(M1)[colnames(M1)!='label']
	d2    <- length(f1)
	indx1 <- d2 +1
	indx2 <- 2*d2
        
	mat1 <- M1[M1$label==label2[1],f1]
	mat2 <- M1[M1$label==label2[2],f1]

	index1 <- c()
	for(i in 1:k){
	   index1 <- rbind(index1,sample(1:dim(mat1)[1],size/2,replace=T))
	}
        
	index2 <- c()
	for(i in 1:k){
	   index2 <- rbind(index2,sample(1:dim(mat2)[1],size/2,replace=T))
	}
	r1  <- t(apply(index1,2,function(v){apply(mat1[v,],2,function(v1){mean(v1,na.rm=T)})}))
	r2  <- t(apply(index2,2,function(v){apply(mat2[v,],2,function(v1){mean(v1,na.rm=T)})}))
	
	s1  <-  cbind(r1,r2)
	k1  <-  data.frame(t(apply(s1,1,function(v){weightAdd(v[1:d2],v[indx1:indx2],weight)})))
	rownames(k1) <- 1:(size/2)
	k1$label <- label2[1]
	
	s1  <-  cbind(r2,r1)
	k2  <-  data.frame(t(apply(s1,1,function(v){weightAdd(v[1:d2],v[indx1:indx2],weight)})))
	rownames(k2) <- (size/2+1):size
	k2$label <- label2[2]

	
	kr <- rbind(k1,k2)
	kr <- kr[,c('label',f1)]
	#print(sprintf("Finished generating %d samples with weight %4.3f",size,weight))
	return(kr)
}


#################### generate mixture of real and Pseudo samples
#### the number of pseudo samples equals the difference between
## parameter size and real sampleNumber 
## size: sample size of each class
generateMixSamples <- function(M1,weight=1,k=2,size=1000){
	if(!'label' %in% colnames(M1)){
		 print('label column is required');
		 return(NULL)
	}
	#print("beginning generating Pseudo samples......")
	sampleTable  <-  table(M1$label)
	label2       <-  names(sampleTable)

	f1    <- colnames(M1)[colnames(M1)!='label']
	d2    <- length(f1)
	indx1 <- d2 +1
	indx2 <- 2*d2
        
	mat1 <- M1[M1$label==label2[1],f1]
	mat2 <- M1[M1$label==label2[2],f1]
	L1 <- size - dim(mat1)[1]
	L2 <- size - dim(mat2)[1]

	index1 <- c()
	index2 <- c()
	for(i in 1:k){
	   index1 <- rbind(index1,sample(1:dim(mat1)[1],L1,replace=T))
		index2 <- rbind(index2,sample(1:dim(mat2)[1],L1,replace=T))
	}
	r1  <- t(apply(index1,2,function(v){apply(mat1[v,],2,function(v1){mean(v1,na.rm=T)})}))
	r2  <- t(apply(index2,2,function(v){apply(mat2[v,],2,function(v1){mean(v1,na.rm=T)})}))
	s1  <-  cbind(r1,r2)
	k1  <-  data.frame(t(apply(s1,1,function(v){weightAdd(v[1:d2],v[indx1:indx2],weight)})))
	k1$label <- label2[1]
	k1[is.na(k1)] <- NA
	####
	
	index1 <- c()
	index2 <- c()
	for(i in 1:k){
	   index1 <- rbind(index1,sample(1:dim(mat2)[1],L2,replace=T))
		index2 <- rbind(index2,sample(1:dim(mat1)[1],L2,replace=T))
	}
	r1  <- t(apply(index1,2,function(v){apply(mat2[v,],2,function(v1){mean(v1,na.rm=T)})}))
	r2  <- t(apply(index2,2,function(v){apply(mat1[v,],2,function(v1){mean(v1,na.rm=T)})}))
	s2  <-  cbind(r1,r2)
	k2  <-  data.frame(t(apply(s2,1,function(v){weightAdd(v[1:d2],v[indx1:indx2],weight)})))
	k2$label <- label2[2]
	k2[is.na(k2)] <- NA

	kr <- rbind(k1,k2)
	kr <- kr[,c('label',f1)]
	kr <- rbind(kr,M1)
	return(kr)
}

Model <- function(M1,measureType='class'){
   mat <- scaleRow(M1[,colnames(M1)!='label'])
   cvfit <- cv.glmnet(as.matrix(mat),as.factor(M1$label),family='binomial',alpha=1,type.measure=measureType)
   cvfit
}
Test <- function(model,mat,predType='class'){
   mat <- scaleRow(M1[,colnames(M1)!='label'])
	pred <- NULL
	if(predType=='class'){
		pred <- data.frame(predict(model,newx = as.matrix(mat),s = model$lambda.1se,type=predType))
		pred$observed <- as.character(M1$label)
		colnames(pred) <- c("predicted","observed")
	}else if(predType=='response'){
      lbls <- levels(as.factor(M1$label))[c(2,1)]
	   pred  <-  data.frame(t(apply(predict(model,newx = as.matrix(mat),s = model$lambda.1se,type="response"),1,function(v){c(v,1-v)})))
	   colnames(pred) <- lbls
		pred$predicted <- apply(pred,1,function(v){names(v)[max(v)==v]})
      pred$observed <- M1$label
	}
	L <- dim(M1)[1]
	acc <- sum(pred$predicted==pred$observed)/L*100
        #print(sprintf(" %d proteins predicted £ºaccruacy=%4.3f",dim(mat)[2], acc))
        #return(sprintf("%s %4.3f",paste0(colnames(mat),collapse=","),acc))
	return(acc)
}
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
		  if(length(selection0) >=3){
				tmp0   <- scaleRow(TM[,selection0])
				cvfit  <- cv.glmnet(as.matrix(tmp0),Y0,family='binomial',alpha=1,type.measure='class')
				cf <- coef(cvfit,s='lambda.1se')
				selection1 <- rownames(cf)[cf[,1]!=0][-1]

            if(length(selection1) >=3){
					tmp1   <- scaleRow(VM[,selection1])
					model  <- cv.glmnet(as.matrix(scaleRow(TM[,selection1])),Y0,family='binomial',alpha=0,type.measure='class')
					pred   <- data.frame(predict(model,newx = as.matrix(tmp1),s = model$lambda.1se,type='class'))
					acc <- as.numeric(sapply(levels(as.factor(Y1)),function(ch){sum(pred[which(Y1==ch),1]==ch)/sum(Y1==ch)}))
					acc <- round(100*c(sum(pred==Y1)/dim(VM)[1],acc),digits=3)
               
					PPV  <- sum(Y1=='B' & pred=='B')/sum(pred=='B')*100
					NPV  <- sum(Y1=='M' & pred=='M')/sum(pred=='M')*100
               strShow <- sprintf("%d selected:%4.3f %4.3f %4.3f-%4.3f %4.3f",length(selection1),acc[1],acc[2],acc[3],PPV,NPV)

					if(!is.null(VM2)){
							Y2    <- VM2$label
							tmp2   <- scaleRow(VM2[,selection1])
							pred1   <- data.frame(predict(model,newx = as.matrix(tmp2),s = model$lambda.1se,type='class'))
							acc1 <- as.numeric(sapply(levels(as.factor(Y2)),function(ch){sum(pred1[which(Y2==ch),1]==ch)/sum(Y2==ch)}))
							acc1 <- round(100*c(sum(pred1==Y2)/dim(VM2)[1],acc1),digits=3)
							if(acc[1]>80 & acc[2]>80 & acc[3]>80 & acc[4]>80 & acc[5] >80 & acc[6] >80 & length(selection1)<=20){ 
								fits[[length(fits)+1]] <- model
								print("Congratrulation! A better model found")
							}else if(length(selection1)<=10 & (acc[3]>=90 | acc[4]>=90)){
                        fits[[length(fits)+1]] <- model
								print("Congratrulation! A schew model found")
							}
					 PPV  <- sum(Y2=='B' & pred1=='B')/sum(pred1=='B')*100
					 NPV  <- sum(Y2=='M' & pred1=='M')/sum(pred1=='M')*100
					 strShow <- paste0(strShow,"|",sprintf("%4.3f %4.3f %4.3f-%4.3f %4.3f",acc1[1],acc1[2],acc1[3],PPV,NPV))
					}
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


if(T){
	rm(list=ls())
	setwd("/work")
	library(glmnet)
	library(caret)
	
   source("src/common.R")
   set.seed(2019314)
        print("loading protein matrix.....")
	df0 <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
	rownames(df0) <- df0$SpecimenID
	df0 <- df0[,-2]

	df0$label[df0$label %in% c('N','M','A')]  <- 'B'
	df0$label[df0$label %in% c('C','P','W')]  <- 'M'
	R0   <- apply(df0,2,function(v){sum(is.na(v))/length(v)*100})
   
	A  <- read.table("data/batchA.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
	Atest <- df0[setdiff(rownames(df0),rownames(A)),]
	B <- read.table("data/batchB.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
	Btest <- df0[setdiff(rownames(df0),rownames(B)),]
	zhe1 <- read.table("data/zhe1_54samples_protMat.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
	zhe1$label <- sapply(rownames(zhe1),function(v){ch=substr(v,1,1);ifelse(ch=='A','B','M')})
	zhe1 <- zhe1[,c(dim(zhe1)[2],1:(dim(zhe1)[2]-1))]

	color3=c(A="red",B="black",C="blue")
	color2=c(M='red',B='blue')

	##weight is the importance of v1
	weightAdd <- function(v1,v2,weight){
		v01 <- weight* v1
		v02 <- (1- weight) * v2
		if(weight>0.5){v02[is.na(v02)] <- 0}
		else {v01[is.na(v01)] <- 0 }
		return(v01+v02)
	}
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

balanceSamples <- function(M1,k=2){
  tmp          <- M1[,colnames(M1)!='label']
  features     <- colnames(M1)[colnames(M1)!='label']
  lblTable     <- table(M1$label)
  smallLabel   <- names(lblTable)[lblTable==min(lblTable)]
  bigLabel     <- names(lblTable)[lblTable==max(lblTable)]
  smallSamples <- M1[M1$label==smallLabel,features]
  size <- abs(lblTable[1]-lblTable[2])
  index <- c()
  for(i in 1:k){  index <- rbind(index,sample(1:dim(smallSamples)[1],size,replace=T))}
  k0 <- data.frame(t(apply(index,2,function(v){apply(smallSamples[v,],2,function(v0){mean(v0,na.rm=T)})})))
  k0$label <- smallLabel
  k0 <- k0[,c('label',features)]
  rbind(M1,k0)
}

rowNorm <- T
scaleRow <- function(M0){
  tmp <- M0;
  if(rowNorm){ tmp[,colnames(tmp)!='label'] <- t(scale(t(tmp[,colnames(tmp)!='label'])))}
  tmp[is.na(tmp)] <- 0
  tmp
}

glmFeatures <- function(M1){
		Y= as.factor(M1$label)
		hitFit <- NULL
		fits <- list()
		mse=100
		flg=T
		selected <- c()
		selection0 <- colnames(M1)[colnames(M1) != 'label']
      while(flg){
				tmp <- scaleRow(M1[,selection0])
				cvfit <- cv.glmnet(as.matrix(tmp),Y,family='binomial',alpha=1,type.measure='class')
				fits[[length(fits)+1]] <- cvfit
				cf <- coef(cvfit, s = 'lambda.1se')
				selection1 <- rownames(cf)[cf[,1]!=0][-1]
				mse1 <- min(cvfit$cvm)
				print(sprintf("%d %4.3f",length(selection1),mse1))
				if(mse1 <= mse){
						mse    <-  mse1;
						hitFit <-  cvfit
				}
				if(length(selection1)>= 3){ ## 3 proteins at least
					selected <- c(selected,sprintf("%s %4.3f",paste0(selection1,collapse=","),mse1))
					flg <- length(selection1) != length(selection0) 
					if(!flg){
                    hitFit <-  cvfit
					     mse <- mse1
					}
					selection0 <- selection1
	         }else{  
					flg=F	 
				}
      }#flg
		cf <- coef(hitFit, s = 'lambda.1se')
		selected <- rownames(cf)[cf[,1]!=0][-1]
		print(sprintf("%d selected %4.3f",length(selected),mse))
		selected
}
############################################################################################################################################################

glmFeatures1 <- function(TM,VM1,VM2=NULL){
     	Y0 <- TM$label
		Y1 <- VM$label
		print(levels(as.factor(Y1)))
      hitFit <- NULL
		fits   <- list()
		acc=0
		flg=T
		selected   <- c()
		selection0 <- colnames(TM)[colnames(TM) != 'label']
      while(flg){
		  
		  if(length(selection0) >=3){
				tmp0   <- scaleRow(TM[,selection0])
				cvfit  <- cv.glmnet(as.matrix(tmp0),Y1,family='binomial',alpha=1,type.measure='class')
				cf <- coef(cvfit,s='lambda.1se')
				selection1 <- rownames(cf)[cf[,1]!=0][-1]

            if(length(selection1) >=2){
					tmp1   <- scaleRow(VM1[,selection1])
					model  <- cv.glmnet(as.matrix(scaleRow(TM[,selection1])),Y1,family='binomial',alpha=0,type.measure='class')
					pred   <- data.frame(predict(model,newx = as.matrix(tmp1),s = model$lambda.1se,type='class'))
					acc <- as.numeric(sapply(levels(as.factor(Y1)),function(ch){sum(pred[which(Y1==ch),1]==ch)/sum(Y1==ch)}))
					acc <- round(100*c(sum(pred==Y1)/dim(VM)[1],acc),digits=3)
				
					if(!is.null(VM2)){
							Y2    <- VM2$label
							tmp2   <- scaleRow(VM2[,selection1])
							pred1   <- data.frame(predict(model,newx = as.matrix(tmp2),s = model$lambda.1se,type='class'))
							acc1 <- as.numeric(sapply(levels(as.factor(Y2)),function(ch){sum(pred1[which(Y2==ch),1]==ch)/sum(Y2==ch)}))
							acc1 <- round(100*c(sum(pred1==Y2)/dim(VM2)[1],acc1),digits=3)
							acc <- c(acc,acc1)
							if(acc[1]>80 & acc[2]>80 & acc[3]>80 & acc[4]>80 & acc[5] >80 & acc[6] >80 & length(selection1)<=20){ 
								fits[[length(fits)+1]] <- model
								print("Congratrulation! A better model found")
							}else if(length(selection1)<=10 & (acc[3]>=90 | acc[4]>=90)){
                        fits[[length(fits)+1]] <- model
								print("Congratrulation! A better model found")
							}
					}
					
					print(sprintf("%3d selected:  %4.3f %4.3f %4.3f | %4.3f %4.3f %4.3f",length(selection1),acc[1],acc[2],acc[3],acc[4],acc[5],acc[6]))
					flg <- length(selection1)!=length(selection0)
				}
				selection0 <- selection1
			}else{flg <- F}
	   }

		
		return(fits)
}

###############################################################################################################################################
Model <- function(M1,measureType='class'){
   mat <- scaleRow(M1[,colnames(M1)!='label'])
   cvfit <- cv.glmnet(as.matrix(mat),as.factor(M1$label),family='binomial',alpha=0,type.measure=measureType)
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


##########################################

#######################################

overlap <- intersect(colnames(zhe1),colnames(A))
R0 <- apply(A[,overlap],2,function(v){sum(is.na(v))/length(v)*100})

k0 <- A[,overlap]
candidates <- list()
seed=0;
mR=0
Flg=T
while(Flg){
	str1=readLines("http://dispatcher:8080/seed/new",warn=FALSE)
	if(str1 != 'NONE'){
		a=strsplit(str1,"\\|")[[1]]
		mR=as.integer(a[2])
		seed=as.integer(a[1])
		set.seed(seed)
                print(sprintf("seed=%d,rate=%d",seed,mR))
		TM <- PseudoSamples(k0[,R0<=mR],weight=0.90,size=1000)
		VM <- PseudoSamples(Atest,weight=0.9)
		tmp <- glmFeatures1(TM,VM,Atest)
		obj <- list("model"=tmp,"rate"=mR,"seed"=seed)
		if(!is.null(tmp) & length(tmp)!=0){
		  candidates[[length(candidates)+1]] <- tmp
		}
	}else{
           Flg=F
        }	
}
if(length(candidates) >0){
   save(candidates,file=paste0("output/",seed,".Rdata")
}

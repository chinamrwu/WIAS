# generating pseudo samples for sample size balance.
# when required sample size less than the parameter size,this function generates pseudo samples,each of which
# is defined by two points P1 and P2. P1 
# averaging two samples randomly picked from the one class,and adding a noise, which is determined by
# random selected two from the samples of all the adversarial classes,the final pseudo sample is defined
#
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
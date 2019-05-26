generateSamples <- function(M0,mixDegree=0.1,size=500){
    M0 <- matA;mixDegree=0.1;size=200
	 tb <- table(M0$label)
	 tmp <- M0[,colnames(M0)!='label']

	 rtv <- c()
	 for(nm in names(tb)){
        index0 <- which(M0$label==nm)
		  index1 <- which(M0$label!=nm)

		  sampleNumber <- size-tb[nm]
        if(size==0){sampleNumber <- length(index0)}
        indexA <- rbind(sample(index0,sampleNumber,replace=T), sample(index0,sampleNumber,replace=T))
		  indexB <- rbind(sample(index1,sampleNumber,replace=T), sample(index1,sampleNumber,replace=T))
        index  <- rbind(indexA,indexB)
        t0 <- apply(index,2,function(idx){
           p1 <- apply(tmp[idx[1:2],],2,function(v){mean(2^v,na.rm=T)})
			  p2 <- apply(tmp[idx[3:4],],2,function(v){mean(2^v,na.rm=T)})
			  log2(p1+mixDegree*(p2-p1))
		  })
		  t0 <- data.frame(t(t0))
		  rownames(t0) <- paste0(nm,1:dim(t0)[1])
	}


}
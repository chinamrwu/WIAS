## continue TPD_ABbatch_20190309.R for test zhe1 54 patient 

zhe0 <- read.table("data/tpdzy_prot_20190127.csv",sep=",",header=T,stringsAsFactors=F)[,-1]
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

########################################################################
zheTest <- zhe1
zheTest$label <- sapply(rownames(zheTest),function(v){ifelse(substr(v,1,1)=='A','B','M')})

zhePredict1 <- list()
for(i in 1:length(allSelected)){ 
  model <- Model(rbind(A[,c('label',allSelected[[i]])],B[,c('label',allSelected[[i]])]))
  if(all(allSelected[[i]] %in% colnames(zheTest))){
    print(i)
    zhePredict1[[length(zhePredict1)+1]] <- Test(model,zheTest[,c('label',allSelected[[i]])],'response')
  }
}

########################



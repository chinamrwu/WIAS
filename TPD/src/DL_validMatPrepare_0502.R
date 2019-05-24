matDL <- read.table('data/prot_TPDDL_190502.txt',sep="\t",header=T,stringsAsFactors=F)[,-1]
matDL <- matDL[grepl('_HUMAN',matDL$prot),]
colnames(matDL)[-1] <- as.character(sapply(colnames(matDL)[-1],function(v){a <- strsplit(v,"_DIA_")[[1]][2];strsplit(a,"_with_")[[1]][1]}))
rownames(matDL) <- as.character(sapply(matDL$prot,function(v){strsplit(v,"\\|")[[1]][2]}))
protIds <- rownames(matDL)
matDL <- matDL[,-1]
matDL <- matDL[,!grepl('pool',colnames(matDL))]
matDL <- 2^matDL
matDL <- data.frame(t(matDL))

sampleIds <- colnames(matDL)[!grepl('repB',colnames(matDL))]

k0 <- sapply(sampleIds,function(v){
   rname <- c(v,paste0(v,'_repB'))
   apply(matDL[rname,],2,function(v0){log2(mean(v0,na.rm=T))})
})

k0[is.na(k0)] <- NA
k0 <- data.frame(t(k0))




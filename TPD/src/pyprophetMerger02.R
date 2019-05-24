remove(list=ls())
library(parallel)
print("Searching target files from pyprophet directory.........")
fnames <- list.files(path = 'output/pyprophet',pattern="*_with_dscore_filtered.csv")
print(sprintf("%d target files found!,beginning to merge them.....",length(fnames)))

coreNumber=12

outputFile <- 'output/TPD_peptide_matrix.txt'
parseSingle <- function(fname){
       df1 <-  read.table(paste0('output/pyprophet/',fname),sep="\t",header=T,stringsAsFactors=F) 
       df1 <- df1[df1$decoy==0 & df1$peak_group_rank==1,]
       df1 <- df1[df1$m_score <0.01,]
       df1 <- df1[grepl('1/sp',df1$ProteinName),c('peptide_group_label','ProteinName','Intensity','filename')]
       return(df1)
}
mergeTwo <- function(fname01,fname02){
    df01 <- parseSingle(fname01)
    df02 <- parseSingle(fname02)
    result <- merge(x=df01[,-4],y=df02[,-4],by=c('peptide_group_label','ProteinName'),all=T)
    clnames <- as.character(sapply(c(fname01,fname02),function(v){
         clname <- strsplit(v,"\\.")[[1]][1]
         clname <- strsplit(clname,"\\/")[[1]]
         clname <- clname[length(clname)]
	 clname
	 }))
    d <- dim(result)[2]
    colnames(result)[(d-1):d] <- clnames
    result
}
##############################################################################################################

count <- 0
index <- seq(2,length(fnames),by=(coreNumber+1))
if(length(fnames)>0){
   df0 <- parseSingle(fnames[1])[,-4]
   colnames(df0)[3] <- strsplit(fnames[1],"_with_")[[1]][1]
   for(i in 1:(length(index)-1)){
      mats <- mclapply(fnames[(index[i]+1):(index[i+1]-1)], parseSingle, mc.cores = coreNumber)
      L <- length(mats)
      for(tmp in mats){
         df0 <- merge(x=df0,y=tmp[,-4],by=c('peptide_group_label','ProteinName'),all=T)
	 clname <- strsplit(as.character(tmp[1,4]),"\\.")[[1]][1]
	 clname <- strsplit(clname,"\\/")[[1]]
	 clname <- clname[length(clname)]
	 colnames(df0)[dim(df0)[2]] <- clname
      }
      count <- count+L
      print(sprintf("%d files have been merged",count)) 
    }
}
dones <- as.character(sapply(colnames(df0)[-c(1:2)],function(v){paste0(v,'_with_dscore_filtered.csv')}))
fnames <- setdiff(fnames,dones)
print(sprintf("Now merging the last %d files",length(fnames)))
if(length(fnames)>0){
  mats <- mclapply(fnames, parseSingle, mc.cores = coreNumber)
  for(tmp in mats){
         df0 <- merge(x=df0,y=tmp[,-4],by=c('peptide_group_label','ProteinName'),all=T)
	 clname <- strsplit(as.character(tmp[1,4]),"\\.")[[1]][1]
	 clname <- strsplit(clname,"\\/")[[1]]
	 clname <- clname[length(clname)]
	 colnames(df0)[dim(df0)[2]] <- clname
  }
}

write.table(df0,file=outputFile,sep="\t",col.names=T,row.names=F,quote=F)
print(sprintf("%d unique peptides found and %d files are merged",dim(df0)[1],(dim(df0)[2]-2)))
print("Work is done! Your are welcome back! ")


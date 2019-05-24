print("Loading data from pyprophet directory.........")
fnames <- list.files(path = "pyprophet",pattern="*_with_dscore_filtered.csv")
df1 <- df1[df1$decoy==0 & df1$peak_group_rank==1,]
df1 <- df1[df1$m_score <0.01,]
df1 <- df1[grepl('1/sp',df1$ProteinName),c('peptide_group_label','ProteinName','Intensity')]
df2 <- df2[df2$decoy==0 & df2$peak_group_rank==1,]
df2 <- df2[df2$m_score <0.01,]
df2 <- df2[grepl('1/sp',df2$ProteinName),c('peptide_group_label','ProteinName','Intensity')]
df0 <- merge(x=df1,y=df2,by=c('peptide_group_label','ProteinName'),all=T)
indx <-2
print(sprintf("There are %d files to be merged......",length(fnames)))
for(fname in fnames[-c(1:2)]){
   tmp <- read.table(paste0("pyprophet/",fname),sep="\t",header=T,stringsAsFactors=F)
   tmp <- tmp[tmp$decoy==0 & tmp$peak_group_rank==1,]
   tmp <- tmp[tmp$m_score <0.01,]
   tmp <- tmp[grepl('1/sp',tmp$ProteinName),c('peptide_group_label','ProteinName','Intensity')]
   df0 <- merge(x=df0,y=tmp,by=c('peptide_group_label','ProteinName'),all=T)
   colnames(df0)[dim(df0)[2]]=fname
   print(sprintf("%d :%s merged",indx,fname))
   indx <- indx+1
}
write.table(df0,file="/work/output/TPD_peptide_matrix.txt",sep="\t",col.names=T,row.names=F,quote=F)


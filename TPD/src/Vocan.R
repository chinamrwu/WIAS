Volcano <- function(df1,outFile=NULL,thresholdFC=1.5,thresholdPValue=0.05){ 
   label <- df1$label
   lbls = unique(df1$label)
   table(df1$label) # A143,C75
   fc <- apply(2^df1[,colnames(df1)!='label'],2,function(x){
       log2( mean(na.omit(x[label==lbls[1]])) /mean(na.omit(x[label==lbls[2]])))
    })
   
   df1[is.na(df1)] <- 0
   pValue <- apply(df1[,colnames(df1)!='label'], 2, function(v) {
       p1 <- t.test(v[label == lbls[1]], v[label == lbls[2]], paired = F, var.equal = F)$p.value
       p1 <-  p.adjust(p1,method="BH")
       p1
    }) 
   name = list(up = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[up],fc = fc[up], p_value = pValue[up],
                            type = rep("Upregulation",sum(up)),stringsAsFactors=F),
            down = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[down],fc = fc[down], p_value = pValue[down],
                              type = rep("Downregulation",sum(down))),stringsAsFactors=F)
   name1 = rbind(name[[1]],name[[2]],stringsAsFactors=F)
   rownames(name1) <- 1:dim(name1)[1]
   name1 <- name1[order(abs(name1$fc),decreasing=T),]
   if( !is.null(outFile) ){
     write.table(name1,file=outFile,sep="\t",col.names=T,row.names=F,quote=F)
   }
   #write.xlsx(name1, "Table_1_TPD_volcano_prot_AC10_fc1.5_190131.xlsx",showNA = T,row.names = F)
   return(name1)
}
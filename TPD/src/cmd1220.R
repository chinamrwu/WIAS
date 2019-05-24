result=c()
for(ch in strsplit("NMACP","")[[1]]){
     tmp <- trainM[trainM$label==ch,-1]
     v1 <- apply(tmp,2,function(v){sum(is.na(v))})
     result <- cbind(result,v1)
     result <- cbind(result,v1*100/dim(tmp)[1])
     result <- cbind(result,dim(tmp)[1])
 }

result <- data.frame(result)
result$protId <- rownames(result)
result <- result[,c(dim(result)[2],1:(dim(result)[2]-1))]
write.csv(result,file="E:/projects/TPD/results/protein_missingValues_stats.csv",col.names=T,row.names=F,quote=F)
colnames(result) <- c("protId",sapply(strsplit('NMACP','')[[1]],function(ch){c
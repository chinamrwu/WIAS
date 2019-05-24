rm(list=ls())
setwd("E:/projects/TPD")

machines <- c('A','B','C')


all3 <- c()
allProts <- c()
for(m in machines){
   df0 <- read.table(sprintf("data/machine%s_20190311.txt",m),sep="\t",header=F,stringsAsFactors=F)
   avgs <- df0$V1[grepl("\\*",df0$V1)]
	avgs <- as.numeric(sapply(avgs,function(v){strsplit(strsplit(v," average acc ")[[1]][2]," with missing rate ")[[1]]}))
	avgs <- data.frame(t(matrix(avgs,nrow=2))[,c(2,1)])
	colnames(avgs) <- c("missingRate","avgAccuracy")
	avgs$machine <- m
	print(avgs)
   df0 <- df0$V1[!grepl("\\*",df0$V1)]
   df0 <- t(sapply(df0,function(v){strsplit(v," ")[[1]]}))
	rownames(df0) <- 1:dim(df0)[1]
	df0 <- data.frame(df0,stringsAsFactors=F)
   colnames(df0) <- c("prots","accuracy","missingRate")
   df0$accuracy <- as.numeric(df0$accuracy)
	df0$missingRate <- as.numeric(df0$missingRate)
	df0 <- df0[order(df0$accuracy,decreasing=T),]
   allFounds <- c()
	for(prot in df0$prots){allFounds <- unique(c(allFounds,as.character(sapply(prot,function(v){strsplit(v,",")[[1]]}))))}
	print(sprintf("Total %d proteins found in machine %s!",length(allFounds),m))
	print(paste0(allFounds,collapse=","))
   allProts <- unique(c(allProts,allFounds))

	df1 =sqldf("SELECT prots,min(accuracy) minAcc,avg(accuracy) avgAcc,max(accuracy) mxAcc,count(*) CNT from df0 group by prots order by avg(accuracy) desc")
   df1$protNumber <- as.integer(sapply(df1$prots,function(v){length(strsplit(v,",")[[1]])}))
	df1$machine <- rep(m,dim(df1)[1])
   all3 <- rbind(all3,df1)
}
print(sprintf("%d proteins found ",length(allProts)))
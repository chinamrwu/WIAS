
###### Analysis of the output of RF
remove(list=ls())
library("party")
library("sqldf")
require(pROC)

rawMat <- read.table("E:/projects/TPD/data/RF_TPDT_expert_trainingM_20181024.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(rawMat) <- rawMat[,1]

testData <- read.table("E:/projects/TPD/data/RF_TPDV_expert_20181024.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
protIds <- intersect(colnames(rawMat)[3:dim(rawMat)[2]],colnames(testData))
trainM <- rawMat[,c("label",protIds)]

missRate <- apply(trainM[,-1],2,function(v){100*sum(v==0)/dim(trainM)[1]})

trainM <- rawMat[,c("label",names(missRate)[missRate<=30])]

getTime <- function(){
  s1 <- as.character(Sys.time())
  a <- strsplit(s1,"-")[[1]]
  b <- strsplit(a[3],":")[[1]]
  c1 <- paste0(strsplit(b[1]," ")[[1]],collapse="_");
  c2 <- c(a[1],a[2],c1,b[2],b[3])
  paste0(c2,collapse="")
}

strTime <- getTime()
uniqLabels =unique(trainM$label)
k <- length(uniqLabels)
nms <- c()
for(i in 1:(k-1)){
   for(j in (i+1):k){
     nms <- c(nms,paste0(uniqLabels[i],uniqLabels[j]))
   }
}

dlt <- 1
bestRFs <- list()
betterRecord=list();
set.seed(1978)

for(lbl in nms){
        print(sprintf("******************************* %s **************************",lbl))
	a <- substr(lbl,1,1);
	b <- substr(lbl,2,2);
	tmp <- trainM[trainM$label==a | trainM$label==b,]
	fullRF <- cforest(formula=as.factor(tmp$label) ~ . ,data=tmp, control =cforest_unbiased(ntree=500,mtry=5))
        imps <- varimp(fullRF,conditional=T)

	imps  <- data.frame(importance(fullRF));
	impScore <- imps$MeanDecreaseAccuracy     # * imps$MeanDecreaseGini

	imps <- imps[order(impScore,decreasing=T),]
	orderedFeatures <- rownames(imps)

	indx <- 1
	bestScore <- 0
	bestForest <- NULL
	while ( indx < 100 ){
		tmpRF <- NULL
		score <- 0
		
		currentFeatures <- c("label",orderedFeatures[1:indx])
		tmp <- trainM[trainM$label==a | trainM$label==b,currentFeatures]

		tmpRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)

		predictions=as.data.frame(tmpRF$votes)
		clss <- colnames(predictions)
		predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
		predictions$observed <- as.character(sapply(rownames(predictions),function(rn){substr(rn,1,1)}))

		RF <- list()
		tmpROC <- roc(ifelse(predictions$observed==clss[1], clss[1], clss[2]), as.numeric(predictions[,clss[1]]))
		RF$ROC <- tmpROC
		RF$label <- lbl
		RF$tree <- tmpRF
		score <- tmpROC$auc

		if(score > bestScore ){
			bestScore <- score
			bestForest <- RF
			bestFeature <- paste0(orderedFeatures[1:indx],collapse="_")
			betterRecord[[length(betterRecord)+1]] <- paste0(bestFeature,"_",lbl,"_",score)
		        print(sprintf("find better forest for %s with score: %f and %d features",lbl,bestScore,indx))
			
		}
      
         indx <- indx+dlt
       }##while()
     bestRFs[[length(bestRFs)+1]] <- bestForest 
   }
    
   ##########################################  draw ROC plot
	#mt <- matrix(1:10,nrow=2,ncol=5)
	pdf(paste0("E:/projects/TPD/results/RF_TPD_expert_ROC_All_v3_",strTime,".pdf"),width=12)
	par(mfrow = c(2, 5),mar=c(0,0,0,0))

	for(obj in bestRFs){
		plot.roc(obj$ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",
		main=obj$label,legacy.axes = TRUE,print.auc.cex=1.2)
	}
	mtext("ROC Plots (expert <= 30% missing values)",outer=T,line=-1)
	dev.off()
   ################# importance proteins 
   protScore=c()
   for(obj in bestRFs){
     impM  <- data.frame(importance(obj$tree));
     impScore <- impM$MeanDecreaseAccuracy    #* impM$MeanDecreaseGini
     impM <- impM[order(impScore,decreasing=T),]
     impM$protId <- rownames(impM)
     impM$score  <- impScore[order(impScore,decreasing=T)]
     impM$label <- obj$label
     protScore <- rbind(protScore,impM[,c("protId","score","label")])
   }
   protScore <- data.frame(protScore,stringsAsFactors=F)
   rownames(protScore) <- 1:dim(protScore)[1]
   protScore$score <- as.numeric(protScore$score)
   geneMapping <- read.table("E:/projects/TPD/data/Gene_Prot_Mapping.txt",sep="\t",header=T,stringsAsFactor=F)
   protScore$symbol <- geneMapping$gene[match(protScore$protId,geneMapping$protId)]
   
   write.table(protScore,file=paste0("E:/projects/TPD/results/geneScores_expert_v3_",strTime,".txt"),sep="\t",quote=F,col.names=T,row.names=F)

   biomarks <- unique(protScore$protId)
   result <- data.frame(t(sapply(biomarks,function(pid){v=protScore$score[protScore$protId==pid];c(mean(v),sum(v))})))
   result$biomark <- biomarks
   colnames(result)<- c("avgScore","totall","protId")
   result <- result[,c("protId","avgScore","totall")]
   result$symbol <- geneMapping$gene[match(result$protId,geneMapping$protId)]
   result <- sqldf("SELECT * FROM result order by avgScore desc")

   write.table(result,file=paste0("E:/projects/TPD/results/biomarks_expert_v3_",strTime,"_.txt"),sep="\t",quote=F,col.names=T,row.names=F)

   ### betterRecord
   betterRecord <- as.character(sapply(betterRecord,function(v){a=strsplit(v,"_")[[1]];paste0(a,collapse="\t")}))
   write.table(betterRecord,file="E:/projects/TPD/results/Trace_expert.txt",col.names=F,quote=F,row.names=T)

   t1 <- c()
   k <- length(betterRecord)

   t1 <- as.character(sapply(betterRecord[[1]],function(v){a=strsplit(v,"_")[[1]];paste0(a,collapse="\t")})) 
   for(i in 2:k){
      v1 <- as.character(sapply(betterRecord[[k-1]],function(v){a=strsplit(v,"_")[[1]]}))
      v2 <- as.character(sapply(betterRecord[[k]]  ,function(v){a=strsplit(v,"_")[[1]]}))
      if(v1[length(v1)-1]==v2[length(v2)-1]){
            s1 <- setdiff(v2[1:(length(v2)-2)],v1[1:(length(v1)-2)]);
	    t1 <- rbind(t1,paste0(s1,collapse="\t"),"\t",v2[(length(v2)-1)],"\t",v2[length(v2)])
      }else{
       t1 <- rbind(t1,as.character(sapply(betterRecord[[k]],function(v){a=strsplit(v,"_")[[1]];paste0(a,collapse="\t")}))) 
      }
   }
   write.table(t1,file="E:/projects/TPD/results/Trace_expert_v3_simple.txt",col.names=F,quote=F,row.names=T,sep="\n")
 

###### Analysis of the output of RF
remove(list=ls())
library("randomForest")
library("sqldf")
require(pROC)

rawMat <- read.table("E:/projects/TPD/data/RF_TPDT_expert_trainingM_20181024.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)

rownames(rawMat) <- rawMat[,1]
trainM <- rawMat[,-1]

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

k=length(uniqLabels)
dlt <- 1
bestRFs <- list()


betterRecord=c();

v1 <- c(4,5,8,6,16,20,7,8,7,37)
names(v1) <- nms
set.seed(1978)
for(lbl in nms){
        print(sprintf("********************************************** %s *******************************",lbl))
	a <- substr(lbl,1,1);
	b <- substr(lbl,2,2);
	tmp <- trainM[trainM$label==a | trainM$label==b,]
	fullRF <- randomForest(formula=as.factor(tmp$label) ~ . ,data=tmp,importance=T,ntree=1000)

	imps  <- data.frame(importance(fullRF));
	impScore <- imps$MeanDecreaseAccuracy * imps$MeanDecreaseGini

	imps <- imps[order(impScore,decreasing=T),]
	orderedFeatures <- rownames(imps)

	indx <- 2
	bestScore <- 0
	bestForest <- NULL
	while ( indx < 100 ){
		tmpRF <- NULL
		score <- 0
		
		if(indx <= v1[lbl]){
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
				bestFeature <- currentFeatures
				betterRecord <- rbind(betterRecord,c(lbl,bestScore,length(bestFeature)))
				print(sprintf("find better tree for %s with score: %f and %d features",lbl,bestScore,indx))
			}
		}
      
         indx <- indx+dlt
       }##while()
     bestRFs[[length(bestRFs)+1]] <- bestForest 
   }

   ##########################################  draw ROC plot
pdf(paste0("E:/projects/TPD/results/RF_TPD_expert_ROC_All_v4_",strTime,".pdf"),width=12)
	par(mfrow = c(2, 5),mar=c(0,0,0,0))

	for(obj in bestRFs){
		plot.roc(obj$ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="best",
		main=obj$label,legacy.axes = TRUE,print.auc.cex=1.2)
	}
	mtext("ROC Plots (openSWATH <= 30% missing values)",outer=T,line=-1)
	dev.off()
   ################# importance proteins 

   protScore=c()
   for(obj in bestRFs){
     print(obj$label)
     impM  <- data.frame(importance(obj$tree));
     impScore <- impM$MeanDecreaseAccuracy * impM$MeanDecreaseGini
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
   geneMapping <- read.table("E:/projects/TPD/data/Gene_Prot_Mapping.txt",sep="\t",header=T,stringsAsFactor=F)
   protScore$symbol <- geneMapping$gene[match(protScore$protId,geneMapping$protId)]
   write.table(protScore,file=paste0("E:/projects/TPD/results/geneScores_expert_v4_",strTime,".txt"),sep="\t",quote=F,col.names=T,row.names=F)

   biomarks <- unique(protScore$protId)
   result <- data.frame(t(sapply(biomarks,function(pid){v=protScore$score[protScore$protId==pid];c(mean(v),sum(v))})))
   result$biomark <- biomarks
   colnames(result)<- c("avgScore","totall","protId")
   result$symbol <-  geneMapping$gene[match(result$protId,geneMapping$protId)]
   result <- result[,c("protId","avgScore","totall")]
   result <- sqldf("SELECT * FROM result order by avgScore desc")
   write.table(result,file=paste0("E:/projects/TPD/results/biomarks_expert_v4_",strTime,".txt"),sep="\t",quote=F,col.names=T,row.names=F)

   
   



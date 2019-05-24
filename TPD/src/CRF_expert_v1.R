
###### Analysis of the output of RF
remove(list=ls())
library("party")
library("sqldf")
require(pROC)
source("F:/OneDrive/Documents/src/NovelVIM.R")

rawMat <- read.table("E:/projects/TPD/data/RF_TPDT_expert_trainingM_20181024.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(rawMat) <- rawMat[,1]

testData <- read.table("E:/projects/TPD/data/RF_TPDV_expert_20181024.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
protIds <- intersect(colnames(rawMat)[3:dim(rawMat)[2]],colnames(testData))
trainM <- rawMat[,c("label",protIds)]

missRate <- apply(trainM[,-1],2,function(v){100*sum(v==0)/dim(trainM)[1]})

trainM <- rawMat[,c("label",names(missRate)[missRate<=30])]
trainM$label <- factor(trainM$label,ordered = TRUE,levels = c("N", "M", "A","C","P"))

getTime <- function(){
  s1 <- as.character(Sys.time())
  a <- strsplit(s1,"-")[[1]]
  b <- strsplit(a[3],":")[[1]]
  c1 <- paste0(strsplit(b[1]," ")[[1]],collapse="_");
  c2 <- c(a[1],a[2],c1,b[2],b[3])
  paste0(c2,collapse="")
}
strTime <- getTime()
set.seed(1978)

ordinalRF <- cforest(label ~ ., data = trainM,control = cforest_unbiased(ntree = 500))

ER_VI <- varimp(ordinalRF) # error rate based variable importance (standard measure)
ER_VI <- ER_VI[order(ER_VI,decreasing=T)]

RPS_VI <- varimpRPS(ordinalRF) # RPS-based variable importance (novel VIM)
RPS_VI <- RPS_VI[order(RPS_VI,decreasing=T)]

MAE_VI <- varimpMAE(ordinalRF) # MAE-based variable importance (novel VIM)
MAE_VI <- MAE_VI[order(MAE_VI,decreasing=T)]

MSE_VI <- varimpMSE(ordinalRF)
MSE_VI <- MSE_VI[order(MSE_VI,decreasing=T)]



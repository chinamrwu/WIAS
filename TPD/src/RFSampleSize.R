remove(list=ls())
library("randomForest")
library("caret")
library("ggplot2")
set.seed(20181119)
source("E:/projects/TPD/src/growRF.R")
options(width=350)

SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
rownames(SWT) <- SWT$patientId
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求

tradeOff <- 5
missingRate <- apply(SWT[,-1],2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
trainM <- SWT[,c('patientId',clnames)]
trainM$label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))
trainM <- trainM[,c("label",clnames)]

missingRate <- apply(SWV,2,function(v){100*sum(is.na(v))/length(v)})
clnames <- names(missingRate)[missingRate<tradeOff]
testM <- SWV[,clnames]
protIds <- intersect(colnames(trainM),colnames(testM))

trainM <- trainM[,c("label",protIds)]
#trainM[,-1] <- t(apply(trainM[,-1],1,function(v){v/mean(v,na.rm=T)}))

testM <- testM[,protIds]
#testM <- data.frame(t(apply(testM,1,function(v){v/mean(v,na.rm=T)})),check.names=F)
trainM[is.na(trainM)] <- 0
testM[is.na(testM)] <- 0
#############################################

LearningCurve <- function(label){
  tmp <- trainM[trainM$label %in% label,]
  tmp$label <- factor(tmp$label)
  learning_curve <- learing_curve_dat(tmp,proportion =(2:10)/10, method = "rf", metric = "ROC", outcome='label',verbose = TRUE,  
  tuneLength = 15, test_prop = 0.2,trControl=trainControl(summaryFunction = twoClassSummary, classProbs = TRUE))
  learning_curve
}


p1 <- LearningCurve(c('N','M'))
p2 <- LearningCurve(c('N','A'))
p3 <- LearningCurve(c('N','C'))
p4 <- LearningCurve(c('N','P'))
p5 <- LearningCurve(c('M','A'))
p6 <- LearningCurve(c('M','C'))
p7 <- LearningCurve(c('M','P'))
p8 <- LearningCurve(c('A','C'))
p9 <- LearningCurve(c('A','P'))
p10 <- LearningCurve(c('C','P'))

showCurve <- function(pfm,strTitle="Learning Curve"){
   ggplot(pfm[pfm$Data=='Testing',], aes(x = Training_Size, y = ROC, color = Data))+
   geom_smooth(method = loess, span = .8)+ scale_x_continuous(expand = c(0, 0))+ 
   scale_y_continuous(expand = c(0, 0))+
   ggtitle(strTitle)+
   scale_colour_manual(values=c("Testing"='aquamarine4'))+
   theme(legend.position="none")+
    theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank())

}



  pdf(file="E:/projects/TPD/results/TPD_SampleSizeEstimation_20181211.pdf",width=15,height=12)
    grid.newpage();
    pushViewport(viewport(layout = grid.layout(4, 3)));

    print(showCurve(p1,'NM'),vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(showCurve(p2,'NA'),vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(showCurve(p3,'NC'),vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
    print(showCurve(p4,'NP'),vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(showCurve(p5,'MA'),vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
    print(showCurve(p6,'MC'),vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
    print(showCurve(p7,'MP'),vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(showCurve(p8,'AC'),vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
    print(showCurve(p9,'AP'),vp=viewport(layout.pos.row = 3, layout.pos.col = 3))
    print(showCurve(p10,'CP'),vp=viewport(layout.pos.row =4, layout.pos.col = 1))
 dev.off()

curves <- list()
curves[[length(curves)+1]] <- p1
curves[[length(curves)+1]] <- p2
curves[[length(curves)+1]] <- p3
curves[[length(curves)+1]] <- p4
curves[[length(curves)+1]] <- p5
curves[[length(curves)+1]] <- p6
curves[[length(curves)+1]] <- p7
curves[[length(curves)+1]] <- p8
curves[[length(curves)+1]] <- p9
curves[[length(curves)+1]] <- p10






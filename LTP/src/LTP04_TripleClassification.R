rm(list=ls())
library(data.table)
library(glmnet)
library(sqldf)

setwd('F:/WIAS/LTP')
source('src/LTP03_MachineLearning.R')
source('src/common.R')
print('Loading data ......')

matA <- read.table('data/matA.txt',sep="\t",header=T,stringsAsFactors=F)
matB <- read.table('data/matB.txt',sep="\t",header=T,stringsAsFactors=F)

features <- intersect(colnames(matA),colnames(matB))
matA <- matA[,features]
matB <- matB[,features]


benign <- c('B','Q','P')
malign <- c('C','S','M','G')
normal <- c('N')

color3 <- c('B'='blue','M'='red','N'='green')

tmpA <- matA;
tmpA$label[tmpA$label %in% benign] <- 'B'
tmpA$label[tmpA$label %in% malign] <- 'M'
tmpA$label[tmpA$label %in% normal] <- 'N'

set.seed(1)
trainingSet <- pseudoMixSamples(tmpA,mixDegree=0.1,size=500)

tmpB <- matB;
tmpB$label[tmpB$label %in% benign] <- 'B'
tmpB$label[tmpB$label %in% malign] <- 'M'
tmpB$label[tmpB$label %in% normal] <- 'N'

set.seed(2)
validSet <- pseudoMixSamples(tmpB,mixDegree=0.30,size=200) 

Y0 <- as.factor(trainingSet$label)


iterateLASSO <- function(DM0){
		Y1 <- as.factor(DM0$label)
		sf1 <- colnames(DM0)[colnames(DM0) != 'label']
		sf2 <- c()
		flg <- T
		result <- c()

		while(flg){
		      sf2 <- sf1
				tmp <- scaleRow(DM0[,c('label',sf2)])
				matT <- as.matrix(tmp[,colnames(tmp)!='label'])
				cvfit  <- cv.glmnet(matT,Y1,family='multinomial',alpha=1,type.measure='class')
				cf     <- coef(cvfit,s='lambda.1se')
				sf2 <- c()
				for(obj in cf){  sf2 <- c(sf2,rownames(obj)[obj[,1]!=0][-1]);}
				sf2 <- unique(sf2)
				print(length(sf2))
				flg <- length(sf1)!=length(sf2)
				sf1 <- sf2
		}
		paste0(sf1,collapse=";")
}

sf <- iterateLASSO(trainingSet)
sf <- strsplit(sf,";")[[1]]
p1 <- drawUMAP(tmpB[,c('label',sf)],color3,rowNormalization=T)


plots <- list()
R0 <- apply(trainingSet,2,function(v){sum(is.na(v))/length(v)*100})
if(F){
	for(miss in seq(100,5,by= -5)){
		trainDS <- trainingSet[,R0 <= miss]
      sf <- iterateLASSO(trainDS);
		sf <- strsplit(sf,";")[[1]]
      plots[[length(plots)+1]] <- drawUMAP(tmpB[,c('label',sf)],color3,rowNormalization=T,strTitle=sprintf("%d %d proteins",miss,length(sf)))
	}
}
rm(list=ls())
library(data.table)
library(glmnet)
library(sqldf)

setwd('E:/WIAS/LTP')
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

trainingSet <- pseudoMixSamples(tmpA,mixDegree=0.1,size=500)

tmpB <- matB;
tmpB$label[tmpB$label %in% benign] <- 'B'
tmpB$label[tmpB$label %in% malign] <- 'M'
tmpB$label[tmpB$label %in% normal] <- 'N'

validSet <- pseudoMixSamples(tmpB,mixDegree=0.30,size=200) 

Y0 <- as.factor(trainingSet$label)

tmp0 <- scaleRow(trainingSet)
matT <- as.matrix(tmp0[,colnames(tmp0)!='label'])
cvfit  <- cv.glmnet(matT,Y0,family='multinomial',alpha=1,type.measure='class')
cf     <- coef(cvfit,s='lambda.min')
sf <- c()
for(obj in cf){  sf <- c(sf,rownames(obj)[obj[,1]!=0][-1]);}
sf <- unique(sf)
p1 <- drawUMAP(tmpB[,c('label',unique(sf))],color3,rowNormalization=T)



tmp1 <- scaleRow(trainingSet[,c('label',sf)])
matT <- as.matrix(tmp1[,colnames(tmp1)!='label'])
cvfit1  <- cv.glmnet(matT,Y0,family='multinomial',alpha=1,type.measure='class')
cf1     <- coef(cvfit1,s='lambda.1se')
sf1 <- c()
for(obj in cf1){  sf1 <- c(sf1,rownames(obj)[obj[,1]!=0][-1]);}
sf1 <- unique(sf1)

predict1 <- predict(cvfit1,newx=matB[,sf1],

p2 <- drawUMAP(tmpB[,c('label',unique(sf1))],color3,rowNormalization=T)


tmp2 <- scaleRow(trainingSet[,c('label',sf1)])
matT <- as.matrix(tmp2[,colnames(tmp2)!='label'])
cvfit2  <- cv.glmnet(matT,Y0,family='multinomial',alpha=1,type.measure='class')
cf2     <- coef(cvfit2,s='lambda.1se')
sf2 <- c()
for(obj in cf2){  sf2 <- c(sf2,rownames(obj)[obj[,1]!=0][-1]);}
sf2 <- unique(sf2)
p3 <- drawUMAP(tmpB[,c('label',unique(sf2))],color3,rowNormalization=T)

tmp3 <- scaleRow(trainingSet[,c('label',sf2)])
matT <- as.matrix(tmp3[,colnames(tmp3)!='label'])
cvfit3  <- cv.glmnet(matT,Y0,family='multinomial',alpha=1,type.measure='class')
cf3     <- coef(cvfit3,s='lambda.1se')
sf3 <- c()
for(obj in cf3){  sf3 <- c(sf3,rownames(obj)[obj[,1]!=0][-1]);}
sf3 <- unique(sf3)
p4 <- drawUMAP(tmpB[,c('label',unique(sf3))],color3,rowNormalization=T)

tmp4 <- scaleRow(trainingSet[,c('label',sf3)])
matT <- as.matrix(tmp4[,colnames(tmp4)!='label'])
cvfit4  <- cv.glmnet(matT,Y0,family='multinomial',alpha=1,type.measure='class')
cf4     <- coef(cvfit4,s='lambda.1se')
sf4 <- c()
for(obj in cf4){  sf4 <- c(sf4,rownames(obj)[obj[,1]!=0][-1]);}
sf4 <- unique(sf4)
p5 <- drawUMAP(tmpB[,c('label',unique(sf4))],color3,rowNormalization=T)

tmp5 <- scaleRow(trainingSet[,c('label',sf4)])
matT <- as.matrix(tmp5[,colnames(tmp5)!='label'])
cvfit5  <- cv.glmnet(matT,Y0,family='multinomial',alpha=1,type.measure='class')
cf5     <- coef(cvfit5,s='lambda.1se')
sf5 <- c()
for(obj in cf5){  sf5 <- c(sf5,rownames(obj)[obj[,1]!=0][-1]);}
sf5 <- unique(sf5)
p6 <- drawUMAP(tmpB[,c('label',unique(sf5))],color3,rowNormalization=T)

iterateLASSO <- function(DM0){
		Y1 <- as.factor(DM0$label)
		sf1 <- colnames(DM0)[colnames(DM0) != 'label']
		sf2 <- c()
		flg <- T
		while(flg){
		      sf2 <- sf1
				tmp <- scaleRow(DM0[,c('label',sf2)])
				matT <- as.matrix(tmp[,colnames(tmp)!='label'])
				cvfit  <- cv.glmnet(matT,Y1,family='multinomial',alpha=1,type.measure='class')
				cf     <- coef(cvfit,s='lambda.1se')
				sf2 <- c()
				for(obj in cf){  sf2 <- c(sf2,rownames(obj)[obj[,1]!=0][-1]);}
				sf2 <- unique(sf2)
				flg <- length(sf1)==length(sf2)
				sf1 <- sf2
		}

}
R0 <- apply(trainingSet,2,function(v){sum(is.na(v))/length(v)*100})
if(F){
for(miss in seq(100,5,by= -5)){
   trainDS <- trainingSet[,R0 <= miss]




}
}
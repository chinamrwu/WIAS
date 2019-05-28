# This source file tests  pseudoMixSamples function

rm(list=ls())
library(data.table)
library(sqldf)
library(umap)
source('src/common.R')
source('src/LTP03_MachineLearning.R')

###
matA <- read.table('data/matA.txt',sep="\t",header=T,stringsAsFactors=F)
matB <- read.table('data/matB.txt',sep="\t",header=T,stringsAsFactors=F)

pseudoA <- pseudoMixSamples(matA,size=30)
pseudoA[is.na(pseudoA)] <- NA
color2 <- c('R'='blue','P'='red','other'='black')

plots <- list();
for(dgr in seq(0.05,0.5,by=0.05)){
   set.seed(1)
	pseudoB <- pseudoMixSamples(matB,mixDegree=dgr,size=60)
	pseudoB[is.na(pseudoB)] <- NA
	tmp1 <- pseudoB[pseudoB$label=='B',]
	tmp1$label[ grepl('#',rownames(tmp1))]  <- 'P'
	tmp1$label[!grepl('#',rownames(tmp1))]  <- 'R'
   plots[[length(plots)+1]] <- drawUMAP(tmp1,color2,strTitle=sprintf("%4.3f",dgr),rowNormalization=T)
}

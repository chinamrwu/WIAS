rm(list=ls())
library(data.table)
library(glmnet)
library(sqldf)

setwd('F:/WIAS/LTP')
source('src/LTP03_MachineLearning.R')
print('Loading data ......')

matA <- read.table('data/matA.txt',sep="\t",header=T,stringsAsFactors=F)
matB <- read.table('data/matB.txt',sep="\t",header=T,stringsAsFactors=F)

benign <- c('B','Q','P')
malign <- c('C','S','M','G')
normal <- c('N')

color3 <- c('B'='blue','M'='red','N'='green')

tmp <- matA;
tmp$label[tmp$label %in% benign] <- 'B'
tmp$label[tmp$label %in% malign] <- 'M'
tmp$label[tmp$label %in% normal] <- 'N'

trainingSet <- generateSamples(tmp)

tmp <- matB;
tmp$label[tmp$label %in% benign] <- 'B'
tmp$label[tmp$label %in% malign] <- 'M'
tmp$label[tmp$label %in% normal] <- 'N'

validSet <- generateSamples(tmp,mixDegree=0.25,size=500) 





rm(list=ls())
library(ggplot2)
library(caret)
library(umap)
library(tsne)
library(glmnet)
library(sqldf)
library(caret)
setwd("E:/projects/TPD")
source("src/common.R")
source("src/MachineLearning.R")
set.seed(190311)
#############################################################################################################
color3=c(A="red",B="black",C="blue")
color2=c(M='red',B='blue')

#################################################################################
print("Loading datasets ......")


A  <- read.table("data/A0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
B  <- read.table("data/B0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
df0 <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(df0) <- df0$SpecimenID
df0 <- df0[,-2]
print("distribution A:");print(table(A$label))
print("distribution B:");print(table(B$label))

A$label[A$label %in% c('N','M','A')]  <- 'B'
A$label[A$label %in% c('C','P','W')]  <- 'M'
B$label[B$label %in% c('N','M','A')]  <- 'B'
B$label[B$label %in% c('C','P','W')]  <- 'M'
df0$label[df0$label %in% c('N','M','A')]  <- 'B'
df0$label[df0$label %in% c('C','P','W')]  <- 'M'

print("distribution A:");print(table(A$label))
print("distribution B:");print(table(B$label))

zhe1 <- read.table("data/zhe1_54samples_protMat.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
zhe1$label <- sapply(rownames(zhe1),function(v){ch=substr(v,1,1);ifelse(ch=='A','B','M')})
zhe1 <- zhe1[,c(dim(zhe1)[2],1:(dim(zhe1)[2]-1))]

TestA <- df0[setdiff(rownames(df0),rownames(A)),]
TestB <- df0[setdiff(rownames(df0),rownames(B)),]
rowNorm <- T

missingA <- apply(A,2,function(v){sum(is.na(v))/length(v)*100})
missingB <- apply(B,2,function(v){sum(is.na(v))/length(v)*100})

##############################################################################################################
seed=1
models <- list()
for(rate in seq(5,20,by=5)){
  for(i in 1:2){
     set.seed(seed)
	  print(sprintf("seed:%d, rate:%d",seed,rate))
	  seed <- seed+1
	  trainA <- PseudoSamples(A[,missingA <= rate])
     modelA <- glmFeatures(trainA,TestA)
	  modelA$seed <- seed
	  modelA$rate <- rate
	  modelA$src <- 'A'
     models[[length(models)+1]] <- modelA
	  print("------------------------------------------------")
	  trainB <- PseudoSamples(B[,missingB <= rate])
     modelB <- glmFeatures(trainB,TestB)
	  modelB$seed <- seed
	  modelB$rate <- rate
	  modelB$src <- 'B'
     models[[length(models)+1]] <- modelB
     print(sprintf("---------------------------------------------------------- %d ---- %d --",rate,i))

	}
}


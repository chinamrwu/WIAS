rm(list=ls())
library(ggplot2)
library(caret)
library(umap)
library(tsne)
library(glmnet)
library(sqldf)
library(caret)
#setwd("E:/projects/TPD")
setwd("/work")
source("src/common.R")
source("src/MachineLearning.R")
set.seed(190311)
#############################################################################################################
color3=c(A="red",B="black",C="blue")
color2=c(M='red',B='blue')

#################################################################################
print("Loading datasets ......")


A0  <- read.table("data/A0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
B0  <- read.table("data/B0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
C0  <- read.table("data/C0.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
df0 <- read.table("data/TPD_prot_matrix_avg_20190304.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(df0) <- df0$SpecimenID
df0 <- df0[,-2]
print("distribution A:");print(table(A0$label))
print("distribution B:");print(table(B0$label))
print("distribution C:");print(table(C0$label))

A0$label[A0$label %in% c('N','M','A')]  <- 'B'
A0$label[A0$label %in% c('C','P','W')]  <- 'M'
B0$label[B0$label %in% c('N','M','A')]  <- 'B'
B0$label[B0$label %in% c('C','P','W')]  <- 'M'
C0$label[C0$label %in% c('N','M','A')]  <- 'B'
C0$label[C0$label %in% c('C','P','W')]  <- 'M'

df0$label[df0$label %in% c('N','M','A')]  <- 'B'
df0$label[df0$label %in% c('C','P','W')]  <- 'M'

print("distribution A:");print(table(A0$label))
print("distribution B:");print(table(B0$label))
print("distribution C:");print(table(B0$label))


zhe1 <- read.table("data/zhe1_54samples_protMat.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
zhe1$label <- sapply(rownames(zhe1),function(v){ch=substr(v,1,1);ifelse(ch=='A','B','M')})
zhe1 <- zhe1[,c(dim(zhe1)[2],1:(dim(zhe1)[2]-1))]

TestA <- df0[setdiff(rownames(df0),rownames(A0)),]
TestB <- df0[setdiff(rownames(df0),rownames(B0)),]
TestC <- df0[setdiff(rownames(df0),rownames(C0)),]

rowNorm <- F

missingA <- apply(A0,2,function(v){sum(is.na(v))/length(v)*100})
missingB <- apply(B0,2,function(v){sum(is.na(v))/length(v)*100})
missingC <- apply(C0,2,function(v){sum(is.na(v))/length(v)*100})
R0       <- apply(df0,2,function(v){sum(is.na(v))/length(v)*100})
##############################################################################################################
if(F){  ## 
	seed=1
	models <- list()
	for(rate in seq(5,20,by=5)){
	  for(i in 1:100){
		  set.seed(seed)
		  print(sprintf("seed:%d, rate:%d",seed,rate))
		  seed <- seed+1
		  trainA <- PseudoSamples(A0[,missingA <= rate])
		  modelA <- glmFeatures(trainA,TestA)
		  if(dim(modelA)[1]>0){
				modelA$seed <- seed
				modelA$rate <- rate
				modelA$src <- 'A'
				models[[length(models)+1]] <- modelA
		  }
		  print("------------------------------------------------")
		  set.seed(seed)
		  trainB <- PseudoSamples(B0[,missingB <= rate])
		  modelB <- glmFeatures(trainB,TestB)
		  if(dim(modelB)[1]>0){
			  modelB$seed <- seed
			  modelB$rate <- rate
			  modelB$src <- 'B'
			  models[[length(models)+1]] <- modelB
		  }
		  print("------------------------------------------------")
			set.seed(seed)
			trainC <- PseudoSamples(C0[,missingB <= rate])
			modelC <- glmFeatures(trainC,TestC)
			if(dim(modelC)[1]>0){
				modelC$seed <- seed
				modelC$rate <- rate
				modelC$src <- 'C'
				models[[length(models)+1]] <- modelC
			}
		  print(sprintf("---------------------------------------------------------- %d ---- %d --",rate,i))
		}
	}
}
if(T){
   models <- list()
	seed=0;
	rate=5
	Flg=T
	mixWeight=0.95
	while(Flg){
	        str1=readLines("http://dispatcher:8080/task/new",warn=FALSE)
			  if(str1!='NONE'){
				  rate=as.integer(str1)
					seed=as.integer(readLines("http://dispatcher:8080/seed/new",warn=FALSE))
					set.seed(seed)
					print(sprintf("rate=%d",rate))
					trainA <- generateMixSamples(A0[,R0 <= rate],weight=mixWeight)
					modelA <- glmIterate(trainA)
					if(dim(modelA)[1]>0){
							modelA$seed <- seed
							modelA$rate <- rate
							modelA$src <- 'A'
							models[[length(models)+1]] <- modelA
					}
					#print("-----------------------------------------")

					seed=as.integer(readLines("http://dispatcher:8080/seed/new",warn=FALSE))
					set.seed(seed)
					trainB <- generateMixSamples(B0[,R0 <= rate],weight=mixWeight)
					modelB <- glmIterate(trainB)
					if(dim(modelB)[1]>0){
							modelB$seed <- seed
							modelB$rate <- rate
							modelB$src <- 'B'
							models[[length(models)+1]] <- modelB
					}
					#print("-----------------------------------------")
					
					seed=as.integer(readLines("http://dispatcher:8080/seed/new",warn=FALSE))
					set.seed(seed)
					trainC <- generateMixSamples(C0[,R0 <= rate],weight=mixWeight)
					modelC <- glmIterate(trainC)
					if(dim(modelC)[1]>0){
							modelC$seed <- seed
							modelC$rate <- rate
							modelC$src <- 'C'
							models[[length(models)+1]] <- modelC
					}
					print("------------------------------------------------------------------------------------------")
			 }else{
				  Flg=F
			 }
	}
	save(models,file=sprintf("/work/output/models_%d_0325.Rdata",seed))
}

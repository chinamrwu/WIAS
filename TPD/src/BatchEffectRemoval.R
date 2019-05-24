rm(list=ls())
library(ggplot2)
library(caret)
library(umap)

setwd("E:/projects/TPD")
source("src/common.R")
set.seed(190311)
pool <- read.csv('data/TPD_SG579_116poolProt_matrix_190304.csv',header=T,stringsAsFactors = F)

pool <- data.frame(t(pool[,-1]))
colnames(pool) <- prots
rownames(pool) <- as.character(sapply(rownames(pool),function(v){a=strsplit(v,"sunyt_TPD_DIA_")[[1]];b=strsplit(a[2],"_pool")[[1]][1];paste0(a[1],"_",b)}))
pool$MS <- machines
pool$label <- machines
pool <- pool[,c('label',prots)]
##########################################################################################################
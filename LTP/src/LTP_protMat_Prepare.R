rm(list=ls())
library(data.table)


setwd('E:/WIAS/LTP')
df0 <- read.csv("data/LTP_protMat_20190413.csv",header = T,stringsAsFactors = F)
tmp <- data.frame(t(df0[,-c(1,2)]))
colnames(tmp) <- df0$prot

sampleId <- as.character(sapply(rownames(tmp),function(v){ a <- strsplit(v,'_with_dscore_')[[1]][1];strsplit(a,'LTP_DIA_')[[1]][2]	}))

rownames(tmp) <- sampleId

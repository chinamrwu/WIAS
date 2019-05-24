remove(list=ls())
library(openxlsx)
library(sqldf)

df0 <- read.xlsx('E:/projects/TPD/data/AllSampleInformation.xlsx',sheet =1,colName=T)
df1 <- read.xlsx('E:/projects/TPD/data/AllSampleInformation.xlsx',sheet =2,colName=T)

df0 <- df0[!(is.na(df0[,1]) | is.na(df0[,2]) | is.na(df0[,3])),]
df0 <- df0[!(grepl('ouse',df0$SpecimenID) | grepl('ool',df0$SpecimenID)),]
df0 <- df0[!(grepl('ool',df0$ID2) | grepl('ous',df0$ID2)),]
#df0 <- df0[!grepl('e',df0$ID2),]

df0$patientId <- sapply(df0$ID2,function(v){substr(v,1,(nchar(v)-1))})
df0$label <- df0$patientId
df0$label[grepl('mC',df0$label)] <- 'C'
df0$label[grepl('wC',df0$label)] <- 'W'
df0$label <- sapply(df0$label,function(v){substr(v,1,1)})
df0$SpecimenID[which(df0$patientId=='P87')] <- '15:PF147_P'
df0$SpecimenID[which(df0$patientId=='A84')] <- '15:PF147_A'
df0 <- sqldf("SELECT * FROM df0 order by groups,patientId")
write.table(df0,file="E:/projects/TPD/data/AllSampleInformation_20190115.txt",sep="\t",col.names = T,row.names = F,quote = F)

     
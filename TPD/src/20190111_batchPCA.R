remove(list=ls())
library("data.table")
library("R.matlab")
library("mice")
options(width=350)
set.seed(20181224)

drawPCA <- function(df0,strTitle){ ## M is a matrix or dataframe, rows are samples and columns are features, rownames are sample names
   M <- df0[,colnames(df0)!='label']
   M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))

   M[is.na(M)] <- 0
   m1 <- prcomp(M,T);
   Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
   Y  <- Y[,c(1,2)]
   
   Y <- data.frame(Y,df0$label);
   colnames(Y) <- c("PC1","PC2","label")

   eigs <- m1$sdev^2
   percentages <- eigs[1:2] / sum(eigs)
   lblColors   <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3',Training='#537e35',Validation='red')
   #lblColors <- c(Training='#537e35',Validation='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank())
   p <- p +  labs(x = sprintf("PC1(%4.2f%%)",percentages[1]*100),
                  y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
	      title =strTitle)
   p <- p +   scale_colour_manual(values=lblColors)
  
   p
}



SWT <- read.table("E:/projects/TPD/data/TPDT_OpenSWATH_Win600_20181115.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
SWV <- read.table("E:/projects/TPD/data/TPDV_openSWATH_avgRepcas_20181108.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
colnames(SWT) <- sapply(colnames(SWT),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
rownames(SWT) <- SWT$patientId
colnames(SWV) <- sapply(colnames(SWV),function(v){strsplit(v,"\\|")[[1]][1]})##random forest 对列名有要求
label <- as.character(sapply(SWT$patientId,function(v){substr(v,1,1)}))

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

pr1 <- c('O75396','P02765','P06744')
pr2 <- c('P07305','P21397','Q13263')
pr3 <- c('O75083','P02765','P10809')
pr4 <- c('P04083','P30044','P10809')
tenprot <- unique(c(pr1,pr2,pr3,pr4))

#######################################
M0 <-  trainM[,protIds]   ## raw training matrix
T0  <-  testM[,protIds]   ## raw validation matrix

M0    <- data.frame(t(apply(M0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
tmean <- apply(M0,2,mean,na.rm=T)
tsd   <- apply(M0,2,sd,na.rm=T)
M0    <- data.frame(apply(M0,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
T0 <- data.frame(t(apply(T0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
T0 <- (T0-tmean)/tsd

M0$label <- label
T0$label <- rep('Validation',dim(T0)[1])

p1 <- drawPCA(rbind(M0,T0),'218protein:training+validation')



M0 <- M0[,tenprot]
T0 <- T0[,tenprot]

###################################
M0 <-  trainM[,tenprot]   ## raw training matrix
T0 <-  testM[,tenprot]   ## raw validation matrix

M0    <- data.frame(t(apply(M0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
tmean <- apply(M0,2,mean,na.rm=T)
tsd   <- apply(M0,2,sd,na.rm=T)
M0    <- data.frame(apply(M0,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)}))
T0 <- data.frame(t(apply(T0,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
T0 <- (T0-tmean)/tsd

M0$label <- label
T0$label <- rep('Validation',dim(T0)[1])

p2 <- drawPCA(rbind(M0,T0),'10 protein:training+validation')

print(p1)
print(p2)

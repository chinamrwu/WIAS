library("ggplot2")


rawMat <- read.table("E:/projects/TPD/data/RF_TPDT_OpenSWATH_trainningM_20181022.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(rawMat) <- rawMat[,1]
testData <- read.table("E:/projects/TPD/data/RF_TPDV_expert_20181024.txt",stringsAsFactors = F,header = T,check.names =F,sep="\t")
protIds <- intersect(colnames(rawMat)[3:dim(rawMat)[2]],colnames(testData))
trainM <- rawMat[,c("label",protIds)]
missRate <- apply(trainM[,-1],2,function(v){100*sum(v==0)/dim(trainM)[1]})
trainM <- rawMat[,c("label",names(missRate)[missRate<=30])]

IDMaps <- read.table("E:/projects/TPD/data/Gene_Prot_Mapping.txt",header=T,sep="\t",stringsAsFactors=F)

lbls <-c("N","M","A","C","P")
protIds <-         c("P04083","P07202","P02751","P17931","P35268","O60504","Q99584","P10909","P46778","Q71UM5") 
protIds <- c(protIds,"P61353","P20700","O75083","P08294","Q96C19","Q16851","Q14764","P12277","P19971","O94875") 

protIds <- c("Q14764","O60504","P10909","O60437","Q13642","P12277","O94875","P02743","P02787","P02765","P04179","Q9HCC0","P07202")
protIds <- c(protIds,"P25705","Q92506","P39023","P52272","P02452","P46940","P61978","P13797","P02751","P04083","P07355","P17931","P01266")




t1 <- sapply(protIds,function(ID){
   
  })
t1 <- data.frame(t1,stringsAsFactors=F)
t1$label <- rownames(t1)



indx <- match(protIds,IDMaps$protId,nomatch=-1)
indx <- indx[which(indx>0)]
protIds <- IDMaps[indx,]


#par(mfrow = c(4, 5),mar=c(0,0,0,0))
pdf("E:/projects/TPD/results/TPD_importantProtein_openSWATH_02.pdf",width=12,height=16)
grid.newpage();
pushViewport(viewport(layout = grid.layout(7, 4)));

indx=0
for(ID in protIds$protId){
  geneName <- protIds$gene[protIds$protId==ID]

   t1 <- data.frame(t(sapply(lbls,function(lbl){
	tmp <- trainM[trainM$label==lbl,ID]
	c(mean(tmp),sd(tmp))
   })))
   colnames(t1) <- c("log2Intensity","SD")
   t1$label <- lbls

  tmp <- t1
 
  p1 <- ggplot(tmp, aes(x=label,y=log2Intensity,fill=label)) +
  geom_bar(position="dodge",stat="identity")+
  geom_errorbar(aes(ymin=log2Intensity-SD, ymax=log2Intensity+SD), width=.2, position=position_dodge(.9))+
  #geom_text(data=tmp,aes(x=label,y=log2Intensity+2,label=sprintf("%0.2f", round(SD, digits = 2))))+
  ggtitle(paste0(ID,"|",geneName))+
  scale_fill_manual("legend", values = c("N" = "black", "M" = "orange", "A" = "blue","C"="brown3","P"="cyan3"))+
  theme(  panel.grid.major = element_blank(),
   	    panel.grid.minor = element_blank(),
   	    panel.border = element_blank(),
   	    panel.background = element_blank(),
	   legend.position='none')
  
  print(p1,vp=viewport(layout.pos.row = indx %/% 4+1, layout.pos.col = indx %% 4+1))
  indx <- indx+1
}
dev.off()



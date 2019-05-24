featureA    <- rownames(coef(models[[1]],s="lambda.1se"))[-1]
featureB    <- rownames(coef(models[[2]],s="lambda.1se"))[-1]
featureC    <- rownames(coef(models[[3]],s="lambda.1se"))[-1]


M0 <- df0[,unique(c('label',featureB))];
ump <- drawUMAP(M0,color2,rowNormalization=T,strTitle=sprintf('UMP:%d features from B',dim(M0)[2]-1))
umpB <- ump + scale_x_continuous(breaks = round(seq(min(ump$data$X), max(ump$data$X), by = 0.5),1))+ 
     scale_y_continuous(breaks = round(seq(min(ump$data$Y), max(ump$data$Y), by = 0.5),1)) 


M0 <- df0[,unique(c('label',featureB,featureA))];
ump <- drawUMAP(M0,color2,rowNormalization=T,strTitle=sprintf('UMP:%d features from A B',dim(M0)[2]-1))
umpBA <- ump + scale_x_continuous(breaks = round(seq(min(ump$data$X), max(ump$data$X), by = 0.5),1))+ 
     scale_y_continuous(breaks = round(seq(min(ump$data$Y), max(ump$data$Y), by = 0.5),1)) 


M0 <- df0[,unique(c('label',featureB,featureC))];
ump <- drawUMAP(M0,color2,rowNormalization=T,strTitle=sprintf('UMP:%d features from B C',dim(M0)[2]-1))
umpBC <- ump + scale_x_continuous(breaks = round(seq(min(ump$data$X), max(ump$data$X), by = 0.5),1))+ 
     scale_y_continuous(breaks = round(seq(min(ump$data$Y), max(ump$data$Y), by = 0.5),1)) 

M0 <- df0[,unique(c('label',featureA,featureC))];
ump <- drawUMAP(M0,color2,rowNormalization=T,strTitle=sprintf('UMP:%d features from A C',dim(M0)[2]-1))
umpAC <- ump + scale_x_continuous(breaks = round(seq(min(ump$data$X), max(ump$data$X), by = 0.5),1))+ 
     scale_y_continuous(breaks = round(seq(min(ump$data$Y), max(ump$data$Y), by = 0.5),1)) 


M0 <- df0[,unique(c('label',featureA,featureB,featureC))];
ump <- drawUMAP(M0,color2,rowNormalization=T,strTitle =sprintf('UMP:%d features from A B C',dim(M0)[2]-1))
umpABC <- ump + scale_x_continuous(breaks = round(seq(min(ump$data$X), max(ump$data$X), by = 0.5),1))+ 
     scale_y_continuous(breaks = round(seq(min(ump$data$Y), max(ump$data$Y), by = 0.5),1)) 


######################################
dat <- umpABC$data
ABC <- c(sum(dat$Y < - 3.2 & dat$label=='M'),sum(dat$Y < -3.2 & dat$label=='B'))

dat <- umpAC$data
AC <- c(sum(dat$Y > 3 & dat$label=='M'),sum(dat$Y > 3 & dat$label=='B'))

dat <- umpBC$data
BC <- c(sum(dat$Y > 2.3 & dat$label=='M'),sum(dat$Y >2.3 & dat$label=='B')) 

dat <- umpBA$data
BA <- c(sum(dat$Y < - 3.0 & dat$label=='M'),sum(dat$Y < -3.0 & dat$label=='B'))

dat <- umpB$data
B0 <- c(sum(dat$Y < - 3.2 & dat$label=='M'),sum(dat$Y < -3.2 & dat$label=='B'))


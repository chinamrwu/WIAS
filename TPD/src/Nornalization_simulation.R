library(MASS)
library(ggplot2)

sampleSize =1000
rds1 <- mvrnorm(sampleSize, mu=c(3,8), Sigma=rbind(c(1, 0),c(0, 1)))
rds2 <- mvrnorm(sampleSize, mu=c(0,5), Sigma=rbind(c(1, 0),c(0, 1)))
ds1 <- rbind(rds1,rds2)
m1   <- apply(rds1,2,mean)
sd1  <- apply(rds1,2,sd)


ds4 <- rbind(ds1,apply(rds1,2,function(v){(v-mean(v))/sd(v)}))
ds4 <- rbind(ds4,apply(rds2,2,function(v){(v-m1)/sd(v)}))
ds4 <- data.frame(ds4)
ds4$label <- 'U'
ds4$label[1:sampleSize] <- 'training'
ds4$label[(sampleSize+1):(2*sampleSize)] <- 'validation'
ds4$label[(2*sampleSize+1):(3*sampleSize)] <- 'feature_normalized_training'
ds4$label[(3*sampleSize+1):(4*sampleSize)] <- 'normalized_validation_usingTrainingMean'
colnames(ds4) <- c('X','Y','label')



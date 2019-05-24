remove(list=ls())
library(randomForest)
library(pROC)

# generate some random data
set.seed(1111)
train <- data.frame(condition = sample(c("mock", "lethal"), replace = T, size = 1000))
 test <- data.frame(condition = sample(c("mock", "lethal"), replace = T, size = 1000))

train$feat01 <- sapply(train$condition, (function(i){ if (i == "mock") { rnorm(n = 1, mean = 0)} else if (i == "lethal") { rnorm(n = 1, mean = 1.5)} else { rnorm(n = 1, mean = -1.5)} }))
train$feat02 <- sapply(train$condition, (function(i){ if (i == "mock") { rnorm(n = 1, mean = 0)} else if (i == "lethal") { rnorm(n = 1, mean = 1.5)} else { rnorm(n = 1, mean = -1.5)} }))
train$feat03 <- sapply(train$condition, (function(i){ if (i == "mock") { rnorm(n = 1, mean = 0)} else if (i == "lethal") { rnorm(n = 1, mean = 1.5)} else { rnorm(n = 1, mean = -1.5)} }))


test$feat01 <- sapply(test$condition, (function(i){ if (i == "mock") { rnorm(n = 1, mean = 0)} else if (i == "lethal") { rnorm(n = 1, mean = 1.5)} else { rnorm(n = 1, mean = -1.5)} }))
test$feat02 <- sapply(test$condition, (function(i){ if (i == "mock") { rnorm(n = 1, mean = 0)} else if (i == "lethal") { rnorm(n = 1, mean = 1.5)} else { rnorm(n = 1, mean = -1.5)} }))
test$feat03 <- sapply(test$condition, (function(i){ if (i == "mock") { rnorm(n = 1, mean = 0)} else if (i == "lethal") { rnorm(n = 1, mean = 1.5)} else { rnorm(n = 1, mean = -1.5)} }))


model <- randomForest(formula = condition ~ ., data = train, ntree = 10, maxnodes= 100, norm.votes = F) 
predictions <- as.data.frame(predict(model, test, type = "prob"))

predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- test$condition

# 1 ROC curve, mock vs non mock
roc.mock <- roc(ifelse(predictions$observed=="mock", "mock", "non-mock"), as.numeric(predictions$mock))
plot(roc.mock, col = "gray60")

# others
roc.lethal <- roc(ifelse(predictions$observed=="lethal", "lethal", "non-lethal"), as.numeric(predictions$mock))
#roc.resist <- roc(ifelse(predictions$observed=="resist", "resist", "non-resist"), as.numeric(predictions$mock))
lines(roc.lethal, col = "blue")
#lines(roc.resist, col = "red")



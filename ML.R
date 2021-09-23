# short introduction
## https://machinelearningmastery.com/compare-the-performance-of-machine-learning-algorithms-in-r/
# caret workflow
## https://www.machinelearningplus.com/machine-learning/caret-package/#1introduction
###install.packages(c('caret', 'skimr', 'RANN', 'randomForest', 'fastAdaboost', 'gbm', 'xgboost', 'caretEnsemble', 'C50', 'earth'))
library(tidyverse)
library(mlbench)
library(caret) # classifier and regression models , 237 models  metric ( ROC )
  #install.packages('e1071', dependencies=TRUE)
library(randomForest)
library(skimr)
library(RANN)# require for knnImpute
library(doMC) 
data("PimaIndiansDiabetes")

# selected variable from mixed model result  
var_train <- unique(trial_p$variable)

#### two datasets : train, test####
# createDataPartition(data$Y,p=0.8,list=F)
# trainData <- data[trainRowNumber,]
# testData  <- data[-trainRowNumber,]
trainData <- df_trial %>%filter(Age %in% c(6,7,8), Group %in% c('HOM_Untreated','WT_Untreated')) %>%  select(Age,Group,var_train)
colnames(trainData) <- make.names(colnames(trainData))
Age<- trainData$Age
Group <- as.factor(trainData$Group)

#### missing values : preProcess()####
anyNA(trainData) #F
anyNA(testData)  #F
# preNA <- preProcess(trainData,method = 'knnImpute')
# trainData <- predict(preNA,newdata = trainData)

#### one-hot encoding #### 
# this for character data , we only 
# dummyModel <- dummyVars(Y~.,data=trainData)
# trainData_mat <- predict(dummies_model, newdata = trainData)
# trainData <- data.frame(trainData_mat)

#### tramsform the data #### 
preP_m <- preProcess(trainData,method = 'range')
trainData <- predict(preP_m,newdata = trainData)
trainData$Age <- Age     # put the age back
trainData$Group <- as.factor(Group) # put the group back

apply(trainData[, 1:10], 2, FUN=function(x){c('min'=min(x), 'max'=max(x))})

#### visualization of feature ####
featurePlot(x=trainData[,3:10],
            y=as.factor(trainData$Group), # y should be factor, or returns NULL
            plot = 'box',
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales=list(x=list(relation='free'),
                        y=list(relation='free')))
featurePlot(x=trainData[,3:10],
            y=as.factor(trainData$Group),
            plot = 'density',
            strip=strip.custom(par.strip.text=list(cex=.5)),
            scales=list(x=list(relation='free'),
                        y=list(relation='free')))
#### feature selection ####
# recursive feature elimination : rfe 
# result : 3 (coupbling variables)
subsets <- c(1:5,10,15,18)
ctrl <- rfeControl(functions = rfFuncs,
                   method='repeatedcv',
                   repeats = 5,
                   verbose = F)
lmProfile <- rfe(x=trainData[,-2],y=as.factor(trainData$Group),sizes = subsets,rfeControl = ctrl)
plot(lmProfile,type=c('g','o'))



#### train model ####

# all ml model supported in caret package (239)
# screen the models which are able to do classification
# exclude the model whose package is not availbale to R
# in the caret package (175)
modelnames <- names(getModelInfo())
model_screen <- list()
lb <- c()
for (i in 1:length(modelnames)) {
  m <- modelnames[[i]]
  for_cla <- unique(modelLookup(modelnames[[1]])[,5]) # forClass

  lt  <- getModelInfo(modelnames[[i]]) # check the libraries needed
  for (j in 1:length(lt)) {
    x <- lt[[j]]$library
    lb <- append(lb,x)
  }
  model_screen[[i]] <- data.frame(model=m,forClass=for_cla)
  
}
model_screen <- do.call(rbind.data.frame,model_screen)
model_screen <- model_screen %>% filter(forClass=='TRUE')
#### installed packages needed ####
lb <- unique(lb)
list.of.packages <- lb
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lb_non <- c('adaptDA', 'CHAID', 'elmNN', 'foba', 'gpls', 'rrlda', 'logicFS', 'FCNN4R', 'mxnet', 'randomGLM', 'vbmp' )


#### exclude the model which's package is not available 
model_screen <- list() 
for (i in 1:length(modelnames)) {
  m <- modelnames[[i]]
  y <- c()
  lt  <- getModelInfo(modelnames[[i]]) # check the libraries needed
  for (j in 1:length(lt)) {
    x <- lt[[j]]$library
    y <- append(y,x)
  }
  if(!any(y %in% lb_non)){
    for_cla <- unique(modelLookup(modelnames[[i]])[,5]) #forClass
    model_screen[[i]] <- data.frame(model=m,forClass=for_cla)
    
    }
}
model_screen <- do.call(rbind.data.frame,model_screen)
model_screen <- model_screen %>% filter(forClass=='TRUE')


#### prepare the test data ####
testData <-  df_trial %>%filter(Age %in% c(6,7,8), Group %in% c('WT_Tanganil_8wk','HOM_Tanganil_8wk')) %>%  select(Age,Group,var_train)
colnames(testData) <- make.names(colnames(testData))
Age<- testData$Age
Group <- testData$Group
Group <- str_replace_all(Group,'_Tanganil_8wk','_Untreated')
preP_m <- preProcess(testData,method = 'range')
testData <- predict(preP_m,newdata = testData)
testData$Age <- Age     # put the age back
testData$Group <- as.factor(Group) # put the group back

#### predict on testData ####

#### confusion matrix ####

#### hyperparameter tuning to optimize the performance ####
# tuneLength : the number of unique values for the 
  # tuning parameters caret will consider
# tuneGrid : 

fitControl <- trainControl(
  method = 'repeatedcv',
  number = 5,
  savePredictions = 'final',       # save prediction for optimal tuning parameter
  classProbs=T,                    # class probabilities
  summaryFunction=twoClassSummary, # results summary func
  allowParallel = T
  )

# model_mars2 = train(Group ~ ., data=trainData, method='earth', tuneLength = 5, metric='ROC', trControl = fitControl)
# model_mars2
# pred_data <- predict(model_mars2,testData[,-2])
# confusionMatrix(reference=testData$Group,data = pred_data,mode='everything')

#### train all the models ####

library(foreach)
library(doParallel)
library(pls)
library(doMC)
cores <- detectCores()
registerDoMC(cores = cores-2)
# cl <- makeCluster(cores[1]-3,type = "PSOCK")
# pls.options(parallel = cl)

# registerDoParallel(cl)

ml_list <- list()
model_selected  <- as.character(model_screen$model)   # the models for class
model_selected <- c('regLogistic','naive_bayes','glm','xgbTree','knn')
modelLookup('logreg')# 'logreg',
for (i in 1:length(model_selected)) {
  tryCatch({    
  ml_list[[i]] <- caret::train(Group~.,data=trainData,method=model_selected[[i]],trControl=fitControl,tuneLength=5,metric='ROC')
  # names(ml_list[i]) <- model_selected[[i]]
  } ,error=function(e){cat('ERROR :',conditionMessage(e),'model :',model_selected[i] )})
}

model_result <- resamples(ml_list[1:5])
summary(model_result)
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(model_result, scales=scales)

# ml_list<- foreach (i = 1:length(model_selected)) %dopar% {
#     caret::train(Group~.,data=trainData,method=model_selected[i],trControl=fitControl,tuneLength=5,metric='ROC')
#   }
stopCluster(cl)


#### select models by tag ####
tag <- read.csv('tag_data.csv')
tag[tag[,'Regression']==1,1]

  results <- resamples(ml_list)
summary(results)

#descriptive statistics
skim_to_wide(trainData)

# colnames(df_trial_t) <- make.names(colnames(df_trial_t))
# 
# control <- trainControl(method = 'repeatedcv',number = 10,repeats = 3)


str(ml_list)
results <- resamples(ml_list)
summary(results)


# box and whisker plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

# density plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))
densityplot(results, scales=scales, pch = "|")

parallelplot(results)

qqplot(ml_list[[1]])

splom(results)
train_rpart <- train(Group~.,data = df_trial_t,method='svmRadial',trControl=control)

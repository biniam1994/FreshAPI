# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()

library("readxl")
library(caret)
library(Rmisc)

source("CreateObj.R")
source("accuracy.R")

# xlsx files
df <- read_excel("chickenBurger_FTIR_raw&sg_v2.xlsx")
dfA <- read_excel("chickenBurger_FTIR_raw&sg_batchA.xlsx")
dfB <- read_excel("chickenBurger_FTIR_raw&sg_batchB.xlsx")

#######################################################################################
# set the row names and remove the first column with the row IDs after this treatment
sampleNames <- c()
for(i in 1:nrow(df[,1])){
  sampleNames[i] <- df[i,1]
}
row.names(df) <- sampleNames
df <- df[,-1]

sampleNames_A <- c() #batch A
for(i in 1:nrow(dfA[,1])){
  sampleNames_A[i] <- dfA[i,1]
}
row.names(dfA) <- sampleNames_A
dfA <- dfA[,-1]

sampleNames_B <- c() #batch B
for(i in 1:nrow(dfB[,1])){
  sampleNames_B[i] <- dfB[i,1]
}
row.names(dfB) <- sampleNames_B
dfB <- dfB[,-1]

######################################################################
#feature extraction (for wavelength between 0 and 1000 nm)
removeID <- c()
features <- colnames(df)
for(i in 6:ncol(df)){
  if(as.numeric(features[i]) <= 1000){
    removeID[length(removeID)+1] <- i
  }
}
df <- df[,-removeID]

#batch A
removeID_A <- c()
featuresA <- colnames(dfA)
for(i in 6:ncol(dfA)){
  if(as.numeric(featuresA[i]) <= 1000){
    removeID_A[length(removeID_A)+1] <- i
  }
}
dfA <- dfA[,-removeID_A]

#batch B
removeID_B <- c()
featuresB <- colnames(dfB)
for(i in 6:ncol(dfB)){
  if(as.numeric(featuresB[i]) <= 1000){
    removeID_B[length(removeID_B)+1] <- i
  }
}
dfB <- dfB[,-removeID_B]

#####################################################################################################
############################# feature selection step 1 ##############################################
############################# remove redundant features #############################################

#ensure the results are repeatable
set.seed(123)
# create specra dataframe
spectra <- df[,6:3000]
#calculate correlation matrix
correlationMatrix <- cor(spectra)
#summarize the correlation matrix
print(correlationMatrix)
# find attributes which are highly correlated (ideally > 0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.75)
# print index of highly correlated attributes
print(highlyCorrelated)
#remove highly correlated attributes
spectra2 <- spectra[,-highlyCorrelated]
# create a dataframe with removed features
df2 <- CreateObj(df[2:5], spectra2)


################################### batch A ###########################################################
#ensure the results are repeatable
set.seed(123)
# create specra dataframe
spectraA <- dfA[,5:2999]
spectraB <- dfB[,5:2999]
#calculate correlation matrix
correlationMatrixA <- cor(spectraA)
#summarize the correlation matrix
print(correlationMatrixA)
# find attributes which are highly correlated (ideally > 0.75)
highlyCorrelatedA <- findCorrelation(correlationMatrixA, cutoff = 0.75)
# print index of highly correlated attributes
print(highlyCorrelatedA)
#remove highly correlated attributes
spectra2A <- spectraA[,-highlyCorrelatedA]
spectra2B <- spectraB[,-highlyCorrelatedA]
# create a dataframe with removed features
dfA2 <- CreateObj(dfA[1:4], spectra2A)
dfB2 <- CreateObj(dfB[1:4], spectra2B)

################################## batch B ###########################################################
#ensure the results are repeatable
set.seed(123)
# create specra dataframe
spectraA_2 <- dfA[,5:ncol(dfA)]
spectraB_2 <- dfB[,5:ncol(dfB)]
#calculate correlation matrix
correlationMatrixA_2 <- cor(spectraB_2)
#summarize the correlation matrix
print(correlationMatrixA_2)
# find attributes which are highly correlated (ideally > 0.75)
highlyCorrelatedA_2 <- findCorrelation(correlationMatrixA_2, cutoff = 0.75)
# print index of highly correlated attributes
print(highlyCorrelatedA_2)
#remove highly correlated attributes
spectra2A_2 <- spectraA_2[,-highlyCorrelatedA_2]
spectra2B_2 <- spectraB_2[,-highlyCorrelatedA_2]
# create a dataframe with removed features
dfA2_B <- CreateObj(dfA[1:4], spectra2A_2)
dfB2_B <- CreateObj(dfB[1:4], spectra2B_2)


#######################################################################################################
############################# feature selection step 2 ################################################
############################# rank features by importance #############################################
#prerate training scheme
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# create df for diffrent bacterial counts separetly
df_TSA_wf <- df2[,-c(2:4)]
df_CFC_wf <- df2[,-c(1,3:4)]
df_STAA_wf <- df2[,-c(1:2,4)]
df_MRS_wf <- df2[,-c(1:3)]

# create df for diffrent bacterial counts separetly batch A based on B
dfA_TSA_wf <- dfA2[,-c(2:4)]
dfA_CFC_wf <- dfA2[,-c(1,3:4)]
dfA_STAA_wf <- dfA2[,-c(1:2,4)]
dfA_MRS_wf <- dfA2[,-c(1:3)]

# create df for diffrent bacterial counts separetly batch B based on A
dfB_TSA_wf <- dfB2[,-c(2:4)]
dfB_CFC_wf <- dfB2[,-c(1,3:4)]
dfB_STAA_wf <- dfB2[,-c(1:2,4)]
dfB_MRS_wf <- dfB2[,-c(1:3)]

# create df for diffrent bacterial counts separetly batch A based on B
dfA_TSA_wf_B <- dfA2_B[,-c(2:4)]
dfA_CFC_wf_B <- dfA2_B[,-c(1,3:4)]

# create df for diffrent bacterial counts separetly batch B based on B
dfB_TSA_wf_B <- dfB2_B[,-c(2:4)]
dfB_CFC_wf_B <- dfB2_B[,-c(1,3:4)]

#######################################################################################################
#################################### TSA ##############################################################
#######################################################################################################

######################### select features for the whole dataset #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_tvc <- train(TSA~., data = df_TSA_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance <- varImp(model_varImp_tvc, scale = FALSE)
#sumarize importance
print(importance)
#plot importance
plot(importance)
# get an index of important variables
feature_ID_TSA <- c(1)
for(i in 1:nrow(importance$importance)){
  if(importance$importance[i,1] >= 0.5){
    feature_ID_TSA[length(feature_ID_TSA)+1] <- i+1
  }
}
# create a dataframe eith selected variables
df_TSA <- df_TSA_wf[,feature_ID_TSA]

######################### select features for the batch A #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_tvc_bA <- train(TSA~., data = dfA_TSA_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_bA <- varImp(model_varImp_tvc_bA, scale = FALSE)
#sumarize importance
print(importance_bA)
#plot importance
plot(importance_bA)
# get an index of important variables
feature_ID_TSA_bA <- c(1)
for(i in 1:nrow(importance_bA$importance)){
  if(importance_bA$importance[i,1] >= 0.5){
    feature_ID_TSA_bA[length(feature_ID_TSA_bA)+1] <- i+1
  }
}
# create a dataframe eith selected variables
dfA_TSA <- dfA_TSA_wf[,feature_ID_TSA_bA]
dfB_TSA <- dfB_TSA_wf[,feature_ID_TSA_bA]

######################### select features for the batch B #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_tvc_bB <- train(TSA~., data = dfB_TSA_wf_B, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_bB <- varImp(model_varImp_tvc_bB, scale = FALSE)
#sumarize importance
print(importance_bB)
#plot importance
plot(importance_bB)
# get an index of important variables
feature_ID_TSA_bB <- c(1)
for(i in 1:nrow(importance_bB$importance)){
  if(importance_bB$importance[i,1] >= 0.5){
    feature_ID_TSA_bB[length(feature_ID_TSA_bB)+1] <- i+1
  }
}
# create a dataframe eith selected variables
dfA_TSA_B <- dfA_TSA_wf[,feature_ID_TSA_bB]
dfB_TSA_B <- dfB_TSA_wf[,feature_ID_TSA_bB]

###############################################################################################################
################################### Prediction for TSA ########################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_tsa1 <- createDataPartition(df_TSA$TSA, p = .7, 
                                       list = FALSE, 
                                       times = 1, groups =3) #caret package

#split the whole dataset
train_tsa1 <- df_TSA[trainIndex_tsa1,]
test_tsa1 <- df_TSA[-trainIndex_tsa1,]



# resembling (method to avoid overfitting)
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)
## repeated ten times

#############################################################################################
############################# feature selection step 3 ################################################
############################# feature selection #############################################

# rfe control (define the control using a random forest selection function)
rctrl1 <- rfeControl(method = "cv",
                     number = 10, #10-fold cv
                     functions = rfFuncs)

######################### select features for the whole dataset #######################################
# run the RFE algorithm
results_TSA <- rfe(TSA ~ ., data = train_tsa1,
                   sizes = c(2:ncol(train_tsa1)), rfeControl = rctrl1)
#sumarize the results
print(results_TSA)
#list the chosen features
predictors_tsa <- predictors(results_TSA)
#plot the results
plot(results_TSA, type = c("g", "o"))
#create a dataframe with selected features
predictors_tsa_id <- c(1)
for(i in 2:ncol(train_tsa1)){
  for(j in 1:length(predictors_tsa)){
    if( as.numeric(colnames(train_tsa1[i]))  == as.numeric(gsub("`", "", (as.name(predictors_tsa[j]))))){
      predictors_tsa_id[length(predictors_tsa_id)+1] <- i
    }
  }
}
train_tsa2 <- train_tsa1[,predictors_tsa_id]
test_tsa2 <- test_tsa1[,predictors_tsa_id]
df2_TSA <- df_TSA[,predictors_tsa_id]

#save selected features in one vector
predictors_tsa_toSave <- c()
len_predictors_tsa <- 0
for(i in predictors_tsa){
  if(as.numeric(lengths(predictors_tsa[i]))== 1){
    len_predictors_tsa <- len_predictors_tsa +1
  }
}
for(i in 1:len_predictors_tsa){
  predictors_tsa_toSave[i] <- as.numeric(gsub("`", "", (as.name(predictors_tsa[i]))))
}
saveRDS(predictors_tsa_toSave, "./models/CB/FTIR/predictors_knn_SVMLM_SVMR_SVMP_LM_TVC.rds") #save model as an rds object

######################### select features for the btch A #######################################
# run the RFE algorithm
results_TSA_bA <- rfe(TSA ~ ., data = dfA_TSA,
                   sizes = c(2:length(feature_ID_TSA_bA)), rfeControl = rctrl1)
#sumarize the results
print(results_TSA_bA)
#list the chosen features
predictors_tsa_bA <- predictors(results_TSA_bA)
#plot the results
plot(results_TSA_bA, type = c("g", "o"))
#create a dataframe with selected features
predictors_tsa_id_bA <- c(1)
for(i in 2:ncol(dfA_TSA)){
  for(j in 1:length(predictors_tsa_bA)){
    if( as.numeric(colnames(dfA_TSA[i]))  == as.numeric(gsub("`", "", (as.name(predictors_tsa_bA[j]))))){
      predictors_tsa_id_bA[length(predictors_tsa_id_bA)+1] <- i
    }
  }
}
dfA2_TSA <- dfA_TSA[,predictors_tsa_id_bA]
dfB2_TSA <- dfB_TSA[,predictors_tsa_id_bA]

######################### select features for the batch B #######################################
# run the RFE algorithm
results_TSA_bB <- rfe(TSA ~ ., data = dfB_TSA_B,
                      sizes = c(2:length(feature_ID_TSA_bB)), rfeControl = rctrl1)
#sumarize the results
print(results_TSA_bB)
#list the chosen features
predictors_tsa_bB <- predictors(results_TSA_bB)
#plot the results
plot(results_TSA_bB, type = c("g", "o"))
#create a dataframe with selected features
predictors_tsa_id_bB <- c(1)
for(i in 2:ncol(dfB_TSA_B)){
  for(j in 1:length(predictors_tsa_bB)){
    if( as.numeric(colnames(dfB_TSA_B[i]))  == as.numeric(gsub("`", "", (as.name(predictors_tsa_bB[j]))))){
      predictors_tsa_id_bB[length(predictors_tsa_id_bB)+1] <- i
    }
  }
}
dfA2_TSA_B <- dfA_TSA_B[,predictors_tsa_id_bB]
dfB2_TSA_B <- dfB_TSA_B[,predictors_tsa_id_bB]

#############################################################################################
#################### split batch A ##########################################################
#split data
set.seed(123) #caret package
#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndexA_tsa1 <- createDataPartition(dfA2_TSA$TSA, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainA_tsa1 <- dfA2_TSA[trainIndexA_tsa1,]
testA_tsa1 <- dfA2_TSA[-trainIndexA_tsa1,]

#################### split batch B ##########################################################
#split data
set.seed(123) #caret package
#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndexB_tsa1 <- createDataPartition(dfB2_TSA_B$TSA, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainB_tsa1 <- dfB2_TSA_B[trainIndexB_tsa1,]
testB_tsa1 <- dfB2_TSA_B[-trainIndexB_tsa1,]

#############################################################################################
################################ RF_TVC #####################################################

##################
# the first method
model_Rf_TSA1 <- train(TSA ~ ., data = train_tsa2,
                       method = "rf", trControl = fitControl)
pred_Rf_TSA1<- predict(model_Rf_TSA1,test_tsa2)
RMSE.Rf_TSA1 <- RMSE(test_tsa2$TSA,pred_Rf_TSA1) #RMSE =0.54

# rf plot predicted vs Observed
plot(pred_Rf_TSA1,test_tsa2$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_TSA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.Rf_TSA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TSA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_Rf_TSA <- train(TSA~., data = train_tsa3, method = "rf", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TSA,predicted.Rf_TSA)
  # put results to the vectors
  rmse.Rf_TSA.vector <- c(rmse.Rf_TSA.vector, rmse.Rf_TSA)
}
# compute the standard deviation
sd.Rf_TSA <- round(sd(rmse.Rf_TSA.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.Rf_TSA <- CI(rmse.Rf_TSA.vector, ci = 0.95) #package: Rmisc #0.5469621 0.5241529 0.5013436 
# rmse mean
rmse.Rf_TSA.mean <- round(mean(rmse.Rf_TSA.vector), digits = 2) #0.52

plot(rmse.Rf_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_TSA.mean,"+/- ",sd.Rf_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_TSA[3], digits = 2),
                  "to",round(ci.Rf_TSA[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_TSA2 <- train(TSA ~ ., data = dfA2_TSA,
                       method = "rf", trControl = fitControl)
pred_Rf_TSA2<- predict(model_Rf_TSA2,dfB2_TSA)
RMSE.Rf_TSA2 <- RMSE(dfB2_TSA$TSA,pred_Rf_TSA2) #RMSE =0.88

# rf plot predicted vs Observed
plot(pred_Rf_TSA2,dfB2_TSA$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ KNN_TVC #####################################################

##################
# the first method
model_knn_TSA1 <- train(TSA ~ ., data = train_tsa2,
                       method = "knn", trControl = fitControl)
pred_knn_TSA1<- predict(model_knn_TSA1,test_tsa2)
RMSE.knn_TSA1 <- RMSE(test_tsa2$TSA,pred_knn_TSA1) #RMSE =0.7
accuracy_knn_TSA <- accuracy(pred_knn_TSA1, test_tsa2$TSA) #86.84
# this is the best k-nn model for TVC
saveRDS(model_knn_TSA1, "./models/CB/FTIR/model_knn_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TSA, "./models/CB/FTIR/tested_knn_TVC.rds") #save tested data as an rds object
saveRDS(pred_knn_TSA1, "./models/CB/FTIR/predicted_knn_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.knn_TSA1, "./models/CB/FTIR/RMSE_knn_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_knn_TSA, "./models/CB/FTIR/accuracy_knn_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_knn_TSA1,test_tsa2$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_TSA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.knn_TSA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TSA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_knn_TSA <- train(TSA~., data = train_tsa3, method = "knn", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.knn_TSA <- predict(model_knn_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.knn_TSA <- RMSE(pred = test_tsa3$TSA,predicted.knn_TSA)
  # put results to the vectors
  rmse.knn_TSA.vector <- c(rmse.knn_TSA.vector, rmse.knn_TSA)
}

# compute the standard deviation
sd.knn_TSA <- round(sd(rmse.knn_TSA.vector), digits = 2) #0.07
# find the 95% confidence intervals parameters
ci.knn_TSA <- CI(rmse.knn_TSA.vector, ci = 0.95) # 0.5484873 0.5286199 0.5087525 
# rmse mean
rmse.knn_TSA.mean <- round(mean(rmse.knn_TSA.vector), digits = 2) #0.53

plot(rmse.knn_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_TSA.mean,"+/- ",sd.knn_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_TSA[3], digits = 2),
                  "to",round(ci.knn_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_knn_TSA2 <- train(TSA ~ ., data = dfA2_TSA,
                       method = "knn", trControl = fitControl)
pred_knn_TSA2<- predict(model_knn_TSA2,dfB2_TSA)
RMSE.knn_TSA2 <- RMSE(dfB2_TSA$TSA,pred_knn_TSA2) #RMSE =0.54

# rf plot predicted vs Observed
plot(pred_knn_TSA2,dfB2_TSA$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#############################################################################################
################################ svmLM_TVC #####################################################

##################
# the first method
model_SVMLM_TSA1 <- train(TSA ~ ., data = train_tsa2,
                        method = "svmLinear", trControl = fitControl)
pred_SVMLM_TSA1<- predict(model_SVMLM_TSA1,test_tsa2)
RMSE.SVMLM_TSA1 <- RMSE(test_tsa2$TSA,pred_SVMLM_TSA1) #RMSE =0.59
accuracy_SVMLM_TSA <- accuracy(pred_SVMLM_TSA1, test_tsa2$TSA)

# this is the best svm linear model for TVC
saveRDS(model_SVMLM_TSA1, "./models/CB/FTIR/model_SVMLM_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TSA, "./models/CB/FTIR/tested_SVMLM_TVC.rds") #save tested data as an rds object
saveRDS(pred_SVMLM_TSA1, "./models/CB/FTIR/predicted_SVMLM_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMLM_TSA1, "./models/CB/FTIR/RMSE_SVMLM_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMLM_TSA, "./models/CB/FTIR/accuracy_SVMLM_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMLM_TSA1,test_tsa2$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("SVMLM_FTIR RMSE:",round(RMSE.SVMLM_TSA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMLM_TSA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TSA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMLM_TSA <- train(TSA~., data = train_tsa3, method = "svmLinear", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMLM_TSA <- predict(model_SVMLM_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.SVMLM_TSA <- RMSE(pred = test_tsa3$TSA,predicted.SVMLM_TSA)
  # put results to the vectors
  rmse.SVMLM_TSA.vector <- c(rmse.SVMLM_TSA.vector, rmse.SVMLM_TSA)
}

# compute the standard deviation
sd.SVMLM_TSA <- round(sd(rmse.SVMLM_TSA.vector), digits = 2) #0.07
# find the 95% confidence intervals parameters
ci.SVMLM_TSA <- CI(rmse.SVMLM_TSA.vector, ci = 0.95) # 0.5392655 0.5201776 0.5010897  
# rmse mean
rmse.SVMLM_TSA.mean <- round(mean(rmse.SVMLM_TSA.vector), digits = 2) #0.53

plot(rmse.SVMLM_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_TSA.mean,"+/- ",sd.SVMLM_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_TSA[3], digits = 2),
                  "to",round(ci.SVMLM_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMLM_TSA2 <- train(TSA ~ ., data = dfA2_TSA,
                        method = "svmLinear", trControl = fitControl)
pred_SVMLM_TSA2<- predict(model_SVMLM_TSA2,dfB2_TSA)
RMSE.SVMLM_TSA2 <- RMSE(dfB2_TSA$TSA,pred_SVMLM_TSA2) #RMSE =0.73

# rf plot predicted vs Observed
plot(pred_SVMLM_TSA2,dfB2_TSA$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMRadial_TVC #####################################################

##################
# the first method
model_SVMR_TSA1 <- train(TSA ~ ., data = train_tsa2,
                          method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA1<- predict(model_SVMR_TSA1,test_tsa2)
RMSE.SVMR_TSA1 <- RMSE(test_tsa2$TSA,pred_SVMR_TSA1) #RMSE =0.47
accuracy_SVMR_TSA <- accuracy(pred_SVMR_TSA1, test_tsa2$TSA)

# this is the best svm Radial model for TVC
saveRDS(model_SVMR_TSA1, "./models/CB/FTIR/model_SVMR_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TSA, "./models/CB/FTIR/tested_SVMR_TVC.rds") #save tested data as an rds object
saveRDS(pred_SVMR_TSA1, "./models/CB/FTIR/predicted_SVMR_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMR_TSA1, "./models/CB/FTIR/RMSE_SVMR_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMR_TSA, "./models/CB/FTIR/accuracy_SVMR_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMR_TSA1,test_tsa2$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_TSA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMR_TSA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TSA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMR_TSA <- train(TSA~., data = train_tsa3, method = "svmRadial", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMR_TSA <- predict(model_SVMR_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.SVMR_TSA <- RMSE(pred = test_tsa3$TSA,predicted.SVMR_TSA)
  # put results to the vectors
  rmse.SVMR_TSA.vector <- c(rmse.SVMR_TSA.vector, rmse.SVMR_TSA)
}

# compute the standard deviation
sd.SVMR_TSA <- round(sd(rmse.SVMR_TSA.vector), digits = 2) #0.06
# find the 95% confidence intervals parameters
ci.SVMR_TSA <- CI(rmse.SVMR_TSA.vector, ci = 0.95) # 0.4415299 0.4249149 0.4083000   
# rmse mean
rmse.SVMR_TSA.mean <- round(mean(rmse.SVMR_TSA.vector), digits = 2) #0.42

plot(rmse.SVMR_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_TSA.mean,"+/- ",sd.SVMR_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_TSA[3], digits = 2),
                  "to",round(ci.SVMR_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMR_TSA2 <- train(TSA ~ ., data = dfA2_TSA,
                          method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA2<- predict(model_SVMR_TSA2,dfB2_TSA)
RMSE.SVMR_TSA2 <- RMSE(dfB2_TSA$TSA,pred_SVMR_TSA2) #RMSE =0.58

# rf plot predicted vs Observed
plot(pred_SVMR_TSA2,dfB2_TSA$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMPoly_TVC #####################################################

##################
# the first method
model_SVMP_TSA1 <- train(TSA ~ ., data = train_tsa2,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_TSA1<- predict(model_SVMP_TSA1,test_tsa2)
RMSE.SVMP_TSA1 <- RMSE(test_tsa2$TSA,pred_SVMP_TSA1) #RMSE =0.49
accuracy_SVMP_TSA <- accuracy(pred_SVMP_TSA1, test_tsa2$TSA) #94.74

# this is the best svm polynomial model for TVC
saveRDS(model_SVMP_TSA1, "./models/CB/FTIR/model_SVMP_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TSA, "./models/CB/FTIR/tested_SVMP_TVC.rds") #save tested data as an rds object
saveRDS(pred_SVMP_TSA1, "./models/CB/FTIR/predicted_SVMP_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMP_TSA1, "./models/CB/FTIR/RMSE_SVMP_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMP_TSA, "./models/CB/FTIR/accuracy_SVMP_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMP_TSA1,test_tsa2$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_TSA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMP_TSA.vector <- c()
for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TSA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMP_TSA <- train(TSA~., data = train_tsa3, method = "svmPoly", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMP_TSA <- predict(model_SVMP_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.SVMP_TSA <- RMSE(pred = test_tsa3$TSA,predicted.SVMP_TSA)
  # put results to the vectors
  rmse.SVMP_TSA.vector <- c(rmse.SVMP_TSA.vector, rmse.SVMP_TSA)
}

# compute the standard deviation
sd.SVMP_TSA <- round(sd(rmse.SVMP_TSA.vector), digits = 2) #0.05
# find the 95% confidence intervals parameters
ci.SVMP_TSA <- CI(rmse.SVMP_TSA.vector, ci = 0.95) # 0.4479760 0.4341542 0.4203323 
# rmse mean
rmse.SVMP_TSA.mean <- round(mean(rmse.SVMP_TSA.vector), digits = 2) #0.43

plot(rmse.SVMP_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_TSA.mean,"+/- ",sd.SVMP_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_TSA[3], digits = 2),
                  "to",round(ci.SVMP_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMP_TSA2 <- train(TSA ~ ., data = dfA2_TSA,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_TSA2<- predict(model_SVMP_TSA2,dfB2_TSA)
RMSE.SVMP_TSA2 <- RMSE(dfB2_TSA$TSA,pred_SVMP_TSA2) #RMSE =0.59

# rf plot predicted vs Observed
plot(pred_SVMP_TSA2,dfB2_TSA$TSA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ LM_TVC #####################################################

##################
# the first method
model_lm_TSA1 <- train(TSA ~ ., data = train_tsa2,
                       method = "lm", trControl = fitControl)
pred_lm_TSA1<- predict(model_lm_TSA1,test_tsa2)
RMSE.lm_TSA1 <- RMSE(test_tsa2$TSA,pred_lm_TSA1) #RMSE =0.43
accuracy_lm_TSA<- accuracy(pred_lm_TSA1, test_tsa2$TSA) #100

# this is the best linear model for TVC
saveRDS(model_lm_TSA1, "./models/CB/FTIR/model_LM_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TSA, "./models/CB/FTIR/tested_LM_TVC.rds") #save tested data as an rds object
saveRDS(pred_lm_TSA1, "./models/CB/FTIR/predicted_LM_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.lm_TSA1, "./models/CB/FTIR/RMSE_LM_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_lm_TSA, "./models/CB/FTIR/accuracy_LM_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_lm_TSA1,test_tsa2$TSA,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.lm_TSA1,digits = 2), 
                "\nAccuracy",accuracy_lm_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.lm_TSA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TSA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_lm_TSA <- train(TSA~., data = train_tsa3, method = "lm", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.lm_TSA <- predict(model_lm_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.lm_TSA <- RMSE(pred = test_tsa3$TSA,predicted.lm_TSA)
  # put results to the vectors
  rmse.lm_TSA.vector <- c(rmse.lm_TSA.vector, rmse.lm_TSA)
}
# compute the standard deviation
sd.lm_TSA <- round(sd(rmse.lm_TSA.vector), digits = 2) #0.05
# find the 95% confidence intervals parameters
ci.lm_TSA <- CI(rmse.lm_TSA.vector, ci = 0.95) #package: Rmisc #0.48, 0.47, 0.46
# rmse mean
rmse.lm_TSA.mean <- round(mean(rmse.lm_TSA.vector), digits = 2) #0.47

plot(rmse.lm_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_TSA.mean,"+/- ",sd.lm_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.lm_TSA[3], digits = 2),
                  "to",round(ci.lm_TSA[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_lm_TSA2 <- train(TSA ~ ., data = dfA2_TSA,
                       method = "lm", trControl = fitControl)
pred_lm_TSA2<- predict(model_lm_TSA2,dfB2_TSA)
RMSE.lm_TSA2 <- RMSE(dfB2_TSA$TSA,pred_lm_TSA2) #RMSE =0.77

# rf plot predicted vs Observed
plot(pred_lm_TSA2,dfB2_TSA$TSA,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.lm_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

###################
# trained based on B, tested based on A
model_lm_TSA3 <- train(TSA ~ ., data = dfB2_TSA_B,
                         method = "lm", trControl = fitControl)
pred_lm_TSA3<- predict(model_lm_TSA3,dfA2_TSA_B)
RMSE.lm_TSA3 <- RMSE(dfA2_TSA_B$TSA,pred_lm_TSA3) 

# rf plot predicted vs Observed
plot(pred_lm_TSA3,dfA2_TSA_B$TSA,xlim= c(0,11), ylim=c(0,11),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",
                round(RMSE.lm_TSA3,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# trained and tested based on batch A
model_lm_TSA4 <- train(TSA ~ ., data = trainA_tsa1,
                         method = "lm", trControl = fitControl)
pred_lm_TSA4<- predict(model_lm_TSA4,testA_tsa1)
RMSE.lm_TSA4 <- RMSE(testA_tsa1$TSA,pred_lm_TSA4) 

# plot predicted vs Observed
plot(pred_lm_TSA4,testA_tsa1$TSA,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",
                round(RMSE.lm_TSA4,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#######################
# trained and tested based on batch B
model_lm_TSA5 <- train(TSA ~ ., data = trainB_tsa1,
                         method = "lm", trControl = fitControl)
pred_lm_TSA5<- predict(model_lm_TSA5,testB_tsa1)
RMSE.lm_TSA5 <- RMSE(testB_tsa1$TSA,pred_lm_TSA5) 

# plot predicted vs Observed
plot(pred_lm_TSA5,testB_tsa1$TSA,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",
                round(RMSE.lm_TSA5,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#######################################################################################################
#######################################################################################################
#################################### CFC ##############################################################
#######################################################################################################

######################### select features for the whole dataset #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_cfc <- train(CFC~., data = df_CFC_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_CFC <- varImp(model_varImp_cfc, scale = FALSE)
#sumarize importance
print(importance_CFC)
#plot importance
plot(importance_CFC)
# get an index of important variables
feature_ID_CFC <- c(1)
for(i in 1:nrow(importance_CFC$importance)){
  if(importance_CFC$importance[i,1] >= 0.5){
    feature_ID_CFC[length(feature_ID_CFC)+1] <- i+1
  }
}
# create a dataframe eith selected variables
df_CFC <- df_CFC_wf[,feature_ID_CFC]

######################### select features for the batch A #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_cfc_bA <- train(CFC~., data = dfA_CFC_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_bA_CFC <- varImp(model_varImp_cfc_bA, scale = FALSE)
#sumarize importance
print(importance_bA_CFC)
#plot importance
plot(importance_bA_CFC)
# get an index of important variables
feature_ID_CFC_bA <- c(1)
for(i in 1:nrow(importance_bA_CFC$importance)){
  if(importance_bA_CFC$importance[i,1] >= 0.5){
    feature_ID_CFC_bA[length(feature_ID_CFC_bA)+1] <- i+1
  }
}
# create a dataframe eith selected variables
dfA_CFC <- dfA_CFC_wf[,feature_ID_CFC_bA]
dfB_CFC <- dfB_CFC_wf[,feature_ID_CFC_bA]

######################### select features for the batch B #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_cfc_bB <- train(CFC~., data = dfB_CFC_wf_B, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_bB <- varImp(model_varImp_cfc_bB, scale = FALSE)
#sumarize importance
print(importance_bB)
#plot importance
plot(importance_bB)
# get an index of important variables
feature_ID_CFC_bB <- c(1)
for(i in 1:nrow(importance_bB$importance)){
  if(importance_bB$importance[i,1] >= 0.5){
    feature_ID_CFC_bB[length(feature_ID_CFC_bB)+1] <- i+1
  }
}
# create a dataframe eith selected variables
dfA_CFC_B <- dfA_CFC_wf_B[,feature_ID_CFC_bB]
dfB_CFC_B <- dfB_CFC_wf_B[,feature_ID_CFC_bB]

###############################################################################################################
################################### Prediction for CFC ########################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_cfc1 <- createDataPartition(df_CFC$CFC, p = .7, 
                                       list = FALSE, 
                                       times = 1, groups =3) #caret package
#split the whole dataset
train_cfc1 <- df_CFC[trainIndex_cfc1,]
test_cfc1 <- df_CFC[-trainIndex_cfc1,]

# resembling (method to avoid overfitting)
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)
## repeated ten times

#############################################################################################
############################# feature selection step 3 ################################################
############################# feature selection #############################################

# rfe control (define the control using a random forest selection function)
rctrl1 <- rfeControl(method = "cv",
                     number = 10, #10-fold cv
                     functions = rfFuncs)

######################### select features for the whole dataset #######################################
# run the RFE algorithm
results_CFC <- rfe(CFC ~ ., data = train_cfc1,
                   sizes = c(2:ncol(train_cfc1)), rfeControl = rctrl1)
#sumarize the results
print(results_CFC)
#list the chosen features
predictors_cfc <- predictors(results_CFC)
#plot the results
plot(results_CFC, type = c("g", "o"))
#create a dataframe with selected features
predictors_cfc_id <- c(1)
for(i in 2:ncol(train_cfc1)){
  for(j in 1:length(predictors_cfc)){
    if( as.numeric(colnames(train_cfc1[i]))  == as.numeric(gsub("`", "", (as.name(predictors_cfc[j]))))){
      predictors_cfc_id[length(predictors_cfc_id)+1] <- i
    }
  }
}
train_cfc2 <- train_cfc1[,predictors_cfc_id]
test_cfc2 <- test_cfc1[,predictors_cfc_id]
df2_CFC <- df_CFC[,predictors_cfc_id]

#save selected features in one vector
predictors_cfc_toSave <- c()
len_predictors_cfc <- 0
for(i in predictors_cfc){
  if(as.numeric(lengths(predictors_cfc[i]))== 1){
    len_predictors_cfc <- len_predictors_cfc +1
  }
}
for(i in 1:len_predictors_cfc){
  predictors_cfc_toSave[i] <- as.numeric(gsub("`", "", (as.name(predictors_cfc[i]))))
}
saveRDS(predictors_cfc_toSave, "./models/CB/FTIR/predictors_SVMLM_SVMR_SVMP_LM_PS.rds") #save model as an rds object


######################### select features for the btch A #######################################
# run the RFE algorithm
results_CFC_bA <- rfe(CFC ~ ., data = dfA_CFC,
                      sizes = c(2:length(feature_ID_CFC_bA)), rfeControl = rctrl1)
#sumarize the results
print(results_CFC_bA)
#list the chosen features
predictors_cfc_bA <- predictors(results_CFC_bA)
#plot the results
plot(results_CFC_bA, type = c("g", "o"))
#create a dataframe with selected features
predictors_cfc_id_bA <- c(1)
for(i in 2:ncol(dfA_CFC)){
  for(j in 1:length(predictors_cfc_bA)){
    if( as.numeric(colnames(dfA_CFC[i]))  == as.numeric(gsub("`", "", (as.name(predictors_cfc_bA[j]))))){
      predictors_cfc_id_bA[length(predictors_cfc_id_bA)+1] <- i
    }
  }
}
dfA2_CFC <- dfA_CFC[,predictors_cfc_id_bA]
dfB2_CFC <- dfB_CFC[,predictors_cfc_id_bA]

######################### select features for the batch B #######################################
# run the RFE algorithm
results_CFC_bB <- rfe(CFC ~ ., data = dfB_CFC_B,
                      sizes = c(2:length(feature_ID_CFC_bB)), rfeControl = rctrl1)
#sumarize the results
print(results_CFC_bB)
#list the chosen features
predictors_cfc_bB <- predictors(results_CFC_bB)
#plot the results
plot(results_CFC_bB, type = c("g", "o"))
#create a dataframe with selected features
predictors_cfc_id_bB <- c(1)
for(i in 2:ncol(dfB_CFC_B)){
  for(j in 1:length(predictors_cfc_bB)){
    if( as.numeric(colnames(dfB_CFC_B[i]))  == as.numeric(gsub("`", "", (as.name(predictors_cfc_bB[j]))))){
      predictors_cfc_id_bB[length(predictors_cfc_id_bB)+1] <- i
    }
  }
}
dfA2_CFC_B <- dfA_CFC_B[,predictors_cfc_id_bB]
dfB2_CFC_B <- dfB_CFC_B[,predictors_cfc_id_bB]

#############################################################################################
#################### split batch A ##########################################################
#split data
set.seed(123) #caret package
#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndexA_cfc1 <- createDataPartition(dfA2_CFC$CFC, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainA_cfc1 <- dfA2_CFC[trainIndexA_cfc1,]
testA_cfc1 <- dfA2_CFC[-trainIndexA_cfc1,]

#################### split batch B ##########################################################
#split data
set.seed(123) #caret package
#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndexB_cfc1 <- createDataPartition(dfB2_CFC_B$CFC, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainB_cfc1 <- dfB2_CFC_B[trainIndexB_cfc1,]
testB_cfc1 <- dfB2_CFC_B[-trainIndexB_cfc1,]

#############################################################################################
################################ RF_CFC #####################################################

##################
# the first method
model_Rf_CFC1 <- train(CFC ~ ., data = train_cfc2,
                       method = "rf", trControl = fitControl)
pred_Rf_CFC1<- predict(model_Rf_CFC1,test_cfc2)
RMSE.Rf_CFC1 <- RMSE(test_cfc2$CFC,pred_Rf_CFC1) #RMSE =0.54

# rf plot predicted vs Observed
plot(pred_Rf_CFC1,test_cfc2$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_CFC1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.Rf_CFC.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_cfc3 <- createDataPartition(df2_CFC$CFC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_Rf_CFC <- train(CFC~., data = train_cfc3, method = "rf", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_cfc3$CFC,predicted.Rf_CFC)
  # put results to the vectors
  rmse.Rf_CFC.vector <- c(rmse.Rf_CFC.vector, rmse.Rf_CFC)
}
# compute the standard deviation
sd.Rf_CFC <- round(sd(rmse.Rf_CFC.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.Rf_CFC <- CI(rmse.Rf_CFC.vector, ci = 0.95) #package: Rmisc #0.6317044 0.6079357 0.5841669
# rmse mean
rmse.Rf_CFC.mean <- round(mean(rmse.Rf_CFC.vector), digits = 2) #0.61

plot(rmse.Rf_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_CFC.mean,"+/- ",sd.Rf_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_CFC[3], digits = 2),
                  "to",round(ci.Rf_CFC[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_CFC2 <- train(CFC ~ ., data = dfA2_CFC,
                       method = "rf", trControl = fitControl)
pred_Rf_CFC2<- predict(model_Rf_CFC2,dfB2_CFC)
RMSE.Rf_CFC2 <- RMSE(dfB2_CFC$CFC,pred_Rf_CFC2) #RMSE =0.93

# rf plot predicted vs Observed
plot(pred_Rf_CFC2,dfB2_CFC$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ KNN_CFC #####################################################

##################
# the first method
model_knn_CFC1 <- train(CFC ~ ., data = train_cfc2,
                        method = "knn", trControl = fitControl)
pred_knn_CFC1<- predict(model_knn_CFC1,test_cfc2)
RMSE.knn_CFC1 <- RMSE(test_cfc2$CFC,pred_knn_CFC1) #RMSE =0.8

# rf plot predicted vs Observed
plot(pred_knn_CFC1,test_cfc2$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_CFC1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.knn_CFC.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_cfc3 <- createDataPartition(df2_CFC$CFC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_knn_CFC <- train(CFC~., data = train_cfc3, method = "knn", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.knn_CFC <- predict(model_knn_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.knn_CFC <- RMSE(pred = test_cfc3$CFC,predicted.knn_CFC)
  # put results to the vectors
  rmse.knn_CFC.vector <- c(rmse.knn_CFC.vector, rmse.knn_CFC)
}

# compute the standard deviation
sd.knn_CFC <- round(sd(rmse.knn_CFC.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.knn_CFC <- CI(rmse.knn_CFC.vector, ci = 0.95) # 0.6776927 0.6518339 0.6259751
# rmse mean
rmse.knn_CFC.mean <- round(mean(rmse.knn_CFC.vector), digits = 2) #0.65

plot(rmse.knn_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_CFC.mean,"+/- ",sd.knn_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_CFC[3], digits = 2),
                  "to",round(ci.knn_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_knn_CFC2 <- train(CFC ~ ., data = dfA2_CFC,
                        method = "knn", trControl = fitControl)
pred_knn_CFC2<- predict(model_knn_CFC2,dfB2_CFC)
RMSE.knn_CFC2 <- RMSE(dfB2_CFC$CFC,pred_knn_CFC2) #RMSE =0.66

# rf plot predicted vs Observed
plot(pred_knn_CFC2,dfB2_CFC$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#############################################################################################
################################ svmLM_CFC #####################################################

##################
# the first method
model_SVMLM_CFC1 <- train(CFC ~ ., data = train_cfc2,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_CFC1<- predict(model_SVMLM_CFC1,test_cfc2)
RMSE.SVMLM_CFC1 <- RMSE(test_cfc2$CFC,pred_SVMLM_CFC1) #RMSE =0.44
accuracy_SVMLM_CFC <-accuracy(pred_SVMLM_CFC1, test_cfc2$CFC) #97.37

# this is the best svmLM model for Pseudomonas
saveRDS(model_SVMLM_CFC1, "./models/CB/FTIR/model_SVMLM_PS.rds") #save model as an rds object
saveRDS(test_cfc2$CFC, "./models/CB/FTIR/tested_SVMLM_PS.rds") #save tested data as an rds object
saveRDS(pred_SVMLM_CFC1, "./models/CB/FTIR/predicted_SVMLM_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMLM_CFC1, "./models/CB/FTIR/RMSE_SVMLM_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMLM_CFC, "./models/CB/FTIR/accuracy_SVMLM_PS.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMLM_CFC1,test_cfc2$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("SVMLM_FTIR RMSE:",round(RMSE.SVMLM_CFC1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMLM_CFC.vector <- c()
for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_cfc3 <- createDataPartition(df2_CFC$CFC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMLM_CFC <- train(CFC~., data = train_cfc3, method = "svmLinear", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMLM_CFC <- predict(model_SVMLM_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.SVMLM_CFC <- RMSE(pred = test_cfc3$CFC,predicted.SVMLM_CFC)
  # put results to the vectors
  rmse.SVMLM_CFC.vector <- c(rmse.SVMLM_CFC.vector, rmse.SVMLM_CFC)
}

# compute the standard deviation
sd.SVMLM_CFC <- round(sd(rmse.SVMLM_CFC.vector), digits = 2) #0.07
# find the 95% confidence intervals parameters
ci.SVMLM_CFC <- CI(rmse.SVMLM_CFC.vector, ci = 0.95) # 0.5686331 0.5479186 0.5272042  
# rmse mean
rmse.SVMLM_CFC.mean <- round(mean(rmse.SVMLM_CFC.vector), digits = 2) #0.55

plot(rmse.SVMLM_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_CFC.mean,"+/- ",sd.SVMLM_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_CFC[3], digits = 2),
                  "to",round(ci.SVMLM_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMLM_CFC2 <- train(CFC ~ ., data = dfA2_CFC,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_CFC2<- predict(model_SVMLM_CFC2,dfB2_CFC)
RMSE.SVMLM_CFC2 <- RMSE(dfB2_CFC$CFC,pred_SVMLM_CFC2) #RMSE =0.7

# rf plot predicted vs Observed
plot(pred_SVMLM_CFC2,dfB2_CFC$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMRadial_CFC #####################################################

##################
# the first method
model_SVMR_CFC1 <- train(CFC ~ ., data = train_cfc2,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_CFC1<- predict(model_SVMR_CFC1,test_cfc2)
RMSE.SVMR_CFC1 <- RMSE(test_cfc2$CFC,pred_SVMR_CFC1) #RMSE =0.4
accuracy_SVMR_CFC <- accuracy(pred_SVMR_CFC1, test_cfc2$CFC) #97.37

# this is the best svmR model for Pseudomonas
saveRDS(model_SVMR_CFC1, "./models/CB/FTIR/model_SVMR_PS.rds") #save model as an rds object
saveRDS(test_cfc2$CFC, "./models/CB/FTIR/tested_SVMR_PS.rds") #save tested data as an rds object
saveRDS(pred_SVMR_CFC1, "./models/CB/FTIR/predicted_SVMR_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMR_CFC1, "./models/CB/FTIR/RMSE_SVMR_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMR_CFC, "./models/CB/FTIR/accuracy_SVMR_PS.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMR_CFC1,test_cfc2$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_CFC1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMR_CFC.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_cfc3 <- createDataPartition(df2_CFC$CFC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMR_CFC <- train(CFC~., data = train_cfc3, method = "svmRadial", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMR_CFC <- predict(model_SVMR_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.SVMR_CFC <- RMSE(pred = test_cfc3$CFC,predicted.SVMR_CFC)
  # put results to the vectors
  rmse.SVMR_CFC.vector <- c(rmse.SVMR_CFC.vector, rmse.SVMR_CFC)
}

# compute the standard deviation
sd.SVMR_CFC <- round(sd(rmse.SVMR_CFC.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.SVMR_CFC <- CI(rmse.SVMR_CFC.vector, ci = 0.95) # 0.5101954 0.4867048 0.4632142   
# rmse mean
rmse.SVMR_CFC.mean <- round(mean(rmse.SVMR_CFC.vector), digits = 2) #0.49

plot(rmse.SVMR_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_CFC.mean,"+/- ",sd.SVMR_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_CFC[3], digits = 2),
                  "to",round(ci.SVMR_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMR_CFC2 <- train(CFC ~ ., data = dfA2_CFC,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_CFC2<- predict(model_SVMR_CFC2,dfB2_CFC)
RMSE.SVMR_CFC2 <- RMSE(dfB2_CFC$CFC,pred_SVMR_CFC2) #RMSE =0.62

# rf plot predicted vs Observed
plot(pred_SVMR_CFC2,dfB2_CFC$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMPoly_CFC #####################################################

##################
# the first method
model_SVMP_CFC1 <- train(CFC ~ ., data = train_cfc2,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC1<- predict(model_SVMP_CFC1,test_cfc2)
RMSE.SVMP_CFC1 <- RMSE(test_cfc2$CFC,pred_SVMP_CFC1) #RMSE =0.4
accuracy_SVMP_CFC <- accuracy(pred_SVMP_CFC1, test_cfc2$CFC) #100

# this is the best svmP model for Pseudomonas
saveRDS(model_SVMP_CFC1, "./models/CB/FTIR/model_SVMP_PS.rds") #save model as an rds object
saveRDS(test_cfc2$CFC, "./models/CB/FTIR/tested_SVMP_PS.rds") #save tested data as an rds object
saveRDS(pred_SVMP_CFC1, "./models/CB/FTIR/predicted_SVMP_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMP_CFC1, "./models/CB/FTIR/RMSE_SVMP_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMP_CFC, "./models/CB/FTIR/accuracy_SVMP_PS.rds") #save accuracy as an rds object


# rf plot predicted vs Observed
plot(pred_SVMP_CFC1,test_cfc2$CFC,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC1,digits = 2), 
                "\nAccuracy",accuracy_SVMP_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMP_CFC.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_cfc3 <- createDataPartition(df2_CFC$CFC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMP_CFC <- train(CFC~., data = train_cfc3, method = "svmPoly", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMP_CFC <- predict(model_SVMP_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.SVMP_CFC <- RMSE(pred = test_cfc3$CFC,predicted.SVMP_CFC)
  # put results to the vectors
  rmse.SVMP_CFC.vector <- c(rmse.SVMP_CFC.vector, rmse.SVMP_CFC)
}

# compute the standard deviation
sd.SVMP_CFC <- round(sd(rmse.SVMP_CFC.vector), digits = 2) #0.06
# find the 95% confidence intervals parameters
ci.SVMP_CFC <- CI(rmse.SVMP_CFC.vector, ci = 0.95) # 0.4415299 0.4249149 0.4083000   
# rmse mean
rmse.SVMP_CFC.mean <- round(mean(rmse.SVMP_CFC.vector), digits = 2) #0.42

plot(rmse.SVMP_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_CFC.mean,"+/- ",sd.SVMP_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_CFC[3], digits = 2),
                  "to",round(ci.SVMP_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMP_CFC2 <- train(CFC ~ ., data = dfA2_CFC,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC2<- predict(model_SVMP_CFC2,dfB2_CFC)
RMSE.SVMP_CFC2 <- RMSE(dfB2_CFC$CFC,pred_SVMP_CFC2) #RMSE =0.61

# rf plot predicted vs Observed
plot(pred_SVMP_CFC2,dfB2_CFC$CFC,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

###################
# trained based on B, tested based on A
model_SVMP_CFC3 <- train(CFC ~ ., data = dfB2_CFC_B,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC3<- predict(model_SVMP_CFC3,dfA2_CFC_B)
RMSE.SVMP_CFC3 <- RMSE(dfA2_CFC_B$CFC,pred_SVMP_CFC3) 

# rf plot predicted vs Observed
plot(pred_SVMP_CFC3,dfA2_CFC_B$CFC,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC3,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# trained and tested based on batch A
model_SVMP_CFC4 <- train(CFC ~ ., data = trainA_cfc1,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC4<- predict(model_SVMP_CFC4,testA_cfc1)
RMSE.SVMP_CFC4 <- RMSE(testA_cfc1$CFC,pred_SVMP_CFC4) 

# plot predicted vs Observed
plot(pred_SVMP_CFC4,testA_cfc1$CFC,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC4,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#######################
# trained and tested based on batch B
model_SVMP_CFC5 <- train(CFC ~ ., data = trainB_cfc1,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC5<- predict(model_SVMP_CFC5,testB_cfc1)
RMSE.SVMP_CFC5 <- RMSE(testB_cfc1$CFC,pred_SVMP_CFC5) 

# plot predicted vs Observed
plot(pred_SVMP_CFC5,testB_cfc1$CFC,xlim= c(0,10), ylim=c(0,10),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC5,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#############################################################################################
################################ LM_CFC #####################################################

##################
# the first method
model_lm_CFC1 <- train(CFC ~ ., data = train_cfc2,
                       method = "lm", trControl = fitControl)
pred_lm_CFC1<- predict(model_lm_CFC1,test_cfc2)
RMSE.lm_CFC1 <- RMSE(test_cfc2$CFC,pred_lm_CFC1) #RMSE =0.47
accuracy_LM_CFC <- accuracy(pred_lm_CFC1, test_cfc2$CFC) #97.37

# this is the best LM model for Pseudomonas
saveRDS(model_lm_CFC1, "./models/CB/FTIR/model_LM_PS.rds") #save model as an rds object
saveRDS(test_cfc2$CFC, "./models/CB/FTIR/tested_LM_PS.rds") #save tested data as an rds object
saveRDS(pred_lm_CFC1, "./models/CB/FTIR/predicted_LM_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.lm_CFC1, "./models/CB/FTIR/RMSE_LM_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_LM_CFC, "./models/CB/FTIR/accuracy_LM_PS.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_lm_CFC1,test_cfc2$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_CFC1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#saveRDS(model_lm_TSA1, "./models/CTF/FTIR/model_LM_TVC.rds") #save model as an rds object

##################
# 50 iterations
rmse.lm_CFC.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_cfc3 <- createDataPartition(df2_CFC$CFC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_lm_CFC <- train(CFC~., data = train_cfc3, method = "lm", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.lm_CFC <- predict(model_lm_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.lm_CFC <- RMSE(pred = test_cfc3$CFC,predicted.lm_CFC)
  # put results to the vectors
  rmse.lm_CFC.vector <- c(rmse.lm_CFC.vector, rmse.lm_CFC)
}
# compute the standard deviation
sd.lm_CFC <- round(sd(rmse.lm_CFC.vector), digits = 2) #0.07
# find the 95% confidence intervals parameters
ci.lm_CFC <- CI(rmse.lm_CFC.vector, ci = 0.95) #package: Rmisc #0.57, 0.55, 0.53
# rmse mean
rmse.lm_CFC.mean <- round(mean(rmse.lm_CFC.vector), digits = 2) #0.55

plot(rmse.lm_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_CFC.mean,"+/- ",sd.lm_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.lm_CFC[3], digits = 2),
                  "to",round(ci.lm_CFC[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_lm_CFC2 <- train(CFC ~ ., data = dfA2_CFC,
                       method = "lm", trControl = fitControl)
pred_lm_CFC2<- predict(model_lm_CFC2,dfB2_CFC)
RMSE.lm_CFC2 <- RMSE(dfB2_CFC$CFC,pred_lm_CFC2) #RMSE =0.74

# rf plot predicted vs Observed
plot(pred_lm_CFC2,dfB2_CFC$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)



######################################################################################################################
###############################################################################################################
################################### Prediction for STAA (Brochothrix thermosphacta)########################################################
###############################################################################################################

######################### select features for the whole dataset #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_staa <- train(STAA~., data = df_STAA_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_STAA <- varImp(model_varImp_staa, scale = FALSE)
#sumarize importance
print(importance_STAA)
#plot importance
plot(importance_STAA)
# get an index of important variables
feature_ID_STAA <- c(1)
for(i in 1:nrow(importance_STAA$importance)){
  if(importance_STAA$importance[i,1] >= 0.5){
    feature_ID_STAA[length(feature_ID_STAA)+1] <- i+1
  }
}
# create a dataframe eith selected variables
df_STAA <- df_STAA_wf[,feature_ID_STAA]

######################### select features for the batch A #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_staa_bA <- train(STAA~., data = dfA_STAA_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_bA_STAA <- varImp(model_varImp_staa_bA, scale = FALSE)
#sumarize importance
print(importance_bA_STAA)
#plot importance
plot(importance_bA_STAA)
# get an index of important variables
feature_ID_STAA_bA <- c(1)
for(i in 1:nrow(importance_bA_STAA$importance)){
  if(importance_bA_STAA$importance[i,1] >= 0.5){
    feature_ID_STAA_bA[length(feature_ID_STAA_bA)+1] <- i+1
  }
}
# create a dataframe eith selected variables
dfA_STAA <- dfA_STAA_wf[,feature_ID_STAA_bA]
dfB_STAA <- dfB_STAA_wf[,feature_ID_STAA_bA]

###############################################################################################################
################################### Prediction for STAA ########################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_staa1 <- createDataPartition(df_STAA$STAA, p = .7, 
                                       list = FALSE, 
                                       times = 1, groups =3) #caret package
#split the whole dataset
train_staa1 <- df_STAA[trainIndex_staa1,]
test_staa1 <- df_STAA[-trainIndex_staa1,]

# resembling (method to avoid overfitting)
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)
## repeated ten times

#############################################################################################
############################# feature selection step 3 ################################################
############################# feature selection #############################################

# rfe control (define the control using a random forest selection function)
rctrl1 <- rfeControl(method = "cv",
                     number = 10, #10-fold cv
                     functions = rfFuncs)

######################### select features for the whole dataset #######################################
# run the RFE algorithm
results_STAA <- rfe(STAA ~ ., data = train_staa1,
                   sizes = c(2:ncol(train_staa1)), rfeControl = rctrl1)
#sumarize the results
print(results_STAA)
#list the chosen features
predictors_staa <- predictors(results_STAA)
#plot the results
plot(results_STAA, type = c("g", "o"))
#create a dataframe with selected features
predictors_staa_id <- c(1)
for(i in 2:ncol(train_staa1)){
  for(j in 1:length(predictors_staa)){
    if( as.numeric(colnames(train_staa1[i]))  == as.numeric(gsub("`", "", (as.name(predictors_staa[j]))))){
      predictors_staa_id[length(predictors_staa_id)+1] <- i
    }
  }
}
train_staa2 <- train_staa1[,predictors_staa_id]
test_staa2 <- test_staa1[,predictors_staa_id]
df2_STAA <- df_STAA[,predictors_staa_id]

#save selected features in one vector
predictors_staa_toSave <- c()
len_predictors_staa <- 0
for(i in predictors_staa){
  if(as.numeric(lengths(predictors_staa[i]))== 1){
    len_predictors_staa <- len_predictors_staa +1
  }
}
for(i in 1:len_predictors_staa){
  predictors_staa_toSave[i] <- as.numeric(gsub("`", "", (as.name(predictors_staa[i]))))
}
saveRDS(predictors_staa_toSave, "./models/CB/FTIR/predictors_knn_SVMLM_SVMP_LM_BRO.rds") #save model as an rds object


######################### select features for the btch A #######################################
# run the RFE algorithm
results_STAA_bA <- rfe(STAA ~ ., data = dfA_STAA,
                      sizes = c(2:length(feature_ID_STAA_bA)), rfeControl = rctrl1)
#sumarize the results
print(results_STAA_bA)
#list the chosen features
predictors_staa_bA <- predictors(results_STAA_bA)
#plot the results
plot(results_STAA_bA, type = c("g", "o"))
#create a dataframe with selected features
predictors_staa_id_bA <- c(1)
for(i in 2:ncol(dfA_STAA)){
  for(j in 1:length(predictors_staa_bA)){
    if( as.numeric(colnames(dfA_STAA[i]))  == as.numeric(gsub("`", "", (as.name(predictors_staa_bA[j]))))){
      predictors_staa_id_bA[length(predictors_staa_id_bA)+1] <- i
    }
  }
}
dfA2_STAA <- dfA_STAA[,predictors_staa_id_bA]
dfB2_STAA <- dfB_STAA[,predictors_staa_id_bA]

#############################################################################################
################################ RF_STAA #####################################################

##################
# the first method
model_Rf_STAA1 <- train(STAA ~ ., data = train_staa2,
                       method = "rf", trControl = fitControl)
pred_Rf_STAA1<- predict(model_Rf_STAA1,test_staa2)
RMSE.Rf_STAA1 <- RMSE(test_staa2$STAA,pred_Rf_STAA1) #RMSE =0.69

# rf plot predicted vs Observed
plot(pred_Rf_STAA1,test_staa2$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_STAA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.Rf_STAA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_staa3 <- createDataPartition(df2_STAA$STAA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_staa3 <- df2_STAA[trainIndex_staa3,]
  test_staa3 <- df2_STAA[-trainIndex_staa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_Rf_STAA <- train(STAA~., data = train_staa3, method = "rf", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_STAA <- predict(model_Rf_STAA,test_staa3)
  
  # calculate the RMSE value
  rmse.Rf_STAA <- RMSE(pred = test_staa3$STAA,predicted.Rf_STAA)
  # put results to the vectors
  rmse.Rf_STAA.vector <- c(rmse.Rf_STAA.vector, rmse.Rf_STAA)
}
# compute the standard deviation
sd.Rf_STAA <- round(sd(rmse.Rf_STAA.vector), digits = 2) #0.1
# find the 95% confidence intervals parameters
ci.Rf_STAA <- CI(rmse.Rf_STAA.vector, ci = 0.95) #package: Rmisc #0.6288063 0.6012709 0.5737354
# rmse mean
rmse.Rf_STAA.mean <- round(mean(rmse.Rf_STAA.vector), digits = 2) #0.6

plot(rmse.Rf_STAA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_STAA.mean,"+/- ",sd.Rf_STAA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_STAA[3], digits = 2),
                  "to",round(ci.Rf_STAA[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_STAA2 <- train(STAA ~ ., data = dfA2_STAA,
                       method = "rf", trControl = fitControl)
pred_Rf_STAA2<- predict(model_Rf_STAA2,dfB2_STAA)
RMSE.Rf_STAA2 <- RMSE(dfB2_STAA$STAA,pred_Rf_STAA2) #RMSE =0.99

# rf plot predicted vs Observed
plot(pred_Rf_STAA2,dfB2_STAA$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_STAA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ KNN_STAA #####################################################

##################
# the first method
model_knn_STAA1 <- train(STAA ~ ., data = train_staa2,
                        method = "knn", trControl = fitControl)
pred_knn_STAA1<- predict(model_knn_STAA1,test_staa2)
RMSE.knn_STAA1 <- RMSE(test_staa2$STAA,pred_knn_STAA1) #RMSE =0.68
accuracy_knn_STAA <- accuracy(pred_knn_STAA1, test_staa2$STAA) #87.47

# this is the best knn model for Brochothrix
saveRDS(model_knn_STAA1, "./models/CB/FTIR/model_knn_BRO.rds") #save model as an rds object
saveRDS(test_staa2$STAA, "./models/CB/FTIR/tested_knn_BRO.rds") #save tested data as an rds object
saveRDS(pred_knn_STAA1, "./models/CB/FTIR/predicted_knn_BRO.rds") #save predicted data as an rds object
saveRDS(RMSE.knn_STAA1, "./models/CB/FTIR/RMSE_knn_BRO.rds") #save RMSE as an rds object
saveRDS(accuracy_knn_STAA, "./models/CB/FTIR/accuracy_knn_BRO.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_knn_STAA1,test_staa2$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_STAA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.knn_STAA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_staa3 <- createDataPartition(df2_STAA$STAA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_staa3 <- df2_STAA[trainIndex_staa3,]
  test_staa3 <- df2_STAA[-trainIndex_staa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_knn_STAA <- train(STAA~., data = train_staa3, method = "knn", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.knn_STAA <- predict(model_knn_STAA,test_staa3)
  
  # calculate the RMSE value
  rmse.knn_STAA <- RMSE(pred = test_staa3$STAA,predicted.knn_STAA)
  # put results to the vectors
  rmse.knn_STAA.vector <- c(rmse.knn_STAA.vector, rmse.knn_STAA)
}

# compute the standard deviation
sd.knn_STAA <- round(sd(rmse.knn_STAA.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.knn_STAA <- CI(rmse.knn_STAA.vector, ci = 0.95) # 0.6203579 0.5981382 0.5759185 
# rmse mean
rmse.knn_STAA.mean <- round(mean(rmse.knn_STAA.vector), digits = 2) #0.6

plot(rmse.knn_STAA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_STAA.mean,"+/- ",sd.knn_STAA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_STAA[3], digits = 2),
                  "to",round(ci.knn_STAA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_knn_STAA2 <- train(STAA ~ ., data = dfA2_STAA,
                        method = "knn", trControl = fitControl)
pred_knn_STAA2<- predict(model_knn_STAA2,dfB2_STAA)
RMSE.knn_STAA2 <- RMSE(dfB2_STAA$STAA,pred_knn_STAA2) #RMSE =0.88

# rf plot predicted vs Observed
plot(pred_knn_STAA2,dfB2_STAA$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_STAA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#############################################################################################
################################ svmLM_STAA #####################################################

##################
# the first method
model_SVMLM_STAA1 <- train(STAA ~ ., data = train_staa2,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_STAA1<- predict(model_SVMLM_STAA1,test_staa2)
RMSE.SVMLM_STAA1 <- RMSE(test_staa2$STAA,pred_SVMLM_STAA1) #RMSE =0.68
accuracy_SVMLM_STAA <- accuracy(pred_SVMLM_STAA1, test_staa2$STAA) #86.84

# this is the best svmLM model for Brochothrix
saveRDS(model_SVMLM_STAA1, "./models/CB/FTIR/model_SVMLM_BRO.rds") #save model as an rds object
saveRDS(test_staa2$STAA, "./models/CB/FTIR/tested_SVMLM_BRO.rds") #save tested data as an rds object
saveRDS(pred_SVMLM_STAA1, "./models/CB/FTIR/predicted_SVMLM_BRO.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMLM_STAA1, "./models/CB/FTIR/RMSE_SVMLM_BRO.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMLM_STAA, "./models/CB/FTIR/accuracy_SVMLM_BRO.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMLM_STAA1,test_staa2$CFC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("SVMLM_FTIR RMSE:",round(RMSE.SVMLM_STAA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMLM_STAA.vector <- c()
for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_staa3 <- createDataPartition(df2_STAA$STAA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_staa3 <- df2_STAA[trainIndex_staa3,]
  test_staa3 <- df2_STAA[-trainIndex_staa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMLM_STAA <- train(STAA~., data = train_staa3, method = "svmLinear", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMLM_STAA <- predict(model_SVMLM_STAA,test_staa3)
  
  # calculate the RMSE value
  rmse.SVMLM_STAA <- RMSE(pred = test_staa3$STAA,predicted.SVMLM_STAA)
  # put results to the vectors
  rmse.SVMLM_STAA.vector <- c(rmse.SVMLM_STAA.vector, rmse.SVMLM_STAA)
}

# compute the standard deviation
sd.SVMLM_STAA <- round(sd(rmse.SVMLM_STAA.vector), digits = 2) #0.06
# find the 95% confidence intervals parameters
ci.SVMLM_STAA <- CI(rmse.SVMLM_STAA.vector, ci = 0.95) # 0.5706869 0.5525666 0.5344464   
# rmse mean
rmse.SVMLM_STAA.mean <- round(mean(rmse.SVMLM_STAA.vector), digits = 2) #0.55

plot(rmse.SVMLM_STAA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_STAA.mean,"+/- ",sd.SVMLM_STAA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_STAA[3], digits = 2),
                  "to",round(ci.SVMLM_STAA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMLM_STAA2 <- train(STAA ~ ., data = dfA2_STAA,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_STAA2<- predict(model_SVMLM_STAA2,dfB2_STAA)
RMSE.SVMLM_STAA2 <- RMSE(dfB2_STAA$STAA,pred_SVMLM_STAA2) #RMSE =1.1

# rf plot predicted vs Observed
plot(pred_SVMLM_STAA2,dfB2_STAA$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_STAA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMRadial_STAA #####################################################

##################
# the first method
model_SVMR_STAA1 <- train(STAA ~ ., data = train_staa2,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_STAA1<- predict(model_SVMR_STAA1,test_staa2)
RMSE.SVMR_STAA1 <- RMSE(test_staa2$STAA,pred_SVMR_STAA1) #RMSE =0.64

# rf plot predicted vs Observed
plot(pred_SVMR_STAA1,test_staa2$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_STAA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMR_STAA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_staa3 <- createDataPartition(df2_STAA$STAA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_staa3 <- df2_STAA[trainIndex_staa3,]
  test_staa3 <- df2_STAA[-trainIndex_staa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMR_STAA <- train(STAA~., data = train_staa3, method = "svmRadial", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMR_STAA <- predict(model_SVMR_STAA,test_staa3)
  
  # calculate the RMSE value
  rmse.SVMR_STAA <- RMSE(pred = test_staa3$STAA,predicted.SVMR_STAA)
  # put results to the vectors
  rmse.SVMR_STAA.vector <- c(rmse.SVMR_STAA.vector, rmse.SVMR_STAA)
}

# compute the standard deviation
sd.SVMR_STAA <- round(sd(rmse.SVMR_STAA.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMR_STAA <- CI(rmse.SVMR_STAA.vector, ci = 0.95) # 0.5353003 0.5107610 0.4862217 
# rmse mean
rmse.SVMR_STAA.mean <- round(mean(rmse.SVMR_STAA.vector), digits = 2) #0.51

plot(rmse.SVMR_STAA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_STAA.mean,"+/- ",sd.SVMR_STAA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_STAA[3], digits = 2),
                  "to",round(ci.SVMR_STAA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMR_STAA2 <- train(STAA ~ ., data = dfA2_STAA,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_STAA2<- predict(model_SVMR_STAA2,dfB2_STAA)
RMSE.SVMR_STAA2 <- RMSE(dfB2_STAA$STAA,pred_SVMR_STAA2) #RMSE =0.66

# rf plot predicted vs Observed
plot(pred_SVMR_STAA2,dfB2_STAA$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_STAA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMPoly_STAA #####################################################

##################
# the first method
model_SVMP_STAA1 <- train(STAA ~ ., data = train_staa2,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_STAA1<- predict(model_SVMP_STAA1,test_staa2)
RMSE.SVMP_STAA1 <- RMSE(test_staa2$STAA,pred_SVMP_STAA1) #RMSE =0.56
accuracy_SVMP_STAA <- accuracy(pred_SVMP_STAA1, test_staa2$STAA) #94.74

# this is the best svmP model for Brochothrix
saveRDS(model_SVMP_STAA1, "./models/CB/FTIR/model_SVMP_BRO.rds") #save model as an rds object
saveRDS(test_staa2$STAA, "./models/CB/FTIR/tested_SVMP_BRO.rds") #save tested data as an rds object
saveRDS(pred_SVMP_STAA1, "./models/CB/FTIR/predicted_SVMP_BRO.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMP_STAA1, "./models/CB/FTIR/RMSE_SVMP_BRO.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMP_STAA, "./models/CB/FTIR/accuracy_SVMP_BRO.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMP_STAA1,test_staa2$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_STAA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMP_STAA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_staa3 <- createDataPartition(df2_STAA$STAA, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_staa3 <- df2_STAA[trainIndex_staa3,]
  test_staa3 <- df2_STAA[-trainIndex_staa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMP_STAA <- train(STAA~., data = train_staa3, method = "svmPoly", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMP_STAA <- predict(model_SVMP_STAA,test_staa3)
  
  # calculate the RMSE value
  rmse.SVMP_STAA <- RMSE(pred = test_staa3$STAA,predicted.SVMP_STAA)
  # put results to the vectors
  rmse.SVMP_STAA.vector <- c(rmse.SVMP_STAA.vector, rmse.SVMP_STAA)
}

# compute the standard deviation
sd.SVMP_STAA <- round(sd(rmse.SVMP_STAA.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMP_STAA <- CI(rmse.SVMP_STAA.vector, ci = 0.95) #0.5666778 0.5403788 0.5140798  
# rmse mean
rmse.SVMP_STAA.mean <- round(mean(rmse.SVMP_STAA.vector), digits = 2) #0.54

plot(rmse.SVMP_STAA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_STAA.mean,"+/- ",sd.SVMP_STAA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_STAA[3], digits = 2),
                  "to",round(ci.SVMP_STAA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMP_STAA2 <- train(STAA ~ ., data = dfA2_STAA,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_STAA2<- predict(model_SVMP_STAA2,dfB2_STAA)
RMSE.SVMP_STAA2 <- RMSE(dfB2_STAA$STAA,pred_SVMP_STAA2) #RMSE =0.66

# rf plot predicted vs Observed
plot(pred_SVMP_STAA2,dfB2_STAA$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_STAA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ LM_STAA #####################################################

##################
# the first method
model_lm_STAA1 <- train(STAA ~ ., data = train_staa2,
                        method = "lm", trControl = fitControl)
pred_lm_STAA1<- predict(model_lm_STAA1,test_staa2)
RMSE.lm_STAA1 <- RMSE(test_staa2$STAA,pred_lm_STAA1) #RMSE =0.67
accuracy_LM_STAA <- accuracy(pred_lm_STAA1, test_staa2$STAA) #86.84

# this is the best LM model for Brochothrix
saveRDS(model_lm_STAA1, "./models/CB/FTIR/model_LM_BRO.rds") #save model as an rds object
saveRDS(test_staa2$STAA, "./models/CB/FTIR/tested_LM_BRO.rds") #save tested data as an rds object
saveRDS(pred_lm_STAA1, "./models/CB/FTIR/predicted_LM_BRO.rds") #save predicted data as an rds object
saveRDS(RMSE.lm_STAA1, "./models/CB/FTIR/RMSE_LM_BRO.rds") #save RMSE as an rds object
saveRDS(accuracy_LM_STAA, "./models/CB/FTIR/accuracy_LM_BRO.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_lm_STAA1,test_staa2$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_STAA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.lm_STAA.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_staa3 <- createDataPartition(df2_STAA$STAA, p = .7, 
                                          list = FALSE, 
                                          times = 1, groups =3) #caret package
  #split the whole dataset
  train_staa3 <- df2_STAA[trainIndex_staa3,]
  test_staa3 <- df2_STAA[-trainIndex_staa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_lm_STAA <- train(STAA~., data = train_staa3, method = "lm", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.lm_STAA <- predict(model_lm_STAA,test_staa3)
  
  # calculate the RMSE value
  rmse.lm_STAA <- RMSE(pred = test_staa3$STAA,predicted.lm_STAA)
  # put results to the vectors
  rmse.lm_STAA.vector <- c(rmse.lm_STAA.vector, rmse.lm_STAA)
}
# compute the standard deviation
sd.lm_STAA <- round(sd(rmse.lm_STAA.vector), digits = 2) #0.07
# find the 95% confidence intervals parameters
ci.lm_STAA <- CI(rmse.lm_STAA.vector, ci = 0.95) #package: Rmisc #0.58, 0.56, 0.54
# rmse mean
rmse.lm_STAA.mean <- round(mean(rmse.lm_STAA.vector), digits = 2) #0.56

plot(rmse.lm_STAA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_STAA.mean,"+/- ",sd.lm_STAA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.lm_STAA[3], digits = 2),
                  "to",round(ci.lm_STAA[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_lm_STAA2 <- train(STAA ~ ., data = dfA2_STAA,
                        method = "lm", trControl = fitControl)
pred_lm_STAA2<- predict(model_lm_STAA2,dfB2_STAA)
RMSE.lm_STAA2 <- RMSE(dfB2_STAA$STAA,pred_lm_STAA2) #RMSE =2.19

# rf plot predicted vs Observed
plot(pred_lm_STAA2,dfB2_STAA$STAA,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_STAA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


################################################################################################################
###############################################################################################################
################################### Prediction for MRS (LAB)########################################################
###############################################################################################################

######################### select features for the whole dataset #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_mrs <- train(MRS~., data = df_MRS_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_MRS <- varImp(model_varImp_mrs, scale = FALSE)
#sumarize importance
print(importance_MRS)
#plot importance
plot(importance_MRS)
# get an index of important variables
feature_ID_MRS <- c(1)
for(i in 1:nrow(importance_MRS$importance)){
  if(importance_MRS$importance[i,1] >= 0.5){
    feature_ID_MRS[length(feature_ID_MRS)+1] <- i+1
  }
}
# create a dataframe eith selected variables
df_MRS <- df_MRS_wf[,feature_ID_MRS]

######################### select features for the batch A #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_mrs_bA <- train(MRS~., data = dfA_MRS_wf, method = "rf", trControl = control, importance = TRUE)
#estimate viable importance
importance_bA_MRS <- varImp(model_varImp_mrs_bA, scale = FALSE)
#sumarize importance
print(importance_bA_MRS)
#plot importance
plot(importance_bA_MRS)
# get an index of important variables
feature_ID_MRS_bA <- c(1)
for(i in 1:nrow(importance_bA_MRS$importance)){
  if(importance_bA_MRS$importance[i,1] >= 0.5){
    feature_ID_MRS_bA[length(feature_ID_MRS_bA)+1] <- i+1
  }
}
# create a dataframe eith selected variables
dfA_MRS <- dfA_MRS_wf[,feature_ID_MRS_bA]
dfB_MRS <- dfB_MRS_wf[,feature_ID_MRS_bA]

###############################################################################################################
################################### Prediction for MRS ########################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_mrs1 <- createDataPartition(df_MRS$MRS, p = .7, 
                                        list = FALSE, 
                                        times = 1, groups =3) #caret package
#split the whole dataset
train_mrs1 <- df_MRS[trainIndex_mrs1,]
test_mrs1 <- df_MRS[-trainIndex_mrs1,]

# resembling (method to avoid overfitting)
fitControl <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)
## repeated ten times

#############################################################################################
############################# feature selection step 3 ################################################
############################# feature selection #############################################

# rfe control (define the control using a random forest selection function)
rctrl1 <- rfeControl(method = "cv",
                     number = 10, #10-fold cv
                     functions = rfFuncs)

######################### select features for the whole dataset #######################################
# run the RFE algorithm
results_MRS <- rfe(MRS ~ ., data = train_mrs1,
                    sizes = c(2:ncol(train_mrs1)), rfeControl = rctrl1)
#sumarize the results
print(results_MRS)
#list the chosen features
predictors_mrs <- predictors(results_MRS)
#plot the results
plot(results_MRS, type = c("g", "o"))
#create a dataframe with selected features
predictors_mrs_id <- c(1)
for(i in 2:ncol(train_mrs1)){
  for(j in 1:length(predictors_mrs)){
    if( as.numeric(colnames(train_mrs1[i]))  == as.numeric(gsub("`", "", (as.name(predictors_mrs[j]))))){
      predictors_mrs_id[length(predictors_mrs_id)+1] <- i
    }
  }
}
train_mrs2 <- train_mrs1[,predictors_mrs_id]
test_mrs2 <- test_mrs1[,predictors_mrs_id]
df2_MRS <- df_MRS[,predictors_mrs_id]

#save selected features in one vector
predictors_mrs_toSave <- c()
len_predictors_mrs <- 0
for(i in predictors_mrs){
  if(as.numeric(lengths(predictors_mrs[i]))== 1){
    len_predictors_mrs <- len_predictors_mrs +1
  }
}
for(i in 1:len_predictors_mrs){
  predictors_mrs_toSave[i] <- as.numeric(gsub("`", "", (as.name(predictors_mrs[i]))))
}
saveRDS(predictors_mrs_toSave, "./models/CB/FTIR/predictors_SVMR_SVMP_LAB.rds") #save model as an rds object


######################### select features for the batch A #######################################
# run the RFE algorithm
results_MRS_bA <- rfe(MRS ~ ., data = dfA_MRS,
                       sizes = c(2:length(feature_ID_MRS_bA)), rfeControl = rctrl1)
#sumarize the results
print(results_MRS_bA)
#list the chosen features
predictors_mrs_bA <- predictors(results_MRS_bA)
#plot the results
plot(results_MRS_bA, type = c("g", "o"))
#create a dataframe with selected features
predictors_mrs_id_bA <- c(1)
for(i in 2:ncol(dfA_MRS)){
  for(j in 1:length(predictors_mrs_bA)){
    if( as.numeric(colnames(dfA_MRS[i]))  == as.numeric(gsub("`", "", (as.name(predictors_mrs_bA[j]))))){
      predictors_mrs_id_bA[length(predictors_mrs_id_bA)+1] <- i
    }
  }
}
dfA2_MRS <- dfA_MRS[,predictors_mrs_id_bA]
dfB2_MRS <- dfB_MRS[,predictors_mrs_id_bA]

#############################################################################################
################################ RF_MRS #####################################################

##################
# the first method
model_Rf_MRS1 <- train(MRS ~ ., data = train_mrs2,
                        method = "rf", trControl = fitControl)
pred_Rf_MRS1<- predict(model_Rf_MRS1,test_mrs2)
RMSE.Rf_MRS1 <- RMSE(test_mrs2$MRS,pred_Rf_MRS1) #RMSE =0.73

# rf plot predicted vs Observed
plot(pred_Rf_MRS1,test_mrs2$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_MRS1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.Rf_MRS.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_mrs3 <- createDataPartition(df2_MRS$MRS, p = .7, 
                                          list = FALSE, 
                                          times = 1, groups =3) #caret package
  #split the whole dataset
  train_mrs3 <- df2_MRS[trainIndex_mrs3,]
  test_mrs3 <- df2_MRS[-trainIndex_mrs3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_Rf_MRS <- train(MRS~., data = train_mrs3, method = "rf", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_MRS <- predict(model_Rf_MRS,test_mrs3)
  
  # calculate the RMSE value
  rmse.Rf_MRS <- RMSE(pred = test_mrs3$MRS,predicted.Rf_MRS)
  # put results to the vectors
  rmse.Rf_MRS.vector <- c(rmse.Rf_MRS.vector, rmse.Rf_MRS)
}
# compute the standard deviation
sd.Rf_MRS <- round(sd(rmse.Rf_MRS.vector), digits = 2) #0.11
# find the 95% confidence intervals parameters
ci.Rf_MRS <- CI(rmse.Rf_MRS.vector, ci = 0.95) #package: Rmisc #0.7163602 0.6862626 0.6561649 
# rmse mean
rmse.Rf_MRS.mean <- round(mean(rmse.Rf_MRS.vector), digits = 2) #0.69

plot(rmse.Rf_MRS.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_MRS.mean,"+/- ",sd.Rf_MRS, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_MRS[3], digits = 2),
                  "to",round(ci.Rf_MRS[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_MRS2 <- train(MRS ~ ., data = dfA2_MRS,
                        method = "rf", trControl = fitControl)
pred_Rf_MRS2<- predict(model_Rf_MRS2,dfB2_MRS)
RMSE.Rf_MRS2 <- RMSE(dfB2_MRS$MRS,pred_Rf_MRS2) #RMSE =1.2

# rf plot predicted vs Observed
plot(pred_Rf_MRS2,dfB2_MRS$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_MRS2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ KNN_MRS #####################################################

##################
# the first method
model_knn_MRS1 <- train(MRS ~ ., data = train_mrs2,
                         method = "knn", trControl = fitControl)
pred_knn_MRS1<- predict(model_knn_MRS1,test_mrs2)
RMSE.knn_MRS1 <- RMSE(test_mrs2$MRS,pred_knn_MRS1) #RMSE =0.81

# rf plot predicted vs Observed
plot(pred_knn_MRS1,test_mrs2$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_MRS1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.knn_MRS.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_mrs3 <- createDataPartition(df2_MRS$MRS, p = .7, 
                                          list = FALSE, 
                                          times = 1, groups =3) #caret package
  #split the whole dataset
  train_mrs3 <- df2_MRS[trainIndex_mrs3,]
  test_mrs3 <- df2_MRS[-trainIndex_mrs3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_knn_MRS <- train(MRS~., data = train_mrs3, method = "knn", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.knn_MRS <- predict(model_knn_MRS,test_mrs3)
  
  # calculate the RMSE value
  rmse.knn_MRS <- RMSE(pred = test_mrs3$MRS,predicted.knn_MRS)
  # put results to the vectors
  rmse.knn_MRS.vector <- c(rmse.knn_MRS.vector, rmse.knn_MRS)
}

# compute the standard deviation
sd.knn_MRS <- round(sd(rmse.knn_MRS.vector), digits = 2) #0.11
# find the 95% confidence intervals parameters
ci.knn_MRS <- CI(rmse.knn_MRS.vector, ci = 0.95) # 0.7911262 0.7599592 0.7287921
# rmse mean
rmse.knn_MRS.mean <- round(mean(rmse.knn_MRS.vector), digits = 2) #0.76

plot(rmse.knn_MRS.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_MRS.mean,"+/- ",sd.knn_MRS, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_MRS[3], digits = 2),
                  "to",round(ci.knn_MRS[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_knn_MRS2 <- train(MRS ~ ., data = dfA2_MRS,
                         method = "knn", trControl = fitControl)
pred_knn_MRS2<- predict(model_knn_MRS2,dfB2_MRS)
RMSE.knn_MRS2 <- RMSE(dfB2_MRS$MRS,pred_knn_MRS2) #RMSE =1.24

# rf plot predicted vs Observed
plot(pred_knn_MRS2,dfB2_MRS$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_MRS2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#############################################################################################
################################ svmLM_MRS #####################################################

##################
# the first method
model_SVMLM_MRS1 <- train(MRS ~ ., data = train_mrs2,
                           method = "svmLinear", trControl = fitControl)
pred_SVMLM_MRS1<- predict(model_SVMLM_MRS1,test_mrs2)
RMSE.SVMLM_MRS1 <- RMSE(test_mrs2$MRS,pred_SVMLM_MRS1) #RMSE =0.91

# rf plot predicted vs Observed
plot(pred_SVMLM_MRS1,test_mrs2$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("SVMLM_FTIR RMSE:",round(RMSE.SVMLM_MRS1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMLM_MRS.vector <- c()
for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_mrs3 <- createDataPartition(df2_MRS$MRS, p = .7, 
                                          list = FALSE, 
                                          times = 1, groups =3) #caret package
  #split the whole dataset
  train_mrs3 <- df2_MRS[trainIndex_mrs3,]
  test_mrs3 <- df2_MRS[-trainIndex_mrs3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMLM_MRS <- train(MRS~., data = train_mrs3, method = "svmLinear", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMLM_MRS <- predict(model_SVMLM_MRS,test_mrs3)
  
  # calculate the RMSE value
  rmse.SVMLM_MRS <- RMSE(pred = test_mrs3$MRS,predicted.SVMLM_MRS)
  # put results to the vectors
  rmse.SVMLM_MRS.vector <- c(rmse.SVMLM_MRS.vector, rmse.SVMLM_MRS)
}

# compute the standard deviation
sd.SVMLM_MRS <- round(sd(rmse.SVMLM_MRS.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMLM_MRS <- CI(rmse.SVMLM_MRS.vector, ci = 0.95) # 0.8536262 0.8277178 0.8018094 
# rmse mean
rmse.SVMLM_MRS.mean <- round(mean(rmse.SVMLM_MRS.vector), digits = 2) #0.83

plot(rmse.SVMLM_MRS.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_MRS.mean,"+/- ",sd.SVMLM_MRS, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_MRS[3], digits = 2),
                  "to",round(ci.SVMLM_MRS[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMLM_MRS2 <- train(MRS ~ ., data = dfA2_MRS,
                           method = "svmLinear", trControl = fitControl)
pred_SVMLM_MRS2<- predict(model_SVMLM_MRS2,dfB2_MRS)
RMSE.SVMLM_MRS2 <- RMSE(dfB2_MRS$MRS,pred_SVMLM_MRS2) #RMSE =1.18

# rf plot predicted vs Observed
plot(pred_SVMLM_MRS2,dfB2_MRS$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_MRS2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMRadial_MRS #####################################################

##################
# the first method
model_SVMR_MRS1 <- train(MRS ~ ., data = train_mrs2,
                          method = "svmRadial", trControl = fitControl)
pred_SVMR_MRS1<- predict(model_SVMR_MRS1,test_mrs2)
RMSE.SVMR_MRS1 <- RMSE(test_mrs2$MRS,pred_SVMR_MRS1) #RMSE =0.63
accuracy_SVMR_MRS <- accuracy(pred_SVMR_MRS1, test_mrs2$MRS) #84.21

# this is the best SVMR model for LAB
saveRDS(model_SVMR_MRS1, "./models/CB/FTIR/model_SVMR_LAB.rds") #save model as an rds object
saveRDS(test_mrs2$MRS, "./models/CB/FTIR/tested_SVMR_LAB.rds") #save tested data as an rds object
saveRDS(pred_SVMR_MRS1, "./models/CB/FTIR/predicted_SVMR_LAB.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMR_MRS1, "./models/CB/FTIR/RMSE_SVMR_LAB.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMR_MRS, "./models/CB/FTIR/accuracy_SVMR_LAB.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMR_MRS1,test_mrs2$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_MRS1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMR_MRS.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_mrs3 <- createDataPartition(df2_MRS$MRS, p = .7, 
                                          list = FALSE, 
                                          times = 1, groups =3) #caret package
  #split the whole dataset
  train_mrs3 <- df2_MRS[trainIndex_mrs3,]
  test_mrs3 <- df2_MRS[-trainIndex_mrs3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMR_MRS <- train(MRS~., data = train_mrs3, method = "svmRadial", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMR_MRS <- predict(model_SVMR_MRS,test_mrs3)
  
  # calculate the RMSE value
  rmse.SVMR_MRS <- RMSE(pred = test_mrs3$MRS,predicted.SVMR_MRS)
  # put results to the vectors
  rmse.SVMR_MRS.vector <- c(rmse.SVMR_MRS.vector, rmse.SVMR_MRS)
}

# compute the standard deviation
sd.SVMR_MRS <- round(sd(rmse.SVMR_MRS.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.SVMR_MRS <- CI(rmse.SVMR_MRS.vector, ci = 0.95) # 0.6625555 0.6388163 0.6150771   
# rmse mean
rmse.SVMR_MRS.mean <- round(mean(rmse.SVMR_MRS.vector), digits = 2) #0.64

plot(rmse.SVMR_MRS.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_MRS.mean,"+/- ",sd.SVMR_MRS, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_MRS[3], digits = 2),
                  "to",round(ci.SVMR_MRS[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMR_MRS2 <- train(MRS ~ ., data = dfA2_MRS,
                          method = "svmRadial", trControl = fitControl)
pred_SVMR_MRS2<- predict(model_SVMR_MRS2,dfB2_MRS)
RMSE.SVMR_MRS2 <- RMSE(dfB2_MRS$MRS,pred_SVMR_MRS2) #RMSE =1.16

# rf plot predicted vs Observed
plot(pred_SVMR_MRS2,dfB2_MRS$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_MRS2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMPoly_MRS #####################################################

##################
# the first method
model_SVMP_MRS1 <- train(MRS ~ ., data = train_mrs2,
                          method = "svmPoly", trControl = fitControl)
pred_SVMP_MRS1<- predict(model_SVMP_MRS1,test_mrs2)
RMSE.SVMP_MRS1 <- RMSE(test_mrs2$MRS,pred_SVMP_MRS1) #RMSE =0.73
accuracy_SVMP_MRS <- accuracy(pred_SVMP_MRS1, test_mrs2$MRS ) #  78.95

# this is the best SVMP model for LAB
saveRDS(model_SVMP_MRS1, "./models/CB/FTIR/model_SVMP_LAB.rds") #save model as an rds object
saveRDS(test_mrs2$MRS, "./models/CB/FTIR/tested_SVMP_LAB.rds") #save tested data as an rds object
saveRDS(pred_SVMP_MRS1, "./models/CB/FTIR/predicted_SVMP_LAB.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMP_MRS1, "./models/CB/FTIR/RMSE_SVMP_LAB.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMP_MRS, "./models/CB/FTIR/accuracy_SVMP_LAB.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMP_MRS1,test_mrs2$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_MRS1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.SVMP_MRS.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_mrs3 <- createDataPartition(df2_MRS$MRS, p = .7, 
                                          list = FALSE, 
                                          times = 1, groups =3) #caret package
  #split the whole dataset
  train_mrs3 <- df2_MRS[trainIndex_mrs3,]
  test_mrs3 <- df2_MRS[-trainIndex_mrs3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMP_MRS <- train(MRS~., data = train_mrs3, method = "svmPoly", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMP_MRS <- predict(model_SVMP_MRS,test_mrs3)
  
  # calculate the RMSE value
  rmse.SVMP_MRS <- RMSE(pred = test_mrs3$MRS,predicted.SVMP_MRS)
  # put results to the vectors
  rmse.SVMP_MRS.vector <- c(rmse.SVMP_MRS.vector, rmse.SVMP_MRS)
}

# compute the standard deviation
sd.SVMP_MRS <- round(sd(rmse.SVMP_MRS.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.SVMP_MRS <- CI(rmse.SVMP_MRS.vector, ci = 0.95) # 0.6949040 0.6721088 0.6493135   
# rmse mean
rmse.SVMP_MRS.mean <- round(mean(rmse.SVMP_MRS.vector), digits = 2) #0.67

plot(rmse.SVMP_MRS.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_MRS.mean,"+/- ",sd.SVMP_MRS, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_MRS[3], digits = 2),
                  "to",round(ci.SVMP_MRS[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMP_MRS2 <- train(MRS ~ ., data = dfA2_MRS,
                          method = "svmPoly", trControl = fitControl)
pred_SVMP_MRS2<- predict(model_SVMP_MRS2,dfB2_MRS)
RMSE.SVMP_MRS2 <- RMSE(dfB2_MRS$MRS,pred_SVMP_MRS2) #RMSE =1.17

# rf plot predicted vs Observed
plot(pred_SVMP_MRS2,dfB2_MRS$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMP_MRS2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ LM_MRS #####################################################

##################
# the first method
model_lm_MRS1 <- train(MRS ~ ., data = train_mrs2,
                       method = "lm", trControl = fitControl)
pred_lm_MRS1<- predict(model_lm_MRS1,test_mrs2)
RMSE.lm_MRS1 <- RMSE(test_mrs2$MRS,pred_lm_MRS1) #RMSE =0.9

# rf plot predicted vs Observed
plot(pred_lm_MRS1,test_mrs2$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_MRS1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# 50 iterations
rmse.lm_MRS.vector <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_mrs3 <- createDataPartition(df2_MRS$MRS, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_mrs3 <- df2_MRS[trainIndex_mrs3,]
  test_mrs3 <- df2_MRS[-trainIndex_mrs3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_lm_MRS <- train(MRS~., data = train_mrs3, method = "lm", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.lm_MRS <- predict(model_lm_MRS,test_mrs3)
  
  # calculate the RMSE value
  rmse.lm_MRS <- RMSE(pred = test_mrs3$MRS,predicted.lm_MRS)
  # put results to the vectors
  rmse.lm_MRS.vector <- c(rmse.lm_MRS.vector, rmse.lm_MRS)
}
# compute the standard deviation
sd.lm_MRS <- round(sd(rmse.lm_MRS.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.lm_MRS <- CI(rmse.lm_MRS.vector, ci = 0.95) #package: Rmisc #0.84, 0.81, 0.79
# rmse mean
rmse.lm_MRS.mean <- round(mean(rmse.lm_MRS.vector), digits = 2) #0.81

plot(rmse.lm_MRS.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_MRS.mean,"+/- ",sd.lm_MRS, 
                  ")\nThere is a 95% likelihood that the range", round(ci.lm_MRS[3], digits = 2),
                  "to",round(ci.lm_MRS[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_lm_MRS2 <- train(MRS ~ ., data = dfA2_MRS,
                       method = "lm", trControl = fitControl)
pred_lm_MRS2<- predict(model_lm_MRS2,dfB2_MRS)
RMSE.lm_MRS2 <- RMSE(dfB2_MRS$MRS,pred_lm_MRS2) #RMSE =3.41

# rf plot predicted vs Observed
plot(pred_lm_MRS2,dfB2_MRS$MRS,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_MRS2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()

library("readxl")
library(caret)
library(Rmisc)

source("CreateObj.R")

# xlsx files
df <- read_excel("chickenThighFillet_FTIR_log2_v3.xlsx")
dfA <- read_excel("chickenThighFillet_FTIR_log2_batchA.xlsx")
dfB <- read_excel("chickenThighFillet_FTIR_log2_batchB.xlsx")

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

############## method to round wavelengths ###########################
###### df
#round colnames
roundedCol_df <- round(as.numeric(colnames(df[,4:length(colnames(df))])))
colnames(df) <- c(colnames(df[,1:3]), roundedCol_df)
df_orgin <- df
data <- df
d <- c()
for(i in 4:(length(colnames(df_orgin))-1)){
  if(as.numeric(colnames(df_orgin[,i])) == as.numeric(colnames(df_orgin[,i+1]))){
    for(j in 1:nrow(df_orgin)){
      d[j] <- mean(as.numeric(df_orgin[j,i]), as.numeric(df_orgin[j,i+1]))
    }
    d <- matrix(d, ncol = 1)
    d <- as.data.frame(d)
    rownames(d) <- rownames(df_orgin)
    colnames(d) <- colnames(df_orgin[,i])
    for(z in 4:(length(colnames(data))-1)){
      if(as.numeric(colnames(data[z])) == as.numeric(colnames(df_orgin[,i]))){
        df <- CreateObj(data[,1:(z-1)], d)
        df <- CreateObj(df, data[,(z+2):length(data)])
        data <- df
      }
    }
    d <- c()
    i = i+1
  }
}

########## dfA
#round colnames
roundedCol_dfA <- round(as.numeric(colnames(dfA[,3:length(colnames(dfA))])))
colnames(dfA) <- c(colnames(dfA[,1:2]), roundedCol_dfA)
dfA_orgin <- dfA
dataA <- dfA
d <- c()
for(i in 3:(length(colnames(dfA_orgin))-1)){
  if(as.numeric(colnames(dfA_orgin[,i])) == as.numeric(colnames(dfA_orgin[,i+1]))){
    for(j in 1:nrow(dfA_orgin)){
      d[j] <- mean(as.numeric(dfA_orgin[j,i]), as.numeric(dfA_orgin[j,i+1]))
    }
    d <- matrix(d, ncol = 1)
    d <- as.data.frame(d)
    rownames(d) <- rownames(dfA_orgin)
    colnames(d) <- colnames(dfA_orgin[,i])
    for(z in 3:(length(colnames(dataA))-1)){
      if(as.numeric(colnames(dataA[z])) == as.numeric(colnames(dfA_orgin[,i]))){
        dfA <- CreateObj(dataA[,1:(z-1)], d)
        dfA <- CreateObj(dfA, dataA[,(z+2):length(dataA)])
        dataA <- dfA
      }
    }
    d <- c()
    i = i+1
  }
}

######## dfB
#round colnames
roundedCol_dfB <- round(as.numeric(colnames(dfB[,4:length(colnames(dfB))])))
colnames(dfB) <- c(colnames(dfB[,1:3]), roundedCol_dfB)
dfB_orgin <- dfB
dataB <- dfB
d <- c()
for(i in 4:(length(colnames(dfB_orgin))-1)){
  if(as.numeric(colnames(dfB_orgin[,i])) == as.numeric(colnames(dfB_orgin[,i+1]))){
    for(j in 1:nrow(dfB_orgin)){
      d[j] <- mean(as.numeric(dfB_orgin[j,i]), as.numeric(dfB_orgin[j,i+1]))
    }
    d <- matrix(d, ncol = 1)
    d <- as.data.frame(d)
    rownames(d) <- rownames(dfB_orgin)
    colnames(d) <- colnames(dfB_orgin[,i])
    for(z in 4:(length(colnames(dataB))-1)){
      if(as.numeric(colnames(dataB[z])) == as.numeric(colnames(dfB_orgin[,i]))){
        dfB <- CreateObj(dataB[,1:(z-1)], d)
        dfB <- CreateObj(dfB, dataB[,(z+2):length(dataB)])
        dataB <- dfB
      }
    }
    d <- c()
    i = i+1
  }
}


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
spectra <- df[,4:ncol(df)]
#calculate correlation matrix
correlationMatrix <- cor(spectra)
#summarize the correlation matrix
print(correlationMatrix)
# find attributes which are highly correlated (ideally > 0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.999)
# print index of highly correlated attributes
print(highlyCorrelated)
#remove highly correlated attributes
spectra2 <- spectra[, - highlyCorrelated]

######################################################
# that data is highly correlated 
# create a dataframe with removed features
df2 <- CreateObj(df[2:3], spectra2)


################################### batch A ###########################################################
#ensure the results are repeatable
set.seed(123)
# create specra dataframe
spectraA <- dfA[,3:ncol(dfA)]
spectraB <- dfB[,3:ncol(dfB)]
#calculate correlation matrix
correlationMatrixA <- cor(spectraA)
#summarize the correlation matrix
print(correlationMatrixA)
# find attributes which are highly correlated (ideally > 0.75)
highlyCorrelatedA <- findCorrelation(correlationMatrixA, cutoff = 0.999)
# print index of highly correlated attributes
print(highlyCorrelatedA)
#remove highly correlated attributes
spectra2A <- spectraA[,-highlyCorrelatedA]
spectra2B <- spectraB[,-highlyCorrelatedA]

######################################################
# that data is highly correlated 
# create a dataframe with removed features
dfA2 <- CreateObj(dfA[1:2], spectra2A)
dfB2 <- CreateObj(dfB[1:2], spectra2B)


#######################################################################################################
############################# feature selection step 2 ################################################
############################# rank features by importance #############################################
#prerate training scheme
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# create df for diffrent bacterial counts separetly
df_TSA_wf <- df2[,-2]
df_CFC_wf <- df2[,-c(1)]

# create df for diffrent bacterial counts separetly batch A
dfA_TSA_wf <- dfA2[,-c(2)]
dfA_CFC_wf <- dfA2[,-c(1)]

# create df for diffrent bacterial counts separetly batch B
dfB_TSA_wf <- dfB2[,-c(2)]
dfB_CFC_wf <- dfB2[,-c(1)]


#######################################################################################################
#################################### TSA ##############################################################
#######################################################################################################

######################### select features for the whole dataset #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_tvc <- train(TVC~., data = df_TSA_wf, method = "rf", trControl = control, importance = TRUE)
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
model_varImp_tvc_bA <- train(TVC~., data = dfA_TSA_wf, method = "rf", trControl = control, importance = TRUE)
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

###############################################################################################################
################################### Prediction for TSA ########################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_tsa1 <- createDataPartition(df_TSA$TVC, p = .7, 
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
results_TSA <- rfe(TVC ~ ., data = train_tsa1,
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

######################### select features for the btch A #######################################
# run the RFE algorithm
results_TSA_bA <- rfe(TVC ~ ., data = dfA_TSA,
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

#############################################################################################
################################ RF_TVC #####################################################

##################
# the first method
model_Rf_TSA1 <- train(TVC ~ ., data = train_tsa2,
                       method = "rf", trControl = fitControl)
pred_Rf_TSA1<- predict(model_Rf_TSA1,test_tsa2)
RMSE.Rf_TSA1 <- RMSE(test_tsa2$TVC,pred_Rf_TSA1) #RMSE =1.36

# rf plot predicted vs Observed
plot(pred_Rf_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_Rf_TSA <- train(TVC~., data = train_tsa3, method = "rf", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TVC,predicted.Rf_TSA)
  # put results to the vectors
  rmse.Rf_TSA.vector <- c(rmse.Rf_TSA.vector, rmse.Rf_TSA)
}
# compute the standard deviation
sd.Rf_TSA <- round(sd(rmse.Rf_TSA.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.Rf_TSA <- CI(rmse.Rf_TSA.vector, ci = 0.95) #package: Rmisc #1.212559 1.185966 1.159373
# rmse mean
rmse.Rf_TSA.mean <- round(mean(rmse.Rf_TSA.vector), digits = 2) #1.19

plot(rmse.Rf_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_TSA.mean,"+/- ",sd.Rf_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_TSA[3], digits = 2),
                  "to",round(ci.Rf_TSA[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_TSA2 <- train(TVC ~ ., data = dfA2_TSA,
                       method = "rf", trControl = fitControl)
pred_Rf_TSA2<- predict(model_Rf_TSA2,dfB2_TSA)
RMSE.Rf_TSA2 <- RMSE(dfB2_TSA$TVC,pred_Rf_TSA2) #RMSE =1.58

# rf plot predicted vs Observed
plot(pred_Rf_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ KNN_TVC #####################################################

##################
# the first method
model_knn_TSA1 <- train(TVC ~ ., data = train_tsa2,
                        method = "knn", trControl = fitControl)
pred_knn_TSA1<- predict(model_knn_TSA1,test_tsa2)
RMSE.knn_TSA1 <- RMSE(test_tsa2$TVC,pred_knn_TSA1) #RMSE =1.32

# rf plot predicted vs Observed
plot(pred_knn_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_knn_TSA <- train(TVC~., data = train_tsa3, method = "knn", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.knn_TSA <- predict(model_knn_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.knn_TSA <- RMSE(pred = test_tsa3$TVC,predicted.knn_TSA)
  # put results to the vectors
  rmse.knn_TSA.vector <- c(rmse.knn_TSA.vector, rmse.knn_TSA)
}

# compute the standard deviation
sd.knn_TSA <- round(sd(rmse.knn_TSA.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.knn_TSA <- CI(rmse.knn_TSA.vector, ci = 0.95) #1.308494 1.281847 1.255200
# rmse mean
rmse.knn_TSA.mean <- round(mean(rmse.knn_TSA.vector), digits = 2) #1.28

plot(rmse.knn_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_TSA.mean,"+/- ",sd.knn_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_TSA[3], digits = 2),
                  "to",round(ci.knn_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_knn_TSA2 <- train(TVC ~ ., data = dfA2_TSA,
                        method = "knn", trControl = fitControl)
pred_knn_TSA2<- predict(model_knn_TSA2,dfB2_TSA)
RMSE.knn_TSA2 <- RMSE(dfB2_TSA$TVC,pred_knn_TSA2) #RMSE =1.46

# rf plot predicted vs Observed
plot(pred_knn_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#############################################################################################
################################ svmLM_TVC #####################################################

##################
# the first method
model_SVMLM_TSA1 <- train(TVC ~ ., data = train_tsa2,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_TSA1<- predict(model_SVMLM_TSA1,test_tsa2)
RMSE.SVMLM_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMLM_TSA1) #RMSE =1.23

# rf plot predicted vs Observed
plot(pred_SVMLM_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMLM_TSA <- train(TVC~., data = train_tsa3, method = "svmLinear", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMLM_TSA <- predict(model_SVMLM_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.SVMLM_TSA <- RMSE(pred = test_tsa3$TVC,predicted.SVMLM_TSA)
  # put results to the vectors
  rmse.SVMLM_TSA.vector <- c(rmse.SVMLM_TSA.vector, rmse.SVMLM_TSA)
}

# compute the standard deviation
sd.SVMLM_TSA <- round(sd(rmse.SVMLM_TSA.vector), digits = 2) #0.08
# find the 95% confidence intervals parameters
ci.SVMLM_TSA <- CI(rmse.SVMLM_TSA.vector, ci = 0.95) # 1.12, 1.1, 1.08 
# rmse mean
rmse.SVMLM_TSA.mean <- round(mean(rmse.SVMLM_TSA.vector), digits = 2) #1.1

plot(rmse.SVMLM_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_TSA.mean,"+/- ",sd.SVMLM_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_TSA[3], digits = 2),
                  "to",round(ci.SVMLM_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMLM_TSA2 <- train(TVC ~ ., data = dfA2_TSA,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_TSA2<- predict(model_SVMLM_TSA2,dfB2_TSA)
RMSE.SVMLM_TSA2 <- RMSE(dfB2_TSA$TVC,pred_SVMLM_TSA2) #RMSE =45.85

# rf plot predicted vs Observed
plot(pred_SVMLM_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMRadial_TVC #####################################################

##################
# the first method
model_SVMR_TSA1 <- train(TVC ~ ., data = train_tsa2,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA1<- predict(model_SVMR_TSA1,test_tsa2)
RMSE.SVMR_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMR_TSA1) #RMSE =1.28

# rf plot predicted vs Observed
plot(pred_SVMR_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMR_TSA <- train(TVC~., data = train_tsa3, method = "svmRadial", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMR_TSA <- predict(model_SVMR_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.SVMR_TSA <- RMSE(pred = test_tsa3$TVC,predicted.SVMR_TSA)
  # put results to the vectors
  rmse.SVMR_TSA.vector <- c(rmse.SVMR_TSA.vector, rmse.SVMR_TSA)
}

# compute the standard deviation
sd.SVMR_TSA <- round(sd(rmse.SVMR_TSA.vector), digits = 2) #0.1
# find the 95% confidence intervals parameters
ci.SVMR_TSA <- CI(rmse.SVMR_TSA.vector, ci = 0.95) #1.21, 1.18, 1.15
# rmse mean
rmse.SVMR_TSA.mean <- round(mean(rmse.SVMR_TSA.vector), digits = 2) #1.18

plot(rmse.SVMR_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_TSA.mean,"+/- ",sd.SVMR_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_TSA[3], digits = 2),
                  "to",round(ci.SVMR_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMR_TSA2 <- train(TVC ~ ., data = dfA2_TSA,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA2<- predict(model_SVMR_TSA2,dfB2_TSA)
RMSE.SVMR_TSA2 <- RMSE(dfB2_TSA$TVC,pred_SVMR_TSA2) #RMSE =1.37

# rf plot predicted vs Observed
plot(pred_SVMR_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMPoly_TVC #####################################################

##################
# the first method
model_SVMP_TSA1 <- train(TVC ~ ., data = train_tsa2,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_TSA1<- predict(model_SVMP_TSA1,test_tsa2)
RMSE.SVMP_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMP_TSA1) #RMSE =1.26

# rf plot predicted vs Observed
plot(pred_SVMP_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMP_TSA <- train(TVC~., data = train_tsa3, method = "svmPoly", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMP_TSA <- predict(model_SVMP_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.SVMP_TSA <- RMSE(pred = test_tsa3$TVC,predicted.SVMP_TSA)
  # put results to the vectors
  rmse.SVMP_TSA.vector <- c(rmse.SVMP_TSA.vector, rmse.SVMP_TSA)
}

# compute the standard deviation
sd.SVMP_TSA <- round(sd(rmse.SVMP_TSA.vector), digits = 2) #0.14
# find the 95% confidence intervals parameters
ci.SVMP_TSA <- CI(rmse.SVMP_TSA.vector, ci = 0.95) # 1.25, 1.2, 1.16
# rmse mean
rmse.SVMP_TSA.mean <- round(mean(rmse.SVMP_TSA.vector), digits = 2) #1.2

plot(rmse.SVMP_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_TSA.mean,"+/- ",sd.SVMP_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_TSA[3], digits = 2),
                  "to",round(ci.SVMP_TSA[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMP_TSA2 <- train(TVC ~ ., data = dfA2_TSA,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_TSA2<- predict(model_SVMP_TSA2,dfB2_TSA)
RMSE.SVMP_TSA2 <- RMSE(dfB2_TSA$TVC,pred_SVMP_TSA2) #RMSE =1.45

# rf plot predicted vs Observed
plot(pred_SVMP_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_TSA2,digits = 2)))
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
model_varImp_cfc <- train(Ps~., data = df_CFC_wf, method = "rf", trControl = control, importance = TRUE)
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
model_varImp_cfc_bA <- train(Ps~., data = dfA_CFC_wf, method = "rf", trControl = control, importance = TRUE)
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

###############################################################################################################
################################### Prediction for CFC ########################################################
###############################################################################################################

#split data
set.seed(123) #caret package

#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndex_cfc1 <- createDataPartition(df_CFC$Ps, p = .7, 
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
results_CFC <- rfe(Ps ~ ., data = train_cfc1,
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

######################### select features for the btch A #######################################
# run the RFE algorithm
results_CFC_bA <- rfe(Ps ~ ., data = dfA_CFC,
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

#############################################################################################
################################ RF_CFC #####################################################

##################
# the first method
model_Rf_CFC1 <- train(Ps ~ ., data = train_cfc2,
                       method = "rf", trControl = fitControl)
pred_Rf_CFC1<- predict(model_Rf_CFC1,test_cfc2)
RMSE.Rf_CFC1 <- RMSE(test_cfc2$Ps,pred_Rf_CFC1) #RMSE =1.36

# rf plot predicted vs Observed
plot(pred_Rf_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_TSA1,digits = 2)))
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
  trainIndex_cfc3 <- createDataPartition(df2_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_Rf_CFC <- train(Ps~., data = train_cfc3, method = "rf", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_cfc3$Ps,predicted.Rf_CFC)
  # put results to the vectors
  rmse.Rf_CFC.vector <- c(rmse.Rf_CFC.vector, rmse.Rf_CFC)
}
# compute the standard deviation
sd.Rf_CFC <- round(sd(rmse.Rf_CFC.vector), digits = 2) #0.1
# find the 95% confidence intervals parameters
ci.Rf_CFC <- CI(rmse.Rf_CFC.vector, ci = 0.95) #package: Rmisc #1.37, 1.34, 1.31
# rmse mean
rmse.Rf_CFC.mean <- round(mean(rmse.Rf_CFC.vector), digits = 2) #1.34

plot(rmse.Rf_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_CFC.mean,"+/- ",sd.Rf_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_CFC[3], digits = 2),
                  "to",round(ci.Rf_CFC[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_CFC2 <- train(Ps ~ ., data = dfA2_CFC,
                       method = "rf", trControl = fitControl)
pred_Rf_CFC2<- predict(model_Rf_CFC2,dfB2_CFC)
RMSE.Rf_CFC2 <- RMSE(dfB2_CFC$Ps,pred_Rf_CFC2) #RMSE =1.55

# rf plot predicted vs Observed
plot(pred_Rf_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ KNN_CFC #####################################################

##################
# the first method
model_knn_CFC1 <- train(Ps ~ ., data = train_cfc2,
                        method = "knn", trControl = fitControl)
pred_knn_CFC1<- predict(model_knn_CFC1,test_cfc2)
RMSE.knn_CFC1 <- RMSE(test_cfc2$Ps,pred_knn_CFC1) #RMSE =1.32

# rf plot predicted vs Observed
plot(pred_knn_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_cfc3 <- createDataPartition(df2_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_knn_CFC <- train(Ps~., data = train_cfc3, method = "knn", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.knn_CFC <- predict(model_knn_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.knn_CFC <- RMSE(pred = test_cfc3$Ps,predicted.knn_CFC)
  # put results to the vectors
  rmse.knn_CFC.vector <- c(rmse.knn_CFC.vector, rmse.knn_CFC)
}

# compute the standard deviation
sd.knn_CFC <- round(sd(rmse.knn_CFC.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.knn_CFC <- CI(rmse.knn_CFC.vector, ci = 0.95) # 1.44, 1.41, 1.39
# rmse mean
rmse.knn_CFC.mean <- round(mean(rmse.knn_CFC.vector), digits = 2) #1.41

plot(rmse.knn_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_CFC.mean,"+/- ",sd.knn_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_CFC[3], digits = 2),
                  "to",round(ci.knn_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_knn_CFC2 <- train(Ps ~ ., data = dfA2_CFC,
                        method = "knn", trControl = fitControl)
pred_knn_CFC2<- predict(model_knn_CFC2,dfB2_CFC)
RMSE.knn_CFC2 <- RMSE(dfB2_CFC$Ps,pred_knn_CFC2) #RMSE =1.54

# rf plot predicted vs Observed
plot(pred_knn_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


#############################################################################################
################################ svmLM_CFC #####################################################

##################
# the first method
model_SVMLM_CFC1 <- train(Ps ~ ., data = train_cfc2,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_CFC1<- predict(model_SVMLM_CFC1,test_cfc2)
RMSE.SVMLM_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMLM_CFC1) #RMSE =1.24

# rf plot predicted vs Observed
plot(pred_SVMLM_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_cfc3 <- createDataPartition(df2_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMLM_CFC <- train(Ps~., data = train_cfc3, method = "svmLinear", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMLM_CFC <- predict(model_SVMLM_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.SVMLM_CFC <- RMSE(pred = test_cfc3$Ps,predicted.SVMLM_CFC)
  # put results to the vectors
  rmse.SVMLM_CFC.vector <- c(rmse.SVMLM_CFC.vector, rmse.SVMLM_CFC)
}

# compute the standard deviation
sd.SVMLM_CFC <- round(sd(rmse.SVMLM_CFC.vector), digits = 2) #0.1
# find the 95% confidence intervals parameters
ci.SVMLM_CFC <- CI(rmse.SVMLM_CFC.vector, ci = 0.95) # 1.1, 1.07, 1.05
# rmse mean
rmse.SVMLM_CFC.mean <- round(mean(rmse.SVMLM_CFC.vector), digits = 2) #1.07

plot(rmse.SVMLM_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_CFC.mean,"+/- ",sd.SVMLM_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_CFC[3], digits = 2),
                  "to",round(ci.SVMLM_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMLM_CFC2 <- train(Ps ~ ., data = dfA2_CFC,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_CFC2<- predict(model_SVMLM_CFC2,dfB2_CFC)
RMSE.SVMLM_CFC2 <- RMSE(dfB2_CFC$Ps,pred_SVMLM_CFC2) #RMSE =39.67

# rf plot predicted vs Observed
plot(pred_SVMLM_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMRadial_CFC #####################################################

##################
# the first method
model_SVMR_CFC1 <- train(Ps ~ ., data = train_cfc2,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_CFC1<- predict(model_SVMR_CFC1,test_cfc2)
RMSE.SVMR_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMR_CFC1) #RMSE =1.28

# rf plot predicted vs Observed
plot(pred_SVMR_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
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
  trainIndex_cfc3 <- createDataPartition(df2_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMR_CFC <- train(Ps~., data = train_cfc3, method = "svmRadial", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMR_CFC <- predict(model_SVMR_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.SVMR_CFC <- RMSE(pred = test_cfc3$Ps,predicted.SVMR_CFC)
  # put results to the vectors
  rmse.SVMR_CFC.vector <- c(rmse.SVMR_CFC.vector, rmse.SVMR_CFC)
}

# compute the standard deviation
sd.SVMR_CFC <- round(sd(rmse.SVMR_CFC.vector), digits = 2) #0.1
# find the 95% confidence intervals parameters
ci.SVMR_CFC <- CI(rmse.SVMR_CFC.vector, ci = 0.95) # 1.33, 1.31, 1.28 
# rmse mean
rmse.SVMR_CFC.mean <- round(mean(rmse.SVMR_CFC.vector), digits = 2) #1.31

plot(rmse.SVMR_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_CFC.mean,"+/- ",sd.SVMR_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_CFC[3], digits = 2),
                  "to",round(ci.SVMR_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMR_CFC2 <- train(Ps ~ ., data = dfA2_CFC,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_CFC2<- predict(model_SVMR_CFC2,dfB2_CFC)
RMSE.SVMR_CFC2 <- RMSE(dfB2_CFC$Ps,pred_SVMR_CFC2) #RMSE =1.52

# rf plot predicted vs Observed
plot(pred_SVMR_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#############################################################################################
################################ SVMPoly_CFC #####################################################

##################
# the first method
model_SVMP_CFC1 <- train(Ps ~ ., data = train_cfc2,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC1<- predict(model_SVMP_CFC1,test_cfc2)
RMSE.SVMP_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMP_CFC1) #RMSE =1.25

# rf plot predicted vs Observed
plot(pred_SVMP_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_CFC1,digits = 2)))
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
  trainIndex_cfc3 <- createDataPartition(df2_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_cfc3 <- df2_CFC[trainIndex_cfc3,]
  test_cfc3 <- df2_CFC[-trainIndex_cfc3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_SVMP_CFC <- train(Ps~., data = train_cfc3, method = "svmPoly", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.SVMP_CFC <- predict(model_SVMP_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.SVMP_CFC <- RMSE(pred = test_cfc3$Ps,predicted.SVMP_CFC)
  # put results to the vectors
  rmse.SVMP_CFC.vector <- c(rmse.SVMP_CFC.vector, rmse.SVMP_CFC)
}

# compute the standard deviation
sd.SVMP_CFC <- round(sd(rmse.SVMP_CFC.vector), digits = 2) #0.1
# find the 95% confidence intervals parameters
ci.SVMP_CFC <- CI(rmse.SVMP_CFC.vector, ci = 0.95) # 1.2, 1.17, 1.15 
# rmse mean
rmse.SVMP_CFC.mean <- round(mean(rmse.SVMP_CFC.vector), digits = 2) #1.17

plot(rmse.SVMP_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_CFC.mean,"+/- ",sd.SVMP_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_CFC[3], digits = 2),
                  "to",round(ci.SVMP_CFC[1], digits = 2),"covers the true error of the model."))


##################
# the second method
# the first method
model_SVMP_CFC2 <- train(Ps ~ ., data = dfA2_CFC,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC2<- predict(model_SVMP_CFC2,dfB2_CFC)
RMSE.SVMP_CFC2 <- RMSE(dfB2_CFC$Ps,pred_SVMP_CFC2) #RMSE =1.55

# rf plot predicted vs Observed
plot(pred_SVMP_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)



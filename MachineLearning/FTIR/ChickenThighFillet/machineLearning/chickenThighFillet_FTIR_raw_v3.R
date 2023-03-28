# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()

library("readxl")
library(caret)
library(Rmisc)

source("CreateObj.R")

# xlsx files
df <- read_excel("chickenThighFillet_FTIR_raw_v3.xlsx")
dfA <- read_excel("chickenThighFillet_FTIR_raw_batchA.xlsx")
dfB <- read_excel("chickenThighFillet_FTIR_raw_batchB.xlsx")

##########################################################################################
#function to compute model accuracy
accuracy <- function(predicted, tested){
  #
  # Fuction to compute the accuracy at ±1 LogCount for regression models
  # predicted = a vector of numeric data with predicted values
  # tested = a vector of numeric data with tested values
  # Hint - the accuracy of the model in percentage is returned
  #
  S<-c()
  for ( i in 1:length(predicted)){
    if(((tested[i])+1)>=predicted[i] & predicted[i]>=((tested[i])-1)){
      S[length(S)+1]<-i
    }
    else{
    }
  }
  round(length(S)/length(tested)*100, digits = 2)
}

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
roundedCol_dfB <- round(as.numeric(colnames(dfB[,3:length(colnames(dfB))])))
colnames(dfB) <- c(colnames(dfB[,1:2]), roundedCol_dfB)
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
    for(z in 3:(length(colnames(dataB))-1)){
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
removeID <- c(1:3)
features <- colnames(df)
for(i in 4:length(features)){
  if(as.numeric(features[i]) >= 1000 && as.numeric(features[i]) < 4000){
    removeID[length(removeID)+1] <- i
  }
}
df <- df[,removeID]

#batch A
removeID_A <- c(1:2)
featuresA <- colnames(dfA)
for(i in 3:ncol(dfA)){
  if(as.numeric(featuresA[i]) >= 1000 && as.numeric(featuresA[i]) < 4000){
    removeID_A[length(removeID_A)+1] <- i
  }
}
dfA <- dfA[,removeID_A]

#batch B
removeID_B <- c(1:2)
featuresB <- colnames(dfB)
for(i in 3:ncol(dfB)){
  if(as.numeric(featuresB[i]) >= 1000 && as.numeric(featuresB[i]) < 4000){
    removeID_B[length(removeID_B)+1] <- i
  }
}  
dfB <- dfB[,removeID_B]

####################################################################################################
############################# choose every 100 feature from spectra ################################

#for the whole dataset
df100 <- df[,seq(4, ncol(df), 100)]
df100 <- CreateObj(df[,2:3], df100)
df100_TSA <- df100[,c(1, 3:ncol(df100))]
df100_CFC <- df100[,c(2:ncol(df100))]

# for batch A
dfA100 <- dfA[,seq(3, ncol(dfA), 100)]
dfA100 <- CreateObj(dfA[,1:2], dfA100)
dfA100_TSA <- dfA100[,c(1, 3:ncol(dfA100))]
dfA100_CFC <- dfA100[,c(2:ncol(dfA100))]  
# for batch B
dfB100 <- dfB[,seq(3, ncol(dfB), 100)]
dfB100 <- CreateObj(dfB[,1:2], dfB100)
dfB100_TSA <- dfB100[,c(1, 3:ncol(dfB100))]
dfB100_CFC <- dfB100[,c(2:ncol(dfB100))]

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

################################### batch B ###########################################################
#ensure the results are repeatable
set.seed(123)
# create specra dataframe
spectraA_2 <- dfA[,3:ncol(dfA)]
spectraB_2 <- dfB[,3:ncol(dfB)]
#calculate correlation matrix
correlationMatrixA_2 <- cor(spectraB_2)
#summarize the correlation matrix
print(correlationMatrixA_2)
# find attributes which are highly correlated (ideally > 0.75)
highlyCorrelatedA_2 <- findCorrelation(correlationMatrixA_2, cutoff = 0.999)
# print index of highly correlated attributes
print(highlyCorrelatedA_2)
#remove highly correlated attributes
spectra2A_2 <- spectraA_2[,-highlyCorrelatedA_2]
spectra2B_2 <- spectraB_2[,-highlyCorrelatedA_2]

######################################################
# that data is highly correlated 
# create a dataframe with removed features
dfA2_B <- CreateObj(dfA[1:2], spectra2A_2)
dfB2_B <- CreateObj(dfB[1:2], spectra2B_2)

#######################################################################################################
############################# feature selection step 2 ################################################
############################# rank features by importance #############################################
#prerate training scheme
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# create df for diffrent bacterial counts separetly
df_TSA_wf <- df2[,-2]
df_CFC_wf <- df2[,-c(1)]

# create df for diffrent bacterial counts separetly batch A based on A
dfA_TSA_wf <- dfA2[,-c(2)]
dfA_CFC_wf <- dfA2[,-c(1)]

# create df for diffrent bacterial counts separetly batch B based on A
dfB_TSA_wf <- dfB2[,-c(2)]
dfB_CFC_wf <- dfB2[,-c(1)]

# create df for diffrent bacterial counts separetly batch A based on B
dfA_TSA_wf_B <- dfA2_B[,-c(2)]
dfA_CFC_wf_B <- dfA2_B[,-c(1)]

# create df for diffrent bacterial counts separetly batch B based on B
dfB_TSA_wf_B <- dfB2_B[,-c(2)]
dfB_CFC_wf_B <- dfB2_B[,-c(1)]


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

######################### select features for the batch B #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_tvc_bB <- train(TVC~., data = dfB_TSA_wf_B, method = "rf", trControl = control, importance = TRUE)
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
dfA_TSA_B <- dfA_TSA_wf_B[,feature_ID_TSA_bB]
dfB_TSA_B <- dfB_TSA_wf_B[,feature_ID_TSA_bB]

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

############################# Randomly split data 100 ################################
trainIndex100_TSA1 <- createDataPartition(df100_TSA$TVC, p = .7, list = FALSE, times = 1, groups =3)
train_TSA100 <- df100_TSA[trainIndex100,]
test_TSA100 <- df100_TSA[-trainIndex100,]




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
saveRDS(predictors_tsa_toSave, "./models/CTF/FTIR/predictors_TVC.rds") #save model as an rds object

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

######################### select features for the batch B #######################################
# run the RFE algorithm
results_TSA_bB <- rfe(TVC ~ ., data = dfB_TSA_B,
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
trainIndexA_tsa1 <- createDataPartition(dfA2_TSA$TVC, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainA_tsa1 <- dfA2_TSA[trainIndexA_tsa1,]
testA_tsa1 <- dfA2_TSA[-trainIndexA_tsa1,]

#################### split batch B ##########################################################
#split data
set.seed(123) #caret package
#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndexB_tsa1 <- createDataPartition(dfB2_TSA_B$TVC, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainB_tsa1 <- dfB2_TSA_B[trainIndexB_tsa1,]
testB_tsa1 <- dfB2_TSA_B[-trainIndexB_tsa1,]

#############################################################################################
################################ RF_TVC #####################################################

##################
# the first method
model_Rf_TSA1 <- train(TVC ~ ., data = train_tsa2,
                       method = "rf", trControl = fitControl)
pred_Rf_TSA1<- predict(model_Rf_TSA1,test_tsa2)
RMSE.Rf_TSA1 <- RMSE(test_tsa2$TVC,pred_Rf_TSA1) #RMSE =1.32
# Random Forest model Accuracy for TVC
accuracy_Rf_TSA <- accuracy(pred_Rf_TSA1, test_tsa2$TVC) #50


# rf plot predicted vs Observed
plot(pred_Rf_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.Rf_TSA1,digits = 2), 
                "\nAccuracy",accuracy_Rf_TSA,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# saveRDS(model_Rf_TSA1, "./models/CTF/FTIR/model_RF_TVC.rds") #save model as an rds object
# saveRDS(test_tsa2$TVC, "./models/CTF/FTIR/tested_RF_TVC.rds") #save tested data as an rds object
# saveRDS(pred_Rf_TSA1, "./models/CTF/FTIR/predicted_RF_TVC.rds") #save predicted data as an rds object
# saveRDS(RMSE.Rf_TSA1, "./models/CTF/FTIR/RMSE_RF_TVC.rds") #save RMSE as an rds object
# saveRDS(accuracy_Rf_TSA, "./models/CTF/FTIR/accuracy_RF_TVC.rds") #save accuracy as an rds object

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
ci.Rf_TSA <- CI(rmse.Rf_TSA.vector, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.Rf_TSA.mean <- round(mean(rmse.Rf_TSA.vector), digits = 2) #1.19

plot(rmse.Rf_TSA.vector, type="l",ylim=c(0,2),xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_TSA.mean,"+/- ",sd.Rf_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_TSA[3], digits = 2),
                  "to",round(ci.Rf_TSA[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_TSA2 <- train(TVC ~ ., data = dfA2_TSA,
                       method = "rf", trControl = fitControl)
pred_Rf_TSA2<- predict(model_Rf_TSA2,dfB2_TSA)
RMSE.Rf_TSA2 <- RMSE(dfB2_TSA$TVC,pred_Rf_TSA2) #RMSE =1.35

# rf plot predicted vs Observed
plot(pred_Rf_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_Rf_TSA100 <- train(TVC ~ ., data = train_TSA100, method = "rf")
pred_Rf_TSA100 <- predict(model_Rf_TSA100, test_TSA100)
RMSE.Rf_TSA100 <- RMSE(test_TSA100$TVC,pred_Rf_TSA100)
accuracy_Rf_TSA100 <- accuracy(pred_Rf_TSA100, test_TSA100$TVC) 

plot(pred_Rf_TSA100,test_TSA100$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.Rf_TSA100,digits = 2), 
                "\nAccuracy",accuracy_Rf_TSA100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.Rf_TSA.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df100_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df100_TSA[trainIndex_tsa3,]
  test_tsa3 <- df100_TSA[-trainIndex_tsa3,]

  # train model
  model_Rf_TSA <- train(TVC~., data = train_tsa3, method = "rf") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TVC,predicted.Rf_TSA)
  # put results to the vectors
  rmse.Rf_TSA.vector100 <- c(rmse.Rf_TSA.vector100, rmse.Rf_TSA)
}
# compute the standard deviation
sd.Rf_TSA100 <- round(sd(rmse.Rf_TSA.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.Rf_TSA100 <- CI(rmse.Rf_TSA.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.Rf_TSA.mean100 <- round(mean(rmse.Rf_TSA.vector100), digits = 2) #1.19

plot(rmse.Rf_TSA.vector100, type="l",ylim = c(0,2),xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_TSA.mean100,"+/- ",sd.Rf_TSA100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_TSA100[3], digits = 2),
                  "to",round(ci.Rf_TSA100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ KNN_TVC #####################################################

##################
# the first method
model_knn_TSA1 <- train(TVC ~ ., data = train_tsa2,
                        method = "knn", trControl = fitControl)
pred_knn_TSA1<- predict(model_knn_TSA1,test_tsa2)
RMSE.knn_TSA1 <- RMSE(test_tsa2$TVC,pred_knn_TSA1) #RMSE =1.37
accuracy_knn_TSA <- accuracy(pred_knn_TSA1, test_tsa2$TVC) #86.84
# this is the best k-nn model for TVC
saveRDS(model_knn_TSA1, "./models/CTF/FTIR/model_knn_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TVC, "./models/CTF/FTIR/tested_knn_TVC.rds") #save tested data as an rds object
saveRDS(pred_knn_TSA1, "./models/CTF/FTIR/predicted_knn_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.knn_TSA1, "./models/CTF/FTIR/RMSE_knn_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_knn_TSA, "./models/CTF/FTIR/accuracy_knn_TVC.rds") #save accuracy as an rds object

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
sd.knn_TSA <- round(sd(rmse.knn_TSA.vector), digits = 2) #0.1
# find the 95% confidence intervals parameters
ci.knn_TSA <- CI(rmse.knn_TSA.vector, ci = 0.95) # 1.310374 1.282562 1.254750
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
RMSE.knn_TSA2 <- RMSE(dfB2_TSA$TVC,pred_knn_TSA2) #RMSE =1.3

# rf plot predicted vs Observed
plot(pred_knn_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("KNN_FTIR RMSE:",round(RMSE.knn_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)


# for df100
model_knn_TSA100 <- train(TVC ~ ., data = train_TSA100, method = "knn")
pred_knn_TSA100 <- predict(model_knn_TSA100, test_TSA100)
RMSE.knn_TSA100 <- RMSE(test_TSA100$TVC,pred_knn_TSA100)
accuracy_knn_TSA100 <- accuracy(pred_knn_TSA100, test_TSA100$TVC) 

plot(pred_knn_TSA100,test_TSA100$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.knn_TSA100,digits = 2), 
                "\nAccuracy",accuracy_knn_TSA100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.knn_TSA.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df100_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df100_TSA[trainIndex_tsa3,]
  test_tsa3 <- df100_TSA[-trainIndex_tsa3,]
  
  # train model
  model_Rf_TSA <- train(TVC~., data = train_tsa3, method = "knn") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TVC,predicted.Rf_TSA)
  # put results to the vectors
  rmse.knn_TSA.vector100 <- c(rmse.knn_TSA.vector100, rmse.Rf_TSA)
}
# compute the standard deviation
sd.knn_TSA100 <- round(sd(rmse.knn_TSA.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.knn_TSA100 <- CI(rmse.knn_TSA.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.knn_TSA.mean100 <- round(mean(rmse.knn_TSA.vector100), digits = 2) #1.19

plot(rmse.knn_TSA.vector100, type="l",ylim=c(0,2),xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_TSA.mean100,"+/- ",sd.knn_TSA100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_TSA100[3], digits = 2),
                  "to",round(ci.knn_TSA100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ svmLM_TVC #####################################################

##################
# the first method
model_SVMLM_TSA1 <- train(TVC ~ ., data = train_tsa2,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_TSA1<- predict(model_SVMLM_TSA1,test_tsa2)
RMSE.SVMLM_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMLM_TSA1) #RMSE =1.25
accuracy_SVMLM_TSA <- accuracy(pred_SVMLM_TSA1, test_tsa2$TVC)

# this is the best svm linear model for TVC
saveRDS(model_SVMLM_TSA1, "./models/CTF/FTIR/model_SVMLM_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TVC, "./models/CTF/FTIR/tested_SVMLM_TVC.rds") #save tested data as an rds object
saveRDS(pred_SVMLM_TSA1, "./models/CTF/FTIR/predicted_SVMLM_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMLM_TSA1, "./models/CTF/FTIR/RMSE_SVMLM_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMLM_TSA, "./models/CTF/FTIR/accuracy_SVMLM_TVC.rds") #save accuracy as an rds object

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
sd.SVMLM_TSA <- round(sd(rmse.SVMLM_TSA.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMLM_TSA <- CI(rmse.SVMLM_TSA.vector, ci = 0.95) #1.111985 1.085658 1.059330
# rmse mean
rmse.SVMLM_TSA.mean <- round(mean(rmse.SVMLM_TSA.vector), digits = 2) #1.09

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
RMSE.SVMLM_TSA2 <- RMSE(dfB2_TSA$TVC,pred_SVMLM_TSA2) #RMSE =1.37

# rf plot predicted vs Observed
plot(pred_SVMLM_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_SVMLM_TSA100 <- train(TVC ~ ., data = train_TSA100, method = "svmLinear")
pred_SVMLM_TSA100 <- predict(model_SVMLM_TSA100, test_TSA100)
RMSE.SVMLM_TSA100 <- RMSE(test_TSA100$TVC,pred_SVMLM_TSA100)
accuracy_SVMLM_TSA100 <- accuracy(pred_SVMLM_TSA100, test_TSA100$TVC) 

plot(pred_SVMLM_TSA100,test_TSA100$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.SVMLM_TSA100,digits = 2), 
                "\nAccuracy",accuracy_SVMLM_TSA100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.SVMLM_TSA.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df100_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df100_TSA[trainIndex_tsa3,]
  test_tsa3 <- df100_TSA[-trainIndex_tsa3,]
  
  # train model
  model_Rf_TSA <- train(TVC~., data = train_tsa3, method = "svmLinear") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TVC,predicted.Rf_TSA)
  # put results to the vectors
  rmse.SVMLM_TSA.vector100 <- c(rmse.SVMLM_TSA.vector100, rmse.Rf_TSA)
}
# compute the standard deviation
sd.SVMLM_TSA100 <- round(sd(rmse.SVMLM_TSA.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMLM_TSA100 <- CI(rmse.SVMLM_TSA.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.SVMLM_TSA.mean100 <- round(mean(rmse.SVMLM_TSA.vector100), digits = 2) #1.19

plot(rmse.SVMLM_TSA.vector100, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_TSA.mean100,"+/- ",sd.SVMLM_TSA100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_TSA100[3], digits = 2),
                  "to",round(ci.SVMLM_TSA100[1], digits = 2),"covers the true error of the model."))


#############################################################################################
################################ SVMRadial_TVC #####################################################

##################
# the first method
model_SVMR_TSA1 <- train(TVC ~ ., data = train_tsa2,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA1<- predict(model_SVMR_TSA1,test_tsa2)
RMSE.SVMR_TSA1 <- RMSE(test_tsa2$TVC,pred_SVMR_TSA1) #RMSE =1.26
accuracy_SVMR_TSA <- accuracy(pred_SVMR_TSA1, test_tsa2$TVC) #56.9

# this is the best svm Radial model for TVC
saveRDS(model_SVMR_TSA1, "./models/CTF/FTIR/model_SVMR_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TVC, "./models/CTF/FTIR/tested_SVMR_TVC.rds") #save tested data as an rds object
saveRDS(pred_SVMR_TSA1, "./models/CTF/FTIR/predicted_SVMR_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMR_TSA1, "./models/CTF/FTIR/RMSE_SVMR_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMR_TSA, "./models/CTF/FTIR/accuracy_SVMR_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMR_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Radial support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.SVMR_TSA1,digits = 2), 
                "\nAccuracy",accuracy_SVMR_TSA,"% - at ±1 LogCount"))
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
ci.SVMR_TSA <- CI(rmse.SVMR_TSA.vector, ci = 0.95) #1.206069 1.177259 1.148448  
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
RMSE.SVMR_TSA2 <- RMSE(dfB2_TSA$TVC,pred_SVMR_TSA2) #RMSE =1.49

# rf plot predicted vs Observed
plot(pred_SVMR_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Radial support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",
                round(RMSE.SVMR_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

###################
# trained based on B, tested based on A
model_SVMR_TSA3 <- train(TVC ~ ., data = dfB2_TSA_B,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA3<- predict(model_SVMR_TSA3,dfA2_TSA_B)
RMSE.SVMR_TSA3 <- RMSE(dfA2_TSA_B$TVC,pred_SVMR_TSA3) 

# rf plot predicted vs Observed
plot(pred_SVMR_TSA3,dfA2_TSA_B$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Radial support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",
                round(RMSE.SVMR_TSA3,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# trained and tested based on batch A
model_SVMR_TSA4 <- train(TVC ~ ., data = trainA_tsa1,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA4<- predict(model_SVMR_TSA4,testA_tsa1)
RMSE.SVMR_TSA4 <- RMSE(testA_tsa1$TVC,pred_SVMR_TSA4) 

# plot predicted vs Observed
plot(pred_SVMR_TSA4,testA_tsa1$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Radial support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",
                round(RMSE.SVMR_TSA4,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#######################
# trained and tested based on batch B
model_SVMR_TSA5 <- train(TVC ~ ., data = trainB_tsa1,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_TSA5<- predict(model_SVMR_TSA5,testB_tsa1)
RMSE.SVMR_TSA5 <- RMSE(testB_tsa1$TVC,pred_SVMR_TSA5) 

# plot predicted vs Observed
plot(pred_SVMR_TSA5,testB_tsa1$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Radial support-vector machine (SVM) model for Total Viable Counts (TVC) \nRMSE:",
                round(RMSE.SVMR_TSA5,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_SVMR_TSA100 <- train(TVC ~ ., data = train_TSA100, method = "svmRadial")
pred_SVMR_TSA100 <- predict(model_SVMR_TSA100, test_TSA100)
RMSE.SVMR_TSA100 <- RMSE(test_TSA100$TVC,pred_SVMR_TSA100)
accuracy_SVMR_TSA100 <- accuracy(pred_SVMR_TSA100, test_TSA100$TVC) 

plot(pred_SVMR_TSA100,test_TSA100$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.SVMR_TSA100,digits = 2), 
                "\nAccuracy",accuracy_SVMR_TSA100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.SVMR_TSA.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df100_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df100_TSA[trainIndex_tsa3,]
  test_tsa3 <- df100_TSA[-trainIndex_tsa3,]
  
  # train model
  model_Rf_TSA <- train(TVC~., data = train_Ttsa3, method = "svmRadial") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TVC,predicted.Rf_TSA)
  # put results to the vectors
  rmse.SVMR_TSA.vector100 <- c(rmse.SVMR_TSA.vector100, rmse.Rf_TSA)
}
# compute the standard deviation
sd.SVMR_TSA100 <- round(sd(rmse.SVMR_TSA.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMR_TSA100 <- CI(rmse.SVMR_TSA.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.SVMR_TSA.mean100 <- round(mean(rmse.SVMR_TSA.vector100), digits = 2) #1.19

plot(rmse.SVMR_TSA.vector100, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_TSA.mean100,"+/- ",sd.SVMR_TSA100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_TSA100[3], digits = 2),
                  "to",round(ci.SVMR_TSA100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ SVMPoly_TVC #####################################################

##################
# the first method
model_SVMP_TSA1 <- train(TVC ~ ., data = train_tsa2,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_TSA1<- predict(model_SVMP_TSA1,test_tsa2)
accuracy_SVMP_TSA <- accuracy(pred_SVMP_TSA1, test_tsa2$TVC) #53.43
RMSE.SVMP_TSA1 <- RMSE(test_tsa2$TVC, pred_SVMP_TSA1)

# this is the best svm polynomial model for TVC
saveRDS(model_SVMP_TSA1, "./models/CTF/FTIR/model_SVMP_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TVC, "./models/CTF/FTIR/tested_SVMP_TVC.rds") #save tested data as an rds object
saveRDS(pred_SVMP_TSA1, "./models/CTF/FTIR/predicted_SVMP_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMP_TSA1, "./models/CTF/FTIR/RMSE_SVMP_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMP_TSA, "./models/CTF/FTIR/accuracy_SVMP_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_SVMP_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_TSA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMP_TSA1, "./models/CTF/FTIR/model_SVMP_TVC.rds") #save model as an rds object

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
sd.SVMP_TSA <- round(sd(rmse.SVMP_TSA.vector), digits = 2) #0.15
# find the 95% confidence intervals parameters
ci.SVMP_TSA <- CI(rmse.SVMP_TSA.vector, ci = 0.95) # 1.234301 1.193029 1.151756
# rmse mean
rmse.SVMP_TSA.mean <- round(mean(rmse.SVMP_TSA.vector), digits = 2) #1.19

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
RMSE.SVMP_TSA2 <- RMSE(dfB2_TSA$TVC,pred_SVMP_TSA2) #RMSE =1.28

# rf plot predicted vs Observed
plot(pred_SVMP_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmPoly_FTIR RMSE:",round(RMSE.SVMP_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_SVMP_TSA100 <- train(TVC ~ ., data = train_TSA100, method = "svmPoly")
pred_SVMP_TSA100 <- predict(model_SVMP_TSA100, test_TSA100)
RMSE.SVMP_TSA100 <- RMSE(test_TSA100$TVC,pred_SVMP_TSA100)
accuracy_SVMP_TSA100 <- accuracy(pred_SVMP_TSA100, test_TSA100$TVC) 

plot(pred_SVMP_TSA100,test_TSA100$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.SVMP_TSA100,digits = 2), 
                "\nAccuracy",accuracy_SVMP_TSA100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.SVMP_TSA.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df100_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df100_TSA[trainIndex_tsa3,]
  test_tsa3 <- df100_TSA[-trainIndex_tsa3,]
  
  # train model
  model_Rf_TSA <- train(TVC~., data = train_tsa3, method = "svmPoly") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TVC,predicted.Rf_TSA)
  # put results to the vectors
  rmse.SVMP_TSA.vector100 <- c(rmse.SVMP_TSA.vector100, rmse.Rf_TSA)
}
# compute the standard deviation
sd.SVMP_TSA100 <- round(sd(rmse.SVMP_TSA.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMP_TSA100 <- CI(rmse.SVMP_TSA.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.SVMP_TSA.mean100 <- round(mean(rmse.SVMP_TSA.vector100), digits = 2) #1.19

plot(rmse.SVMP_TSA.vector100, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_TSA.mean100,"+/- ",sd.SVMP_TSA100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_TSA100[3], digits = 2),
                  "to",round(ci.SVMP_TSA100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ LM_TVC #####################################################

##################
# the first method
model_lm_TSA1 <- train(TVC ~ ., data = train_tsa2,
                       method = "lm", trControl = fitControl)
pred_lm_TSA1<- predict(model_lm_TSA1,test_tsa2)
RMSE.lm_TSA1 <- RMSE(test_tsa2$TVC,pred_lm_TSA1) #RMSE =1.4
accuracy_lm_TSA<- accuracy(pred_lm_TSA1, test_tsa2$TVC) #100

# this is the best linear model for TVC
saveRDS(model_lm_TSA1, "./models/CTF/FTIR/model_LM_TVC.rds") #save model as an rds object
saveRDS(test_tsa2$TVC, "./models/CTF/FTIR/tested_LM_TVC.rds") #save tested data as an rds object
saveRDS(pred_lm_TSA1, "./models/CTF/FTIR/predicted_LM_TVC.rds") #save predicted data as an rds object
saveRDS(RMSE.lm_TSA1, "./models/CTF/FTIR/RMSE_LM_TVC.rds") #save RMSE as an rds object
saveRDS(accuracy_lm_TSA, "./models/CTF/FTIR/accuracy_LM_TVC.rds") #save accuracy as an rds object

# rf plot predicted vs Observed
plot(pred_lm_TSA1,test_tsa2$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_TSA1,digits = 2)))
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
  trainIndex_tsa3 <- createDataPartition(df2_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df2_TSA[trainIndex_tsa3,]
  test_tsa3 <- df2_TSA[-trainIndex_tsa3,]
  # resembling (method to avoid overfitting)
  fitControl <- trainControl(method = "cv", number = 10)
  # train model
  model_lm_TSA <- train(TVC~., data = train_tsa3, method = "lm", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.lm_TSA <- predict(model_lm_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.lm_TSA <- RMSE(pred = test_tsa3$TVC,predicted.lm_TSA)
  # put results to the vectors
  rmse.lm_TSA.vector <- c(rmse.lm_TSA.vector, rmse.lm_TSA)
}
# compute the standard deviation
sd.lm_TSA <- round(sd(rmse.lm_TSA.vector), digits = 2) #0.14
# find the 95% confidence intervals parameters
ci.lm_TSA <- CI(rmse.lm_TSA.vector, ci = 0.95) #package: Rmisc #1.31, 1.27, 1.23 
# rmse mean
rmse.lm_TSA.mean <- round(mean(rmse.lm_TSA.vector), digits = 2) #1.27

plot(rmse.lm_TSA.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_TSA.mean,"+/- ",sd.lm_TSA, 
                  ")\nThere is a 95% likelihood that the range", round(ci.lm_TSA[3], digits = 2),
                  "to",round(ci.lm_TSA[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_lm_TSA2 <- train(TVC ~ ., data = dfA2_TSA,
                       method = "lm", trControl = fitControl)
pred_lm_TSA2<- predict(model_lm_TSA2,dfB2_TSA)
RMSE.lm_TSA2 <- RMSE(dfB2_TSA$TVC,pred_lm_TSA2) #RMSE =1.39

# rf plot predicted vs Observed
plot(pred_lm_TSA2,dfB2_TSA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_TSA2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_lm_TSA100 <- train(TVC ~ ., data = train_TSA100, method = "lm")
pred_lm_TSA100 <- predict(model_lm_TSA100, test_TSA100)
RMSE.lm_TSA100 <- RMSE(test_TSA100$TVC,pred_lm_TSA100)
accuracy_lm_TSA100 <- accuracy(pred_lm_TSA100, test_TSA100$TVC) 

plot(pred_lm_TSA100,test_TSA100$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.lm_TSA100,digits = 2), 
                "\nAccuracy",accuracy_lm_TSA100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.lm_TSA.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_tsa3 <- createDataPartition(df100_TSA$TVC, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_tsa3 <- df100_TSA[trainIndex_tsa3,]
  test_tsa3 <- df100_TSA[-trainIndex_tsa3,]
  
  # train model
  model_Rf_TSA <- train(TVC~., data = train_tsa3, method = "lm") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_TSA <- predict(model_Rf_TSA,test_tsa3)
  
  # calculate the RMSE value
  rmse.Rf_TSA <- RMSE(pred = test_tsa3$TVC,predicted.Rf_TSA)
  # put results to the vectors
  rmse.lm_TSA.vector100 <- c(rmse.lm_TSA.vector100, rmse.Rf_TSA)
}
# compute the standard deviation
sd.lm_TSA100 <- round(sd(rmse.lm_TSA.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.lm_TSA100 <- CI(rmse.lm_TSA.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.lm_TSA.mean100 <- round(mean(rmse.lm_TSA.vector100), digits = 2) #1.19

plot(rmse.lm_TSA.vector100, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_TSA.mean100,"+/- ",sd.lm_TSA100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.lm_TSA100[3], digits = 2),
                  "to",round(ci.lm_TSA100[1], digits = 2),"covers the true error of the model."))



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

######################### select features for the batch B #######################################
#ensure results are repeatable
set.seed(123)
#train the model
model_varImp_cfc_bB <- train(Ps~., data = dfB_CFC_wf_B, method = "rf", trControl = control, importance = TRUE)
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
trainIndex_cfc1 <- createDataPartition(df_CFC$Ps, p = .7, 
                                       list = FALSE, 
                                       times = 1, groups =3) #caret package
#split the whole dataset
train_cfc1 <- df_CFC[trainIndex_cfc1,]
test_cfc1 <- df_CFC[-trainIndex_cfc1,]

############################# Randomly split data 100 ################################
trainIndex100_CFC <- createDataPartition(df100_CFC$Ps, p = .7, list = FALSE, times = 1, groups =3)
train_CFC100 <- df100_CFC[trainIndex100_CFC,]
test_CFC100 <- df100_CFC[-trainIndex100_CFC,]

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
saveRDS(predictors_cfc_toSave, "./models/CTF/FTIR/predictors_PS.rds") #save model as an rds object

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

######################### select features for the batch B #######################################
# run the RFE algorithm
results_CFC_bB <- rfe(Ps ~ ., data = dfB_CFC_B,
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
trainIndexA_cfc1 <- createDataPartition(dfA2_CFC$Ps, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainA_cfc1 <- dfA2_CFC[trainIndexA_cfc1,]
testA_cfc1 <- dfA2_CFC[-trainIndexA_cfc1,]

#################### split batch B ##########################################################
#split data
set.seed(123) #caret package
#data split into 70% for training and 30% for testing, based on df (whole data set)
trainIndexB_cfc1 <- createDataPartition(dfB2_CFC_B$Ps, p = .7, list = FALSE, times = 1, groups =3) #caret package
#split the whole dataset
trainB_cfc1 <- dfB2_CFC_B[trainIndexB_cfc1,]
testB_cfc1 <- dfB2_CFC_B[-trainIndexB_cfc1,]

#############################################################################################
################################ RF_CFC #####################################################

##################
# the first method
model_Rf_CFC1 <- train(Ps ~ ., data = train_cfc2,
                       method = "rf", trControl = fitControl)
pred_Rf_CFC1<- predict(model_Rf_CFC1,test_cfc2)
RMSE.Rf_CFC1 <- RMSE(test_cfc2$Ps,pred_Rf_CFC1) #RMSE =1.32
# Random Forest model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_Rf_CFC <- accuracy(pred_Rf_CFC1, test_cfc2$Ps)

# this is the best rf model for Pseudomonas
saveRDS(model_Rf_CFC1, "./models/CTF/FTIR/model_RF_PS.rds") #save model as an rds object
saveRDS(test_cfc2$Ps, "./models/CTF/FTIR/tested_RF_PS.rds") #save tested data as an rds object
saveRDS(pred_Rf_CFC1, "./models/CTF/FTIR/predicted_RF_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.Rf_CFC1, "./models/CTF/FTIR/RMSE_RF_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_Rf_CFC, "./models/CTF/FTIR/accuracy_RF_PS.rds") #save accuracy as an rds object

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
ci.Rf_CFC <- CI(rmse.Rf_CFC.vector, ci = 0.95) #package: Rmisc #1.367338 1.337521 1.307704 
# rmse mean
rmse.Rf_CFC.mean <- round(mean(rmse.Rf_CFC.vector), digits = 2) #1.34

plot(rmse.Rf_CFC.vector, type="b",ylim=c(0,2),xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_CFC.mean,"+/- ",sd.Rf_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_CFC[3], digits = 2),
                  "to",round(ci.Rf_CFC[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_Rf_CFC2 <- train(Ps ~ ., data = dfA2_CFC,
                       method = "rf", trControl = fitControl)
pred_Rf_CFC2<- predict(model_Rf_CFC2,dfB2_CFC)
RMSE.Rf_CFC2 <- RMSE(dfB2_CFC$Ps,pred_Rf_CFC2) #RMSE =1.47

# rf plot predicted vs Observed
plot(pred_Rf_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("RF_FTIR RMSE:",round(RMSE.Rf_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_Rf_CFC100 <- train(Ps ~ ., data = train_CFC100, method = "rf")
pred_Rf_CFC100 <- predict(model_Rf_CFC100, test_CFC100)
RMSE.Rf_CFC100 <- RMSE(test_CFC100$Ps,pred_Rf_CFC100)
accuracy_Rf_CFC100 <- accuracy(pred_Rf_CFC100, test_CFC100$Ps) 

plot(pred_Rf_CFC100,test_CFC100$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.Rf_CFC100,digits = 2), 
                "\nAccuracy",accuracy_Rf_CFC100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.Rf_CFC.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_CFC3 <- createDataPartition(df100_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_CFC3 <- df100_CFC[trainIndex_CFC3,]
  test_CFC3 <- df100_CFC[-trainIndex_CFC3,]
  
  # train model
  model_Rf_CFC <- train(Ps~., data = train_CFC3, method = "rf") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_CFC100)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_CFC3$Ps,predicted.Rf_TSA)
  # put results to the vectors
  rmse.Rf_CFC.vector100 <- c(rmse.Rf_CFC.vector100, rmse.Rf_CFC)
}
# compute the standard deviation
sd.Rf_CFC100 <- round(sd(rmse.Rf_CFC.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.Rf_CFC100 <- CI(rmse.Rf_CFC.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.Rf_CFC.mean100 <- round(mean(rmse.Rf_CFC.vector100), digits = 2) #1.19

plot(rmse.Rf_CFC.vector100, type="b",ylim=c(0,2), xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.Rf_CFC.mean100,"+/- ",sd.Rf_CFC100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_CFC100[3], digits = 2),
                  "to",round(ci.Rf_CFC100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ KNN_CFC #####################################################

##################
# the first method
model_knn_CFC1 <- train(Ps ~ ., data = train_cfc2,
                        method = "knn", trControl = fitControl)
pred_knn_CFC1<- predict(model_knn_CFC1,test_cfc2)
RMSE.knn_CFC1 <- RMSE(test_cfc2$Ps,pred_knn_CFC1) #RMSE =1.33
# k-nearest neighbours (k-NN) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_knn_CFC <- accuracy(pred_knn_CFC1, test_cfc2$Ps)

# this is the best k-nn model for Pseudomonas
saveRDS(model_knn_CFC1, "./models/CTF/FTIR/model_knn_PS.rds") #save model as an rds object
saveRDS(test_cfc2$Ps, "./models/CTF/FTIR/tested_knn_PS.rds") #save tested data as an rds object
saveRDS(pred_knn_CFC1, "./models/CTF/FTIR/predicted_knn_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.knn_CFC1, "./models/CTF/FTIR/RMSE_knn_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_knn_CFC, "./models/CTF/FTIR/accuracy_knn_PS.rds") #save accuracy as an rds object

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
ci.knn_CFC <- CI(rmse.knn_CFC.vector, ci = 0.95) # 1.443843 1.417400 1.390957
# rmse mean
rmse.knn_CFC.mean <- round(mean(rmse.knn_CFC.vector), digits = 2) #1.42

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

# for df100
model_knn_CFC100 <- train(Ps ~ ., data = train_CFC100, method = "knn")
pred_knn_CFC100 <- predict(model_knn_CFC100, test_CFC100)
RMSE.knn_CFC100 <- RMSE(test_CFC100$Ps,pred_knn_CFC100)
accuracy_knn_CFC100 <- accuracy(pred_knn_CFC100, test_CFC100$Ps) 

plot(pred_knn_CFC100,test_CFC100$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.knn_CFC100,digits = 2), 
                "\nAccuracy",accuracy_knn_CFC100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.knn_CFC.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_CFC3 <- createDataPartition(df100_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_CFC3 <- df100_CFC[trainIndex_CFC3,]
  test_CFC3 <- df100_CFC[-trainIndex_CFC3,]
  
  # train model
  model_Rf_CFC <- train(Ps~., data = train_CFC3, method = "knn") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_CFC3)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_CFC3$Ps,predicted.Rf_TSA)
  # put results to the vectors
  rmse.knn_CFC.vector100 <- c(rmse.knn_CFC.vector100, rmse.Rf_CFC)
}
# compute the standard deviation
sd.knn_CFC100 <- round(sd(rmse.knn_CFC.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.knn_CFC100 <- CI(rmse.knn_CFC.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.knn_CFC.mean100 <- round(mean(rmse.knn_CFC.vector100), digits = 2) #1.19

plot(rmse.knn_CFC.vector100, type="l",ylim=c(0,2), xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.knn_CFC.mean100,"+/- ",sd.knn_CFC100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.knn_CFC100[3], digits = 2),
                  "to",round(ci.knn_CFC100[1], digits = 2),"covers the true error of the model."))


#############################################################################################
################################ svmLM_CFC #####################################################

##################
# the first method
model_SVMLM_CFC1 <- train(Ps ~ ., data = train_cfc2,
                          method = "svmLinear", trControl = fitControl)
pred_SVMLM_CFC1<- predict(model_SVMLM_CFC1,test_cfc2)
RMSE.SVMLM_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMLM_CFC1) #RMSE =1.24
# Linear support-vector machine (SVM) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_SVMLM_CFC <- accuracy(pred_SVMLM_CFC1, test_cfc2$Ps)

# rf plot predicted vs Observed
plot(pred_SVMLM_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("SVMLM_FTIR RMSE:",round(RMSE.SVMLM_CFC1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMLM_CFC1, "./models/CTF/FTIR/model_SVMLM_PS.rds") #save model as an rds object
saveRDS(test_cfc2$Ps, "./models/CTF/FTIR/tested_SVMLM_PS.rds") #save tested data as an rds object
saveRDS(pred_SVMLM_CFC1, "./models/CTF/FTIR/predicted_SVMLM_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMLM_CFC1, "./models/CTF/FTIR/RMSE_SVMLM_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMLM_CFC, "./models/CTF/FTIR/accuracy_SVMLM_PS.rds") #save accuracy as an rds object

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
sd.SVMLM_CFC <- round(sd(rmse.SVMLM_CFC.vector), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMLM_CFC <- CI(rmse.SVMLM_CFC.vector, ci = 0.95) # 1.09, 1.06, 1.04 
# rmse mean
rmse.SVMLM_CFC.mean <- round(mean(rmse.SVMLM_CFC.vector), digits = 2) #1.06

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
RMSE.SVMLM_CFC2 <- RMSE(dfB2_CFC$Ps,pred_SVMLM_CFC2) #RMSE =1.69

# rf plot predicted vs Observed
plot(pred_SVMLM_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmLinear_FTIR RMSE:",round(RMSE.SVMLM_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_SVMLM_CFC100 <- train(Ps ~ ., data = train_CFC100, method = "svmLinear")
pred_SVMLM_CFC100 <- predict(model_SVMLM_CFC100, test_CFC100)
RMSE.SVMLM_CFC100 <- RMSE(test_CFC100$Ps,pred_SVMLM_CFC100)
accuracy_SVMLM_CFC100 <- accuracy(pred_SVMLM_CFC100, test_CFC100$Ps) 

plot(pred_SVMLM_CFC100,test_CFC100$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.SVMLM_CFC100,digits = 2), 
                "\nAccuracy",accuracy_SVMLM_CFC100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.SVMLM_CFC.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_CFC3 <- createDataPartition(df100_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_CFC3 <- df100_CFC[trainIndex_CFC3,]
  test_CFC3 <- df100_CFC[-trainIndex_CFC3,]
  
  # train model
  model_Rf_CFC <- train(Ps~., data = train_CFC3, method = "svmLinear") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_CFC3)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_CFC3$Ps,predicted.Rf_TSA)
  # put results to the vectors
  rmse.SVMLM_CFC.vector100 <- c(rmse.SVMLM_CFC.vector100, rmse.Rf_CFC)
}
# compute the standard deviation
sd.SVMLM_CFC100 <- round(sd(rmse.SVMLM_CFC.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMLM_CFC100 <- CI(rmse.SVMLM_CFC.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.SVMLM_CFC.mean100 <- round(mean(rmse.SVMLM_CFC.vector100), digits = 2) #1.19

plot(rmse.SVMLM_CFC.vector100, type="l",ylim=c(0,2), xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMLM_CFC.mean100,"+/- ",sd.SVMLM_CFC100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMLM_CFC100[3], digits = 2),
                  "to",round(ci.SVMLM_CFC100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ SVMRadial_CFC #####################################################

##################
# the first method
model_SVMR_CFC1 <- train(Ps ~ ., data = train_cfc2,
                         method = "svmRadial", trControl = fitControl)
pred_SVMR_CFC1<- predict(model_SVMR_CFC1,test_cfc2)
RMSE.SVMR_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMR_CFC1) #RMSE =1.29
# Radial support-vector machine (SVM) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_SVMR_CFC <- accuracy(pred_SVMR_CFC1, test_cfc2$Ps)

# rf plot predicted vs Observed
plot(pred_SVMR_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_CFC1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMR_CFC1, "./models/CTF/FTIR/model_SVMR_PS.rds") #save model as an rds object
saveRDS(test_cfc2$Ps, "./models/CTF/FTIR/tested_SVMR_PS.rds") #save tested data as an rds object
saveRDS(pred_SVMR_CFC1, "./models/CTF/FTIR/predicted_SVMR_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMR_CFC1, "./models/CTF/FTIR/RMSE_SVMR_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMR_CFC, "./models/CTF/FTIR/accuracy_SVMR_PS.rds") #save accuracy as an rds object

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
ci.SVMR_CFC <- CI(rmse.SVMR_CFC.vector, ci = 0.95) # 1.34, 1.31, 1.28 
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
RMSE.SVMR_CFC2 <- RMSE(dfB2_CFC$Ps,pred_SVMR_CFC2) #RMSE =1.5

# rf plot predicted vs Observed
plot(pred_SVMR_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("svmRadial_FTIR RMSE:",round(RMSE.SVMR_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_SVMR_CFC100 <- train(Ps ~ ., data = train_CFC100, method = "svmRadial")
pred_SVMR_CFC100 <- predict(model_SVMR_CFC100, test_CFC100)
RMSE.SVMR_CFC100 <- RMSE(test_CFC100$Ps,pred_SVMR_CFC100)
accuracy_SVMR_CFC100 <- accuracy(pred_SVMR_CFC100, test_CFC100$Ps) 

plot(pred_SVMR_CFC100,test_CFC100$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.SVMR_CFC100,digits = 2), 
                "\nAccuracy",accuracy_SVMR_CFC100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.SVMR_CFC.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_CFC3 <- createDataPartition(df100_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_CFC3 <- df100_CFC[trainIndex_CFC3,]
  test_CFC3 <- df100_CFC[-trainIndex_CFC3,]
  
  # train model
  model_Rf_CFC <- train(Ps~., data = train_CFC3, method = "svmRadial") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_CFC3)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_CFC3$Ps,predicted.Rf_TSA)
  # put results to the vectors
  rmse.SVMR_CFC.vector100 <- c(rmse.SVMR_CFC.vector100, rmse.Rf_CFC)
}
# compute the standard deviation
sd.SVMR_CFC100 <- round(sd(rmse.SVMR_CFC.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMR_CFC100 <- CI(rmse.SVMR_CFC.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.SVMR_CFC.mean100 <- round(mean(rmse.SVMR_CFC.vector100), digits = 2) #1.19

plot(rmse.SVMR_CFC.vector100, type="l",ylim=c(0,2), xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMR_CFC.mean100,"+/- ",sd.SVMR_CFC100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMR_CFC100[3], digits = 2),
                  "to",round(ci.SVMR_CFC100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ SVMPoly_CFC #####################################################

##################
# the first method
model_SVMP_CFC1 <- train(Ps ~ ., data = train_cfc2,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC1<- predict(model_SVMP_CFC1,test_cfc2)
RMSE.SVMP_CFC1 <- RMSE(test_cfc2$Ps,pred_SVMP_CFC1) #RMSE =1.26
# Polynomial support-vector machine (SVM) model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_SVMP_CFC <- accuracy(pred_SVMP_CFC1, test_cfc2$Ps)

# rf plot predicted vs Observed
plot(pred_SVMP_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC1,digits = 2), 
                "\nAccuracy",accuracy_SVMP_CFC,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_SVMP_CFC1, "./models/CTF/FTIR/model_SVMP_PS.rds") #save model as an rds object
saveRDS(test_cfc2$Ps, "./models/CTF/FTIR/tested_SVMP_PS.rds") #save tested data as an rds object
saveRDS(pred_SVMP_CFC1, "./models/CTF/FTIR/predicted_SVMP_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.SVMP_CFC1, "./models/CTF/FTIR/RMSE_SVMP_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_SVMP_CFC, "./models/CTF/FTIR/accuracy_SVMP_PS.rds") #save accuracy as an rds object

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
sd.SVMP_CFC <- round(sd(rmse.SVMP_CFC.vector), digits = 2) #0.2
# find the 95% confidence intervals parameters
ci.SVMP_CFC <- CI(rmse.SVMP_CFC.vector, ci = 0.95) # 1.25, 1.20, 1.14  
# rmse mean
rmse.SVMP_CFC.mean <- round(mean(rmse.SVMP_CFC.vector), digits = 2) #1.2

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
RMSE.SVMP_CFC2 <- RMSE(dfB2_CFC$Ps,pred_SVMP_CFC2) #RMSE =1.53

# rf plot predicted vs Observed
plot(pred_SVMP_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

###################
# trained based on B, tested based on A
model_SVMP_CFC3 <- train(Ps ~ ., data = dfB2_CFC_B,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC3<- predict(model_SVMP_CFC3,dfA2_CFC_B)
RMSE.SVMP_CFC3 <- RMSE(dfA2_CFC_B$Ps,pred_SVMP_CFC3) 

# rf plot predicted vs Observed
plot(pred_SVMP_CFC3,dfA2_CFC_B$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC3,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##################
# trained and tested based on batch A
model_SVMP_CFC4 <- train(Ps ~ ., data = trainA_cfc1,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC4<- predict(model_SVMP_CFC4,testA_cfc1)
RMSE.SVMP_CFC4 <- RMSE(testA_cfc1$Ps,pred_SVMP_CFC4) 

# plot predicted vs Observed
plot(pred_SVMP_CFC4,testA_cfc1$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC4,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

#######################
# trained and tested based on batch B
model_SVMP_CFC5 <- train(Ps ~ ., data = trainB_cfc1,
                         method = "svmPoly", trControl = fitControl)
pred_SVMP_CFC5<- predict(model_SVMP_CFC5,testB_cfc1)
RMSE.SVMP_CFC5 <- RMSE(testB_cfc1$Ps,pred_SVMP_CFC5) 

# plot predicted vs Observed
plot(pred_SVMP_CFC5,testB_cfc1$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 Pseudomonas spp. counts/g",
     ylab="Actual log10 Pseudomonas spp. counts/g",col = "blue", 
     main=paste("Polynomial support-vector machine (SVM) model for Pseudomonas spp. counts \nRMSE:",
                round(RMSE.SVMP_CFC5,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_SVMP_CFC100 <- train(Ps ~ ., data = train_CFC100, method = "svmPoly")
pred_SVMP_CFC100 <- predict(model_SVMP_CFC100, test_CFC100)
RMSE.SVMP_CFC100 <- RMSE(test_CFC100$Ps,pred_SVMP_CFC100)
accuracy_SVMP_CFC100 <- accuracy(pred_SVMP_CFC100, test_CFC100$Ps) 

plot(pred_SVMP_CFC100,test_CFC100$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.SVMP_CFC100,digits = 2), 
                "\nAccuracy",accuracy_SVMP_CFC100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.SVMP_CFC.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_CFC3 <- createDataPartition(df100_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_CFC3 <- df100_CFC[trainIndex_CFC3,]
  test_CFC3 <- df100_CFC[-trainIndex_CFC3,]
  
  # train model
  model_Rf_CFC <- train(Ps~., data = train_CFC3, method = "svmPoly") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_CFC3)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_CFC3$Ps,predicted.Rf_TSA)
  # put results to the vectors
  rmse.SVMP_CFC.vector100 <- c(rmse.SVMP_CFC.vector100, rmse.Rf_CFC)
}
# compute the standard deviation
sd.SVMP_CFC100 <- round(sd(rmse.SVMP_CFC.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.SVMP_CFC100 <- CI(rmse.SVMP_CFC.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.SVMP_CFC.mean100 <- round(mean(rmse.SVMP_CFC.vector100), digits = 2) #1.19

plot(rmse.SVMP_CFC.vector100, type="l",ylim=c(0,2), xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.SVMP_CFC.mean100,"+/- ",sd.SVMP_CFC100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.SVMP_CFC100[3], digits = 2),
                  "to",round(ci.SVMP_CFC100[1], digits = 2),"covers the true error of the model."))

#############################################################################################
################################ LM_CFC #####################################################

##################
# the first method
model_lm_CFC1 <- train(Ps ~ ., data = train_cfc2,
                       method = "lm", trControl = fitControl)
pred_lm_CFC1<- predict(model_lm_CFC1,test_cfc2)
RMSE.lm_CFC1 <- RMSE(test_cfc2$Ps,pred_lm_CFC1) #RMSE =1.4
# Linear model Accuracy for Pseudomonas spp. colony forming units (CFU) counts
accuracy_lm_CFC <- accuracy(pred_lm_CFC1, test_cfc2$Ps)

# rf plot predicted vs Observed
plot(pred_lm_CFC1,test_cfc2$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_TSA1,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(model_lm_CFC1, "./models/CTF/FTIR/model_LM_PS.rds") #save model as an rds object
saveRDS(test_cfc2$Ps, "./models/CTF/FTIR/tested_LM_PS.rds") #save tested data as an rds object
saveRDS(pred_lm_CFC1, "./models/CTF/FTIR/predicted_LM_PS.rds") #save predicted data as an rds object
saveRDS(RMSE.lm_CFC1, "./models/CTF/FTIR/RMSE_LM_PS.rds") #save RMSE as an rds object
saveRDS(accuracy_lm_CFC, "./models/CTF/FTIR/accuracy_LM_PS.rds") #save accuracy as an rds object

##################
# 50 iterations
rmse.lm_CFC.vector <- c()

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
  model_lm_CFC <- train(Ps~., data = train_cfc3, method = "lm", trControl = fitControl) #package: caret
  # predicted bacterial counts for the test set
  predicted.lm_CFC <- predict(model_lm_CFC,test_cfc3)
  
  # calculate the RMSE value
  rmse.lm_CFC <- RMSE(pred = test_cfc3$Ps,predicted.lm_CFC)
  # put results to the vectors
  rmse.lm_CFC.vector <- c(rmse.lm_CFC.vector, rmse.lm_CFC)
}
# compute the standard deviation
sd.lm_CFC <- round(sd(rmse.lm_CFC.vector), digits = 2) #0.53
# find the 95% confidence intervals parameters
ci.lm_CFC <- CI(rmse.lm_CFC.vector, ci = 0.95) #package: Rmisc #2.77, 2.62, 2.47
# rmse mean
rmse.lm_CFC.mean <- round(mean(rmse.lm_CFC.vector), digits = 2) #2.62

plot(rmse.lm_CFC.vector, type="b",xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_CFC.mean,"+/- ",sd.lm_CFC, 
                  ")\nThere is a 95% likelihood that the range", round(ci.Rf_CFC[3], digits = 2),
                  "to",round(ci.lm_CFC[1], digits = 2),"covers the true error of the model."))

##################
# the second method
# the first method
model_lm_CFC2 <- train(Ps ~ ., data = dfA2_CFC,
                       method = "lm", trControl = fitControl)
pred_lm_CFC2<- predict(model_lm_CFC2,dfB2_CFC)
RMSE.lm_CFC2 <- RMSE(dfB2_CFC$Ps,pred_lm_CFC2) #RMSE =1.48

# rf plot predicted vs Observed
plot(pred_lm_CFC2,dfB2_CFC$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("LM_FTIR RMSE:",round(RMSE.lm_CFC2,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# for df100
model_lm_CFC100 <- train(Ps ~ ., data = train_CFC100, method = "lm")
pred_lm_CFC100 <- predict(model_lm_CFC100, test_CFC100)
RMSE.lm_CFC100 <- RMSE(test_CFC100$Ps,pred_lm_CFC100)
accuracy_lm_CFC100 <- accuracy(pred_lm_CFC100, test_CFC100$Ps) 

plot(pred_lm_CFC100,test_CFC100$Ps,xlim= c(0,9), ylim=c(0,9),xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forest model for Total Viable Counts \nRMSE:",round(RMSE.lm_CFC100,digits = 2), 
                "\nAccuracy",accuracy_lm_CFC100,"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

# 500 iterations for df100
rmse.lm_CFC.vector100 <- c()

for(i in 1:50){
  # split the data into a training set (70% of observations) and a test set(30%) 
  set.seed(i)
  #data split into 70% for training and 30% for testing, based on df (whole data set)
  trainIndex_CFC3 <- createDataPartition(df100_CFC$Ps, p = .7, 
                                         list = FALSE, 
                                         times = 1, groups =3) #caret package
  #split the whole dataset
  train_CFC3 <- df100_CFC[trainIndex_CFC3,]
  test_CFC3 <- df100_CFC[-trainIndex_CFC3,]
  
  # train model
  model_Rf_CFC <- train(Ps~., data = train_CFC3, method = "lm") #package: caret
  # predicted bacterial counts for the test set
  predicted.Rf_CFC <- predict(model_Rf_CFC,test_CFC3)
  
  # calculate the RMSE value
  rmse.Rf_CFC <- RMSE(pred = test_CFC3$Ps,predicted.Rf_TSA)
  # put results to the vectors
  rmse.lm_CFC.vector100 <- c(rmse.lm_CFC.vector100, rmse.Rf_CFC)
}
# compute the standard deviation
sd.lm_CFC100 <- round(sd(rmse.lm_CFC.vector100), digits = 2) #0.09
# find the 95% confidence intervals parameters
ci.lm_CFC100 <- CI(rmse.lm_CFC.vector100, ci = 0.95) #package: Rmisc #1.215049 1.188691 1.162332 
# rmse mean
rmse.lm_CFC.mean100 <- round(mean(rmse.lm_CFC.vector100), digits = 2) #1.19

plot(rmse.lm_CFC.vector100, type="l",ylim=c(0,2), xlab="iteration",ylab="RMSE", 
     main = paste("RMSE for 50 iterations \n( RMSE.mean = ",rmse.lm_CFC.mean100,"+/- ",sd.lm_CFC100, 
                  ")\nThere is a 95% likelihood that the range", round(ci.lm_CFC100[3], digits = 2),
                  "to",round(ci.lm_CFC100[1], digits = 2),"covers the true error of the model."))

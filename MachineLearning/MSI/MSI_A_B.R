##############################
#                            #
#Group_project MSI Raw Data  #
#   Based on batches         #
#   Models trained on A      #
#   Tested on B              #
#   Machine Learning         #
##############################
# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()


library(caret)
library(kernlab)
library(openxlsx)
source("CreateObj.R")
source("Percentage_Accuracy.R")
# Read datasets
MSI   <- read.xlsx("MSIdata_CT_R1_2_updated.xlsx", sheet=1, rowName=T)
#create two datasets based on batches 
bacteria_A<- MSI[1:100,37:38] ## Bacterial count for batch A
bacteria_B<- MSI[101:196,37:38]## Bacterial count for batch A
MSI_A <- MSI[1:100,1:18] ### Bath A
MSI_B <- MSI[101:196,1:18] ### Bath B

MSI_raw_A<- CreateObj(MSI_A, bacteria_A)##Merging bacterial counts with spectra of batch A
MSI_raw_B<- CreateObj(MSI_B, bacteria_B)##Merging bacterial counts with spectra of batch B

MSI_TVC_A <- MSI_raw_A[,-(20)]###remove pseudomnas A
MSI_PS_A <- MSI_raw_A[,-19]###remove TVC A

MSI_TVC_B <- MSI_raw_B[,-(20)]###remove pseudomnas B
MSI_PS_B <- MSI_raw_B[,-19]###remove TVC B

###############################################################################################################################################################################################
##################################################################TVC_Regression###########################################################################################################
# Split training and test set for TVC count based on batch A
set.seed(123)

train.index <- createDataPartition(MSI_TVC_A$TVC, p = 1,list = FALSE, times = 1)
trainNA <- MSI_TVC_A[train.index,]

train.index <- createDataPartition(MSI_TVC_B$TVC, p = 1,list = FALSE, times = 1)
testNB <- MSI_TVC_B[train.index,]
###Control function
control <- trainControl(method="repeatedcv", number=10, repeats=3)

##################################################################################
################### Linear Model for TVC count by MSI_data   #####################
##################################################################################

lm.model.fit.TVC <- train(TVC ~., method= 'lm', data=trainNA, trControl= control)##Training model
predicted.lm <- predict(lm.model.fit.TVC, testNB[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.lm.TVC <- RMSE(testNB$TVC, predicted.lm)##RMSE calculation
##Linear model Accuracy for TVC
percentage_lm_TVC<- Percentage(predicted.lm,testNB$TVC)
# Linear Model plot predicted values vs Actual/Observed values
plot(predicted.lm, testNB$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Total Viable Counts(TVC) \nRMSE::",round(RMSE.lm.TVC,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
##################################################################################
################ k-nearest neighbours for TVC count by  MSI_data  ################
##################################################################################

knn.model.fit.TVC <- train(TVC ~ ., method='knn', data=trainNA, trcontrol = control, tuneGrid=expand.grid(k=1:20))##Training model
predicted.knn <- predict(knn.model.fit.TVC,testNB[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.knn <- RMSE(testNB$TVC,predicted.knn)##RMSE calculation
##K-nearest neighbours model Accuracy for TVC
percentage_knn_TVC<- Percentage(predicted.knn,testNB$TVC)
# kNN plot predicted values vs Actual/Observed values
plot(predicted.knn,testNB$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)",col = "blue", 
     main=paste("k-nearest neighbours Model for TVC - RMSE:",round(RMSE.knn,digits = 2)," Precentage :",round(percentage_knn_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

###################################################################################
################### Random Forests for TVC count by MSI_data   ####################
###################################################################################

RF.model.fit.TVC <- train(TVC ~ ., method='rf', trcontrol = control, trainNA,tuneGrid=expand.grid(mtry=1:20))##Training model
predicted.rf <- predict(RF.model.fit.TVC, testNB[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.rf <- RMSE(testNB$TVC,predicted.rf)##RMSE calculation
##randomForest model Accuracy for TVC
percentage_rf_TVC<- Percentage(predicted.rf,testNB$TVC)
## RandomForest plot predicted values vs Actual/Observed values
plot(predicted.rf, testNB$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Random Forests Model for TVC - RMSE:",round(RMSE.rf,digits = 2)," Precentage :",round(percentage_rf_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVMPloy for TVC count by MSI_data          ####################
###################################################################################

SVM.model.fit.TVC <- train(TVC ~ ., method='svmPoly', data= trainNA,tuneLength = 5)##Training model
predicted.SVMP <- predict(SVM.model.fit.TVC,testNB[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMP <- RMSE(testNB$TVC, predicted.SVMP)##calculating RMSE
##SVM Poly model Accuracy for TVC
percentage_SVMPoly_TVC<- Percentage(predicted.SVMP, testNB$TVC)
#SVMPloy plot predicted values vs Actual/Observed values
plot(predicted.SVMP,testNB$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynomial SVM Model for TVC - RMSE:",round(RMSE.SVMP,digits = 2),"Precentage :",round(percentage_SVMPoly_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
#########################SVMRadial for TVC count by MSI_data#######################
###################################################################################


regressor = train(TVC ~ .,
                data = trainNA,
                method = 'svmRadial',
                trcontrol = control)##Training model
y_pred = predict(regressor, testNB[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMR.TVC <- RMSE(testNB$TVC, y_pred)##calculating RMSE
##SVM Poly model Accuracy for TVC
percentage_SVMRadial_TVC<- Percentage(y_pred, testNB$TVC)
#SVMRadial plot predicted values vs Actual/Observed values
plot(y_pred,testNB$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for TVC RMSE:",round(RMSE.SVMR.TVC,digits = 2),"Precentage :",round(percentage_SVMRadial_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
############################################################################################
###################################SVM linear for TVC count MSI data########################
############################################################################################
SVMLM.model.fit.TVC <- train(TVC ~ ., method='svmLinear', data= trainNA)##Training model
predicted.SVMLM.TVC <- predict(SVMLM.model.fit.TVC,testNB[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMLM.TVC <- RMSE(testNB$TVC, predicted.SVMLM.TVC)##RMSE calculation
##SVM Linear model Accuracy for TVC
percentage_SVMLM_TVC<- Percentage(predicted.SVMLM.TVC, testNB$TVC)
# SVM linear plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.TVC,xlim= c(0,9), ylim=c(0,9),testNB$TVC,xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear support-vector machine(SVM) Model for Total Viable Counts(TVC) \nRMSE:", round(RMSE.SVMLM.TVC,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
####################################################################################################################################################
###########################################Pseudomonas MSI_data#####################################################################################

set.seed(123)

train.index <- createDataPartition(MSI_PS_A$Pseudomonas.spp., p = 1,list = FALSE, times = 1)
train.index1 <- createDataPartition(MSI_PS_B$Pseudomonas.spp., p = 1,list = FALSE, times = 1)
trainNa <- MSI_PS_A[train.index,]
testNb <- MSI_PS_B[train.index1,]
##################################################################################
################### Linear Model for TVC count by MSI_data  #####################
##################################################################################

lm.model.fit.ps <- train(Pseudomonas.spp. ~ ., method= 'lm', data=trainNa, trcontrol= control) ##Training model
predicted.lm.ps <- predict(lm.model.fit.ps,testNb[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.lm.ps <- RMSE(testNb$Pseudomonas.spp., predicted.lm.ps)##RMSE calculation
##Linear model Accuracy for Pseudomonas
percentage_LM_Pseudomonas<- Percentage(predicted.lm.ps, testNb$Pseudomonas.spp.)
##Linear Model plot predicted values vs Actual/Observed values
plot(predicted.lm.ps, testNb$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Pseudomonas \nRMSE:", round(RMSE.lm.ps,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### Random Forests for TVC count by MSI_data   ####################
###################################################################################

RF.model.fit.Ps <- train(Pseudomonas.spp. ~ ., method='rf', trcontrol = control, trainNa)##Training model
predicted.rf.ps <- predict(RF.model.fit.Ps, testNb[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.rf.ps <- RMSE(testNb$Pseudomonas.spp.,predicted.rf.ps)##RMSE calculation
##RandomForest model Accuracy for Pseudomonas
percentage_rf_Pseudomonas<- Percentage(predicted.rf.ps, testNb$Pseudomonas.spp.)
# RandomForest plot predicted values vs Actual/Observed values
plot(predicted.rf.ps,testNb$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Random Forest Model for Pseudomonas -  RMSE:",round(RMSE.rf.ps,digits = 2),"Precentage :",round(percentage_rf_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVMPloy for Pseudomonas count by MSI_data          ############
###################################################################################

SVM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmPoly', data= trainNa)
predicted.SVMP.PS <- predict(SVM.model.fit.PS,testNb[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMP.ps <- RMSE(testNb$Pseudomonas.spp., predicted.SVMP.PS)##RMSE calculation
##SVM Ploy model Accuracy for Pseudomonas
percentage_SVMP_Pseudomonas<- Percentage(predicted.SVMP.PS, testNb$Pseudomonas.spp.)
# SVM Poly plot predicted values vs Actual/Observed values
plot(predicted.SVMP.PS,testNb$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynomial SVM Model for Pseudomonas -  RMSE:",round(RMSE.SVMP.ps,digits = 2),"Precentage :",round(percentage_SVMP_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVM linear for Pseudomonas count by MSI_data       ############
###################################################################################

SVMLM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmLinear', data= trainNa)##Training model
predicted.SVMLM.PS <- predict(SVMLM.model.fit.PS,testNb[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMLM.ps <- RMSE(testNb$Pseudomonas.spp., predicted.SVMLM.PS)##RMSE calculation
##SVM Linear model Accuracy for Pseudomonas
percentage_SVMLM_Pseudomonas<- Percentage(predicted.SVMLM.PS, testNb$Pseudomonas.spp.)
# SVM linear plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.PS,testNb$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear support-vector machine(SVM) Model for Pseudomonas.spp \nRMSE:",round(RMSE.SVMLM.ps,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
####################################################################################
################ k-nearest neighbours for Pseudomonas count by  MSI_data ###########
####################################################################################

knn.model.fit.ps <- train(Pseudomonas.spp. ~ ., method='knn', data=trainNa, trcontrol = control, tuneGrid=expand.grid(k=1:20))##Training model
predicted.knn.ps <- predict(knn.model.fit.ps, testNb[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.knn.ps <- RMSE(testNb$Pseudomonas.spp.,predicted.knn.ps)##RMSE calculation
## k-nearest neighbours model Accuracy for Pseudomonas
percentage_knn_Pseudomonas<- Percentage(predicted.knn.ps, testNb$Pseudomonas.spp.)
# K-nearest neighbours plot predicted values vs Actual/Observed values
plot(predicted.knn.ps, testNb$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("k-nearest neighbours Model for Pseudomonas -  RMSE:",round(RMSE.knn.ps,digits = 2),"Precentage :",round(percentage_SVMLM_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
#########################SVM Radial for PS count by MSI_data ######################
###################################################################################


regressor.ps = train(Pseudomonas.spp. ~ .,
                   data = trainNa,
                   method = 'svmRadial',
                   trcontrol = control)##Training model
y_pred.ps = predict(regressor.ps, testNb[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMR.ps <- RMSE(testNb$Pseudomonas.spp., y_pred.ps)##RMSE calculation
## k-nearest neighbours model Accuracy for Pseudomonas
percentage_SVMR_Pseudomonas<- Percentage(y_pred.ps, testNb$Pseudomonas.spp.)
#SVM Radial plot predicted values vs Actual/Observed values
plot(y_pred.ps, testNb$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for Pseudomonas -  RMSE:",round(RMSE.SVMR.ps,digits = 2),"Precentage :",round(percentage_SVMR_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

##############################
#                            #
#Group_project MSI Raw Data  #
#   Based on batches         #
#   Models trained on B      #
#   Tested on A              #
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

train.indexS <- createDataPartition(MSI_TVC_B$TVC, p = 1,list = FALSE, times = 1)
trainNB <- MSI_TVC_B[train.indexS,]

train.index <- createDataPartition(MSI_TVC_A$TVC, p = 1,list = FALSE, times = 1)
testNA <- MSI_TVC_A[train.index,]
###Control function
control <- trainControl(method="repeatedcv", number=10, repeats=3)

##################################################################################
################### Linear Model for TVC count by MSI_data   #####################
##################################################################################

lm.model.fit.TVC <- train(TVC ~., method= 'lm', data=trainNB, trControl= control)##Training model
predicted.lm <- predict(lm.model.fit.TVC, testNA[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.lm.TVC <- RMSE(testNA$TVC, predicted.lm)##RMSE calculation
##Linear model Accuracy for TVC
percentage_lm_TVC<- Percentage(predicted.lm,testNA$TVC)
# Linear Model plot predicted values vs Actual/Observed values
plot(predicted.lm, testNA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.lm.TVC,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
##################################################################################
################ k-nearest neighbours for TVC count by  MSI_data  ################
##################################################################################

knn.model.fit.TVC <- train(TVC ~ ., method='knn', data=trainNB, trcontrol = control, tuneGrid=expand.grid(k=1:20))##Training model
predicted.knn <- predict(knn.model.fit.TVC,testNA[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.knn <- RMSE(testNA$TVC,predicted.knn)##RMSE calculation
##K-nearest neighbours model Accuracy for TVC
percentage_knn_TVC<- Percentage(predicted.knn,testNA$TVC)
# kNN plot predicted values vs Actual/Observed values
plot(predicted.knn,testNA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)",col = "blue", 
     main=paste("k-nearest neighbours Model for TVC \nRMSE:",round(RMSE.knn,digits = 2)," Accuracy :",round(percentage_knn_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

###################################################################################
################### Random Forests for TVC count by MSI_data   ####################
###################################################################################

RF.model.fit.TVC <- train(TVC ~ ., method='rf', trcontrol = control, trainNB,tuneGrid=expand.grid(mtry=1:20))##Training model
predicted.rf <- predict(RF.model.fit.TVC, testNA[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.rf <- RMSE(testNA$TVC,predicted.rf)##RMSE calculation
##randomForest model Accuracy for TVC
percentage_rf_TVC<- Percentage(predicted.rf,testNA$TVC)
## RandomForest plot predicted values vs Actual/Observed values
plot(predicted.rf, testNA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Random Forests Model for TVC \nRMSE:",round(RMSE.rf,digits = 2)," Accuracy :",round(percentage_rf_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVMPloy for TVC count by MSI_data          ####################
###################################################################################

SVM.model.fit.TVC <- train(TVC ~ ., method='svmPoly', data= trainNB,tuneLength = 5)##Training model
predicted.SVMP <- predict(SVM.model.fit.TVC,testNA[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMP <- RMSE(testNA$TVC, predicted.SVMP)##calculating RMSE
##SVM Poly model Accuracy for TVC
percentage_SVMPoly_TVC<- Percentage(predicted.SVMP, testNA$TVC)
#SVMPloy plot predicted values vs Actual/Observed values
plot(predicted.SVMP,testNA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynomial SVM Model for TVC \nRMSE:",round(RMSE.SVMP,digits = 2),"Accuracy :",round(percentage_SVMPoly_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
#########################SVMRadial for TVC count by MSI_data#######################
###################################################################################


regressor = train(TVC ~ .,
                data = trainNB,
                method = 'svmRadial',
                trcontrol = control)##Training model
y_pred = predict(regressor, testNA[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMR.TVC <- RMSE(testNA$TVC, y_pred)##calculating RMSE
##SVM Poly model Accuracy for TVC
percentage_SVMRadial_TVC<- Percentage(y_pred, testNA$TVC)
#SVMRadial plot predicted values vs Actual/Observed values
plot(y_pred,testNA$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for TVC \nRMSE:",round(RMSE.SVMR.TVC,digits = 2),"Accuracy :",round(percentage_SVMRadial_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
############################################################################################
###################################SVM linear for TVC count MSI data########################
############################################################################################
SVMLM.model.fit.TVC <- train(TVC ~ ., method='svmLinear', data= trainNB)##Training model
predicted.SVMLM.TVC <- predict(SVMLM.model.fit.TVC,testNA[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMLM.TVC <- RMSE(testNA$TVC, predicted.SVMLM.TVC)##RMSE calculation
##SVM Linear model Accuracy for TVC
percentage_SVMLM_TVC<- Percentage(predicted.SVMLM.TVC, testNA$TVC)
# SVM linear plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.TVC,xlim= c(0,9), ylim=c(0,9),testNA$TVC,xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear SVM Model for TVC \nRMSE:",round(RMSE.SVMLM.TVC,digits = 2),"Accuracy :",round(percentage_SVMLM_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
####################################################################################################################################################
###########################################Pseudomonas MSI_data#####################################################################################

set.seed(123)

train.index <- createDataPartition(MSI_PS_B$Pseudomonas.spp., p = 1,list = FALSE, times = 1)
train.index1 <- createDataPartition(MSI_PS_A$Pseudomonas.spp., p = 1,list = FALSE, times = 1)
trainNb <- MSI_PS_B[train.index,]
testNa <- MSI_PS_A[train.index1,]
##################################################################################
################### Linear Model for TVC count by MSI_data  #####################
##################################################################################

lm.model.fit.ps <- train(Pseudomonas.spp. ~ ., method= 'lm', data=trainNb, trcontrol= control) ##Training model
predicted.lm.ps <- predict(lm.model.fit.ps,testNa[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.lm.ps <- RMSE(testNa$Pseudomonas.spp., predicted.lm.ps)##RMSE calculation
##Linear model Accuracy for Pseudomonas
percentage_LM_Pseudomonas<- Percentage(predicted.lm.ps, testNa$Pseudomonas.spp.)
##Linear Model plot predicted values vs Actual/Observed values
plot(predicted.lm.ps, testNa$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Pseudomonas.spp \nRMSE:", round(RMSE.lm.ps,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### Random Forests for TVC count by MSI_data   ####################
###################################################################################

RF.model.fit.Ps <- train(Pseudomonas.spp. ~ ., method='rf', trcontrol = control, trainNb)##Training model
predicted.rf.ps <- predict(RF.model.fit.Ps, testNa[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.rf.ps <- RMSE(testNa$Pseudomonas.spp.,predicted.rf.ps)##RMSE calculation
##RandomForest model Accuracy for Pseudomonas
percentage_rf_Pseudomonas<- Percentage(predicted.rf.ps, testNa$Pseudomonas.spp.)
# RandomForest plot predicted values vs Actual/Observed values
plot(predicted.rf.ps,testNa$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Random Forest Model for Pseudomonas \nRMSE:",round(RMSE.rf.ps,digits = 2),"Accuracy :",round(percentage_rf_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVMPloy for Pseudomonas count by MSI_data          ############
###################################################################################

SVM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmPoly', data= trainNb)
predicted.SVMP.PS <- predict(SVM.model.fit.PS,testNa[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMP.ps <- RMSE(testNa$Pseudomonas.spp., predicted.SVMP.PS)##RMSE calculation
##SVM Ploy model Accuracy for Pseudomonas
percentage_SVMP_Pseudomonas<- Percentage(predicted.SVMP.PS, testNa$Pseudomonas.spp.)
# SVM Poly plot predicted values vs Actual/Observed values
plot(predicted.SVMP.PS,testNa$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynomial SVM Model for Pseudomonas \nRMSE:",round(RMSE.SVMP.ps,digits = 2),"Accuracy :",round(percentage_SVMP_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVM linear for Pseudomonas count by MSI_data       ############
###################################################################################

SVMLM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmLinear', data= trainNb)##Training model
predicted.SVMLM.PS <- predict(SVMLM.model.fit.PS,testNa[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMLM.ps <- RMSE(testNa$Pseudomonas.spp., predicted.SVMLM.PS)##RMSE calculation
##SVM Linear model Accuracy for Pseudomonas
percentage_SVMLM_Pseudomonas<- Percentage(predicted.SVMLM.PS, testNa$Pseudomonas.spp.)
# SVM linear plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.PS,testNa$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("linear SVM Model for Pseudomonas \nRMSE:",round(RMSE.SVMLM.ps,digits = 2),"Accuracy :",round(percentage_SVMLM_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
####################################################################################
################ k-nearest neighbours for Pseudomonas count by  MSI_data ###########
####################################################################################

knn.model.fit.ps <- train(Pseudomonas.spp. ~ ., method='knn', data=trainNb, trcontrol = control, tuneGrid=expand.grid(k=1:20))##Training model
predicted.knn.ps <- predict(knn.model.fit.ps, testNa[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.knn.ps <- RMSE(testNa$Pseudomonas.spp.,predicted.knn.ps)##RMSE calculation
## k-nearest neighbours model Accuracy for Pseudomonas
percentage_knn_Pseudomonas<- Percentage(predicted.knn.ps, testNa$Pseudomonas.spp.)
# K-nearest neighbours plot predicted values vs Actual/Observed values
plot(predicted.knn.ps, testNa$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("k-nearest neighbours Model for Pseudomonas \nRMSE:",round(RMSE.knn.ps,digits = 2),"Accuracy :",round(percentage_SVMLM_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
#########################SVM Radial for PS count by MSI_data ######################
###################################################################################


regressor.ps = train(Pseudomonas.spp. ~ .,
                   data = trainNb,
                   method = 'svmRadial',
                   trcontrol = control)##Training model
y_pred.ps = predict(regressor.ps, testNa[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMR.ps <- RMSE(testNa$Pseudomonas.spp., y_pred.ps)##RMSE calculation
## k-nearest neighbours model Accuracy for Pseudomonas
percentage_SVMR_Pseudomonas<- Percentage(y_pred.ps, testNa$Pseudomonas.spp.)
#SVM Radial plot predicted values vs Actual/Observed values
plot(y_pred.ps, testNa$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for Pseudomonas \nRMSE:",round(RMSE.SVMR.ps,digits = 2),"Accuracy :",round(percentage_SVMR_Pseudomonas, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

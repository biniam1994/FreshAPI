##############################
#                            #
# Group_project MSI Data     #
#   Machine Learning         #
#    Mixed batches           #
##############################
# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()


library(caret)
library(kernlab)
library(openxlsx)
source("CreateObj.R")
source("Percentage_Accuracy.r")
# Read datasets
MSI   <- read.xlsx("MSIdata_CT_R1_2_updated.xlsx", sheet=1, rowName=T)
bacteria<- MSI[,37:38]##Extracting last 2 bacterial column to variable 

MSI <- MSI[,1:18] ##Extracting all the 18 Mean spectra in variable
MSI_raw<- CreateObj(MSI, bacteria) ##Merging spectra with Bactria column by rownames
MSI_TVC <- MSI_raw[,-(20)]###remove Pseudomonas
MSI_Ps <- MSI_raw[,-(19)]###remove TVC

# Split training and test set for TVC count
set.seed(123)##set the seed to have reproducable results

train.index <- createDataPartition(MSI_TVC$TVC, p = .7,list = FALSE, groups = 3, times = 1)
trainN <- MSI_TVC[train.index,]
testN <- MSI_TVC[-train.index,]

# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
##################################################################################
################### Linear Model for TVC count by MSI_data  #####################

lm.model.fit.TVC <- train(TVC ~., method= 'lm', data=trainN, trControl= control)##Training model
predicted.lm <- predict(lm.model.fit.TVC, testN[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.lm.TVC <- RMSE(testN$TVC, predicted.lm)##Calculating RMSE
##Linear model Accuracy for TVC
percentage_lm_TVC<-Percentage(predicted.lm,testN$TVC)
# LM plot predicted values vs Actual/Observed values
plot(predicted.lm, testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)", ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.lm.TVC,digits = 2)," \nAccuracy :",round(percentage_lm_TVC, digits = 2),"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(lm.model.fit.TVC, "Model.MSI_LM_TVC.rds")##save model as rds object
saveRDS(percentage_lm_TVC, "Accuracy.MSI_LM_TVC.rds")##save model Accuracy
saveRDS(predicted.lm, "Prediction.MSI_LM_TVC.rds")##save model Prediction
saveRDS(RMSE.lm.TVC, "RMSE.MSI_LM_TVC.rds")##save model RMSE
saveRDS(testN$TVC, "Test_MSI_TVC.rds")##save Test in rds

##################################################################################
################ k-nearest neighbours for TVC count by  MSI_data ################

knn.model.fit.TVC <- train(TVC ~ ., method='knn', data=trainN, trcontrol = control, tuneGrid=expand.grid(k=1:20))##Training model
predicted.knn <- predict(knn.model.fit.TVC,testN[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.knn <- RMSE(testN$TVC,predicted.knn)##Calculating RMSE
##KNN model Accuracy for TVC
percentage_knn_TVC<- Percentage(predicted.knn, testN$TVC)

# kNN plot predicted values vs Actual/Observed values
plot(predicted.knn,testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)",col = "blue", 
     main=paste("k-nearest neighbours Model for TVC - RMSE:",round(RMSE.knn,digits = 2)," Accuracy :",round(percentage_knn_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(knn.model.fit.TVC, "Model.MSI_KNN_TVC.rds")##save model as rds object
saveRDS(percentage_knn_TVC, "Accuracy.MSI_knn_TVC.rds")##save model Accuracy
saveRDS(predicted.knn, "Prediction.MSI_knn_TVC.rds")##save model Prediction
saveRDS(RMSE.knn, "RMSE.MSI_knn_TVC.rds")##save model RMSE
###################################################################################
################### Random Forests for TVC count by MSI_data  #####################

RF.model.fit.TVC <- train(TVC ~ ., method='rf', trcontrol = control, trainN,tuneGrid=expand.grid(mtry=1:20))##Training model
predicted.rf <- predict(RF.model.fit.TVC, testN[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.rf <- RMSE(testN$TVC,predicted.rf)##calculating RMSE

## RandomForest model Accuracy for TVC
percentage_RF_TVC<- Percentage(predicted.rf,testN$TVC)
# RandomForest plot predicted values vs Actual/Observed values
plot(predicted.rf, testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Random Forests Model for TVC - RMSE:",round(RMSE.rf,digits = 2)," Accuracy :",round(percentage_RF_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(RF.model.fit.TVC, "Model.MSI_RF_TVC.rds")##save model as rds object
saveRDS(percentage_RF_TVC, "Accuracy.MSI_RF_TVC.rds")##save model Accuracy
saveRDS(predicted.rf, "Prediction.MSI_RF_TVC.rds")##save model Prediction
saveRDS(RMSE.rf, "RMSE.MSI_RF_TVC.rds")##save model RMSE
###################################################################################
################### SVMPloy for TVC count by MSI_data         #####################

SVM.model.fit.TVC <- train(TVC ~ ., method='svmPoly', data= trainN,tuneLength = 5,trcontrol = control)##model training
predicted.SVMP <- predict(SVM.model.fit.TVC,testN[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMP <- RMSE(testN$TVC,predicted.SVMP)##calculating RMSE

## SVM_Poly model Accuracy for TVC
percentage_SVMP_TVC<- Percentage(predicted.SVMP,testN$TVC)
#SVMPloy plot predicted values vs Actual/Observed values
plot(predicted.SVMP,testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynomial SVM Model for TVC - RMSE:",round(RMSE.SVMP,digits = 2)," Accuracy :",round(percentage_SVMP_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(SVM.model.fit.TVC, "Model.MSI_SVMPoly_TVC.rds")##save model as rds object
saveRDS(percentage_SVMP_TVC, "Accuracy.MSI_SVMPoly_TVC.rds")##save model Accuracy
saveRDS(predicted.SVMP, "Prediction.MSI_SVMPoly_TVC.rds")##save model Prediction
saveRDS(RMSE.SVMP, "RMSE.MSI_SVMPoly_TVC.rds")##save model RMSE

###################################################################################
#########################SVM_Radial for TVC count by MSI_data######################
###################################################################################
regressor = train(TVC ~ .,
                data = trainN,
                method= 'svmRadial',
                trcontrol = control)##Training model
y_pred = predict(regressor, testN[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMR.TVC <- RMSE(testN$TVC, y_pred)##calculating RMSE

## SVM_radial model Accuracy for TVC
percentage_SVMRAD_TVC<- Percentage(y_pred, testN$TVC)
#SVMRadial plot predicted values vs Actual/Observed values
plot(y_pred,testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for TVC - RMSE:",round(RMSE.SVMR.TVC,digits = 2)," Accuracy :",round(percentage_SVMRAD_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(regressor, "Model.MSI_SVMRadial_TVC.rds")##save model as rds object
saveRDS(percentage_SVMRAD_TVC, "Accuracy.MSI_SVMRadial_TVC.rds")##save model Accuracy
saveRDS(y_pred, "Prediction.MSI_SVMRadial_TVC.rds")##save model Prediction
saveRDS(RMSE.SVMR.TVC, "RMSE.MSI_SVMRadial_TVC.rds")##save model RMSE

#####################################################################################
###################################Svm_Linear for TVC count by MSI_data##############

SVMLM.model.fit.TVC <- train(TVC ~ ., method='svmLinear', data= trainN, tuneLength = 4,trcontrol = control)##Training model
predicted.SVMLM.TVC <- predict(SVMLM.model.fit.TVC,testN[,-19])## Predicting bacterial count TVC using the trained model.
RMSE.SVMLM.TVC <- RMSE(testN$TVC, predicted.SVMLM.TVC)##calculating RMSE

## SVM_linear model Accuracy for TVC
percentage_SVMLM_TVC<-Percentage(predicted.SVMLM.TVC, testN$TVC)
#SVMLinear plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.TVC,xlim= c(0,9), ylim=c(0,9),testN$TVC,xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear SVM Model for TVC - RMSE:",round(RMSE.SVMLM.TVC,digits = 2), "Accuracy:",round(percentage_SVMLM_TVC,digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(SVMLM.model.fit.TVC, "Model.MSI_SVMlinear_TVC.rds")##save model as rds object
saveRDS(percentage_SVMLM_TVC, "Accuracy.MSI_SVMlinear_TVC.rds")##save model Accuracy
saveRDS(predicted.SVMLM.TVC, "Prediction.MSI_SVMlinear_TVC.rds")##save model Prediction
saveRDS(RMSE.SVMLM.TVC, "RMSE.MSI_SVMlinear_TVC.rds")##save model RMSE
###########################################Pseudomonas MSI_data#####################################################################################
#####################################################################################################################################################
set.seed(123)
train.index <- createDataPartition(MSI_Ps$Pseudomonas.spp., p = .7,list = FALSE,groups = 3, times = 1)
trainNP <- MSI_Ps[train.index,]
testNP <- MSI_Ps[-train.index,]

####################################################################################
################### Linear Model for Pseudomonas count by MSI_data  ################

lm.model.fit.ps <- train(Pseudomonas.spp. ~ ., method= 'lm', data=trainNP, trcontrol= control)##Training model 
predicted.lm.ps <- predict(lm.model.fit.ps,testNP)## Predicting bacterial count Pseudomonas using the trained model.
RMSE.lm.ps <- RMSE(testNP$Pseudomonas.spp., predicted.lm.ps)##RMSE calculation
##Linear model Accuracy for Pseudomonas
percentage_lm_PS<-Percentage(predicted.lm.ps, testNP$Pseudomonas.spp.)
#plot predicted values vs Actual/Observed values
plot(predicted.lm.ps, testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Pseudomonas.spp counts \nRMSE:", round(RMSE.lm.ps,digits = 2)," \nAccuracy:",round(percentage_lm_PS, digits = 2),"% - at ±1 LogCount"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(lm.model.fit.ps, "Model.MSI_LM_Pseudomonas.rds")##save model as rds object
saveRDS(percentage_lm_PS, "Accuracy.MSI_LM_Pseudomonas.rds")##save model Accuracy
saveRDS(predicted.lm.ps, "Prediction.MSI_LM_Pseudomonas.rds")##save model Prediction
saveRDS(RMSE.lm.ps, "RMSE.MSI_LM_Pseudomonas.rds")##save model RMSE
saveRDS(testNP$Pseudomonas.spp., "Test_MSI_Pseudomonas.rds")##save Test in rds
####################################################################################
################### Random Forests for Pseudomonas count by MSI_data  ##############

RF.model.fit.Ps <- train(Pseudomonas.spp. ~ ., method='rf', trcontrol = control, trainNP)##Training model
predicted.rf.ps <- predict(RF.model.fit.Ps, testNP[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.rf.ps <- RMSE(testNP$Pseudomonas.spp.,predicted.rf.ps)##RMSE calculation

##RandomForest model Accuracy for Pseudomonas
percentage_RF_PS<- Percentage(predicted.rf.ps, testNP$Pseudomonas.spp.)
# RF plot predicted values vs Actual/Observed values
plot(predicted.rf.ps,testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("RandomForest Model for Pseudomonas.spp - RMSE:",round(RMSE.rf.ps,digits = 2)," Accuracy :",round(percentage_RF_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(RF.model.fit.Ps, "Model.MSI_RF_Pseudomonas.rds")##save model as rds object
saveRDS(percentage_RF_PS, "Accuracy.MSI_RF_Pseudomonas.rds")##save model Accuracy
saveRDS(predicted.rf.ps, "Prediction.MSI_RF_Pseudomonas.rds")##save model Prediction
saveRDS(RMSE.rf.ps, "RMSE.MSI_RF_Pseudomonas.rds")##save model RMSE
###################################################################################
################### SVM_Ploy for Pseudomonas count by MSI_data ####################


SVM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmPoly', data= trainNP, trcontrol = control)##Training model
predicted.SVMP.PS <- predict(SVM.model.fit.PS,testNP[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMP.ps <- RMSE(testNP$Pseudomonas.spp., predicted.SVMP.PS)##RMSE calculation

##SVM_Poly model Accuracy for Pseudomonas
percentage_SVMP_PS<- Percentage(predicted.SVMP.PS, testNP$Pseudomonas.spp.)
# SVMPoly plot predicted values vs Actual/Observed values
plot(predicted.SVMP.PS,testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynomial SVM Model for Pseudomonas.spp - RMSE:",round(RMSE.SVMP.ps,digits = 2)," Accuracy :",round(percentage_SVMP_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(SVM.model.fit.PS, "Model.MSI_SVMPoly_Pseudomonas.rds")##save model as rds object
saveRDS(percentage_SVMP_PS, "Accuracy.MSI_SVMPoly_Pseudomonas.rds")##save model Accuracy
saveRDS(predicted.SVMP.PS, "Prediction.MSI_SVMPoly_Pseudomonas.rds")##save model Prediction
saveRDS(RMSE.SVMP.ps, "RMSE.MSI_SVMPoly_Pseudomonas.rds")##save model RMSE
###################################################################################
################### SVM_Linear for Pseudomonas count by MSI_data ##################


SVMLM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmLinear', data= trainNP, trcontrol = control)##Training model
predicted.SVMLM.PS <- predict(SVMLM.model.fit.PS,testNP[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.SVMLM.ps <- RMSE(testNP$Pseudomonas.spp., predicted.SVMLM.PS)##RMSE calculation

##SVM_Linear model Accuracy for Pseudomonas

percentage_SVMLM_PS<- Percentage(predicted.SVMLM.PS,testNP$Pseudomonas.spp.)

# SVMLinear plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.PS,testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)",
     main=paste("Linear SVM Model for Pseudomonas.spp - RMSE:",round(RMSE.SVMLM.ps,digits = 2)," Accuracy :",round(percentage_SVMLM_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(SVMLM.model.fit.PS, "Model.MSI_SVMlinear_Pseudomonas.rds")##save model as rds object
saveRDS(percentage_SVMLM_PS, "Accuracy.MSI_SVMlinear_Pseudomonas.rds")##save model Accuracy
saveRDS(predicted.SVMLM.PS, "Prediction.MSI_SVMlinear_Pseudomonas.rds")##save model Prediction
saveRDS(RMSE.SVMLM.ps, "RMSE.MSI_SVMlinear_Pseudomonas.rds")##save model RMSE
##################################################################################
################ k-nearest neighbours for Pseudomonas count by  MSI_data #########


knn.model.fit.ps <- train(Pseudomonas.spp. ~ ., method='knn', data=trainNP, trcontrol = control, tuneGrid=expand.grid(k=1:40))##Training model
predicted.knn.ps <- predict(knn.model.fit.ps, testNP[,-19])## Predicting bacterial count Pseudomonas using the trained model.
RMSE.knn.ps <- RMSE(testNP$Pseudomonas.spp.,predicted.knn.ps)##RMSE calculation

##KNN model Accuracy for Pseudomonas
percentage_knn_PS<- Percentage(predicted.knn.ps,testNP$Pseudomonas.spp.)
# KNN plot predicted values vs Actual/Observed values
plot(predicted.knn.ps, testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("K-nearest neighbours Model for Pseudomonas.spp - RMSE:",round(RMSE.knn.ps,digits = 2)," Accuracy :",round(percentage_knn_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(knn.model.fit.ps, "Model.MSI_knn_Pseudomonas.rds")##save model as rds object
saveRDS(percentage_knn_PS, "Accuracy.MSI_knn_Pseudomonas.rds")##save model Accuracy
saveRDS(predicted.knn.ps, "Prediction.MSI_knn_Pseudomonas.rds")##save model Prediction
saveRDS(RMSE.knn.ps, "RMSE.MSI_knn_Pseudomonas.rds")##save model RMSE
#####################################################################################
#########################SVM_Radial for Pseudomonas count by MSI_data################


regressor.ps = train(Pseudomonas.spp. ~ .,
                   data = trainNP,
                   method = 'svmRadial',
                   trcontrol = control)# Predicting bacterial count Pseudomonas using the trained model.
y_pred.ps = predict(regressor.ps, testNP[,-19])##Training model
RMSE.SVMR.ps <- RMSE(testNP$Pseudomonas.spp., y_pred.ps)##RMSE calculation

##SVM_radial model Accuracy for Pseudomonas
percentage_SVMRAD_PS<- Percentage(y_pred.ps,testNP$Pseudomonas.spp.)
# SVMRadial plot predicted values vs Actual/Observed values
plot(y_pred.ps, testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for Pseudomonas.spp - RMSE:",round(RMSE.SVMR.ps,digits = 2)," Accuracy :",round(percentage_SVMRAD_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

saveRDS(regressor.ps, "Model.MSI_SVMRadial_Pseudomonas.rds")##save model as rds object
saveRDS(percentage_SVMRAD_PS, "Accuracy.MSI_SVMRadial_Pseudomonas.rds")##save model Accuracy
saveRDS(y_pred.ps, "Prediction.MSI_SVMRadial_Pseudomonas.rds")##save model Prediction
saveRDS(RMSE.SVMR.ps, "RMSE.MSI_SVMRadial_Pseudomonas.rds")##save model RMSE
#################################################################################################################################################

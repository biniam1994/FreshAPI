##############################
#                            #
#   Group_project MSI Data  #
#   Machine Learning         #
##############################
# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()



library(caret)
library(kernlab)
library(gmodels)
library(openxlsx)
source("CreateObj.R")
source("Percentage_Accuracy.r")
# Read datasets
MSI   <- read.xlsx("MSIdata_CT_R1_2_updated.xlsx", sheet=1, rowName=T)
#create two datasets based on batches 
bacteria_B<- MSI[100:196,37:38] ## Bacterial count for batch A
MSI_B <- MSI[101:196,1:18] ### Bath A


MSI_raw_B<- CreateObj(MSI_B, bacteria_B)##Merging bacterial counts with spectra of batch A

MSI_TVC_B <- MSI_raw_B[,-(20)]###remove pseudomnas A
MSI_PS_B <- MSI_raw_B[,-19]###remove TVC A

###############################################################################################################################################################################################
##################################################################TVC_Regression_RAW###########################################################################################################
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# Split training and test set for TVC count
set.seed(125)
train.index <- createDataPartition(MSI_TVC_B$TVC, p = .7,list = FALSE, groups=3, times = 50)
##################################################################################
################### Linear Model for TVC count by MSI_data   #####################
##################################################################################
RMSELMT<-c()
for (i in 1:50) {
  trainN <- MSI_TVC_B[train.index[,i],]
  testN <- MSI_TVC_B[-train.index[,i],]
  lm.model.fit.TVC <- train(TVC ~., method= 'lm', data=trainN, trControl= control)##train
  predicted.lm <- predict(lm.model.fit.TVC, testN[,-19])##Prediction
  RMSE.lm.TVC <- RMSE(testN$TVC, predicted.lm)##RMSE calculation
  RMSELMT<-c(RMSELMT, RMSE.lm.TVC)
}
lm_TVCM<-mean(RMSELMT)##calculating mean
iteration <- as.array(c(1:50))
lm_TVCS<-sd(RMSELMT)##calculating SD
LM_TVC_95<-ci(RMSELMT, confidence = 0.95)
plot(iteration, RMSELMT,ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("RMSE for 50 iterations \n( RMSE.mean = ",
                                                                                                      round(lm_TVCM, digits = 3)," +/- ",
                                                                                                      round(lm_TVCS,digits = 3), " )\nThere is a 95% likelihood that the range",
                                                                                                      round(LM_TVC_95[2], digits = 3), "to", 
                                                                                                      round(LM_TVC_95[3], digits = 3), "covers the true error of the model."))##Plotting RMSE

##Linear model Accuracy for TVC
percentage_lm_TVC<-Percentage(predicted.lm,testN$TVC)
# LM plot predicted values vs Actual/Observed values
plot(predicted.lm, testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Total Viable Counts (TVC) \nRMSE:",round(RMSE.lm.TVC,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
##################################################################################
################ k-nearest neighbours for TVC count by  MSI_data  ################
##################################################################################
RMSEKNNT<-c()
for (i in 1:50) {
  trainN <- MSI_TVC_B[train.index[,i],]
  testN <- MSI_TVC_B[-train.index[,i],]
  knn.model.fit.TVC <- train(TVC ~ ., method='knn', data=trainN, trcontrol = control, tuneGrid=expand.grid(k=1:20))
  predicted.knn <- predict(knn.model.fit.TVC,testN[,-19])
  RMSE.knn <- RMSE(testN$TVC, predicted.knn)
  RMSEKNNT<-c(RMSE.knn, RMSEKNNT)
}
KNN_TVCM<-mean(RMSEKNNT)##calculating mean
KNN_TVCS<-sd(RMSEKNNT)##calculating SD
Knn_TVC_95<-ci(RMSEKNNT, confidence = 0.95)
plot(iteration, RMSEKNNT,ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("Knn Model(MSI_TVC) Mean RMSE:",
                                                                                                       round(KNN_TVCM, digits = 3), "CI 95%:", 
                                                                                                       round(Knn_TVC_95[2], digits = 3), "-", 
                                                                                                       round(Knn_TVC_95[3], digits = 3), "SD:",
                                                                                                       round(KNN_TVCS,digits = 3)))##Plotting RMSE

##KNN model Accuracy for TVC
percentage_knn_TVC<- Percentage(predicted.knn, testN$TVC)
# kNN plot predicted values vs Actual/Observed values
plot(predicted.knn,testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)",col = "blue", 
     main=paste("k-nearest neighbours Model for TVC - RMSE:",round(RMSE.knn,digits = 2)," Accuracy :",round(percentage_knn_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### Random Forests for TVC count by MSI_data   ####################
###################################################################################
RMSERFT<-c()
for (i in 1:50) {
  trainN <- MSI_TVC_B[train.index[,i],]
  testN <- MSI_TVC_B[-train.index[,i],]
  RF.model.fit.TVC <- train(TVC ~ ., method='rf', trcontrol = control, trainN,tuneGrid=expand.grid(mtry=1:30))
  predicted.rf <- predict(RF.model.fit.TVC, testN[,-19])
  RMSE.rf <- RMSE(testN$TVC,predicted.rf)
  RMSERFT<-c(RMSERFT, RMSE.rf)
}
RF_TVCM<-mean(RMSERFT)##calculating mean
RF_TVCS<-sd(RMSERFT)##calculating SD
RF_TVC_95<-ci(RMSERFT, confidence = 0.95)
plot(iteration, RMSERFT, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("RandomForest Model(MSI_TVC) Mean RMSE:",
                                                                                                       round(RF_TVCM, digits = 3), "CI 95%:", 
                                                                                                       round(RF_TVC_95[2], digits = 3), "-", 
                                                                                                       round(RF_TVC_95[3], digits = 3), "SD:",
                                                                                                       round(RF_TVCS,digits = 3)))##Plotting RMSE
## RandomForest model Accuracy for TVC
percentage_RF_TVC<- Percentage(predicted.rf,testN$TVC)
# RF plot predicted values vs Actual/Observed values
plot(predicted.rf, testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Random Forests Model for TVC - RMSE:",round(RMSE.rf,digits = 2)," Accuracy :",round(percentage_RF_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVM_Ploy for TVC count by MSI_data          ###################
###################################################################################
RMSESVMPLOYT<-c()
for (i in 1:50) {
  trainN <- MSI_TVC_B[train.index[,i],]
  testN <- MSI_TVC_B[-train.index[,i],]
  SVM.model.fit.TVC <- train(TVC ~ ., method='svmPoly', data= trainN, trcontrol=control, tuneLength = 5)
  predicted.SVM <- predict(SVM.model.fit.TVC,testN[,-19])
  RMSE.SVMP <- RMSE(testN$TVC,predicted.SVM)##calculating RMSE
  RMSESVMPLOYT<-c(RMSESVMPLOYT, RMSE.SVMP)
}
SPoly_TVCM<-mean(RMSESVMPLOYT)##calculating mean
SPoly_TVCS<-sd(RMSESVMPLOYT)##calculating SD
SP_TVC_95<-ci(RMSESVMPLOYT, confidence = 0.95)
plot(iteration, RMSESVMPLOYT, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("SVM_Poly Model(MSI_TVC) Mean RMSE:",
                                                                                                            round(SPoly_TVCM, digits = 3), "CI 95%:", 
                                                                                                            round(SP_TVC_95[2], digits = 3), "-", 
                                                                                                            round(SP_TVC_95[3], digits = 3), "SD:",
                                                                                                            round(SPoly_TVCS,digits = 3)))##Plotting RMSE

## SVM_Poly model Accuracy for TVC
percentage_SVMP_TVC<- Percentage(predicted.SVMP,testN$TVC)
#SVMPloy plot predicted values vs Actual/Observed values
plot(predicted.SVMP,testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynoimal SVM Model for TVC - RMSE:",round(RMSE.SVMP,digits = 2)," Accuracy :", round(percentage_SVMP_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
#########################SVM_Radial for TVC count by MSI_data #####################
###################################################################################
RMSESVMRAD<-c()
for (i in 1:50) {
  trainN <- MSI_TVC_B[train.index[,i],]
  testN <- MSI_TVC_B[-train.index[,i],]
  regressor = train(TVC ~ .,
                  data = trainN,
                  method = 'svmRadial',
                  trcontrol = control)
  y_pred = predict(regressor, testN[,-19])
  RMSE.SVMR.TVC <- RMSE(testN$TVC, y_pred)##calculating RMSE
  RMSESVMRAD<-c(RMSESVMRAD,RMSE.SVMR.TVC)
}
SRAD_TVCM<-mean(RMSESVMRAD)##calculating mean
SRAD_TVCS<-sd(RMSESVMRAD)##calculating SD
SR_TVC_95<-ci(RMSESVMRAD, confidence = 0.95)
plot(iteration, RMSESVMRAD, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("SVM_Radial Model(MSI_TVC) Mean RMSE:",
                                                                                                          round(SRAD_TVCM, digits = 3), "CI 95%:", 
                                                                                                          round(SR_TVC_95[2], digits = 3), "-", 
                                                                                                          round(SR_TVC_95[3], digits = 3), "SD:",
                                                                                                          round(SRAD_TVCS,digits = 3)))##Plotting RMSE
## SVM_radial model Accuracy for TVC
percentage_SVMRAD_TVC<- Percentage(y_pred, testN$TVC)
#SVM_Radial plot predicted values vs Actual/Observed values
plot(y_pred,testN$TVC,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for TVC - RMSE:",round(RMSE.SVMR.TVC,digits = 2)," Accuracy :", round(percentage_SVMRAD_TVC, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###########################################################################################
###################################Svm_Linear for TVC count by MSI_data####################
###########################################################################################
RMSESVMLM<-c()
for (i in 1:50){
  trainN <- MSI_TVC_B[train.index[,i],]
  testN <- MSI_TVC_B[-train.index[,i],]
  SVMLM.model.fit.TVC <- train(TVC ~ ., method='svmLinear', data= trainN, trcontrol=control)
  predicted.SVMLM.TVC <- predict(SVMLM.model.fit.TVC,testN[,-19])
  ##plot predicted vs Observed with cl-95%
  RMSE.SVMLM.TVC <- RMSE(testN$TVC, predicted.SVMLM.TVC)##RMSE calculation
  RMSESVMLM<- c(RMSESVMLM,RMSE.SVMLM.TVC)
}
SLM_TVCM<-mean(RMSESVMLM)
iteration <- as.array(c(1:50))
SLM_TVCS<-sd(RMSESVMLM)
SLM_TVC_95<-ci(RMSESVMLM, confidence = 0.95)
plot(iteration, RMSESVMLM, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("SVM_Linear Model(MSI_TVC) Mean RMSE:",
                                                                                                         round(SLM_TVCM, digits = 3), "CI 95%:", 
                                                                                                         round(SLM_TVC_95[2], digits = 3), "-", 
                                                                                                         round(SLM_TVC_95[3], digits = 3), "SD:",
                                                                                                         round(SLM_TVCS,digits = 3)))##Plotting RMSE

## SVM_linear model Accuracy for TVC
percentage_SVMLM_TVC<-Percentage(predicted.SVMLM.TVC, testN$TVC)
#plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.TVC,xlim= c(0,9), ylim=c(0,9),testN$TVC,xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear support-vector machine(SVM) Model for Total Viable Counts(TVC) \nRMSE:",round(RMSE.SVMLM.TVC,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
#####################################################################################################################################################
###########################################Pseudomonas MSI_data #####################################################################################
#####################################################################################################################################################

set.seed(125)
train.indexs <- createDataPartition(MSI_PS_B$Pseudomonas.spp., p = .7, list = FALSE,groups = 3, times = 50)

##################################################################################
################### Linear Model for TVC count by MSI_data   #####################
##################################################################################

RMSELMP<-c()
for (i in 1:50){
  trainNP <- MSI_PS_B[train.indexs[,i],]
  testNP <- MSI_PS_B[-train.indexs[,i],]
  lm.model.fit.ps <- train(Pseudomonas.spp. ~ ., method= 'lm', data=trainNP, trcontrol= control) 
  predicted.lm.ps <- predict(lm.model.fit.ps,testNP[,-19])
  RMSE.lm.ps <- RMSE(testNP$Pseudomonas.spp., predicted.lm.ps)##RMSE calculation
  RMSELMP<- c(RMSELMP, RMSE.lm.ps)
}
LM_PM<-mean(RMSELMP)
iteration <- as.array(c(1:50))
LM_PS<-sd(RMSELMP)
LM_PS_95<-ci(RMSELMP, confidence = 0.95)

plot(iteration, RMSELMP, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("RMSE for 50 iterations \n( RMSE.mean = ",
                                                                                                       round(LM_PM, digits = 3), " +/- ",
                                                                                                       round(LM_PS,digits = 3)," )\nThere is a 95% likelihood that the range",
                                                                                                       round(LM_PS_95[2], digits = 3), "to", 
                                                                                                       round(LM_PS_95[3], digits = 3), "covers the true error of the model."))##Plotting RMSE

##Linear model Accuracy for Pseudomonas
percentage_lm_PS<-Percentage(predicted.lm.ps, testNP$Pseudomonas.spp.)
#lm plot predicted values vs Actual/Observed values
plot(predicted.lm.ps, testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear Model for Pseudomonas.spp \nRMSE:", round(RMSE.lm.ps,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### Random Forests for TVC count by MSI_data   ####################
###################################################################################
RMSERFP<-c()

for (i in 1:50){
  trainNP <- MSI_PS_B[train.indexs[,i],]
  testNP <- MSI_PS_B[-train.indexs[,i],]
  RF.model.fit.Ps <- train(Pseudomonas.spp. ~ ., method='rf', trcontrol = control, trainNP)
  predicted.rf.ps <- predict(RF.model.fit.Ps, testNP[,-19])
  RMSE.rf.ps <- RMSE(testNP$Pseudomonas.spp.,predicted.rf.ps)
  RMSERFP<-c(RMSERFP,RMSE.rf.ps)
}
RF_PM<-mean(RMSERFP)
iteration <- as.array(c(1:50))
RF_PS<-sd(RMSERFP)
RF_PS_95<-ci(RMSERFP, confidence = 0.95)
plot(iteration, RMSERFP, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("RandomForest Model(MSI_PS) Mean RMSE:",
                                                                                                       round(RF_PM, digits = 3), "CI 95%:", 
                                                                                                       round(RF_PS_95[2], digits = 3), "-", 
                                                                                                       round(RF_PS_95[3], digits = 3), "SD:",
                                                                                                       round(RF_PS,digits = 3)))##Plotting RMSE

##RandomForest model Accuracy for Pseudomonas
percentage_RF_PS<- Percentage(predicted.rf.ps, testNP$Pseudomonas.spp.)
# RF plot predicted values vs Actual/Observed values
plot(predicted.rf.ps,testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Random Forest Model for Pseudomonas.spp. - RMSE:",round(RMSE.rf.ps,digits = 2)," Accuracy :", round(percentage_RF_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVMPloy for Pseudomonas count by MSI_data          ############
###################################################################################
RMSESVPP<-c()
#train.index <- createDataPartition(MSI_Ps$Pseudomonas.spp., p = .7,list = FALSE, times = 50)
for (i in 1:50){
  trainNP <- MSI_PS_B[train.indexs[,i],]
  testNP <- MSI_PS_B[-train.indexs[,i],]
  SVM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmPoly', data= trainNP, trcontrol=control, tuneLength = 5)
  predicted.SVM.PS <- predict(SVM.model.fit.PS,testNP[,-19])
  RMSE.SVMP.ps <- RMSE(testNP$Pseudomonas.spp., predicted.SVM.PS)##RMSE calculation
  RMSESVPP<- c(RMSESVPP, RMSE.SVMP.ps)
}
SP_PM<-mean(RMSESVPP)
iteration <- as.array(c(1:50))
SP_PS<-sd(RMSESVPP) 
SP_PS_95<-ci(RMSESVPP, confidence = 0.95)
plot(iteration, RMSESVPP, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("SVM_Poly Model(MSI_PS) Mean RMSE:",
                                                                                                        round(SP_PM, digits = 3), "CI 95%:", 
                                                                                                        round(SP_PS_95[2], digits = 3), "-", 
                                                                                                        round(SP_PS_95[3], digits = 3), "SD:",
                                                                                                        round(SP_PS,digits = 3)))##Plotting RMSE

##SVM_Poly model Accuracy for Pseudomonas
percentage_SVMP_PS<- Percentage(predicted.SVMP.PS, testNP$Pseudomonas.spp.)
# SVMPoly plot predicted values vs Actual/Observed values
plot(predicted.SVMP.PS,testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Polynomial support-vector machine Model for Pseudomonas.spp \nRMSE:",round(RMSE.SVMP.ps,digits = 2)))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
################### SVM Linear for Pseudomonas count by MS_data        ############
###################################################################################
RMSESLMP<-c()

for (i in 1:50){
  trainNP <- MSI_PS_B[train.indexs[,i],]
  testNP <- MSI_PS_B[-train.indexs[,i],]
  SVMLM.model.fit.PS <- train(Pseudomonas.spp. ~ ., method='svmLinear', data= trainNP)
  predicted.SVMLM.PS <- predict(SVMLM.model.fit.PS,testNP[,-19])
  RMSE.SVMLM.ps <- RMSE(testNP$Pseudomonas.spp., predicted.SVMLM.PS)##RMSE calculation
  RMSESLMP<-c(RMSESLMP, RMSE.SVMLM.ps)
}
LM_PSM<-mean(RMSESLMP)
iteration <- as.array(c(1:50))
LM_PSS<-sd(RMSESLMP)
SLM_PS_95<-ci(RMSESLMP, confidence = 0.95)

plot(iteration, RMSESLMP, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("SVM_Linear Model(MSI_PS) Mean RMSE:",
                                                                                                        round(LM_PSM, digits = 3), "CI 95%:", 
                                                                                                        round(SLM_PS_95[2], digits = 3), "-", 
                                                                                                        round(SLM_PS_95[3], digits = 3), "SD:",
                                                                                                        round(LM_PSS,digits = 3)))##Plotting RMSE
##SVM_Linear model Accuracy for Pseudomonas
percentage_SVMLM_PS<- Percentage(predicted.SVMLM.PS,testNP$Pseudomonas.spp.)
##SVMLinear plot predicted values vs Actual/Observed values
plot(predicted.SVMLM.PS,testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Linear SVM Model for Pseudomonas.spp. - RMSE:",round(RMSE.SVMLM.ps,digits = 2)," Accuracy :", round(percentage_SVMLM_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
##################################################################################
################ k-nearest neighbours for PS count by  FTIR_data #################
##################################################################################
RMSESKNNP<-c()

for (i in 1:50){
  trainNP <- MSI_PS_B[train.indexs[,i],]
  testNP <- MSI_PS_B[-train.indexs[,i],]
  knn.model.fit.ps <- train(Pseudomonas.spp. ~ ., method='knn', data=trainNP, trcontrol = control, tuneGrid=expand.grid(k=1:40))
  predicted.knn.ps <- predict(knn.model.fit.ps, testNP[,-19])
  RMSE.knn.ps <- RMSE(testNP$Pseudomonas.spp.,predicted.knn.ps)##RMSE calculation
  RMSESKNNP<- c(RMSESKNNP, RMSE.knn.ps)
}
KNN_PSM<-mean(RMSESKNNP)
iteration <- as.array(c(1:50))
KNN_PSS<-sd(RMSESKNNP)
Knn_PS_95<-ci(RMSESKNNP, confidence = 0.95)
plot(iteration, RMSESKNNP, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("Knn Model(MSI_PS) Mean RMSE:",
                                                                                                         round(KNN_PSM, digits = 3), "CI 95%:", 
                                                                                                         round(Knn_PS_95[2], digits = 3), "-", 
                                                                                                         round(Knn_PS_95[3], digits = 3), "SD:",
                                                                                                         round(KNN_PSS,digits = 3)))##Plotting RMSE

##KNN model Accuracy for Pseudomonas
percentage_knn_PS<- Percentage(predicted.knn.ps,testNP$Pseudomonas.spp.)
# KNN plot predicted values vs Actual/Observed values
plot(predicted.knn.ps, testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("k-nearest neighbours Model for Pseudomonas.spp. - RMSE:",round(RMSE.knn.ps,digits = 2)," Accuracy :",round(percentage_knn_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)
###################################################################################
#########################SVMRadial for PS count by FTIR_data#######################
###################################################################################
RMSESRADP<-c()

for (i in 1:50){
  trainNP <- MSI_PS_B[train.indexs[,i],]
  testNP <- MSI_PS_B[-train.indexs[,i],]
  regressor.ps = train(Pseudomonas.spp. ~ .,
                     data = trainNP,
                     method = 'svmRadial',
                     trcontrol = control)
  y_pred.ps = predict(regressor.ps, testNP[,-19])
  RMSE.SVMR.ps <- RMSE(testNP$Pseudomonas.spp., y_pred.ps)
  RMSESRADP<- c(RMSESRADP, RMSE.SVMR.ps)
}
RAD_PSM<-mean(RMSESRADP)
iteration <- as.array(c(1:50))
RAD_PSS<-sd(RMSESRADP)
SR_PS_95<-ci(RMSESRADP, confidence = 0.95)

plot(iteration, RMSESRADP, ylim=c(0,2), xlim=c(0,50), xlab="Iteration", ylab="RMSE",type='l', main=paste("SVM_Radial Model(MSI_PS) Mean RMSE:",
                                                                                                         round(RAD_PSM, digits = 3), "CI 95%:", 
                                                                                                         round(SR_PS_95[2], digits = 3), "-", 
                                                                                                         round(SR_PS_95[3], digits = 3), "SD:",
                                                                                                         round(RAD_PSS,digits = 3)))##Plotting RMSE


##SVM_radial model Accuracy for Pseudomonas
percentage_SVMRAD_PS<- Percentage(y_pred.ps,testNP$Pseudomonas.spp.)
# SVMRadial plot predicted values vs Actual/Observed values
plot(y_pred.ps, testNP$Pseudomonas.spp.,xlim= c(0,9), ylim=c(0,9),xlab="Estimated microbial Population (log CFU/g)",ylab="Observed microbial Population (log CFU/g)", col="blue",
     main=paste("Radial SVM Model for Pseudomonas.spp. - RMSE:",round(RMSE.SVMR.ps,digits = 2)," Accuracy :",round(percentage_SVMRAD_PS, digits = 2),"%"))
abline(a=0,b=1)
abline(a=-1,b=1)
abline(a=1,b=1)

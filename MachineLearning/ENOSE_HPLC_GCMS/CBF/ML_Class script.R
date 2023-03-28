library(openxlsx)
library(MLmetrics)
require(caret)
require(randomForest)
library(caret)
library(gmodels)
library(ggpubr)
require(pROC)
require(HH)
require(klaR)
library(factoextra)
#install.packages("haven", version='2.3.1')
#library(DMwR)##smote
source("ML_function.R")

##enose NA and replicate treatment
data<- openxlsx::read.xlsx("enose_breastAUA.xlsx", rowNames = T, check.names = F,sep.names=" ")##enose raw data
data<-NA_tre_re(data,3)##NA imputing
write.xlsx(data,"enose_N.xlsx",  col.names = TRUE, row.names = TRUE)##save data 
data<-repli_tre(data,3) ##triplicate averaging
write.xlsx(data,"enose_R.xlsx",  col.names = TRUE, row.names = TRUE)##save data 

##GCMS NA and replicate treatment
data<- openxlsx::read.xlsx("GCMS_breastAUA.xlsx", sheet=2, rowNames = T, check.names = F, sep.names=" ")##GCMS raw data
data<-NA_tre_re(data,3)##NA imputing
write.xlsx(data,"GCMS_N.xlsx",  col.names = TRUE, row.names = TRUE)##save data 
data<-repli_tre(data,3) ##triplicate averaging
write.xlsx(data,"GCMS_R.xlsx",  col.names = TRUE, row.names = TRUE)##save data 

##HPLC NA and replicate treatment
data<- openxlsx::read.xlsx("HPLC_breastAUA.xlsx", sheet=1, rowNames = T, check.names = F, sep.names=" ")##enose raw data
data<-NA_tre_re(data,3)##NA imputing
write.xlsx(data,"HPLC_N.xlsx",  col.names = TRUE, row.names = TRUE)##save data 
data<-repli_tre(data,3)##triplicate averaging
write.xlsx(data,"HPLC_R.xlsx",  col.names = TRUE, row.names = TRUE)##save data
######

##rounding sensory scores or odor after replicate and NA treatment
##enose
data<- openxlsx::read.xlsx("enose_R.xlsx", rowNames = T, check.names = F,sep.names=" ")
data$Odor<-round(data$Odor)
r.names<-paste0(row.names(data), "_R1")##adding batch identification
row.names(data)<-r.names
write.xlsx(data,"enose_RR.xlsx",  col.names = TRUE, row.names = TRUE)
##HPLC
data<- openxlsx::read.xlsx("GCMS_R.xlsx", rowNames = T, check.names = F,sep.names=" ")
data$Odor<-round(data$Odor)
r.names<-paste0(row.names(data), "_R1")##adding batch identification
row.names(data)<-r.names
write.xlsx(data,"GCMS_RR.xlsx",  col.names = TRUE, row.names = TRUE)
##GCMS
data<- openxlsx::read.xlsx("HPLC_R.xlsx", rowNames = T, check.names = F,sep.names=" ")
data$Odor<-round(data$Odor)
r.names<-paste0(row.names(data), "_R1")##adding batch identification
row.names(data)<-r.names
write.xlsx(data,"HPLC_RR.xlsx",  col.names = TRUE, row.names = TRUE)
######
##rounding sensory scores or odor after NA treatment only
##enose
data<- openxlsx::read.xlsx("enose_N.xlsx", rowNames = T, check.names = F,sep.names=" ")
data$Odor<-round(data$Odor)
r.names<-paste0(row.names(data), "_R1")##adding batch identification
row.names(data)<-r.names
write.xlsx(data,"enose_NR.xlsx",  col.names = TRUE, row.names = TRUE)
##HPLC
data<- openxlsx::read.xlsx("GCMS_N.xlsx", rowNames = T, check.names = F,sep.names=" ")
data$Odor<-round(data$Odor)
r.names<-paste0(row.names(data), "_R1")##adding batch identification
row.names(data)<-r.names
write.xlsx(data,"GCMS_NR.xlsx",  col.names = TRUE, row.names = TRUE)
##GCMS
data<- openxlsx::read.xlsx("HPLC_N.xlsx", rowNames = T, check.names = F,sep.names=" ")
data$Odor<-round(data$Odor)
r.names<-paste0(row.names(data), "_R1")##adding batch identification
row.names(data)<-r.names
write.xlsx(data,"HPLC_NR.xlsx",  col.names = TRUE, row.names = TRUE)
########################################################################## PCA plots ################################################################

##Enose PCA before replicate aggregating
data<- openxlsx::read.xlsx("enose_NR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
pca12(data, "Enose")
pca13(data, "Enose")
##Enose PCA After replicate aggregating
data<- openxlsx::read.xlsx("enose_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
pca12(data, "Enose")
pca13(data, "Enose")
##HPLC PCA before replicate aggregating
data<- openxlsx::read.xlsx("HPLC_NR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
pca12(data, "HPLC")
pca13(data, "HPLC")
##HPLC PCA After replicate aggregating
data<- openxlsx::read.xlsx("HPLC_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
pca12(data, "HPLC")
pca13(data, "HPLC")
##GCMS PCA before replicate aggregating
data<- openxlsx::read.xlsx("GCMS_NR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
pca12(data, "GCMS")
pca13(data, "GCMS")
##GCMS PCA After replicate aggregating
data<- openxlsx::read.xlsx("GCMS_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
pca12(data, "GCMS")
pca13(data, "GCMS")
########################################################################## Classification Machine Learning################################################################
source("ML_function.R")
#####################################################################Enose###################################################
###Enose before replicate aggregation
##load data enose
data<- openxlsx::read.xlsx("enose_NR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
data$Odor<-as.factor(data$Odor)
##split data
set.seed(1234)
trainIndex <- createDataPartition(as.factor(data$Odor), p = .7, 
                                  list = FALSE, group = 3, 
                                  times = 100)

##training Classification
set.seed(124)
sL<-classification(data,"svmLinear","Enose")
knn<-classification(data,"knn","Enose")
rf<-classification(data,"rf","Enose")
sP<-classification(data,"svmPoly","Enose")
sR<-classification(data,"svmRadial","Enose")
nb<-classification(data,"nb","Enose")
#saveRDS(sL$model, "SVMLinearEnose_odor.rds")
#######################################################################################
###Enose After replicate aggregation
##load data enose
data<- openxlsx::read.xlsx("enose_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
data$Odor<-as.factor(data$Odor)
##split data
set.seed(1234)
trainIndex <- createDataPartition(as.factor(data$Odor), p = .7, 
                                  list = FALSE, group = 3, 
                                  times = 100)

##training Classification
set.seed(124)
sL<-classification(data,"svmLinear","Enose")
knn<-classification(data,"knn","Enose")
rf<-classification(data,"rf","Enose")
sP<-classification(data,"svmPoly","Enose")
sR<-classification(data,"svmRadial","Enose")
nb<-classification(data,"nb","Enose")
##save winner model in rds 
saveRDS(sL$model, "SVMLinearEnose_odor.rds")
#############################################################################HPLC###########################################
###HPLC
###HPLC before replicate aggregation
##load data HPLC
data<- openxlsx::read.xlsx("HPLC_NR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
data$Odor<-as.factor(data$Odor)
##split data
set.seed(1235)
trainIndex <- createDataPartition(as.factor(data$Odor), p = .7, 
                                  list = FALSE, group = 3, 
                                  times = 100)

##training Classification
set.seed(125)
sL<-classification(data,"svmLinear","HPLC")
knn<-classification(data,"knn","HPLC")
rf<-classification(data,"rf","HPLC")
sP<-classification(data,"svmPoly","HPLC")
sR<-classification(data,"svmRadial","HPLC")
nb<-classification(data,"nb","HPLC")
saveRDS(rf$model, "HPLC.CBF.rds")
##confusion matrix
rf$confusionMatrix$table
#######################################################################################
###HPLC After replicate aggregation
#load data HPLC
data<- openxlsx::read.xlsx("HPLC_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
data$Odor<-as.factor(data$Odor)
##split data
set.seed(1235)
trainIndex <- createDataPartition(as.factor(data$Odor), p = .7, 
                                  list = FALSE, group = 3, 
                                  times = 100)

##training Classification
set.seed(125)
sL<-classification(data,"svmLinear","HPLC")
knn<-classification(data,"knn","HPLC")
rf<-classification(data,"rf","HPLC")
sP<-classification(data,"svmPoly","HPLC")
sR<-classification(data,"svmRadial","HPLC")
nb<-classification(data,"nb","HPLC")
saveRDS(sL$model, "SVMLinearHPLC_odor.rds")

#################################################################GCMS#######################################################
###GCMS Before replicate aggregation
##load data GCMS
data<- openxlsx::read.xlsx("GCMS_NR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
data$Odor<-as.factor(data$Odor)
##split data
set.seed(1235)
trainIndex <- createDataPartition(as.factor(data$Odor), p = .7, 
                                  list = FALSE, group = 3, 
                                  times = 100)

##training Classification
set.seed(125)
sL<-classification(data,"svmLinear","GCMS")
knn<-classification(data,"knn","GCMS")
rf<-classification(data,"rf","GCMS")
sP<-classification(data,"svmPoly","GCMS")
sR<-classification(data,"svmRadial","GCMS")
nb<-classification(data,"nb","GCMS")
##confusion matrix
sR$confusionMatrix$table
saveRDS(sR$model, "BGCMS.CBF.rds")
#######################################################################################
###HPLC After replicate aggregation
##load data GCMS
data<- openxlsx::read.xlsx("GCMS_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data<-data <- subset(data, select = -c(TVC))
data$Odor<-as.factor(data$Odor)
##split data
set.seed(1235)
trainIndex <- createDataPartition(as.factor(data$Odor), p = .7, 
                                  list = FALSE, group = 3, 
                                  times = 100)

##training Classification
set.seed(125)
sL<-classification(data,"svmLinear","GCMS")
knn<-classification(data,"knn","GCMS")
rf<-classification(data,"rf","GCMS")
sP<-classification(data,"svmPoly","GCMS")
sR<-classification(data,"svmRadial","GCMS")
nb<-classification(data,"nb","GCMS")
#saveRDS(sL$model, "SVMLinearHPLC_odor.rds")


########################################################################## Regression Machine Learning################################################################
##Enose
##Load data
data<- openxlsx::read.xlsx("enose_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data <- subset(data, select = -c(Odor))
#prepare training scheme
control <- trainControl(method="repeatedcv", number=5, repeats=1)

# Split training and test
set.seed(123)
train.index <- createDataPartition(data$TVC, p = .7,list = FALSE, groups=3, times = 100)

##training regression
set.seed(125)
sL<-reg.run(data,"svmLinear","Enose")
knn<-reg.run(data,"knn","Enose")
rf<-reg.run(data,"rf","Enose")
sP<-reg.run(data,"svmPoly","Enose")
sR<-reg.run(data,"svmRadial","Enose")
lm<-reg.run(data,"lm","Enose")

##GCMS
##Load data
data<- openxlsx::read.xlsx("GCMS_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data <- subset(data, select = -c(Odor))
#prepare training scheme
control <- trainControl(method="repeatedcv", number=5, repeats=1)

# Split training and test
set.seed(123)
train.index <- createDataPartition(data$TVC, p = .7,list = FALSE, groups=3, times = 100)

##training regression
set.seed(125)
sL<-reg.run(data,"svmLinear","GCMS")
knn<-reg.run(data,"knn","GCMS")
rf<-reg.run(data,"rf","GCMS")
sP<-reg.run(data,"svmPoly","GCMS")
sR<-reg.run(data,"svmRadial","GCMS")
lm<-reg.run(data,"lm","GCMS")
#saveRDS(model, "SVMLinearHPLC_TVC.rds")

##HPLC
##Load data
data<- openxlsx::read.xlsx("HPLC_RR.xlsx", rowNames = T, check.names = F,sep.names=" ")
data <- subset(data, select = -c(Odor))
#prepare training scheme
control <- trainControl(method="repeatedcv", number=5, repeats=1)

# Split training and test
set.seed(123)
train.index <- createDataPartition(data$TVC, p = .7,list = FALSE, groups=3, times = 100)

##training regression
set.seed(125)
sL<-reg.run(data,"svmLinear","HPLC")
knn<-reg.run(data,"knn","HPLC")
rf<-reg.run(data,"rf","HPLC")
sP<-reg.run(data,"svmPoly","HPLC")
sR<-reg.run(data,"svmRadial","HPLC")
lm<-reg.run(data,"lm","HPLC")
#saveRDS(sL$model, "SVMLinearHPLC_odor.rds")

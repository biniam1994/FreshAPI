###################################################################################################################
###################################################################################################################
#Biniam  thesis functions

library(ggplot2)
library(ggpubr)
library(assertthat)
###NA function########################################################################
library("readxl")
library("openxlsx")
library("zoo")
NA.rm.im<-function(data, Rm){
  ##Remove rows
  
  if (Rm=="yes") {
    data<-data[complete.cases(data), ]
  }
  else{
    
  }
  ##Impute Missing
  
  data[] <- lapply(data, na.aggregate)##impute NA with mean
  
  
}

######function to impute NA's based on replicates############
#############################################################
NA_tre_re<-function(data, replicates){
  n<-as.numeric(replicates)
  df_ls <- lapply(seq(1, nrow(data), by=n), function(i) 
    data[i: pmin((i+(n-1)), nrow(data)),])
  length(df_ls)
  
  
  da<-c()
  ##imputing NA's with mean of replicate
  for (i in 1:length(df_ls)) {
    dk<-(df_ls[i])
    dk<-as.data.frame(dk, check.names = F)
    dk[] <- lapply(dk, na.aggregate)
    dk<-list(dk)
    da<-c(da, dk)
    
  }
  
  ##combinig back all the subset(replicate) list of dataframe to one dataframe
  db<-rbind()
  ##imputing NA's with mean of replicate
  for (i in 1:length(da)) {
    dc<-(da[i])
    dc<-as.data.frame(dc, check.names = F)
    db<-rbind(db, dc)
    
  }
  data<-db
}
#####################################################################
####X character removal function from col names in dataframe#########
destroyX = function(data) {
  df = data
  for (col in c(1:ncol(df))){ #for each column in dataframe
    if (startsWith(colnames(df)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(df)[col] <- substr(colnames(df)[col], 2, 100) #get rid of 'X'
    }
  }
  assign(deparse(substitute(data)), df, inherits = TRUE) #assign corrected data to original name
}
#####################################################################

#########bacterial type selector for ML training
Bacteria_sel<-function(data, bacteria, platform){
  if(platform=="MSI"){
    if(bacteria=="TVC"){
      cols<-c("TVC","Mean.01", "Mean.02", "Mean.03", "Mean.04", "Mean.05", "Mean.06",
              "Mean.07", "Mean.08", "Mean.09", "Mean.10", "Mean.11", "Mean.12",
              "Mean.13", "Mean.14", "Mean.15", "Mean.16", "Mean.17", "Mean.18")
      data<-data[,which(names(data) %in% c(cols))]
    }else if(bacteria=="Ps"){
      cols<-c("Ps","Mean.01", "Mean.02", "Mean.03", "Mean.04", "Mean.05", "Mean.06",
              "Mean.07", "Mean.08", "Mean.09", "Mean.10", "Mean.11", "Mean.12",
              "Mean.13", "Mean.14", "Mean.15", "Mean.16", "Mean.17", "Mean.18")
      data<-data[,cols]
    }
    else if(bacteria=="LAB"){
      cols<-c("LAB","Mean.01", "Mean.02", "Mean.03", "Mean.04", "Mean.05", "Mean.06",
              "Mean.07", "Mean.08", "Mean.09", "Mean.10", "Mean.11", "Mean.12",
              "Mean.13", "Mean.14", "Mean.15", "Mean.16", "Mean.17", "Mean.18")
      data<-data[,cols]
    }
    else if(bacteria=="Bth"){
      cols<-c("Bth","Mean.01", "Mean.02", "Mean.03", "Mean.04", "Mean.05", "Mean.06",
              "Mean.07", "Mean.08", "Mean.09", "Mean.10", "Mean.11", "Mean.12",
              "Mean.13", "Mean.14", "Mean.15", "Mean.16", "Mean.17", "Mean.18")
      data<-data[,cols]
    }
    else if(bacteria=="Etb"){
      cols<-c("Etb","Mean.01", "Mean.02", "Mean.03", "Mean.04", "Mean.05", "Mean.06",
              "Mean.07", "Mean.08", "Mean.09", "Mean.10", "Mean.11", "Mean.12",
              "Mean.13", "Mean.14", "Mean.15", "Mean.16", "Mean.17", "Mean.18")
      data<-data[,cols]
    }
  }else{
    if (bacteria=="TVC") {
      data<-data[,-which(names(data) %in% c("Ps","Bth","LAB","Odor"))]
    }
    else if (bacteria=="Ps") {
      data<-data[,-which(names(data) %in% c("TVC","Odor","Bth","LAB"))]
    }
    else if (bacteria=="Bth"){
      data<-data[,-which(names(data) %in% c("TVC","Odor","Ps","LAB"))]
    }
    else if (bacteria=="LAB"){
      data<-data[,-which(names(data) %in% c("TVC","Odor","Ps","Bth"))]
    }
    else if (bacteria=="Odor"){
      data<-data[,-which(names(data) %in% c("TVC","Bth","Ps","LAB"))]
    }
  }
}
########################################################################################################
####Rep function########################################################################################
library("readxl")
library("openxlsx")
library("zoo")##imputing NA
#source("CreateObj.R")

#################### 1. CALCULATE MEANS BASED ON REPLICATES FOR FTIR Data only##########################
replicate_tre<-function(data, products, replicates){
  #Hint
  #data = FTIR data frame that contain spectra in col and sample in row
  #products = Type of product(chicken meat) eg.(Liver= for checkenLiver product, Burger= for checkenBurger....ect)
  #n= number of replicates that will be aggregated e.(2 microbiological replicates and 3 spectral/technical replicates that means n=2*3) 
  n=replicates
  df <- aggregate(data[,2:ncol(data)],list(rep(1:(nrow(data[,2:ncol(data)])%/%n+1),each=n,len=nrow(data[,2:ncol(data)]))),mean)[-1];
  #labels for samples created based on mean from 6 replicates
  samplesM <- c()
  timeStorM <- c()
  
  samplesV <- c()
  for (i in 1:nrow(data)){
    samplesV[i] <- data[i, 1]
  }
  
  timeStor <- c()
  for (i in 1:nrow(data)){
    timeStor[i] <- data[i, 2]
  }
  
  nth_element <- function(vector, starting_position, n) { 
    vector[seq(starting_position, length(vector), n)] 
  }
  samplesM <- nth_element(samplesV, 1, n)##get nth sample name 
  
  timeStorM <- nth_element(timeStor, 1, n)##get nth time storage
  #############################################Preparing new rowNames#################################################
  #Liver=product
  products=products
  product<-c()
  for(i in 1:length(samplesM)){
    product[i] <- substring(products, 1, 1)
  }
  
  #temp 
  temp <- c()
  for(i in 1:length(samplesM)){
    #temp[i] <- sub(".*", "", samplesM[i])
    temp[i]<-substr(samplesM[i], 1,2)
    temp[i]<-gsub("[^0-9\\D]", "", temp[i])##omit non-numeric
  }
  #sample batch extracting.
  
  b <- c()#storing last chara in vector
  for(i in 1:length(samplesM)){
    b[i]<-sub("^.*(.)$", "\\1", samplesM[i])##get the last character from the rownames to identify batch
    
  }
  s<-letters
  K<-replace(b, b=="A" | b=="B" | b=="C" | b=="D" | b=="E" | b=="F","R1")
  K<-replace(K, K==2, "R2")##Replace batch number to Letters eg.(R1,R2,R3,R4, R5, R6)
  K<-replace(K, K==3, "R3")
  K<-replace(K, K==4, "R4")
  K<-replace(K, K==5, "R5")
  K<-replace(K, K==6, "R6")
  # Merge temp + product + batch (A, B, c...) + time
  sample_name<- paste0(temp, product,timeStorM, "_", K)
  sample_name<- as.matrix(sample_name)
  Renamed_df<- CreateObj(sample_name, df)
  row.names(Renamed_df) <- Renamed_df$V1
  Renamed_df<- Renamed_df[,-c(1:2)]
}
#####replicate treatmeant for HPLC,eNose, GCMS
repli_tre<-function(data, replicates){
  n<-as.numeric(replicates)
  df <- aggregate(data[,1:ncol(data)],list(rep(1:(nrow(data[,1:ncol(data)])%/%n+1),each=n,len=nrow(data[,1:ncol(data)]))),mean)[-1];
  
  samplesM<-(row.names(data))
  nth_element <- function(vector, starting_position, n) { 
    vector[seq(starting_position, length(vector), n)] 
  }
  samplesV <- nth_element(samplesM, 1, n)
  data<-CreateObj(samplesV, df)
  row.names(data)<-data[,1]
  data<-data[,-1]
}
####PCA1 function for FTIR platfrom########################################################################################

library("openxlsx")
library("dplyr")##used to shift columns
library(factoextra)
#library("mixOmics")
#source("CreateObj.R")

## read xlsx files
#data <- read.xlsx("56LiverData.xlsx", sheet=1, rowNames = F)
pca3.function<-function(data){
  data$Batch <- ifelse(grepl("A", data[,1], ignore.case = T), "A", 
                       ifelse(grepl("B", data[,1], ignore.case = T), "B",
                              ifelse(grepl("C", data[,1], ignore.case = T), "C", 
                                     ifelse(grepl("D", data[,1], ignore.case = T), "D","" ))))##check For batch identity from rowname
  row.names(data)<-data[,1]
  data<-data[,-1]
  data <- data %>% dplyr::select(Batch, everything())
  Batches<-data$Batch
  #PCA and plots
  pca <- prcomp(data[,-c(1:5)],center=T, scale=T)
  fviz_pca_ind(pca,
               axes = c(1,2),
               geom = c("point", "text"),
               col.ind = Batches,
               addEllipses = TRUE,
               ellipse.level = 0.95,
               legend.title = "Batches",
               title = "PCA PLOTS"
  )
  
}
####PCA2 function for FTIR platfrom########################################################################################

library("openxlsx")
library("dplyr")##used to shift columns
library(factoextra)
#library("mixOmics")
#source("CreateObj.R")

## read xlsx files

pca4.function<-function(data){
  data$Batch <- ifelse(grepl("A", data[,1], ignore.case = T), "A", 
                       ifelse(grepl("B", data[,1], ignore.case = T), "B",
                              ifelse(grepl("C", data[,1], ignore.case = T), "C", 
                                     ifelse(grepl("D", data[,1], ignore.case = T), "D","" ))))##check For batch identity from rowname
  row.names(data)<-data[,1]
  data<-data[,-1]
  data <- data %>% dplyr::select(Batch, everything())
  Batches<-data$Batch
  #PCA and plots
  pca <- prcomp(data[,-c(1:5)],center=T, scale=T)
  fviz_pca_ind(pca,
               axes = c(1,3),
               geom = c("point", "text"),
               col.ind = Batches,
               addEllipses = TRUE,
               ellipse.level = 0.95,
               legend.title = "Batches",
               title = "PCA PLOTS"
  )
  
}
####PCA3 function for MSI/FTIR/Enose/GCMS/HPLC platfrom########################################################################################

library("openxlsx")
library("dplyr")##used to shift columns
library(factoextra)
#library("mixOmics")
#source("CreateObj.R")

pca1.function<-function(data){
  data$Batch <- ifelse(grepl("R1", row.names(data), ignore.case = T), "A", 
                       ifelse(grepl("R2", row.names(data), ignore.case = T), "B",
                              ifelse(grepl("R3", row.names(data), ignore.case = T), "C", 
                                     ifelse(grepl("R4", row.names(data), ignore.case = T), "D","" ))))##check For batch identity from rowname
  
  data <- data %>% dplyr::select(Batch, everything())
  Batches<-data$Batch
  data<-data[,-which(names(data) %in% c("TVC","Ps","LAB","Bth", "Batch"))]
  #PCA and plots
  pca <- prcomp(data,center=T, scale=T)
  fviz_pca_ind(pca,
               axes = c(1,2),
               geom = c("point", "text"),
               col.ind = Batches,
               addEllipses = TRUE,
               ellipse.level = 0.95,
               legend.title = "Batches",
               title = "PCA PLOTS"
  )
  
}
####PCA4 function for MSI/FTIR/Enose/GCMS/HPLC platfrom########################################################################################

library("openxlsx")
library("dplyr")##used to shift columns
library(factoextra)
#library("mixOmics")
#source("CreateObj.R")

pca2.function<-function(data){
  data$Batch <- ifelse(grepl("R1", row.names(data), ignore.case = T), "A", 
                       ifelse(grepl("R2", row.names(data), ignore.case = T), "B",
                              ifelse(grepl("R3", row.names(data), ignore.case = T), "C", 
                                     ifelse(grepl("R4", row.names(data), ignore.case = T), "D","" ))))##check For batch identity from rowname
  
  data <- data %>% dplyr::select(Batch, everything())
  Batches<-data$Batch
  data<-data[,-which(names(data) %in% c("TVC","Ps","Bth","LAB", "Batch"))]
  #PCA and plots
  pca <- prcomp(data,center=T, scale=T)
  fviz_pca_ind(pca,
               axes = c(1,3),
               geom = c("point", "text"),
               col.ind = Batches,
               addEllipses = TRUE,
               ellipse.level = 0.95,
               legend.title = "Batches",
               title = "PCA PLOTS"
  )
  
}

####ML function for mixed batch training########################################################################################
####ML function for mixed batch training########################################################################################
library(gmodels)
library("randomForest")
library(caret)
library(kernlab)
library(ggpubr)
####
#set.seed(123)
regression.run.pra <- function (data,platform, method,bacteria, iteration, pra) {
  #Hint
  #data= data frame with whole samples containg target variable(bacteria counts)
  #method= machine learning algorithms eg("rf", "knn","svmPoly", "svmRadial","svmLinear", "lm")
  #platform= analytical platform eg("FTIR","MSI", "HPLC", "Enose")
  #bacteria= is bacterial type i.e in the data frame that we want to predict eg(TSA, MRS, CFC, STAA, PS)
  #iteration = is number of iteration you want to run eg(50 or 100)
  # prepare training scheme
  iteration<-as.numeric(iteration)
  
  control <- trainControl(method="repeatedcv", number=5, repeats=3)
  
  set.seed(123)
  if(bacteria=="TVC"){
    train.index <- createDataPartition(data$TVC, p = .7,list = FALSE, groups=3, times = iteration)
  }else if(bacteria=="LAB"){
    train.index <- createDataPartition(data$LAB, p = .7,list = FALSE, groups=3, times = iteration)
  }else if(bacteria=="Bth"){
    train.index <- createDataPartition(data$Bth, p = .7,list = FALSE, groups=3, times = iteration)
  }else if(bacteria=="Ps"){
    train.index <- createDataPartition(data$Ps, p = .7,list = FALSE, groups=3, times = iteration)
  }
  RMSEC<-c()
  for (i in 1:iteration) {
    trainN <- data[train.index[,i],]
    testN <- data[-train.index[,i],]
    
    if(method=="rf"){
      if(bacteria=="TVC"){
        if(nrow(data)>=50){
        model <<- train(TVC ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
        }else{
          model <<- train(TVC ~ ., method= method, data=trainN, tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
          }
      }
      else if (bacteria=="Ps") {
        if(nrow(data)>=50){
        model <<- train(Ps ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry),ntree=pra$ntree )##train
        }else{
          model <<- train(Ps ~ ., method= method, data=trainN, tuneGrid=expand.grid(mtry=pra$mtry),ntree=pra$ntree )##train
          
        }
      }
      else if (bacteria=="LAB") {
        if(nrow(data)>=50){
        model <<- train(LAB ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
        }else{
          model <<- train(LAB ~ ., method= method, data=trainN,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
          
        }
      }
      else if (bacteria=="Bth") {
        if(nrow(data)>=50){
        model <<- train(Bth ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
        }else{
          model <<- train(Bth ~ ., method= method, data=trainN,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
        }
      }
      
      
    }
    else if(method=="knn"){
      if(bacteria=="TVC"){
        if(nrow(data)>=50){
        model <<- train(TVC ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
        }else{
          model <<- train(TVC ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train
        }
      }
      else if(bacteria=="Ps"){
        if(nrow(data)>=50){
        model <<- train(Ps ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
        }else{
          model <<- train(Ps ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train
        }
      }
      else if(bacteria=="LAB"){
        if(nrow(data)>=50){
        model <<- train(LAB ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
        }else{
          model <<- train(LAB ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train
        }
      }
      else if(bacteria=="Bth"){
        if(nrow(data)>=50){
        model <<- train(Bth ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
        }else{
          model <<- train(Bth ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train 
        }
      }
      
    }
    else if(method=="svmLinear"){
      if(bacteria=="TVC"){
        if(nrow(data)>=50){
        model <<- train(TVC ~., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(C=pra$C))##train
        }else{
          model <<- train(TVC ~., method= method, data=trainN, tuneGrid=expand.grid(C=pra$C))##train
        }
      }
      else if(bacteria=="Ps"){
        if(nrow(data)>=50){
        model <<- train(Ps ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(C=pra$C))##train
        }else{
          model <<- train(Ps ~., method= method, data=trainN,tuneGrid=expand.grid(C=pra$C))##train
        }
      }
      else if(bacteria=="LAB"){
        if(nrow(data)>=50){
        model <<- train(LAB ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(C=pra$C))##train
        }else{
          model <<- train(LAB ~., method= method, data=trainN,tuneGrid=expand.grid(C=pra$C))##train
        }
      }
      else if(bacteria=="Bth"){
        if(nrow(data)>=50){
        model <<- train(Bth ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(C=pra$C))##train
        }else{
          model <<- train(Bth ~., method= method, data=trainN,tuneGrid=expand.grid(C=pra$C))##train
          }
      }
      
    }
    else if(method=="svmRadial"){
      if(bacteria=="TVC"){
        if(nrow(data)>=50){
        model <<- train(TVC ~., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }else{
          model <<- train(TVC ~., method= method, data=trainN,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }
      }
      else if(bacteria=="Ps"){
        if(nrow(data)>=50){
        model <<- train(Ps ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }else{
          model <<- train(Ps ~., method= method, data=trainN,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }
      }
      else if(bacteria=="Bth"){
        if(nrow(data)>=50){
        model <<- train(Bth ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }else{
          model <<- train(Bth ~., method= method, data=trainN, tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }
      }
      else if(bacteria=="LAB"){
        if(nrow(data)>=50){
        model <- train(LAB ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }else{
          model <- train(LAB ~., method= method, data=trainN, tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
        }
      }
      
    }
    else if(method=="svmPoly"){
      if(bacteria=="TVC"){
        if(nrow(data)>=50){
        model <<- train(TVC ~., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }else{
          model <<- train(TVC ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }
      }
      else if(bacteria=="Ps"){
        if(nrow(data)>=50){
        model <<- train(Ps ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }else{
          model <<- train(Ps ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }
      }
      else if(bacteria=="LAB"){
        if(nrow(data)>=50){
        model <<- train(LAB ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }else{
          model <<- train(LAB ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }
      }
      else if(bacteria=="Bth"){
        if(nrow(data)>=50){
        model <<- train(Bth ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }else{
          model <<- train(Bth ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
        }
      }
      
    }
    predicted <- predict(model, testN)##Prediction
    
    if (bacteria=="TVC") {
      RMSE <- RMSE(predicted, testN$TVC)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }else if (bacteria=="Ps") {
      RMSE <- RMSE(predicted, testN$Ps)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }else if (bacteria=="LAB") {
      RMSE <- RMSE(predicted, testN$LAB)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }else if (bacteria=="Bth") {
      RMSE <- RMSE(predicted, testN$Bth)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }
  
  }
  
  Mean.RMSE<-mean(RMSEC)##calculating mean
  SD.RMSE<-sd(RMSEC)##calculating SD
  Cl <-gmodels::ci(RMSEC, confidence = 0.95)
  iterations <- as.array(c(1:iteration))
  ##combining both predicted and actual bacterial counts
  df<-cbind(testN[,bacteria], predicted)
  df<-as.data.frame(df)
  colnames(df)=c("Actual", "Predicted")
  dm<<-df
  ##combinig iteration and RMSE
  df1<-cbind(iterations, RMSEC)
  df1<-as.data.frame(df1)
  colnames(df1)=c("Iteration", "RMSE")
  ##plot of RMSE over iteration
  st<-ggscatter(df1,x="Iteration", y="RMSE", ylim=c(0,2), xlim=c(0,iteration), xlab="Iteration", ylab="RMSE", main=paste(method," Model Mean RMSE:",
                                                                                                                         round(Mean.RMSE, digits = 3), "CI 95%:", 
                                                                                                                         round(Cl[2], digits = 3), "-", 
                                                                                                                         round(Cl[3], digits = 3), "SD:",
                                                                                                                         round(SD.RMSE,digits = 3)))##Plotting RMSE
  st<<-st + geom_line( color="black") + geom_point(shape=1, color="black", fill="#69b3a2", size=1)                                                                                         
  ## plot predicted values vs Actual/Observed values
  if (bacteria=="TVC") {
    Accuracy<- Percentage(predicted,testN$TVC)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "TVC", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
  }
  else if (bacteria=="Ps") {
    Accuracy<- Percentage(predicted,testN$Ps)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Ps", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="LAB") {
    Accuracy<- Percentage(predicted,testN$LAB)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "LAB", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="Bth") {
    Accuracy<- Percentage(predicted,testN$Bth)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Bth", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
  }
  
}
###Train Function on Batch A
set.seed(123)
regression.batchA.pra<- function (data, platform, method,Batch, bacteria,pra) {
  #Hint 
  # prepare training scheme
  control <- trainControl(method="repeatedcv", number=5, repeats=3)
  #train.index <- createDataPartition(data$bacteria, p = .7,list = FALSE, groups=3, times = iteration)
  ##separate data into Batches
  if (platform=="MSI" | platform=="FTIR" | platform=="HPLC" | platform=="GCMS" | platform=="eNose") {
    Batchtr<-Batch$Batchtrain
    Batchts<-Batch$Batchtest
    Batchtrain<-data[unlist(sapply(Batchtr, grep, row.names(data), USE.NAMES = F)), ]
    Batchtest<-data[unlist(sapply(Batchts, grep, row.names(data), USE.NAMES = F)), ]
  
  }
  ##
  set.seed(123)
  trainN <- Batchtrain
  testN <- Batchtest
  ##model training
  if(method=="rf"){
    if(bacteria=="TVC"){
      if(nrow(data)>=50){
        model <<- train(TVC ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
      }else{
        model <<- train(TVC ~ ., method= method, data=trainN, tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
      }
    }
    else if (bacteria=="Ps") {
      if(nrow(data)>=50){
        model <<- train(Ps ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry),ntree=pra$ntree )##train
      }else{
        model <<- train(Ps ~ ., method= method, data=trainN, tuneGrid=expand.grid(mtry=pra$mtry),ntree=pra$ntree )##train
        
      }
    }
    else if (bacteria=="LAB") {
      if(nrow(data)>=50){
        model <<- train(LAB ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
      }else{
        model <<- train(LAB ~ ., method= method, data=trainN,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
        
      }
    }
    else if (bacteria=="Bth") {
      if(nrow(data)>=50){
        model <<- train(Bth ~ ., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
      }else{
        model <<- train(Bth ~ ., method= method, data=trainN,tuneGrid=expand.grid(mtry=pra$mtry), ntree=pra$ntree )##train
      }
    }
    
    
  }
  else if(method=="knn"){
    if(bacteria=="TVC"){
      if(nrow(data)>=50){
        model <<- train(TVC ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
      }else{
        model <<- train(TVC ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train
      }
    }
    else if(bacteria=="Ps"){
      if(nrow(data)>=50){
        model <<- train(Ps ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
      }else{
        model <<- train(Ps ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train
      }
    }
    else if(bacteria=="LAB"){
      if(nrow(data)>=50){
        model <<- train(LAB ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
      }else{
        model <<- train(LAB ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train
      }
    }
    else if(bacteria=="Bth"){
      if(nrow(data)>=50){
        model <<- train(Bth ~ ., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(k=pra$k))##train
      }else{
        model <<- train(Bth ~ ., method= method, data=trainN, tuneGrid=expand.grid(k=pra$k))##train 
      }
    }
    
  }
  else if(method=="svmLinear"){
    if(bacteria=="TVC"){
      if(nrow(data)>=50){
        model <<- train(TVC ~., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(C=pra$C))##train
      }else{
        model <<- train(TVC ~., method= method, data=trainN, tuneGrid=expand.grid(C=pra$C))##train
      }
    }
    else if(bacteria=="Ps"){
      if(nrow(data)>=50){
        model <<- train(Ps ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(C=pra$C))##train
      }else{
        model <<- train(Ps ~., method= method, data=trainN,tuneGrid=expand.grid(C=pra$C))##train
      }
    }
    else if(bacteria=="LAB"){
      if(nrow(data)>=50){
        model <<- train(LAB ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(C=pra$C))##train
      }else{
        model <<- train(LAB ~., method= method, data=trainN,tuneGrid=expand.grid(C=pra$C))##train
      }
    }
    else if(bacteria=="Bth"){
      if(nrow(data)>=50){
        model <<- train(Bth ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(C=pra$C))##train
      }else{
        model <<- train(Bth ~., method= method, data=trainN,tuneGrid=expand.grid(C=pra$C))##train
      }
    }
    
  }
  else if(method=="svmRadial"){
    if(bacteria=="TVC"){
      if(nrow(data)>=50){
        model <<- train(TVC ~., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }else{
        model <<- train(TVC ~., method= method, data=trainN,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }
    }
    else if(bacteria=="Ps"){
      if(nrow(data)>=50){
        model <<- train(Ps ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }else{
        model <<- train(Ps ~., method= method, data=trainN,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }
    }
    else if(bacteria=="Bth"){
      if(nrow(data)>=50){
        model <<- train(Bth ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }else{
        model <<- train(Bth ~., method= method, data=trainN, tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }
    }
    else if(bacteria=="LAB"){
      if(nrow(data)>=50){
        model <<- train(LAB ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }else{
        model <<- train(LAB ~., method= method, data=trainN, tuneGrid=expand.grid(sigma=pra$S,C=pra$C))##train
      }
    }
    
  }
  else if(method=="svmPoly"){
    if(bacteria=="TVC"){
      if(nrow(data)>=50){
        model <<- train(TVC ~., method= method, data=trainN, trControl= control, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }else{
        model <<- train(TVC ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }
    }
    else if(bacteria=="Ps"){
      if(nrow(data)>=50){
        model <<- train(Ps ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }else{
        model <<- train(Ps ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }
    }
    else if(bacteria=="LAB"){
      if(nrow(data)>=50){
        model <<- train(LAB ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }else{
        model <<- train(LAB ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }
    }
    else if(bacteria=="Bth"){
      if(nrow(data)>=50){
        model <<- train(Bth ~., method= method, data=trainN, trControl= control,tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }else{
        model <<- train(Bth ~., method= method, data=trainN, tuneGrid=expand.grid(degree=pra$D, C=pra$C, scale=pra$S))##train
      }
    }
    
  }
  
  predicted <- predict(model, testN)##Prediction
  if (bacteria=="TVC") {
    RMSE <- RMSE(predicted, testN$TVC)##RMSE calculation
  }else if (bacteria=="Ps") {
    RMSE <- RMSE(predicted, testN$Ps)##RMSE calculation
  }else if (bacteria=="LAB") {
    RMSE <- RMSE(predicted, testN$LAB)##RMSE calculation
  }else if (bacteria=="Bth") {
    RMSE <- RMSE(predicted, testN$Bth)##RMSE calculation
  }
  ## model Accuracy +- 1 log count
  if (bacteria=="TVC") {
    Accuracy<- Percentage(predicted,testN$TVC)
  }else if (bacteria=="Ps") {
    Accuracy<- Percentage(predicted,testN$Ps)
  }else if (bacteria=="LAB") {
    Accuracy<- Percentage(predicted,testN$LAB)
  }else if (bacteria=="Bth") {
    Accuracy<- Percentage(predicted,testN$Bth)
  }
  
  ##combining both predicted and actual bacterial counts
  df<-cbind(testN[,bacteria], predicted)
  df<-as.data.frame(df)
  colnames(df)=c("Actual", "Predicted")
  dm<<-df
  ## plot predicted values vs Actual/Observed values
  if (bacteria=="TVC") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "TVC", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="Ps") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Ps", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="LAB") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "LAB", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="Bth") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Bth", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  
}

###########################################################
# Model prediction Accuracy calculating function.
#The accuracy calculation function work based on +- one log bacterial count of predicted values from the actual values.

# x = is predicted values using the trained model
# y = is the actual test data 
Percentage<-function(x, y){
  SVMLM_TVC<-c()
  for ( i in 1:length(x)){
    if((x[i])<=((y[i])+1) & x[i]>=((y[i])-1)){
      G<-i
      SVMLM_TVC<-c(SVMLM_TVC, G)
      
    }
    else{
      
    }
    
  }
  percentage_SVMLM_TVC<- length(SVMLM_TVC)/length(y)*100
}
###########################################
#####multiple object saving to RDS object
ObjSave <- function(..., folder) {
  objects <- list(...)
  object_names <- sapply(substitute(list(...))[-1], deparse)
  sapply(1:length(objects), function(i) {
    filename = paste0(folder, "/", object_names[i], ".rds")
    saveRDS(objects[i], filename)
  })
}
############################################################
####adding csv serializer content-type function#############
serializer_csv <- function(){
  function(val, req, res, errorHandler){
    tryCatch({
      res$setHeader("Content-Type", "text/plain")
      res$setHeader("Content-Disposition", 'attachment; filename="xxx.csv"')
      res$body <- paste0(val, collapse="\n")
      return(res$toResponse())
    }, error=function(e){
      errorHandler(req, res, e)
    })
  }
}

if (is.null(plumber::plumber$parent_env$.globals$serializers$csv)){
  plumber::addSerializer("csv", serializer_csv)
}

###########################################################
##############Prediction of bacterial counts ##############
predict_B <- function(data,platform, product, bacteria){
  prediction<- 0
  platform<-platform
  product<- product
  bacteria<-bacteria
  model_name<-paste0(platform,".",product,".",bacteria, ".rds")
  wd <- getwd()
  f <- paste0(wd,"/","models","/", model_name)
  if (see_if(has_extension(model_name, 'rds')) == TRUE){
    model<- readRDS(f)
    if(length(model)==1){##extract if model in list
    model<-model[[1]]
    }else if(length(model)!=1){
      model<-model
    }
  }else {
    print("Error in reading model")
  }
  prediction <- stats::predict(model, data) #predict 
}
###########################################################
##############Prediction of sensory scores   ##############
predict_C <- function(data,platform, product){
  prediction<- 0
  platform<-platform
  product<- product
  model_name<-paste0(platform,".",product, ".rds")
  wd <- getwd()
  f <- paste0(wd,"/","models","/", model_name)
  if (see_if(has_extension(model_name, 'rds')) == TRUE){
    model<- readRDS(f)
    if(length(model)==1){##extract if model in list
      model<-model[[1]]
    }else if(length(model)!=1){
      model<-model
    }
  }else {
    print("Error in reading model")
  }
  predicts <- stats::predict(model, data) #predict
  s <- as.numeric(predicts)
  s<-replace(s, s==1, "Fresh")
  s<-replace(s, s==2, "Semi-Fresh")
  prediction<-replace(s, s==3, "Spoiled")
}
####combining datasets old and new  for the same platform and product
#####################################################################

##Batch combinig  Function For MSI data
Batch_comb_M<- function(fname1, fname2, path){
  wd <- getwd()
  ##first dataset
  f1 <- paste0(wd,"/", path, "/",fname1)
  if (see_if(has_extension(fname1, 'csv')) == TRUE){
    data1<- utils::read.csv(f1, header = T, sep= ",", row.names = 1, check.names = F)
  }
  else if (see_if(has_extension(fname1, 'xlsx')) ==TRUE) {
    data1<- openxlsx::read.xlsx(f1, rowNames =T, check.names = F,sep.names=" ")
  }
  ##get number of Batches in the first dataset
  g<-sub(".*R", "", row.names(data1))
  g<-as.numeric(g)
  e<-max(g)
  ##second dataset
  f2 <- paste0(wd,"/", path, "/", fname2)
  if (see_if(has_extension(fname2, 'csv')) == TRUE){
    data2<- utils::read.csv(f2, check.names = F)
  }
  else if (see_if(has_extension(fname2, 'xlsx')) ==TRUE) {
    data2<- openxlsx::read.xlsx(f2,check.names = F,sep.names=" ")
  }
  colnames(data2)[1] <- "R_names" 
  ##rename rownames with new batch id
  data2$N<-sub(".*R", "", data2$R_names)
  data2$N<-(as.numeric(data2$N) + e) 
  data2$N<-paste("R", data2$N, sep = "" )
  data2$R_names<-gsub('.{2}$', "", data2$R_names)
  data2$R_names<-paste(data2$R_names,data2$N)
  row.names(data2)<- data2$R_names
  data2<-data2[,-which(names(data2) %in% c("R_names","N"))]
  ##combine both data set
  data3<-rbind(data1,data2)
  
}
###
##Batch combinig  Function For Enose, HPLC, GCMS
##Hint
##fname1 is the file name for existing  dataset and fname2 is the file name for new data
##fname2 must have only one batch as always that is the case
##
Batch_comb_EGH<- function(fname1, fname2, path){
  wd <- getwd()
  ##first dataset
  f1 <- paste0(wd,"/", path, "/",fname1)
  if (see_if(has_extension(fname1, 'csv')) == TRUE){
    data1<- utils::read.csv(f1, header = T, sep= ",", row.names = 1)
  }
  else if (see_if(has_extension(fname1, 'xlsx')) ==TRUE) {
    data1<- openxlsx::read.xlsx(f1, rowNames =T, check.names = F,sep.names=" ")
  }
  ##get number of Batches in the first dataset
  g<-sub(".*R", "", row.names(data1))
  g<-as.numeric(g)
  e<-max(g)
  ##second dataset
  f2 <- paste0(wd,"/", path, "/", fname2)
  if (see_if(has_extension(fname2, 'csv')) == TRUE){
    data2<- utils::read.csv(f2,header = T, sep= ",", row.names = 1)
  }
  else if (see_if(has_extension(fname2, 'xlsx')) ==TRUE) {
    data2<- openxlsx::read.xlsx(f2, rowNames =T, check.names = F,sep.names=" ")
  }
  
  data2$Batch<-paste0("_R", e+1) 
  row.names(data2)<-paste0(row.names(data2), data2$Batch)
  data2<-data2[,-which(names(data2) %in% c("Batch"))]
  ##combine both data set
  data3<-rbind(data1,data2)
  
}
################################################################################################################
###############################################################03/07/2020#######################################
#Parameter handling function for ML algorithms

######randomForest#####################################
#######################################################

rf_pram<-function(pram){
  prm<-pram
  p<-strsplit(prm,split=';', fixed=TRUE)
  ####
  z<-strsplit(p[[1]], " ")
  nam<-c("mtry", "ntree")
  names(z)<-nam
  ####
  ##remove ntree/mtry from string
  z$mtry<-sub(".*~", "", z$mtry)
  z$ntree<-sub(".*~", "", z$ntree)
  nt<-as.numeric(z$ntree)
  
  if (grepl(",",z$mtry[1] , fixed = TRUE)==TRUE){
    z$mtry<-strsplit(z$mtry,split=',', fixed=TRUE)
    z$mtry[[1]]<-as.numeric(z$mtry[[1]])
    mt<-c(z$mtry[[1]])
  }else if(grepl(":", z$mtry[1], fixed = TRUE)==TRUE){
    z$mtry[1]<-strsplit(z$mtry[1],split=':', fixed=TRUE)
    z$mtry[[1]]<-as.numeric(z$mtry[[1]])
    ps<-c(z$mtry[[1]])
    z$mtry$a<-ps[[1]]
    z$mtry$b<-ps[[2]]
    mt<-c(z$mtry$a:z$mtry$b)
  }
  ####
  
  df <- list(mt,nt)
  
}

###################Knn###################################################

knn_pram<-function(pram){
  prm<-pram
  p<-strsplit(prm,split=';', fixed=TRUE)
  ####
  z<-strsplit(p[[1]], " ")
  nam<-c("k")
  names(z)<-nam
  ####
  ##remove K from string
  z$k<-sub(".*~", "", z$k)
  ##check for number separater
  if (grepl(",",z$k[1] , fixed = TRUE)==TRUE){
    z$k<-strsplit(z$k,split=',', fixed=TRUE)
    z$k[[1]]<-as.numeric(z$k[[1]])
    nk<-c(z$k[[1]])
  }else if(grepl(":", z$k[1], fixed = TRUE)==TRUE){
    z$k[1]<-strsplit(z$k[1],split=':', fixed=TRUE)
    z$k[[1]]<-as.numeric(z$k[[1]])
    ps<-c(z$k[[1]])
    z$k$a<-ps[[1]]
    z$k$b<-ps[[2]]
    nk<-c(z$k$a:z$k$b)
  }else{
    nk<-as.numeric(z$k)
  }
  ####
  
  df <- list(nk)
  
}

###########
###########svmLinear############################################################

svl_pram<-function(pram){
  prm<-pram
  p<-strsplit(prm,split=';', fixed=TRUE)
  ####
  z<-strsplit(p[[1]], " ")
  nam<-c("C")
  names(z)<-nam
  ####
  ##remove ntree/mtry from string
  z$C<-sub(".*~", "", z$C)
  ##check for number separater
  if (grepl(",",z$C[1] , fixed = TRUE)==TRUE){
    z$C<-strsplit(z$C,split=',', fixed=TRUE)
    z$C[[1]]<-as.numeric(z$C[[1]])
    nc<-c(z$C[[1]])
  }else if(grepl(":", z$C[1], fixed = TRUE)==TRUE){
    z$C[1]<-strsplit(z$C[1],split=':', fixed=TRUE)
    z$C[[1]]<-as.numeric(z$C[[1]])
    ps<-c(z$C[[1]])
    z$C$a<-ps[[1]]
    z$C$b<-ps[[2]]
    nc<-c(z$C$a:z$C$b)
  }else{
    nc<-as.numeric(z$C)
  }
  ####
  
  df <- list(nc)
  
}

######
###########svmRadial############################################################

svr_pram<-function(pram){
  prm<-pram
  p<-strsplit(prm,split=';', fixed=TRUE)
  ####
  z<-strsplit(p[[1]], " ")
  nam<-c("C","S")
  names(z)<-nam
  ####
  ##remove c/s from string
  z$C<-sub(".*~", "", z$C)
  z$S<-sub(".*~", "", z$S)
  ##check for number separater for c
  if (grepl(",",z$C[1] , fixed = TRUE)==TRUE){
    z$C<-strsplit(z$C,split=',', fixed=TRUE)
    z$C[[1]]<-as.numeric(z$C[[1]])
    nc<-c(z$C[[1]])
  }else if(grepl(":", z$C[1], fixed = TRUE)==TRUE){
    z$C[1]<-strsplit(z$C[1],split=':', fixed=TRUE)
    z$C[[1]]<-as.numeric(z$C[[1]])
    ps<-c(z$C[[1]])
    z$C$a<-ps[[1]]
    z$C$b<-ps[[2]]
    nc<-c(z$C$a:z$C$b)
  }else{
    nc<-as.numeric(z$C)
  }
  ######check for number separater for sigma
  if (grepl(",",z$S[1] , fixed = TRUE)==TRUE){
    z$S<-strsplit(z$S,split=',', fixed=TRUE)
    z$S[[1]]<-as.numeric(z$S[[1]])
    ns<-c(z$S[[1]])
  }else if(grepl(":", z$S[1], fixed = TRUE)==TRUE){
    z$S[1]<-strsplit(z$S[1],split=':', fixed=TRUE)
    z$S[[1]]<-as.numeric(z$S[[1]])
    psc<-c(z$S[[1]])
    z$S$a<-psc[[1]]
    z$S$b<-psc[[2]]
    ns<-c(z$S$a:z$S$b)
  }else{
    ns<-as.numeric(z$S)
  }
  ####
  df <- list(nc,ns)
  
}

############################################################################
###svmpoly##

svp_pram<-function(pram){
  prm<-pram
  p<-strsplit(prm,split=';', fixed=TRUE)
  ####
  z<-strsplit(p[[1]], " ")
  nam<-c("C","S","D")
  names(z)<-nam
  ####
  ##remove c/s from string
  z$C<-sub(".*~", "", z$C)
  z$S<-sub(".*~", "", z$S)
  z$D<-sub(".*~", "", z$D)
  ##check for number separater for c
  if (grepl(",",z$C[1] , fixed = TRUE)==TRUE){
    z$C<-strsplit(z$C,split=',', fixed=TRUE)
    z$C[[1]]<-as.numeric(z$C[[1]])
    nc<-c(z$C[[1]])
  }else if(grepl(":", z$C[1], fixed = TRUE)==TRUE){
    z$C[1]<-strsplit(z$C[1],split=':', fixed=TRUE)
    z$C[[1]]<-as.numeric(z$C[[1]])
    ps<-c(z$C[[1]])
    z$C$a<-ps[[1]]
    z$C$b<-ps[[2]]
    nc<-c(z$C$a:z$C$b)
  }else{
    nc<-as.numeric(z$C)
  }
  ######check for number separater for scale
  if (grepl(",",z$S[1] , fixed = TRUE)==TRUE){
    z$S<-strsplit(z$S,split=',', fixed=TRUE)
    z$S[[1]]<-as.numeric(z$S[[1]])
    ns<-c(z$S[[1]])
  }else if(grepl(":", z$S[1], fixed = TRUE)==TRUE){
    z$S[1]<-strsplit(z$S[1],split=':', fixed=TRUE)
    z$S[[1]]<-as.numeric(z$S[[1]])
    psc<-c(z$S[[1]])
    z$S$a<-psc[[1]]
    z$S$b<-psc[[2]]
    ns<-c(z$S$a:z$S$b)
  }else{
    ns<-as.numeric(z$S)
  }
  ####
  ######check for number separater for degree
  if (grepl(",",z$D[1] , fixed = TRUE)==TRUE){
    z$D<-strsplit(z$D,split=',', fixed=TRUE)
    z$D[[1]]<-as.numeric(z$D[[1]])
    nd<-c(z$D[[1]])
  }else if(grepl(":", z$D[1], fixed = TRUE)==TRUE){
    z$D[1]<-strsplit(z$D[1],split=':', fixed=TRUE)
    z$D[[1]]<-as.numeric(z$D[[1]])
    pscd<-c(z$D[[1]])
    z$D$a<-pscd[[1]]
    z$D$b<-pscd[[2]]
    nd<-c(z$D$a:z$D$b)
  }else{
    nd<-as.numeric(z$D)
  }
  ####
  df <- list(nc,ns, nd)
  
}
########################################################
####ML function for mixed batch training on with random search hyperparameters#####################################################################################
library(gmodels)
library("randomForest")
library(caret)
library(kernlab)
library(ggpubr)
#######
set.seed(123)
regression.run <- function (data,platform,method,bacteria, iteration, tunelength) {
  #Hint
  #data= data frame with whole samples containg target variable(bacteria counts)
  #method= machine learning algorithms eg("rf", "knn","svmPoly", "svmRadial","svmLinear", "lm")
  #platform= analytical platform eg("FTIR","MSI", "HPLC", "Enose")
  #bacteria= is bacterial type i.e in the data frame that we want to predict eg(TSA, MRS, CFC, STAA, PS)
  #iteration = is number of iteration you want to run eg(50 or 100)
  #prepare training scheme
  #
  tunelength<-as.numeric(tunelength)
  iteration<-as.numeric(iteration)
  control <- trainControl(method="repeatedcv", number=5, repeats=1, search = "random")
  set.seed(123)
  if(bacteria=="TVC"){
    train.index <- createDataPartition(data$TVC, p = .7,list = FALSE, groups=3, times = iteration)
  }else if(bacteria=="LAB"){
    train.index <- createDataPartition(data$LAB, p = .7,list = FALSE, groups=3, times = iteration)
  }else if(bacteria=="Bth"){
    train.index <- createDataPartition(data$Bth, p = .7,list = FALSE, groups=3, times = iteration)
  }else if(bacteria=="Ps"){
    train.index <- createDataPartition(data$Ps, p = .7,list = FALSE, groups=3, times = iteration)
  }
  RMSEC<-c()
  for (i in 1:iteration) {
    trainN <- data[train.index[,i],]
    testN <- data[-train.index[,i],]
    
   
      if(bacteria=="TVC"){
        if(nrow(data)>=50){
        model <<- train(TVC ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
        }else{
          model <<- train(TVC ~ ., method= method, data=trainN, tuneLength=tunelength )##train
        }
      }
      else if (bacteria=="Ps") {
        if(nrow(data)>=50){
          model <<- train(Ps ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
        }else{
          model <<- train(Ps ~ ., method= method, data=trainN, tuneLength=tunelength )##train
        }
      }
      else if (bacteria=="LAB") {
        if(nrow(data)>=50){
          model <<- train(LAB ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
        }else{
          model <<- train(LAB ~ ., method= method, data=trainN, tuneLength=tunelength )##train
        }
      }
      else if (bacteria=="Bth") {
        if(nrow(data)>=50){
          model <<- train(Bth ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
        }else{
          model <<- train(Bth ~ ., method= method, data=trainN, tuneLength=tunelength )##train
        }
      }
      
      ##
    
    predicted <- predict(model, testN)##Prediction
    
    if (bacteria=="TVC") {
      RMSE <- RMSE(predicted, testN$TVC)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }else if (bacteria=="Ps") {
      RMSE <- RMSE(predicted, testN$Ps)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }else if (bacteria=="LAB") {
      RMSE <- RMSE(predicted, testN$LAB)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }else if (bacteria=="Bth") {
      RMSE <- RMSE(predicted, testN$Bth)##RMSE calculation
      RMSEC<-c(RMSEC, RMSE)
    }
  }
  
  Mean.RMSE<-mean(RMSEC)##calculating mean
  SD.RMSE<-sd(RMSEC)##calculating SD
  Cl <-gmodels::ci(RMSEC, confidence = 0.95)
  iterations <- as.array(c(1:iteration))
  ##combining both predicted and actual bacterial counts
  df<-cbind(testN[,bacteria], predicted)
  df<-as.data.frame(df)
  colnames(df)=c("Actual", "Predicted")
  dm<<-df
  ##combinig iteration and RMSE
  df1<-cbind(iterations, RMSEC)
  df1<-as.data.frame(df1)
  colnames(df1)=c("Iteration", "RMSE")
  ##plot of RMSE over iteration
  st<-ggscatter(df1,x="Iteration", y="RMSE", ylim=c(0,2), xlim=c(0,iteration), xlab="Iteration", ylab="RMSE", main=paste(method," Model Mean RMSE:",
                                                                                                                         round(Mean.RMSE, digits = 3), "CI 95%:", 
                                                                                                                         round(Cl[2], digits = 3), "-", 
                                                                                                                         round(Cl[3], digits = 3), "SD:",
                                                                                                                         round(SD.RMSE,digits = 3)))##Plotting RMSE
  st<<-st + geom_line( color="black") + geom_point()                                                                                         
  ## plot predicted values vs Actual/Observed values
  if (bacteria=="TVC") {
    Accuracy<- Percentage(predicted,testN$TVC)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "TVC", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
  }
  else if (bacteria=="Ps") {
    Accuracy<- Percentage(predicted,testN$Ps)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Ps", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="LAB") {
    Accuracy<- Percentage(predicted,testN$LAB)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "LAB", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="Bth") {
    Accuracy<- Percentage(predicted,testN$Bth)
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Bth", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
  }
  
}
###Train Function on Batches######################################################
set.seed(123)
regression.batchA <- function (data,platform, method,Batch,bacteria, tunelength) {
  #Hint 
  # prepare training scheme
  control <- trainControl(method="repeatedcv", number=5, repeats=3, search = 'random')
 
  ##separate data into Batches
  if (platform=="MSI" | platform=="FTIR" | platform=="HPLC" | platform=="GCMS" | platform=="eNose") {
    Batchtr<-Batch$Batchtrain
    Batchts<-Batch$Batchtest
    Batchtrain<-data[unlist(sapply(Batchtr, grep, row.names(data), USE.NAMES = F)), ]
    Batchtest<-data[unlist(sapply(Batchts, grep, row.names(data), USE.NAMES = F)), ]
    
  }
  ##
  tunelength<-as.numeric(tunelength)
  
  set.seed(123)
  trainN <- Batchtrain
  testN <- Batchtest
  
  if(bacteria=="TVC"){
    if(nrow(data)>=50){
      model <<- train(TVC ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
    }else{
      model <<- train(TVC ~ ., method= method, data=trainN, tuneLength=tunelength )##train
    }
  }
  else if (bacteria=="Ps") {
    if(nrow(data)>=50){
      model <<- train(Ps ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
    }else{
      model <<- train(Ps ~ ., method= method, data=trainN, tuneLength=tunelength )##train
    }
  }
  else if (bacteria=="LAB") {
    if(nrow(data)>=50){
      model <<- train(LAB ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
    }else{
      model <<- train(LAB ~ ., method= method, data=trainN, tuneLength=tunelength )##train
    }
  }
  else if (bacteria=="Bth") {
    if(nrow(data)>=50){
      model <<- train(Bth ~ ., method= method, data=trainN, trControl= control, tuneLength=tunelength )##train
    }else{
      model <<- train(Bth ~ ., method= method, data=trainN, tuneLength=tunelength )##train
    }
  }
  
  
  predicted <- predict(model, testN)##Prediction
  
  
  if (bacteria=="TVC") {
    RMSE <- RMSE(predicted, testN$TVC)##RMSE calculation
  }else if (bacteria=="Ps") {
    RMSE <- RMSE(predicted, testN$Ps)##RMSE calculation
  }else if (bacteria=="LAB") {
    RMSE <- RMSE(predicted, testN$LAB)##RMSE calculation
  }else if (bacteria=="Bth") {
    RMSE <- RMSE(predicted, testN$Bth)##RMSE calculation
  }
  ## model Accuracy +- 1 log count
  if (bacteria=="TVC") {
    Accuracy<- Percentage(predicted,testN$TVC)
  }else if (bacteria=="Ps") {
    Accuracy<- Percentage(predicted,testN$Ps)
  }else if (bacteria=="LAB") {
    Accuracy<- Percentage(predicted,testN$LAB)
  }else if (bacteria=="Bth") {
    Accuracy<- Percentage(predicted,testN$Bth)
  }
  ##combining both predicted and actual bacterial counts
  df<-cbind(testN[,bacteria], predicted)
  df<-as.data.frame(df)
  colnames(df)=c("Actual", "Predicted")
  dm<<-df
  ## plot predicted values vs Actual/Observed values
  if (bacteria=="TVC") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "TVC", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="Ps") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Ps", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="LAB") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "LAB", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  else if (bacteria=="Bth") {
    sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "Bth", "- RMSE:",round(RMSE,digits = 2)," Accuracy :",round(Accuracy, digits = 2),"%"),
                  xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                  xlim= c(0,11), ylim=c(0,11)
    )
    sp<<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 11, color="blue")+ theme_bw() +theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect()
    ) + theme(axis.line = element_line(color = 'black'))
    
  }
  
}

######replicate checking function
rep_check<-function(data, replicates){
  n<-as.numeric(replicates)
  df_ls <- lapply(seq(1, nrow(data), by=n), function(i) 
    data[i: pmin((i+(n-1)), nrow(data)),])
  
  da<-c()
  ##converting list of replicates as separet df
  for (i in 1:length(df_ls)) {
    dk<-(df_ls[i])
    dk<-as.data.frame(dk)
    dk<-list(dk)
    da<-c(da, dk)
    
  }
  
  ##matching sampleID of all replicates of one sample.
  db<-c()
  for (i in 1:length(da)) {
    dc<-as.data.frame(da[i])
    #dc<-as.data.frame(dc)
    R<-row.names(dc)
    m<-adist(R)
    q<-sum(m)
    db<-c(db, q)
    
  }
  num<-(sum(db)/length(da))##sum all the matches and divide by total samples
  if(n==3){##when replicate are three
    if(num==2*n){##checking condition to know all replicates grouped proprely for each sample.
      resl<-"yes"
    }else{
      resl<-"NO"
    }
  }else if(n==6){##when replicate are six
    if(num==7*n){##checking condition to know all replicates grouped proprely for each sample.
      resl<-"yes"
    }else{
      resl<-"NO"
    }
  }else{
    resl<-"null"
  }
}
#############################################################################################
##Function that rounds numeric wavelength for FTIR Data######################################

roundWavelengths <- function(df){
  #
  # function to round column names/wavelengs 
  #
  # df = data frame where column names are numeric
  #
  # Hint - function works just for numeric column names
  #      - CreateObj function is required
  #
  roundedCol_df <- round(as.numeric(colnames(df))) #round column names/wavelengths
  colnames(df) <- roundedCol_df #set rounded wavelengths as column names
  df_orgin <- df 
  data <- df
  d <- c()
  for (i in 1:(ncol(df_orgin)-1)){
    if(as.numeric(colnames(df_orgin[i])) == as.numeric(colnames(df_orgin[i+1]))){
      for(j in 1:nrow(df_orgin)){ #calculate means based on columns of which wavelengths are the same after rounding
        d[j] <- mean(as.numeric(df_orgin[j,i]), as.numeric(df_orgin[j,i+1]))
      }
      d <- matrix(d, ncol = 1)
      d <- as.data.frame(d)
      rownames(d) <- rownames(df_orgin)
      colnames(d) <- colnames(df_orgin[i])
      for(z in 1:(length(colnames(data))-1)){ #create a dataframe where wavelengths are not repeated
        if(colnames(data[z]) == colnames(df_orgin[i])){
          df <- CreateObj(data[,1:(z-1)], d)
          df <- CreateObj(df, data[,(z+2):length(data)])
          data <- df
        }
      }
      d <- c()
      i = i+1
    }
  }
  dataFrame <- df
  
}
#############################################################################################################
#####Functin merges two datasets according rownames##########################################################
CreateObj <- function(data1, data2){
  #
  # Data pretreatment function to create combined data sets
  #
  # data1 = the first data set to combining
  # data2 = the second data set to combining
  #
  # Hint - each data set has to have the same number of rows
  #
  #combine all rows from the both dataset 
  merged <- merge(data1, data2, by = 'row.names')
  rownames(merged) = merged[,1]
  # remove the row names column, which was added th the merged dataset during mergining (additional one)
  as.data.frame(merged[,-1])
}

##Function to extract batches detaile for ML training  
##Batches<-"R1,R2,R3;R4,R5"

Batch_pram<-function(pram){
  prm<-pram
  p<-strsplit(prm,split=';', fixed=TRUE)
  ####
  z<-strsplit(p[[1]], " ")
  nam<-c("Batchtrain", "Batchtest")
  names(z)<-nam
  ####
  if (grepl(",",z$Batchtrain[1] , fixed = TRUE)==TRUE){
    z$Batchtrain<-strsplit(z$Batchtrain,split=',', fixed=TRUE)
    mtr<-c(z$Batchtrain[[1]])
  }else{
    mtr<-c(z$Batchtrain)
  }
  ##
  if (grepl(",",z$Batchtest[1] , fixed = TRUE)==TRUE){
    z$Batchtest<-strsplit(z$Batchtest,split=',', fixed=TRUE)
    mts<-c(z$Batchtest[[1]])
  }else{
    mts<-c(z$Batchtest)
  }
  ####
  Batchtt <- list(mtr,mts)
}



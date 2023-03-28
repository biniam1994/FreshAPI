library(MLmetrics)
require(caret)
require(randomForest)
library(caret)
library(gmodels)
library(ggpubr)
require(pROC)
require(HH)
library(dplyr)
####regression function outputing plots actual vs predicted, accuracy over iteration

reg.run<-function(data, method, machine){
  set.seed(124)
  RMSEC<-c()
  for (i in 1:100) {
    trainN <- data[train.index[,i],]
    testN <- data[-train.index[,i],]
    if(method=="rf"){
      set.seed(12)
      model <<- train(TVC ~ ., method='rf', data=trainN, trControl= control, tuneGrid = expand.grid(mtry=1:20), ntree=150)
    }else if(method=="knn"){
      set.seed(13)
      model <<- train(TVC ~ ., method='knn', data=trainN, trControl= control, tuneGrid = expand.grid(k=1:20))
    }else if(method=="svmLinear"){
      set.seed(14)
      model <<- train(TVC ~ ., method='svmLinear', data=trainN, trControl= control, tuneGrid = expand.grid(C=c(0.1,0.3,0.01,0.5,0.6,0.7,0.8,0.9,1,2,3)))
    }else if(method=="svmRadial"){
      set.seed(15)
      model <<- train(TVC ~ ., method='svmRadial', data=trainN, trControl= control, tuneLength=10)
    }else if(method=="svmPoly"){
    
        model <<- train(TVC ~ ., method='svmPoly', data=trainN, trControl= control, tuneLength=10)
     
    }else if(method=="lm"){
      set.seed(16)
      model <<- train(TVC ~ ., method='lm', data=trainN, trControl= control, tuneLength=10)
    }
   
    predicted <- predict(model, testN)##Prediction
    RMSE <- RMSE(testN$TVC, predicted)##RMSE calculation
    RMSEC<-c(RMSEC, RMSE)
  }
  model.name<-paste0(method,"_R_",machine,".rds")
  saveRDS(model, model.name)
  Mean.RMSE<-mean(RMSEC)##calculating mean
  SD.RMSE<-sd(RMSEC)##calculating SD
  Cl <-gmodels::ci(RMSEC, confidence = 0.95)
  iterations <- as.array(c(1:100))
  ##
  ##combining both predicted and actual bacterial counts
  df<-cbind(testN[,"TVC"], predicted)
  df<-as.data.frame(df)
  colnames(df)=c("Actual", "Predicted")
  dm<-df
  ##combinig iteration and RMSE
  df1<-cbind(iterations, RMSEC)
  df1<-as.data.frame(df1)
  colnames(df1)=c("Iteration", "RMSE")
  
  ##plot of RMSE over iteration
  st<-ggscatter(df1,x="Iteration", y="RMSE", ylim=c(0,4), xlim=c(0,100), xlab="Iteration", ylab="RMSE", main=paste(method," Model Mean RMSE:",
                                                                                                                  round(Mean.RMSE, digits = 3), "CI 95%:", 
                                                                                                                  round(Cl[2], digits = 3), "-", 
                                                                                                                  round(Cl[3], digits = 3), "SD:",
                                                                                                                  round(SD.RMSE,digits = 3)))##Plotting RMSE
  st<-st + geom_line( color="black") + geom_point(shape=1, color="black", fill="#69b3a2", size=1)
  plotname<-paste0(method,"_R_",machine,"_","T.png")
  ggsave(plotname, st,width = 9, height = 6, dpi = 300, units = "in", device='png')
  ##Linear model Accuracy for TVC
  Accuracy<-Percentage(predicted, testN$TVC)
  # LM plot predicted values vs Actual/Observed values
  sp<-ggscatter(df, x="Actual", y="Predicted", main=paste(method, "Model for", "TVC", "- RMSE:",round(RMSE,digits = 4)," Accuracy :",round(Accuracy, digits = 2),"%"),
                xlab="Observed microbial Population (log CFU/g)", ylab="Estimated microbial Population (log CFU/g)",
                xlim= c(0,10), ylim=c(0,10)
  )
  sp<-sp + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = 1, slope = 1, color="blue") + geom_abline(intercept = -1, slope = 1, color="blue")+ stat_cor( method = "pearson",label.x = 0,label.y = 8, color="blue")+ theme_bw() +theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect()
  ) + theme(axis.line = element_line(color = 'black'))
  plotname1<-paste0(method,"_R_",machine,"_","AP.png")
  ggsave(plotname1,sp,width = 9, height = 6, dpi = 300, units = "in", device='png')
  
}

##
#classification function outputting plots assesing accuracy over 100 iterations
#ROC plot
#Input is data and model descriptor
#eg. "svmLinear", "rf" for random forest or "knn" for k nearest neighbouring.
#machine is analytical platform or method
#Example use classification(data,"rf","HPLC")
classification<- function (data, method, machine) {
  accuracy<<- c()
  misclassifications <- c()
  
  for (i in 1:100){
    
    trainSet <-  data[trainIndex[,i],]
    testSet <- data[-trainIndex[,i],]
    
    
    trainCl <- trainSet[,ncol(trainSet)]
    trainCl <- as.factor(trainCl)
    testCl <- testSet[,ncol(trainSet)]
    testCl <<- as.factor(testCl)
    
    trainSet <- trainSet[,1:(ncol(trainSet)-1)]
    testSet <- testSet[,1:(ncol(testSet)-1)]
    
    #arrange resampling with replacement for imbalanced classes
    ctrl <- trainControl(method = "boot", number = 20,
                         classProbs = F,
                         summaryFunction = multiClassSummary,
                         ## oversampling/undersampling/combined:
                         sampling = "up") 
    
    set.seed(123)
    if(method=="rf"){
      fit<- train(trainSet,as.factor(trainCl), 
                  method = method,
                  metric= "Accuracy",
                  trControl = ctrl, tuneGrid = expand.grid(mtry=1:20), ntree=100)
    }else if(method=="svmLinear"){
      if(machine=="GCMS"){
        set.seed(123)
        fit<- train(trainSet,as.factor(trainCl), 
                    method = method,
                    metric= "Accuracy",
                    trControl = ctrl, tuneGrid = expand.grid(C=c(0.1,0.3,0.01,0.5,0.6,0.7,0.8,0.9,1,2,3))) 
      }else{
      set.seed(123)
      fit<- train(trainSet,as.factor(trainCl), 
                  method = method,
                  metric= "Accuracy",
                  trControl = ctrl, tuneLength = 10)
      }#,tuneGrid = expand.grid(C=c(0.1,0.3,0.01,0.5,0.6,0.7,0.8,0.9,1,2,3)))
    }else if(method=="knn"){
      set.seed(123)
      fit<- train(trainSet,as.factor(trainCl), 
                  method = method,
                  metric= "Accuracy",
                  trControl = ctrl, tuneGrid = expand.grid(k=c(1:10)))
    }else if(method=="svmPoly"){
      set.seed(123)
      fit<- train(trainSet,as.factor(trainCl), 
                  method = method,
                  metric= "Accuracy",
                  trControl = ctrl, tuneLength = 10)
    }else if(method=="svmRadial"){
      set.seed(123)
      fit<- train(trainSet,as.factor(trainCl), 
                  method = method,
                  metric= "Accuracy",
                  trControl = ctrl, tuneLength = 10)
    }else if(method=="nb"){
      set.seed(123)
      grid <- data.frame(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0,0.5,1.0))
      fit<- train(trainSet,as.factor(trainCl), 
                  method = method,
                  metric= "Accuracy",
                  trControl = ctrl, tuneGrid = grid)
    }else{
      print("Error Occured: Please select method: knn | rf | svmLinear | svmRadial | svmPoly")
    }
    predicts<<- predict(fit,testSet)
    confus<- confusionMatrix(as.factor(predicts), as.factor(testCl), positive = "3")
    accuracy<<- c(confus$overall[1],accuracy)
    accuracy.mean<<-round((mean(accuracy)), digits = 2)
    keys <- unique(testCl)
    for(key in keys) {
      if (i == 1) { misclassifications[[key]] <- 0 }
      falseNegatives <- which((testCl == key) & (testCl != predicts))
      misclassifications[[key]] <- misclassifications[[key]] + length(falseNegatives)
    }
  }
  ##save model
  model.name<-paste0(method,"_C_",machine,".rds")
  saveRDS(fit, model.name)
  ##calculate SD and Ci 95
  
  SD<-sd(accuracy)##calculating SD
  CI95<-gmodels::ci(accuracy, confidence = 0.95)
  x11()
  plot(accuracy, type="l", ylab = "Accuracy", xlab = "Iterations",ylim = c(0,1),
       main = paste0("Accuracy over 100 iterations : \n Accuracy.mean = ",
                     round(accuracy.mean, digits = 3),"; SD ",
                     round(SD,digits = 3), " ; 95% Confidence Interval: ",
                     round(CI95[2], digits = 3), " - ", 
                     round(CI95[3], digits = 3), "\n "))##Plotting Accuracy))
  #plotname1<-paste0(method,"_C_",machine,"_","Acc.png")
  #ggsave(plotname1,Accu,width = 9, height = 6, dpi = 300, units = "in", device='png')

  ##
  miscl1<-as.data.frame(setNames(list(unname(misclassifications["1"])),'Fresh'))
  miscl2<-as.data.frame(setNames(list(unname(misclassifications["2"])),'SemiFresh'))
  miscl3<-as.data.frame(setNames(list(unname(misclassifications["3"])),'Spoiled'))
  ##
  names(miscl1)[1]<-paste("Fresh")
  names(miscl2)[1]<-paste("SemiFresh")
  names(miscl3)[1]<-paste("Spoiled")
  ##
  misclass<-as.matrix(cbind(miscl1,miscl2, miscl3))
  misclass1<-misclass
  D<-sum(misclass1)
  misclass1[1]<-round((misclass1[1]/D*100), digits=0)
  misclass1[2]<-round((misclass1[2]/D*100), digits=0)
  misclass1[3]<-round((misclass1[3]/D*100), digits=0)
  x11()
  
  barplot(misclass1, ylim=c(1,100),ylab="Precentage of misclassifications over 100 iterations",xlab=" ", 
          main=paste(method," model",": misclassified samples per class"), cex.names=1,
          legend = c(paste0("Fresh : ",misclass1[1], "%"), paste0("SemiFresh : ",misclass1[2], "%"), paste0("Spoiled : ",misclass1[3], "%")))
  #plotname2<-paste0(method,"_C_",machine,"_","bar.png")
  #ggsave(plotname2,bar,width = 9, height = 6, dpi = 300, units = "in", device='png')
  
  ####### ROC curve  model  #######
  Roc <<- pROC::multiclass.roc(as.numeric(predicts),as.numeric(testCl), levels=c(1,2,3))
  Roc_e <- Roc$rocs
  A<<-Roc$auc
  B<<-round(as.numeric(A), digits = 3)
  #Aua<-as.character(auc(Roc))
  x11()
  plot.roc(main=paste0(method, "  ROC Plot for ",machine, " Sensory Scores; Multi-class Auc : ",B, " Accuacy: ", round(accuracy[100], digits = 2)),(Roc_e[[1]]))
  #sapply(1:length(Roc_e),function(i) lines.roc(Roc_e[[i]],col=i))
  lines(Roc_e[[2]],col="red")
  lines(Roc_e[[3]],col="green")
  legend("bottomright",y=NULL, legend = c("Class 1-2","Class 1-3","Class 2-3"), c("black","red","green"))
  #plotname3<-paste0(method,"_C_",machine,"_","roc.png")
  #ggsave(plotname3,ro,width = 9, height = 6, dpi = 300, units = "in", device='png')
  return(
    list(
      confusionMatrix = confus,
      misclassifications = misclassifications,
      model = fit
    )
  )
  
  
}

##
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
    dk<-as.data.frame(dk)
    dk[] <- lapply(dk, na.aggregate)
    dk<-list(dk)
    da<-c(da, dk)
    
  }
  
  ##combinig back all the subset(replicate) list of dataframe to one dataframe
  db<-rbind()
  ##imputing NA's with mean of replicate
  for (i in 1:length(da)) {
    dc<-(da[i])
    dc<-as.data.frame(dc)
    db<-rbind(db, dc)
    
  }
  data<-db
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

##pca plot function
pca12<-function(data, machine){
  data$class <- ifelse(grepl("1", data$Odor, ignore.case = T), "Fresh", 
                       ifelse(grepl("2", data$Odor, ignore.case = T), "SemiFresh",
                              ifelse(grepl("3",data$Odor, ignore.case = T), "Spoiled","")))##check For batch identity from rowname
  
  data <- subset(data, select = -c(Odor))
  data1 <- subset(data, select = -c(class))
  pca <- prcomp(data1,center=T, scale=T)
  fviz_pca_ind(pca,
               axes = c(1,2),
               geom  = c("point"),
               col.ind  = data$class,
               palette=c("green","blue", "red"),
               addEllipses = TRUE,
               ellipse.level = 0.95,
               legend.title = "class",
               title = paste0("PCA PLOTS of ", "machine")
  )
}
##
pca13<-function(data, machine){
  data$class <- ifelse(grepl("1", data$Odor, ignore.case = T), "Fresh", 
                       ifelse(grepl("2", data$Odor, ignore.case = T), "SemiFresh",
                              ifelse(grepl("3",data$Odor, ignore.case = T), "Spoiled","")))##check For batch identity from rowname
  
data <- subset(data, select = -c(Odor))
data1 <- subset(data, select = -c(class))
pca <- prcomp(data1,center=T, scale=T)
fviz_pca_ind(pca,
             axes = c(1,3),
             geom  = c("point"),
             col.ind  = data$class,
             palette=c("green", "blue", "red"),
             addEllipses = TRUE,
             ellipse.level = 0.95,
             legend.title = "class",
             title = paste0("PCA PLOTS of ", "machine")
)

}
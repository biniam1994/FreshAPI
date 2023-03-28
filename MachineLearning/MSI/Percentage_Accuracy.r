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
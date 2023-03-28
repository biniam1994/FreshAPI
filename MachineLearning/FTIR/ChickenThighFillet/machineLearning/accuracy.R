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

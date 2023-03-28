CreateObj <- function(data1, data2){
  #
  # Data pretreatment function to create combined data sets
  # data1 = the first data set to combining
  # data2 = the second data set to combining
  # Hint - each data set has to have the same number of rows
  #
  #combine all rows from the both dataset 
  merged <- merge(data1, data2, by = 'row.names')
  rownames(merged) = merged[,1]
  # remove the row names column, which was added th the merged dataset during mergining (additional one)
  as.data.frame(merged[,-1])
}
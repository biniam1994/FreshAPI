#* @apiTitle Food Spoilage API
#* @apiDescription Predicting Bacterial counts or Sebsory score to asses food spoilage level and Automating regression Machine Learning 


source("functions.R")


##########Prediction of Bacterial counts using existing model
#* parse JSON
#* @post /Predict
function(req, platform, product, bacteria){
  #Hint:-
  #   req = analytical platform machine request (a single row entry)
  #   platform = FTIR | MSI | eNose | HPLC | GCMS
  #   product = CTF | CBT | CBG
  #   bacteria = TVC | Ps | Bth | LAB | null 
  #   ML = classification or regression
  data <- tryCatch(jsonlite::fromJSON(req$postBody),
                   error = function(e) NULL)
  if (is.null(data)) {
    res$status <- 400
    return(list(error = "No data provided"))
  }
  #data <- jsonlite::fromJSON(req$postBody)
  if(bacteria=="TVC" | bacteria=="Bth"  | bacteria=="Ps"  | bacteria=="LAB"){
    if(product=="CBG" | product=="CTF" | product=="CBF"){
    if(platform=="FTIR"){
    # round if required for FT-IR wavelength
      row_colnames<- colnames(data)
      row_col_integer <- 0
      for(i in 1:length(row_colnames)){
        if(as.numeric(row_colnames[i])==round(as.numeric(row_colnames[i]), 0)){
          row_col_integer <- row_col_integer+1
      }
      else{
        row_col_integer <- row_col_integer
      }
    }
    if(length(row_colnames) == row_col_integer){
      data<-data  
    } else{
      data <- data
      data <- roundWavelengths(data)
    }
  }else if(platform=="HPLC" | platform=="eNose" | platform=="GCMS" | platform=="MSI"){
    data<-data
  }else{
    print("Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI ")
  }
  #######
  if(platform=="HPLC" | platform=="eNose" | platform=="GCMS" | platform=="MSI" | platform=="FTIR"){ 
  if(nrow(data) == 1){
    prediction <- predict_B(data,platform, product, bacteria)
    return(paste0("predicted ", bacteria," counts: ", round(prediction,digits = 3)))
  }
  else if(nrow(data) > 1){
    for (i in 1:length(data)) {
      names<-row.names(data)
      prediction <- predict_B(data,platform, product, bacteria)
      return(paste0( names, ";",round(prediction,digits = 3)))
    } 
  }
  }else{
    print("Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI ")
  }
      
}else{
    print("Error occured: Please make sure product type is correct (CBG or CFT or CBF)")
  }
  
}else{
    print("Error occured: Please make sure bacteria type is correct (TVC or Ps or LAB or Bth )")
  }
  
}
##########Classify meat samples into sensory score(Fresh, Semi-Fresh, Spoiled)
#* parse JSON
#* @post /classify
function(req, platform, product){
  #Hint:-
  #   req = analytical platform machine request (a single row entry)
  #   platform = FTIR | MSI | eNose | HPLC | GCMS
  #   product = CTF | CBT | CBG
  #   
  
  data <- tryCatch(jsonlite::fromJSON(req$postBody),
                   error = function(e) NULL)
  if (is.null(data)) {
    res$status <- 400
    return(list(error = "No data provided"))
  }
  #data <- jsonlite::fromJSON(req$postBody)
  if(product=="CBG" | product=="CTF" | product=="CBF"){
    if(platform=="FTIR"){
      # round if required for FT-IR wavelength
      row_colnames<- colnames(data)
      row_col_integer <- 0
      for(i in 1:length(row_colnames)){
        if(as.numeric(row_colnames[i])==round(as.numeric(row_colnames[i]), 0)){
          row_col_integer <- row_col_integer+1
        }
        else{
          row_col_integer <- row_col_integer
        }
      }
      if(length(row_colnames) == row_col_integer){
        data<-data  
      } else{
        data <- data
        data <- roundWavelengths(data)
      }
    }else if(platform=="HPLC" | platform=="eNose" | platform=="GCMS" | platform=="MSI"){
      data<-data
    }
  #######
  if(platform=="HPLC" | platform=="eNose" | platform=="GCMS" | platform=="MSI" | platform=="FTIR"){  
  if(nrow(data) == 1){
    prediction <- predict_C(data,platform, product)
    return(paste0("classified as: ", prediction))
  }
  else if(nrow(data) > 1){
    for (i in 1:length(data)) {
      names<-row.names(data)
      prediction <- predict_C(data,platform, product)
      return(paste0( names, "; ",prediction))
    } 
   }
  }else{
    print("Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI ")
  }
  }else{
    print("Error occured: Please make sure product type is correct (CBG or CFT or CBF)")
  }
}
#########################Endpoint used to check missing values(NA) in dataset#############################
#* parse XLSX or CSV
#* @post /row_index_NA
function(filename, path){
  #Hint:-
  #   filename = is a filename of the data eg.(enose_breastAUA.xlsx or FTIR_data.csv)
  #   path = is the folder name without the root directory eg.()
  #   
  ## read xlsx or csv files
  
  wd <- getwd()
  f <- paste0(wd,"/", path, "/", filename)
  if (see_if(has_extension(filename, 'csv')) == TRUE){
    data<- utils::read.csv(f, header = T, sep= ",", row.names = 1, check.names = F)
  }
  else if (see_if(has_extension(filename, 'xlsx')) ==TRUE) {
    data<- openxlsx::read.xlsx(f, rowNames =T, check.names = F,sep.names=" ")
  }
  
  else{
    print("Error occured: Please check your file it should contain .xlsx or .csv ")
  }
  M <- 0
  K <- 0
  # Get rows with NA and number of NA on the rows
    S<-rowSums(is.na(data)); S[S>0]##Number of NA's sum by row in data
    M<-S[S>0]
    K<-which(rowSums(is.na(data))>0)##get row index with NA
  if(length(K)==0){
    print("There is no missing value(NA) found")
  }else{
    return(paste0(" Indexs of Row with NA:  ", K, "  Sum of NA's in rows of data:  ", M))
  } 
} 

###############EndPoint used to treat missing valuses(NA)################################################
#* @serializer csv
#* @get /rm-im-NA
function(platform, filename,Rm, path, replicates){
  #
  # Hint:-
  # replicates = number of replicate that you have in the data eg.(1 | 2 | 3 | .....) 
  # platform = is the analytical platform used for data generation eg.(FTIR, eNose, MSI, HPLC, GCMS)
  # filename = is a filename of the data eg.(enose_breastAUA.xlsx or FTIR_data.csv)
  # path = is the folder name without the root directory eg.(Demo)
  # Rm= is wheather to remove them or impute eg.(yes or no)
  # 
  #
  n<-as.numeric(replicates)
  fname<-filename
  wd <- getwd()
  f <- paste0(wd,"/", path, "/", fname)
  
  if(platform == "FTIR" | platform =="MSI" | platform =="HPLC"  | platform =="GCMS"  | platform =="eNose"){
    if(n==1){##if data doesn't have replicate NA will be imputed using whole samples column mean 
      if (see_if(has_extension(fname, 'csv')) == TRUE){
        data<- utils::read.csv(f, header = TRUE, sep = ",", check.names = F)
      }
      else if (see_if(has_extension(fname, 'xlsx')) ==TRUE) {
        data<- openxlsx::read.xlsx(f, check.names = F,sep.names=" ")
        
      }else{
        print("Error Occured: Please check your file contains .csv or .xlsx extension")
      }
      data <- NA.rm.im(data,Rm)
      data <- as.data.frame(data,check.names = F)
      row.names(data)<-data[,1]
      data<-data[,-1]
    }else if(n>1){## if data contains replicates then NA will be imputed using replicate column mean
      if (see_if(has_extension(fname, 'csv')) == TRUE){
        data<- utils::read.csv(f, header = TRUE, sep = ",", row.names = 1,check.names = F)
      }
      else if (see_if(has_extension(fname, 'xlsx')) ==TRUE) {
        data<- openxlsx::read.xlsx(f, rowName=T, check.names = F,sep.names=" ")
        
      }else{
        print("Error Occured: Please check your file contains .csv or .xlsx extension")
      }
      data<-NA_tre_re(data,n)##there is no optiion  to omit NA here
    }else{
      print("Error Occured: Please check your file contains appropraite replicate number")
    }
    df <- data.frame(CHAR = letters, NUM = rnorm(length(letters)), stringsAsFactors = F)
    csv_file <- tempfile(fileext = ".csv")
    on.exit(unlink(csv_file), add = TRUE)
    write.csv(data, file = csv_file, row.names = T)
    readBin(csv_file,'raw',n = file.info(csv_file)$size)
    readLines(csv_file)
  }else{
    print("Error Occured: please make sure you choose correct platform: FTIR | MSI | HPLC  | GCMS  | eNose")
  }
  
 
} 

###############EndPoint used to treat replicate #######################################################################

#* @serializer csv
#* @get /replicates_tre

function(products, replicates, platform, filename, path){
  #
  # Function to give an user repond for REST API request 
  # Hint:-
  #   products = example can be chicken Burger then you pass Burger only as product
  #   platform = platform = is the analytical platform used for data generation eg.(FTIR, eNose, MSI, HPLC, GCMS)
  #   replicates = number of replicate that you have in the data eg.(1 | 2 | 3 | .....)
  #   filename  = is a filename of the data eg.(enose_breastAUA.xlsx or FTIR_data.csv)
  #   path  = is the folder name without the root directory eg.(Demo)
 
  fname<-filename
  wd <- getwd()
  f <- paste0(wd,"/", path, "/", fname)
  products<- products
  replicates<- as.numeric(replicates)
  
  if(platform=="FTIR"){
    if (see_if(has_extension(fname, 'xlsx')) == TRUE){
    data<- openxlsx::read.xlsx(f, check.names = F,sep.names=" ")
    }
    else if (see_if(has_extension(fname, 'csv')) ==TRUE) {
    
      data<- utils::read.csv(f, header = T, sep = ",", check.names = F)
    }
    else{
     print("Error ocurred : Filename doesn't have .csv or .xlsx extention")
    }
    data <- replicate_tre(data, products, replicates)##this function will aggrigate and rename rownames
    data <- as.data.frame(data, check.names = F)
  }
  
  
  else if(platform=="HPLC" | platform=="GCMS"  | platform=="eNose"){
    if (see_if(has_extension(fname, 'xlsx')) == TRUE){
      data<- openxlsx::read.xlsx(f, rowNames=T, check.names = F,sep.names=" ")
    }
    else if (see_if(has_extension(fname, 'csv')) ==TRUE) {
      
      data<- utils::read.csv(f, header = T, sep = ",", row.names = 1, check.names = F)
    }
    else{
      print("Error ocurred : Filename doesn't have .csv or .xlsx extention")
    }
    
    data <- repli_tre(data, replicates)## will aggrigate only no rename of rowsnames for enose,GCMS, enose
    data <- as.data.frame(data, check.names = F)
  }
  

  df <- data.frame(CHAR = letters, NUM = rnorm(length(letters)), stringsAsFactors = F)
  csv_file <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_file), add = TRUE)
  write.csv(data, file = csv_file, row.names = TRUE)
  readBin(csv_file,'raw',n = file.info(csv_file)$size)
  readLines(csv_file)
}
###############EndPoint used to plot PCA for based on batches only #######################################################################


#* @get /PC
#* @serializer contentType list(type='image/png')
PC = function(PC, platform, filename, path){
  #Hint:-
  #   PC = is the component you want to plot eg.(PC one vs two or Pc one vs three)
  #   platform = is the analytical platform used for data generation eg.(FTIR, eNose, MSI, HPLC, GCMS)
  #   filename = is a filename of the data eg.(enose_breastAUA.xlsx or FTIR_data.csv)
  #   path = is the folder name without the root directory eg.(Demo)
  
  fname<-filename
  wd <- getwd()
  f <- paste0(wd,"/", path, "/", fname)
  if (PC=="PC12") {
    if (platform=="FTIR" | platform=="MSI" | platform=="HPLC" | platform=="GCMS" | platform=="eNose") {
      if (see_if(has_extension(fname, 'csv')) == TRUE){
        data<- utils::read.csv(f,header = T, sep = ",", row.names = 1, check.names = F)
      }
      else if (see_if(has_extension(fname, 'xlsx')) ==TRUE) {
        data<- openxlsx::read.xlsx(f, rowName=T, check.names = F,sep.names=" ")
      }
      else{
        print("Filename doesn't have .csv or .xlsx extention")
      }
      plot<-pca1.function(data)
      file <- tempfile(fileext = ".png")
      on.exit(unlink(file), add = TRUE)
      ggsave(file,plot,width = 9, height = 6, dpi = 300, units = "in", device='png')
      readBin(file,'raw',n = file.info(file)$size)
    }else{
      print("Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI ")
    }
    
  }else if (PC=="PC13") {
    if (platform=="FTIR" | platform=="MSI" | platform=="HPLC" | platform=="GCMS" | platform=="eNose") {
      if (see_if(has_extension(fname, 'csv')) == TRUE){
        data<- utils::read.csv(f, header = T, sep = ",", row.names = 1, check.names = F)
      }
      else if (see_if(has_extension(fname, 'xlsx')) ==TRUE) {
        data<- openxlsx::read.xlsx(f, rowName=T, check.names = F,sep.names=" ")
      }
      else{
        print("Error Occured: Filename doesn't have .csv or .xlsx extention")
      }
      plot<-pca2.function(data)
      file <- tempfile(fileext = ".png")
      on.exit(unlink(file), add = TRUE)
      ggsave(file,plot,width = 9, height = 6, dpi = 300, units = "in", device='png')
      readBin(file,'raw',n = file.info(file)$size)
    }else{
      print("Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI ")
    }
    
  }else{
    print("Error occured: Please check PC argument has to be PC13 or PC12")
  }
}

###############EndPoint used to train machine learning model#####################################################

#* @serializer contentType list(type="application/octet-stream")
#* @get /MLearning
rest_rds = function(platform, bacteria, Batches, method, iteration, filename, path, pram, tunelength){
  #Hint:-
  #   platform = is the analytical platform used for data generation eg.(FTIR, eNose, MSI, HPLC, GCMS)
  #   bacteria = is bacterial type eg. (TVC, Ps, Bth, LAB)
  #   Batches = train in mixed batch or one batch and test with another batch
  #   method = is machine learning algorith (rf, svmLinear, knn, svmPoly, svmRadial)
  #   iteration = number of iteration to train(1,....100)
  #   filename = is a filename of the data eg.(enose_breastAUA.xlsx or FTIR_data.csv)
  #   path = is the folder name without the root directory eg.(Demo)
  #   pram = tuning hyperparameters eg.(k=1:20 or C=0.1,1,3,2 or ntree=100, mtry=1:10, null....etc)
  #   tunelength = to what lenght to tune using randomserch hyperparameters
  #
  Batches<-Batches
  iteration<-as.numeric(iteration)
  tunelength<-as.numeric(tunelength)
  
  
   if(bacteria=="TVC" | bacteria=="Ps" | bacteria=="Bth" | bacteria=="LAB"){
      if(method=="rf" | method=="knn" | method=="svmLinear" | method=="svmPoly" | method=="svmRadial"){
        if(iteration>=1){
          if(tunelength>=0){
  
  #pram<-"mtry-1:10;ntree-300"
  pram<-pram
  ####pram handling
  if(tunelength == 0){
    if(method=="rf"){
      pra<-rf_pram(pram)
      nam<-c("mtry", "ntree")
      names(pra)<-nam
    }else if(method=="knn"){
      pra<-knn_pram(pram)
      namk<-c("k")
      names(pra)<-namk
    }else if(method=="svmLinear"){
      pra<-svl_pram(pram)
      namk<-c("C")
      names(pra)<-namk
    }else if(method=="svmPoly"){
      pra<-svp_pram(pram)
      namk<-c("C","S","D")
      names(pra)<-namk
    }else if(method=="svmRadial"){
    pra<-svr_pram(pram)
    namk<-c("C","S")
    names(pra)<-namk
    }
  }
  
  if(Batches =="mixed"){
   Batches<-Batches
  }else{
    Batch<-Batch_pram(Batches)
    namB<-c("Batchtrain", "Batchtest")
    names(Batch)<-namB
  }
  
  ################################
  fname<-filename
  wd <- getwd()
  f <- paste0(wd,"/", path, "/", fname)
  folder <- paste0(wd,"/", path)
  if (platform=="FTIR" | platform=="MSI" | platform=="HPLC" | platform=="GCMS" | platform=="eNose"){
    if (see_if(has_extension(fname, 'csv')) == TRUE){
      data<- utils ::read.csv(f, header = TRUE, sep = ",",  row.names=1, check.names = F)
    }
    else if (see_if(has_extension(fname, 'xlsx')) ==TRUE) {
      data<- openxlsx::read.xlsx(f, rowNames = T, check.names = F,sep.names=" ")
    }else{
      print("Error Occured: Filename doesn't have .csv or .xlsx extention")
    }
    
    data<- destroyX(data)
    data<-Bacteria_sel(data, bacteria, platform)
    if(Batches=="mixed"){
      if (tunelength != 0){
        regression.run(data,platform, method, bacteria, iteration,tunelength)
      }else{
        regression.run.pra(data,platform, method, bacteria, iteration,pra)
      }
    }else if(Batches !="mixed"){
      if(tunelength != 0){
        regression.batchA(data,platform, method, Batch, bacteria, tunelength)
      }else{
        regression.batchA.pra(data,platform, method,Batch, bacteria,pra)
      }
    }
    
    tfile = tempfile(fileext = ".rds")
    on.exit(unlink(tfile), add = TRUE)
    if(Batches =="mixed"){
      ObjSave(st, sp, dm, model, folder = folder)
    }else{
      ObjSave(sp, dm, model, folder = folder)
    }
    readBin(tfile, "raw", n = file.info(tfile)$size)
    
  }else {
    print("Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI ")
  }
          }else{
            print("Error occured: Please make sure tunelength is zero or greater")
          }
        }else{
          print("Error occured: Please make sure iteration is 1 to 100 values only")
        }
      }else{
        print("Error occured: Please choose corret machine learning algorithm: svmRadial | svmPoly | svmLinear | knn | rf")
      }
   }else{
     print("Error occured: Please make sure bacteria type is correct (TVC or Ps or LAB or Bth )")
   }
}

##rendering model plots from saved rds object
#* @get /plots
#* @serializer contentType list(type='image/png')
plots = function(model_plot, path, filename){
  #Hint:-
  #   model_plot = is to specifiy whether you are ploting actual(Pre_Act) vs predicted or RMSE overiteration(Iter)
  #   path = is the folder name without the root directory eg.(Demo)
  #   filename = is a filename of the plot generated from MLearning endpoint eg.(sp.rds or st.rds)
  if(filename=="sp.rds" | filename=="st.rds"){
  if (model_plot=="Pre_Act") {
    wd <- getwd()
    f <- paste0(wd,"/", path, "/", filename)
    sp<- base::readRDS(f)
    file = 'Actual_predicted_PLOT.png'
    ggsave(file,sp[[1]],width = 9, height = 6, dpi = 300, units = "in", device='png')
    readBin(file,'raw',n = file.info(file)$size)
  }
  else if (model_plot=="Iter") {
    wd <- getwd()
    f <- paste0(wd,"/", path, "/", filename)
    st<- base::readRDS(f)
    file = 'Iteration_PLOT.png'
    ggsave(file,st[[1]],width = 9, height = 6, dpi = 300, units = "in", device='png')
    readBin(file,'raw',n = file.info(file)$size)
  }
  else{
    print("Error occured: Please choose correct model_plot argument: Iter or Pre_Act")
  }
  }else{
    print("Error occured: Please choose correct filename: sp.rds or st.rds")
  }
}
###########Endpoint to serialze csv file from rds object for Actual vs predicted bacterial counts

#* @serializer csv
#* @get /Predicted_Actual
function(platform, filename, path, round) {
  #Hint:-
  #   platform = is the analytical platform used for data generation eg.(FTIR, eNose, MSI, HPLC, GCMS)
  #   filename = is a filename of the data eg. ( FTIR_data.csv or FTIR.xlsx) or rds object generated from MLearinig
  #   path  = is the folder name without the root directory eg.(Demo)
  wd <- getwd()
  f <- paste0(wd,"/", path, "/", filename)
  
  if (platform=="FTIR" && round=="yes") {
    if((see_if(has_extension(filename, 'csv')) == TRUE) | (see_if(has_extension(filename, 'xlsx')) == TRUE)){
      if (see_if(has_extension(filename, 'csv')) == TRUE){
        data<- utils ::read.csv(f, header = TRUE, sep = ",",  row.names=1, check.names = F)
      }else if (see_if(has_extension(filename, 'xlsx')) ==TRUE) {
        data<- openxlsx::read.xlsx(f, rowNames = T, check.names = F,sep.names=" ")
      }else{
        print("Error Occured: please choose correct file name with .xlsx | .csv extention")
      } 
    ##take only wavelength FTIR by removing other column
    col_nume <-function(x) !grepl("[^0-9]\\D", x)
    col_na<-c()
    ba<-c()
    for (i in colnames(data)) {
      if(col_nume(colnames(data[i]))==T){
        s<-colnames(data[i])
        col_na<-c(col_na, s)
      }else{
        ba<-c(ba, colnames(data[i]))
      }
    }
    data_wavelength<-data[,c(col_na)]
    data_bac<-data[,c(ba)]
     
    # round FTIR data wavelength or columns
      row_colnames<- colnames(data_wavelength)
      row_col_integer <- 0
      for(i in 1:length(row_colnames)){
        if(as.numeric(row_colnames[i])==round(as.numeric(row_colnames[i]), 0)){
          row_col_integer <- row_col_integer+1
        }
        else{
          row_col_integer <- row_col_integer
        }
      }
      if(length(row_colnames) == row_col_integer){
        data_wavelength<-data_wavelength  
      } else{
        data_wavelength <- data_wavelength
        data_wavelength <- roundWavelengths(data_wavelength)
      }
      ##combine the wavelength with other columns which was removed
      data<-cbind(data_bac,data_wavelength) 
      
    }else if (see_if(has_extension(filename, 'rds')) ==TRUE) {
      data<- base::readRDS(f)
      data<-data[[1]]
    }else{
      print("Error Occured: please choose correct file name with .xlsx | .csv extention")
    } 
     
  }else if(platform=="FTIR"){
      if (see_if(has_extension(filename, 'rds')) ==TRUE) {
        data<- base::readRDS(f)
        data<-data[[1]]
      }else{
        print("Error Occured: please choose correct file name with .rds extention")
      }
  }else if(platform=="HPLC" | platform=="eNose" | platform=="GCMS" | platform=="MSI"){
      if (see_if(has_extension(filename, 'rds')) ==TRUE) {
        data<- base::readRDS(f)
        data<-data[[1]]
      }else{
        print("Error Occured: please choose correct file name with .rds extention")
      }
    
  }else{
      print("Error occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI")
    }
##serialze data to csv file to be sent to client
df <- data.frame(CHAR = letters, NUM = rnorm(length(letters)), stringsAsFactors = F)
csv_file <- tempfile(fileext = ".csv")
on.exit(unlink(csv_file), add = TRUE)
write.csv(data, file = csv_file, row.names=T)
readBin(csv_file,'raw',n = file.info(csv_file)$size)
readLines(csv_file)

}


#####Combining existing dataset with new dataset for the same bacterial type, analytical and product type######

#* @serializer csv
#* @get /Combine_dataset
function(platform, filename1, filename2, path) {
  #Hint:-
  #   platform  = is the analytical platform used for data generation eg.(FTIR, eNose, MSI, HPLC, GCMS)
  #   filename1 = is a filename of the first data eg. ( FTIR1.csv or FTIR1.xlsx)
  #   filename2 = is a filename of the second data eg. ( FTIR1.csv or FTIR1.xlsx)
  #   path = is the folder name without the root directory eg.(Demo)
  fname1<-filename1
  fname2<-filename2
 
  if (platform=="MSI" | platform=="FTIR" | platform=="HPLC" | platform=="eNose" | platform=="GCMS") {
    data<-Batch_comb_M(fname1, fname2, path)
  }else{
    print("Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI ")
    
  }
  
  df <- data.frame(CHAR = letters, NUM = rnorm(length(letters)), stringsAsFactors = F)
  csv_file <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_file), add = TRUE)
  write.csv(data, file = csv_file, row.names=T)
  readBin(csv_file,'raw',n = file.info(csv_file)$size)
  readLines(csv_file)
} 

##########replicate arrangement checking 

#* @post /check_rep
function(filename, replicates, path){
  #Hint:-
  #  filename = is a filename of the  data eg. ( FTIR1.csv or FTIR1.xlsx)
  #  replicates = number of replicate in data (3,6)
  #  path = is the folder name without the root directory eg.(Demo)
  if (is.null(c(filename,path, replicates))) {
    res$status <- 400
    return(list(error = "Warrning: Please check Arguement is missing"))
  }
  n<-as.numeric(replicates)
  wd <- getwd()
  f <- paste0(wd,"/", path, "/", filename)
  
  if (see_if(has_extension(filename, 'csv')) == TRUE){
      data<- utils ::read.csv(f, header = TRUE, sep = ",",  row.names=1, check.names = F)
    }
  else if (see_if(has_extension(filename, 'xlsx')) ==TRUE) {
      data<- openxlsx::read.xlsx(f, rowNames = T, check.names = F, sep.names=" ")
    }
  resl<-rep_check(data,n)
 
  if((nrow(data)%%n!=0)==TRUE){
    msg1<-paste0("Warning message: Please check some of your samples doesn't contain replicate of :", n)
    print(msg1)
  }else if(resl=="NO"){
    m<-c(1:n)
    msg2<-paste0("Warning message: Please check your samples replicates aren't underneath of each other ", 1, ":", n,  " or some of your samples doesn't contain replicate of ", n)
    print(msg2)
  }else if(resl=="yes"){
    print("Replicates of Data meets the requirement!")
  }else if(resl=="null"){
    print("Warning message: unable to check replicates")
  }

}


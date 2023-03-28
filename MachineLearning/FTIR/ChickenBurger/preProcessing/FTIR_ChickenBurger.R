rm(list = ls())
graphics.off()


# Loading
library("readxl")
library("mixOmics")
library(prospectr) #savitzkyGolay
library (hyperSpec) #spectra plot
library(factoextra)
library(ggbiplot)
library(stringr) #sub and str_c
library(ggstatsplot) #to create a boxplot that labels the outliers 
library(writexl) #write xlsx 

# read R Code From A File (access generic functions)
source("CreateObj.R")

# read xlsx files
data <- read_excel("Chicken Burger_counts_ ftir spectra.xlsx")


# Extract the names of samples which will use to build prediction models
samples <- data[,1:2]
#marge columns(sampleID_time)
samples$sample <- as.character(interaction(samples,sep="_"))
samples <- samples[,-1]
samples <- samples[,-1]
# extract spectra
spectra <- data[,7:3357]
# extract bacterial counts
bacterialCounts <- data[,3:6]
TSA <- data[,3]
CFC <- data[,4]
STAA <- data[,5]
MRS <- data[,6]

####################################################################################################
#################### 1. CALCULATE MEANS BASED ON REPLICATES ########################################
####################################################################################################
# number of rows per group
n=6
df <- aggregate(data[,3:3357],list(rep(1:(nrow(data[,3:3357])%/%n+1),each=n,len=nrow(data[,3:3357]))),mean)[-1];
spectraM <- df[,5:3355]
#write a spectra as rds object
#saveRDS(spectraM, "./models/CB/FTIR/spectra_CB.rds")
bacterialCountsM <- df[,1:4]
TSA_M <- df[,1]
CFC_M <- df[,2]
STAA_M <- df[,3]
MRS_M <- df[,4]

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
samplesM <- nth_element(samplesV, 1, 6)

timeStorM <- nth_element(timeStor, 1, 6)


#temp 
t <- c()
for(i in 1:length(samplesM)){
  t[i] <- sub("B.*", "", samplesM[i])
}
#batch
b <- c()
for(i in 1:length(samplesM)){
  
  b[i] <- sub("\\d*B|DB", "", samplesM[i])
  b[i] <- substr(b[i], 1, 1)

}
for(i in 1:11){
  b[length(b)-i +1] <- "B"
}


############################################################################################
######################### 2. OUTLIERS HANDLING #############################################
############################################################################################
######################### a. BACTERIAL COUNTS ##############################################
############################################################################################

# create a labeled bacterial counts data frame
bacterialCountsM <- CreateObj(samplesNamneM, bacterialCountsM)
dfM <- CreateObj(bacterialCountsM, spectraM)
rownames(dfM) <- dfM[,1]
dfM <- dfM[,-1]

#################### seperate batches on raw data ##########################################
#create a dataframe with mean from batch A and batch B
batchA <- c()
for(i in 1:length(b)){
  if(b[i] == "A"){
    batchA[length(batchA)+1] <- i
  }
}
dfA <- dfM[batchA,]
dfA <- dfA[,-1]
dfB <- dfM[-batchA,]
dfB <- dfB[,-1]

#sample name for data with removed replicates
# temp + B=burger + batch (A or B) + time
samplesNamneM <- c()
for(i in 1:length(t)){
  samplesNamneM[i] <- str_c(t[i], "B",b[i] ,timeStorM[i], "_", i)
}
samplesNamneM <- matrix(samplesNamneM, length(samplesNamneM), 1)

#boxplot to identify outliers
boxplot(dfM[,1:4], 
        xlab = "medium",
        ylab = "bacterial counts") 
title("BACTERIAL COUNTS - BOXPLOT")
#get id of the outliers 
bc_outNum <- 0
id_outNum <- 0
for (i in 1: nrow(dfM)){
  if(dfM$STAA[i] < 4.15){
    bc_outNum = bc_outNum + 1
    id_outNum[length(id_outNum)+1] <- i
  }
}
id_outNum <- id_outNum[-1]
row.names(dfM[9:10,]) #DBA0_106 (3.440407) & DBA16_107 (4.008182)


#############################################################################################
##################### b. SPECTRA ############################################################
#############################################################################################

#split data into batches
batchNum <- c()
dfM_names <- row.names(dfM)
for(i in 1:length(dfM_names)){
  batchNum[i] <- sub("\\d*B|DB", "", dfM_names[i])
  batchNum[i] <- substr(batchNum[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM <- prcomp(dfM[,5:3355]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum, addEllipses = TRUE, ellipse.level = 0.95) 
# PC2 vs PC3
fviz_pca_ind(pca_dfM3, axes = c(2,3), geom = c("point", "text"), 
             col.ind = batchNum, addEllipses = TRUE) 

# PC1 vs PC3
fviz_pca_ind(pca_dfM3, axes = c(1,3), geom = c("point", "text"), 
             col.ind = batchNum, addEllipses = TRUE) 

#####################################################################################################
################ 3. NORMALISATION (spectra) #########################################################
#####################################################################################################

#range for data frame with means and outliers
min_dfM <- c() #vector with minimum values
max_dfM <- c() #vector with maximum values
for (i in 1:ncol(dfM[5:3355])){
  min_dfM[i] <- min(dfM[,i+4])
  max_dfM[i] <- max(dfM[,i+4])
}
min_dfM <- as.matrix(min_dfM)
max_dfM <- as.matrix(max_dfM)
num_dfM <- c(1: length(min_dfM)) #vector with sample numbers
num_dfM <- as.matrix(num_dfM)

#plot ranges
range_dfM <- CreateObj(min_dfM, max_dfM) #function to marge to datasets (it's on the github in the separate file)
range_dfM <- CreateObj(num_dfM, range_dfM)
colnames(range_dfM) <- c("num", "minimum", "maximum")
x11()
ggplot(range_dfM, aes(range_dfM[,1])) + geom_line(aes(y=range_dfM[,2]), colour="red")+
  labs( title="SPECTRUM RANGE")+
  geom_line(aes(y=range_dfM[,3]), colour="blue") +
  xlab("sample") + ylab('spectra') 


##############################################
# transform spectra
spectraT_dfM4 <- log2(dfM[,5:3355] + 1)
dfM5 <- CreateObj(dfM[,1:4], spectraT_dfM4) #df with transformed specrum values
# range
min_dfM5 <- c()
max_dfM5 <- c()
for (i in 1:ncol(dfM5[5:3355])){
  min_dfM5[i] <- min(dfM5[,i+4])
  max_dfM5[i] <- max(dfM5[,i+4])
}
min_dfM5 <- as.matrix(min_dfM5)
max_dfM5 <- as.matrix(max_dfM5)
num_dfM5 <- c(1: length(min_dfM5))
num_dfM5 <- as.matrix(num_dfM5)

#plot ranges
range_dfM5 <- CreateObj(min_dfM5, max_dfM5)
range_dfM5 <- CreateObj(num_dfM5, range_dfM5)
x11()
ggplot(range_dfM5, aes(range_dfM5[,1])) + geom_line(aes(y=range_dfM5[,2]), colour="red", show.legend = T)+ 
                                                   labs( title="SPECTRUM RANGE \n(log2 transformed)")+
  geom_line(aes(y=range_dfM5[,3]), colour="blue", show.legend = T) +
  xlab("sample") + ylab('spectra') 

#pca
#split data into batches
batchNum5 <- c()
dfM_names5 <- row.names(dfM5)
for(i in 1:length(dfM_names5)){
  batchNum5[i] <- sub("\\d*B|DB", "", dfM_names5[i])
  batchNum5[i] <- substr(batchNum5[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM5 <- prcomp(dfM5[,5:3355]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM5, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum5, addEllipses = TRUE, ellipse.level = 0.95) 

#####################################################################################################
################ 3. NORMALISATION (spectra) #########################################################
#####################################################################################################

sg_spectra6 <- savitzkyGolay(dfM5[,5:3355], w=11, p=3, m=1)
dfM6 <- CreateObj(dfM5[,1:4], sg_spectra6) #df with transformed specrum values
sg_spectra7 <- savitzkyGolay(dfM5[,5:3355], w=11, p=3, m=0)
dfM7 <- CreateObj(dfM5[,1:4], sg_spectra7) 

#hyperspec based on raw spectra
waveLength <- c(as.numeric(colnames(dfM[,5:3355])))
spc <- new ("hyperSpec", spc = dfM[,5:3355], wavelength = waveLength, data = dfM[,5:3355])
plotspc(spc)
#hyperspec based on svn raw spectra
waveLength8 <- c(as.numeric(colnames(dfM8[,5:3355])))
spc8 <- new ("hyperSpec", spc = dfM8[,5:3355], wavelength = waveLength8, data = dfM8[,5:3355])
plotspc(spc8)
#hyperspec based on log2 transformed spectra
waveLength5 <- c(as.numeric(colnames(dfM5[,5:3355])))
spc5 <- new ("hyperSpec", spc = dfM5[,5:3355], wavelength = waveLength5, data = dfM5[,5:3345])
plotspc(spc5)
#hyperspec based on log2 transformed spectra + SG ( w=11, p=3, m=1)
waveLength6 <- c(as.numeric(colnames(dfM6[,5:3345])))
spc6 <- new ("hyperSpec", spc = dfM6[,5:3345], wavelength = waveLength6, data = dfM6[,5:3345])
plotspc(spc6)
#hyperspec based on log2 transformed spectra + SG ( w=11, p=3, m=0)
waveLength7 <- c(as.numeric(colnames(dfM7[,5:3345])))
spc7 <- new ("hyperSpec", spc = dfM7[,5:3345], wavelength = waveLength7, data = dfM7[,5:3345])
plotspc(spc7)


###################################################################################################3
########################## NORMALISATION (after log2 + snv)
X1<-dfM[,5:3355]
Xt1<-t(X1)
Xt1_snv<-scale(Xt1,center=TRUE,scale=TRUE)
wavelengths<-as.numeric(colnames(dfM[,5:3355]))

matplot(wavelengths,(Xt1_snv),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
matplot(wavelengths,(Xt1),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
dfM8 <- t(Xt1_snv)
dfM8 <- CreateObj(dfM5[, 1:4], dfM8)
#pca
#split data into batches
batchNum8 <- c()
dfM_names8 <- row.names(dfM8)
for(i in 1:length(dfM_names8)){
  batchNum8[i] <- sub("\\d*B|DB", "", dfM_names8[i])
  batchNum8[i] <- substr(batchNum8[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM8 <- prcomp(dfM8[,5:3355]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM8, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum8, addEllipses = TRUE, ellipse.level = 0.95) 

# for log transformed data
X<-dfM5[,5:3355]
Xt<-t(X)
Xt_snv<-scale(Xt,center=TRUE,scale=TRUE)
waveLength5<-as.numeric(colnames(dfM5[,5:3355]))
x11()
matplot(waveLength5,(Xt_snv),lty=1,pch=21,
          xlab="data_points(nm)",ylab="log(1/R)")
matplot(waveLength5,(Xt),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
dfM7 <- t(Xt_snv)
dfM7 <- CreateObj(dfM5[, 1:4], dfM7)
#pca
#split data into batches
batchNum7 <- c()
dfM_names7 <- row.names(dfM7)
for(i in 1:length(dfM_names7)){
  batchNum7[i] <- sub("\\d*B|DB", "", dfM_names7[i])
  batchNum7[i] <- substr(batchNum7[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM7 <- prcomp(dfM7[,5:3355]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM7, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum7, addEllipses = TRUE, ellipse.level = 0.95) 

######################## batch normalisation #############################################################
#BATCH A
XA<-dfA[,5:3355]
XtA<-t(XA)
XtA_snv<-scale(XtA,center=TRUE,scale=TRUE)
wavelengths<-as.numeric(colnames(dfA[,5:3355]))

matplot(wavelengths,(XtA_snv),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
matplot(wavelengths,(XtA),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
dfA_svn <- t(XtA_snv)
dfA_svn <- CreateObj(dfA[, 1:4], dfA_svn)

#BATCH B
XB<-dfB[,5:3355]
XtB<-t(XB)
XtB_snv<-scale(XtB,center=TRUE,scale=TRUE)
wavelengths<-as.numeric(colnames(dfB[,5:3355]))

matplot(wavelengths,(XtB_snv),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
matplot(wavelengths,(XtB),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
dfB_svn <- t(XtB_snv)
dfB_svn <- CreateObj(dfB[, 1:4], dfB_svn)

##########################################################################################################
################## Savitzky-Golay ######################################################
#based on raw data
sg_spectra9 <- savitzkyGolay(dfM[,5:3355], w=11, p=3, m=1)
dfM9 <- CreateObj(dfM[,1:4], sg_spectra9) #df with transformed specrum values
sg_spectra10 <- savitzkyGolay(dfM[,5:3355], w=11, p=3, m=0)
dfM10 <- CreateObj(dfM[,1:4], sg_spectra10) 

waveLength9 <- c(as.numeric(colnames(dfM9[,5:3345])))
spc9 <- new ("hyperSpec", spc = dfM9[,5:3345], wavelength = waveLength9, data = dfM9[,5:3345])
plotspc(spc9)

waveLength10 <- c(as.numeric(colnames(dfM10[,5:3345])))
spc10 <- new ("hyperSpec", spc = dfM10[,5:3345], wavelength = waveLength10, data = dfM10[,5:3345])
plotspc(spc10)

#pca
#split data into batches
batchNum9 <- c()
dfM_names9 <- row.names(dfM9)
for(i in 1:length(dfM_names9)){
  batchNum9[i] <- sub("\\d*B|DB", "", dfM_names9[i])
  batchNum9[i] <- substr(batchNum9[i], 1, 1)
}
batchNum10 <- c()
dfM_names10 <- row.names(dfM10)
for(i in 1:length(dfM_names10)){
  batchNum10[i] <- sub("\\d*B|DB", "", dfM_names10[i])
  batchNum10[i] <- substr(batchNum10[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM9 <- prcomp(dfM9[,5:3345]) 
pca_dfM10 <- prcomp(dfM10[,5:3345]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM9, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum9, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM10, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum10, addEllipses = TRUE, ellipse.level = 0.95) 

##################################################################################
#based on snv data
sg_spectra11 <- savitzkyGolay(dfM8[,5:3355], w=11, p=3, m=1)
dfM11 <- CreateObj(dfM[,1:4], sg_spectra11) #df with transformed specrum values
sg_spectra12 <- savitzkyGolay(dfM8[,5:3355], w=11, p=3, m=0)
dfM12 <- CreateObj(dfM[,1:4], sg_spectra12) 

waveLength11 <- c(as.numeric(colnames(dfM11[,5:3345])))
spc11 <- new ("hyperSpec", spc = dfM11[,5:3345], wavelength = waveLength11, data = dfM11[,5:3345])
plotspc(spc11)

waveLength12 <- c(as.numeric(colnames(dfM12[,5:3345])))
spc12 <- new ("hyperSpec", spc = dfM12[,5:3345], wavelength = waveLength12, data = dfM12[,5:3345])
plotspc(spc12)

#pca
#split data into batches
batchNum11 <- c()
dfM_names11 <- row.names(dfM11)
for(i in 1:length(dfM_names11)){
  batchNum11[i] <- sub("\\d*B|DB", "", dfM_names11[i])
  batchNum11[i] <- substr(batchNum11[i], 1, 1)
}
batchNum12 <- c()
dfM_names12 <- row.names(dfM12)
for(i in 1:length(dfM_names12)){
  batchNum12[i] <- sub("\\d*B|DB", "", dfM_names12[i])
  batchNum12[i] <- substr(batchNum12[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM11 <- prcomp(dfM11[,5:3345]) 
pca_dfM12 <- prcomp(dfM12[,5:3345]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM11, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum11, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM12, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum12, addEllipses = TRUE, ellipse.level = 0.95)

##############################################################################################
#based on log2 transformed data
sg_spectra13 <- savitzkyGolay(dfM5[,5:3355], w=11, p=3, m=1)
dfM13 <- CreateObj(dfM5[,1:4], sg_spectra13) #df with transformed specrum values
sg_spectra14 <- savitzkyGolay(dfM5[,5:3355], w=11, p=3, m=0)
dfM14 <- CreateObj(dfM5[,1:4], sg_spectra14) 

waveLength13 <- c(as.numeric(colnames(dfM13[,5:3345])))
spc13 <- new ("hyperSpec", spc = dfM13[,5:3345], wavelength = waveLength13, data = dfM13[,5:3345])
plotspc(spc13)

waveLength14 <- c(as.numeric(colnames(dfM14[,5:3345])))
spc14 <- new ("hyperSpec", spc = dfM14[,5:3345], wavelength = waveLength14, data = dfM14[,5:3345])
plotspc(spc14)

#pca
#split data into batches
batchNum13 <- c()
dfM_names13 <- row.names(dfM13)
for(i in 1:length(dfM_names13)){
  batchNum13[i] <- sub("\\d*B|DB", "", dfM_names13[i])
  batchNum13[i] <- substr(batchNum13[i], 1, 1)
}
batchNum14 <- c()
dfM_names14 <- row.names(dfM14)
for(i in 1:length(dfM_names14)){
  batchNum14[i] <- sub("\\d*B|DB", "", dfM_names14[i])
  batchNum14[i] <- substr(batchNum14[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM13 <- prcomp(dfM13[,5:3345]) 
pca_dfM14 <- prcomp(dfM14[,5:3345]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM13, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum13, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM14, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum14, addEllipses = TRUE, ellipse.level = 0.95) 

##############################################################################################
#based on log2 transformed + snv normalised data
sg_spectra15 <- savitzkyGolay(dfM7[,5:3355], w=11, p=3, m=1)
dfM15 <- CreateObj(dfM7[,1:4], sg_spectra15) #df with transformed specrum values
sg_spectra16 <- savitzkyGolay(dfM7[,5:3355], w=11, p=3, m=0)
dfM16 <- CreateObj(dfM7[,1:4], sg_spectra16) 

waveLength15 <- c(as.numeric(colnames(dfM15[,5:3345])))
spc15 <- new ("hyperSpec", spc = dfM15[,5:3345], wavelength = waveLength15, data = dfM15[,5:3345])
plotspc(spc15)

waveLength16 <- c(as.numeric(colnames(dfM16[,5:3345])))
spc16 <- new ("hyperSpec", spc = dfM16[,5:3345], wavelength = waveLength16, data = dfM16[,5:3345])
plotspc(spc16)

#pca
#split data into batches
batchNum15 <- c()
dfM_names15 <- row.names(dfM15)
for(i in 1:length(dfM_names15)){
  batchNum15[i] <- sub("\\d*B|DB", "", dfM_names15[i])
  batchNum15[i] <- substr(batchNum15[i], 1, 1)
}
batchNum16 <- c()
dfM_names16 <- row.names(dfM16)
for(i in 1:length(dfM_names16)){
  batchNum16[i] <- sub("\\d*B|DB", "", dfM_names16[i])
  batchNum16[i] <- substr(batchNum16[i], 1, 1)
}
#pca based on spectra (without removing any outliers)
pca_dfM15 <- prcomp(dfM15[,5:3345]) 
pca_dfM16 <- prcomp(dfM16[,5:3345]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM15, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum15, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM16, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum16, addEllipses = TRUE, ellipse.level = 0.95)

#################################################################################################
######################## batch SG #############################################################
### raw
sg_spectra17 <- savitzkyGolay(dfA[,5:3355], w=11, p=3, m=1)
dfA_raw_sg <- CreateObj(dfA[,1:4], sg_spectra17) #batch A
sg_spectra18 <- savitzkyGolay(dfB[,5:3355], w=11, p=3, m=1)
dfB_raw_sg <- CreateObj(dfB[,1:4], sg_spectra18) # batch B

### normalised
sg_spectra19 <- savitzkyGolay(dfA_svn[,5:3355], w=11, p=3, m=1)
dfA_svn_sg <- CreateObj(dfA_svn[,1:4], sg_spectra19) #batch A
sg_spectra20 <- savitzkyGolay(dfB_svn[,5:3355], w=11, p=3, m=1)
dfB_svn_sg <- CreateObj(dfB_svn[,1:4], sg_spectra20) # batch B

###########################################################################################################
############## identify sample 12BA12_27
rowNam_dfM5 <- rownames(dfM5)
for(i in 1: length(rowNam_dfM5)){
  if(rowNam_dfM5[i]== "12BA12_27"){
    tsa_12BA12_27 <- dfM5$TSA[i]
    cfc_12BA12_27 <- dfM5$CFC[i]
    staa_12BA12_27 <- dfM5$STAA[i]
    mrs_12BA12_27 <- dfM5$MRS[i]
  }
}


###################################################################
########## eliminate the outlier from the raw dataset (12BA12_27)
dfM17 <- dfM
for (i in 1:nrow(dfM17)){
  if(rownames(dfM17[i,])=="12BA12_27"){
    dfM17 <- dfM17[-i,]
  }
}
#pca
#split data into batches
batchNum17 <- c()
dfM_names17 <- row.names(dfM17)
for(i in 1:length(dfM_names17)){
  batchNum17[i] <- sub("\\d*B|DB", "", dfM_names17[i])
  batchNum17[i] <- substr(batchNum17[i], 1, 1)
}

#pca based on spectra (without removing any outliers)
pca_dfM17 <- prcomp(dfM17[,5:3355]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM17, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum17, addEllipses = TRUE, ellipse.level = 0.95) 

# eliminate the outlier from the raw dataset (12BA12_27) in batches
out <- c()
for(i in 1:nrow(dfB)){
  if(row.names(dfB[i,])=="12BA12_27"){
    out[length(out)+1] <- i 
  }
}
#12BA12_27 is in the batch B
dfA_woo <- dfA
dfB_woo <- dfB[-out,]


###########################################################################################################
####################### write an xlsx with pre-processed data ###############################################
###########################################################################################################

# raw data
batchNum1df <- data.frame(batchNum)
row.names(batchNum1df) <- dfM_names
dfMb <- CreateObj(batchNum1df, dfM)
dfM_samplesNames <- row.names(dfM) 
dfM_samplesNames <- as.matrix(dfM_samplesNames )
row.names(dfM_samplesNames) <- dfM_samplesNames[,1]
dfM_withSampleNames <- CreateObj(dfM_samplesNames, dfMb)
write_xlsx(dfM_withSampleNames, path = "chickenBurger_FTIR_raw_v2.xlsx", col_names = TRUE)

# raw data where one outlier was removed version 2
batchNum17df <- data.frame(batchNum17)
row.names(batchNum17df) <- dfM_names17
dfM17b <- CreateObj(batchNum17df, dfM17)
dfM17_samplesNames <- row.names(dfM17b) 
dfM17_samplesNames <- as.matrix(dfM17_samplesNames )
row.names(dfM17_samplesNames) <- dfM17_samplesNames[,1]
dfM17_withSampleNames <- CreateObj(dfM17_samplesNames, dfM17b)
write_xlsx(dfM17_withSampleNames, path = "chickenBurger_FTIR_raw_woo_v2.xlsx", col_names = TRUE)

# data smoothed by Savitzky-Golay filter version 2
batchNum9df <- data.frame(batchNum9)
row.names(batchNum17df) <- dfM_names9
dfM9b <- CreateObj(batchNum17df, dfM9)
dfM9_samplesNames <- row.names(dfM9) 
dfM9_samplesNames <- as.matrix(dfM9_samplesNames )
row.names(dfM9_samplesNames) <- dfM9_samplesNames[,1]
dfM9_withSampleNames <- CreateObj(dfM9_samplesNames, dfM9b)
write_xlsx(dfM9_withSampleNames, path = "chickenBurger_FTIR_raw&sg_v2.xlsx", col_names = TRUE)

# normatised data smoothed using Savitzky-Golay filter version 2
batchNum11df <- data.frame(batchNum11)
row.names(batchNum11df) <- dfM_names11
dfM11b <- CreateObj(batchNum11df, dfM11)
dfM11_samplesNames <- row.names(dfM11) 
dfM11_samplesNames <- as.matrix(dfM11_samplesNames )
row.names(dfM11_samplesNames) <- dfM11_samplesNames[,1]
dfM11_withSampleNames <- CreateObj(dfM11_samplesNames, dfM11b)
write_xlsx(dfM11_withSampleNames, path = "chickenBurger_FTIR_snv&sg_v2.xlsx", col_names = TRUE)

######################### write batches ##################################################################
### raw data
# batch A
dfA_samplesNames <- row.names(dfA) 
dfA_samplesNames <- as.matrix(dfA_samplesNames )
row.names(dfA_samplesNames) <- dfA_samplesNames[,1]
dfA_withSampleNames <- CreateObj(dfA_samplesNames, dfA)
write_xlsx(dfA_withSampleNames, path = "chickenBurger_FTIR_raw_batchA.xlsx", col_names = TRUE)
# batch B
dfB_samplesNames <- row.names(dfB) 
dfB_samplesNames <- as.matrix(dfB_samplesNames )
row.names(dfB_samplesNames) <- dfB_samplesNames[,1]
dfB_withSampleNames <- CreateObj(dfB_samplesNames, dfB)
write_xlsx(dfB_withSampleNames, path = "chickenBurger_FTIR_raw_batchB.xlsx", col_names = TRUE)

### raw data where one outlier was removed
# batch A
dfA_woo_samplesNames <- row.names(dfA_woo) 
dfA_woo_samplesNames <- as.matrix(dfA_woo_samplesNames )
row.names(dfA_woo_samplesNames) <- dfA_woo_samplesNames[,1]
dfA_woo_withSampleNames <- CreateObj(dfA_woo_samplesNames, dfA_woo)
write_xlsx(dfA_woo_withSampleNames, path = "chickenBurger_FTIR_raw_woo_batchA.xlsx", col_names = TRUE)
# batch B
dfB_woo_samplesNames <- row.names(dfB_woo) 
dfB_woo_samplesNames <- as.matrix(dfB_woo_samplesNames )
row.names(dfB_woo_samplesNames) <- dfB_woo_samplesNames[,1]
dfB_woo_withSampleNames <- CreateObj(dfB_woo_samplesNames, dfB_woo)
write_xlsx(dfB_woo_withSampleNames, path = "chickenBurger_FTIR_raw_woo_batchB.xlsx", col_names = TRUE)

### data smoothed by Savitzky-Golay filter 
# batch A
dfA_raw_sg_samplesNames <- row.names(dfA_raw_sg) 
dfA_raw_sg_samplesNames <- as.matrix(dfA_raw_sg_samplesNames )
row.names(dfA_raw_sg_samplesNames) <- dfA_raw_sg_samplesNames[,1]
dfA_raw_sg_withSampleNames <- CreateObj(dfA_raw_sg_samplesNames, dfA_raw_sg)
write_xlsx(dfA_raw_sg_withSampleNames, path = "chickenBurger_FTIR_raw&sg_batchA.xlsx", col_names = TRUE)
# batch B
dfB_raw_sg_samplesNames <- row.names(dfB_raw_sg) 
dfB_raw_sg_samplesNames <- as.matrix(dfB_raw_sg_samplesNames )
row.names(dfB_raw_sg_samplesNames) <- dfB_raw_sg_samplesNames[,1]
dfB_raw_sg_withSampleNames <- CreateObj(dfB_raw_sg_samplesNames, dfB_raw_sg)
write_xlsx(dfB_raw_sg_withSampleNames, path = "chickenBurger_FTIR_raw&sg_batchB.xlsx", col_names = TRUE)

### normatised data smoothed using Savitzky-Golay filter
# batch A
dfA_svn_sg_samplesNames <- row.names(dfA_svn_sg) 
dfA_svn_sg_samplesNames <- as.matrix(dfA_svn_sg_samplesNames )
row.names(dfA_svn_sg_samplesNames) <- dfA_svn_sg_samplesNames[,1]
dfA_svn_sg_withSampleNames <- CreateObj(dfA_svn_sg_samplesNames, dfA)
write_xlsx(dfA_svn_sg_withSampleNames, path = "chickenBurger_FTIR_snv&sg_batchA.xlsx", col_names = TRUE)
# batch B
dfB_svn_sg_samplesNames <- row.names(dfB_svn_sg) 
dfB_svn_sg_samplesNames <- as.matrix(dfB_svn_sg_samplesNames )
row.names(dfB_svn_sg_samplesNames) <- dfB_svn_sg_samplesNames[,1]
dfB_svn_sg_withSampleNames <- CreateObj(dfB_svn_sg_samplesNames, dfB_svn_sg)
write_xlsx(dfB_svn_sg_withSampleNames, path = "chickenBurger_FTIR_snv&sg_batchB.xlsx", col_names = TRUE)



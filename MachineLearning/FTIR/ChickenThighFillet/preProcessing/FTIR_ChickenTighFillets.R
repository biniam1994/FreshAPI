# raw and log2 transformed data to train model
###############################################################################################3

rm(list = ls())
graphics.off()

# Loading
library("readxl") #read xlsx file
library("mixOmics")
library(factoextra)
library(prospectr) #SavitzkyGolay
library (hyperSpec) #spectra plot
library(graphics)
library(writexl) #write xlsx file

# read R Code From A File (access generic functions)
source("CreateObj.R")

# xlsx files
data <- read_excel("FTIRdata_CT_R1_2 _updated.xlsx")

# Extract the names of samples which will use to build prediction models
Samples <- data[,1]
# extract spectra
spectra <- data[,2:3737]
#write a spectra as rds object
saveRDS(spectra, "./models/CTF/FTIR/spectra_CTF.rds")
# extract bacterial counts
bacterialCounts <- data[,3738:3739]
TVC <- data[,3738]
Pseudomonas <- data[,3739]


############################################################################################
######################### 1. OUTLIERS HANDLING #############################################
############################################################################################
######################### a. BACTERIAL COUNTS ##############################################
############################################################################################

# create a labeled bacterial counts data frame
bacterialCountsM <- CreateObj(Samples, bacterialCounts)
dfM <- CreateObj(bacterialCountsM, spectra)
rownames(dfM) <- dfM[,1]
dfM <- dfM[,-1]

#boxplot to identify outliers
boxplot(dfM[,1:2], 
        xlab = "bacteria",
        ylab = "bacterial counts") 
title("BACTERIAL COUNTS - BOXPLOT")


#############################################################################################
##################### b. SPECTRA ############################################################
#############################################################################################

#split data into batches
batchNum <- c()
for(i in 1:nrow(Samples)){
  batchNum[i] <- sub(".*R", "", Samples[i,1])
}

# remove data with plateau
removeID <- c()
features <- colnames(dfM)
for(i in 3:ncol(dfM)){
  if(as.numeric(features[i]) <= 800){
    removeID[length(removeID)+1] <- i
  }
}
dfM2 <- dfM[,-removeID]

#split data into batches
batchNum <- c()
dfM2_names <- row.names(dfM2)
for(i in 1:length(dfM2_names)){
  if(as.numeric(sub(".*R", "", dfM2_names[i])) == 1){
    batchNum[i] <- "A"
  } else {
    batchNum[i] <- "B"
  }
}

#################### seperate batches on raw data ##########################################
#create a dataframe with mean from batch A and batch B
batchA <- c()
for(i in 1:length(batchNum)){
  if(batchNum[i] == "A"){
    batchA[length(batchA)+1] <- i
  }
}
dfA <- dfM2[batchA,]
dfB <- dfM2[-batchA,]


#pca based on spectra
pca_dfM2 <- prcomp(dfM2[,3:3337]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM2, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum, addEllipses = TRUE) 
# PC2 vs PC3
fviz_pca_ind(pca_dfM2, axes = c(2,3), geom = c("point", "text"), 
             col.ind = batchNum, addEllipses = TRUE) 

# PC1 vs PC3
fviz_pca_ind(pca_dfM2, axes = c(1,3), geom = c("point", "text"), 
             col.ind = batchNum, addEllipses = TRUE) 

##############################################
# transform spectra
spectraT_dfM4 <- log2(dfM2[,3:3322]+3 )
dfM5 <- CreateObj(dfM2[,1:2], spectraT_dfM4) #df with transformed specrum values
#split data into batches
batchNum5 <- c()
dfM5_names <- row.names(dfM5)
for(i in 1:length(dfM5_names)){
  if(as.numeric(sub(".*R", "", dfM5_names[i])) == 1){
    batchNum5[i] <- "A"
  } else {
    batchNum5[i] <- "B"
  }
}

################################### log2 transformed batches #############################
dfA_log2 <- log2(dfA[,3:3322]+3)
dfA_log2 <- CreateObj(dfA[,1:2], dfA_log2)
dfB_log2 <- log2(dfB[,3:3322]+2)
dfB_log2 <- CreateObj(dfB[,1:2], dfB_log2)

#pca based on spectra
pca_dfM5 <- prcomp(dfM5[ ,3:3337]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM5, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum5, addEllipses = TRUE) 

########################## NORMALISATION (after log2 + snv)
#### raw data
X1<-dfM2[,3:3337]
Xt1<-t(X1)
Xt1_snv<-scale(Xt1,center=TRUE,scale=TRUE)
wavelengths<-as.numeric(colnames(dfM2[,3:3337]))

matplot(wavelengths,(Xt1_snv),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
matplot(wavelengths,(Xt1),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
dfM8 <- t(Xt1_snv)
dfM8 <- CreateObj(dfM2[, 1:2], dfM8)
#pca
#split data into batches
batchNum8 <- c()
dfM_names8 <- row.names(dfM8)
for(i in 1:length(dfM_names8)){
  batchNum8[i] <- sub(".*R", "", dfM_names8[i])
}
#pca based on spectra (without removing any outliers)
pca_dfM8 <- prcomp(dfM8[,3:3337]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM8, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum8, addEllipses = TRUE, ellipse.level = 0.95) 


##############################
# for log transformed data
X<-dfM5[,3:3337]
Xt<-t(X)
Xt_snv<-scale(Xt,center=TRUE,scale=TRUE)
waveLength5<-as.numeric(colnames(dfM5[,3:3337]))
x11()
matplot(waveLength5,(Xt_snv),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
matplot(waveLength5,(Xt),lty=1,pch=21,
        xlab="data_points(nm)",ylab="log(1/R)")
dfM7 <- t(Xt_snv)
dfM7 <- CreateObj(dfM5[, 1:2], dfM7)
#pca
#split data into batches
batchNum7 <- c()
dfM_names7 <- row.names(dfM7)
for(i in 1:length(dfM_names7)){
  batchNum7[i] <- sub(".*R", "", dfM_names7[i])
}
#pca based on spectra (without removing any outliers)
pca_dfM7 <- prcomp(dfM7[,3:3337]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM7, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum7, addEllipses = TRUE, ellipse.level = 0.95) 

#hyperspec based on raw spectra
waveLength <- c(as.numeric(colnames(dfM2[,3:3337])))
spc <- new ("hyperSpec", spc = dfM2[,3:3337], wavelength = waveLength, data = dfM2[,3:3337])
plotspc(spc)
#hyperspec based on svn raw spectra
waveLength8 <- c(as.numeric(colnames(dfM8[,3:3337])))
spc8 <- new ("hyperSpec", spc = dfM8[,3:3337], wavelength = waveLength8, data = dfM8[,3:3337])
plotspc(spc8)
#hyperspec based on log2 transformed spectra
waveLength5 <- c(as.numeric(colnames(dfM5[,3:3337])))
spc5 <- new ("hyperSpec", spc = dfM5[,3:3337], wavelength = waveLength5, data = dfM5[,3:3337])
plotspc(spc5)
#hyperspec based on log2 transformed spectra + snv
waveLength7 <- c(as.numeric(colnames(dfM7[,3:3337])))
spc7 <- new ("hyperSpec", spc = dfM7[,3:3337], wavelength = waveLength7, data = dfM7[,3:3337])
plotspc(spc7)


##########################################################################################################
################## Savitzky-Golay
######################################################
#based on raw data
sg_spectra9 <- savitzkyGolay(dfM2[,3:3337], w=11, p=3, m=1)
dfM9 <- CreateObj(dfM2[,1:2], sg_spectra9) #df with transformed specrum values
sg_spectra10 <- savitzkyGolay(dfM2[,3:3337], w=11, p=3, m=0)
dfM10 <- CreateObj(dfM2[,1:2], sg_spectra10) 

waveLength9 <- c(as.numeric(colnames(dfM9[,3:3327])))
spc9 <- new ("hyperSpec", spc = dfM9[,3:3327], wavelength = waveLength9, data = dfM9[,3:3327])
plotspc(spc9)

waveLength10 <- c(as.numeric(colnames(dfM10[,3:3327])))
spc10 <- new ("hyperSpec", spc = dfM10[,3:3327], wavelength = waveLength10, data = dfM10[,3:3327])
plotspc(spc10)

#pca
#split data into batches
batchNum9 <- c()
dfM_names9 <- row.names(dfM9)
for(i in 1:length(dfM_names9)){
  batchNum9[i] <- sub(".*R", "", dfM_names9[i])
}
batchNum10 <- c()
dfM_names10 <- row.names(dfM10)
for(i in 1:length(dfM_names10)){
  batchNum10[i] <- sub(".*R", "", dfM_names10[i])
}
#pca based on spectra (without removing any outliers)
pca_dfM9 <- prcomp(dfM9[,3:3327]) 
pca_dfM10 <- prcomp(dfM10[,3:3327]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM9, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum9, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM10, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum10, addEllipses = TRUE, ellipse.level = 0.95) 

######################################################
#based on log2 transformed data
sg_spectra11 <- savitzkyGolay(dfM5[,3:3337], w=11, p=3, m=1)
dfM11 <- CreateObj(dfM5[,1:2], sg_spectra11) #df with transformed specrum values
sg_spectra12 <- savitzkyGolay(dfM5[,3:3337], w=11, p=3, m=0)
dfM12 <- CreateObj(dfM5[,1:2], sg_spectra12) 

waveLength11 <- c(as.numeric(colnames(dfM11[,3:3327])))
spc11 <- new ("hyperSpec", spc = dfM11[,3:3327], wavelength = waveLength11, data = dfM11[,3:3327])
plotspc(spc11)

waveLength12 <- c(as.numeric(colnames(dfM12[,3:3327])))
spc12 <- new ("hyperSpec", spc = dfM12[,3:3327], wavelength = waveLength12, data = dfM12[,3:3327])
plotspc(spc12)

#pca
#split data into batches
batchNum11 <- c()
dfM_names11 <- row.names(dfM11)
for(i in 1:length(dfM_names11)){
  batchNum11[i] <- sub(".*R", "", dfM_names11[i])
}
batchNum12 <- c()
dfM_names12 <- row.names(dfM12)
for(i in 1:length(dfM_names12)){
  batchNum12[i] <- sub(".*R", "", dfM_names12[i])
}
#pca based on spectra (without removing any outliers)
pca_dfM11 <- prcomp(dfM11[,3:3327]) 
pca_dfM12 <- prcomp(dfM12[,3:3327]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM11, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum11, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM12, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum12, addEllipses = TRUE, ellipse.level = 0.95) 

######################################################
#based on log2 transformed data + svn
sg_spectra13 <- savitzkyGolay(dfM7[,3:3337], w=11, p=3, m=1)
dfM13 <- CreateObj(dfM7[,1:2], sg_spectra13) #df with transformed specrum values
sg_spectra14 <- savitzkyGolay(dfM7[,3:3337], w=11, p=3, m=0)
dfM14 <- CreateObj(dfM7[,1:2], sg_spectra14) 

waveLength13 <- c(as.numeric(colnames(dfM13[,3:3327])))
spc13 <- new ("hyperSpec", spc = dfM13[,3:3327], wavelength = waveLength13, data = dfM13[,3:3327])
plotspc(spc13)

waveLength14 <- c(as.numeric(colnames(dfM14[,3:3327])))
spc14 <- new ("hyperSpec", spc = dfM14[,3:3327], wavelength = waveLength14, data = dfM14[,3:3327])
plotspc(spc14)

#pca
#split data into batches
batchNum13 <- c()
dfM_names13 <- row.names(dfM13)
for(i in 1:length(dfM_names13)){
  batchNum13[i] <- sub(".*R", "", dfM_names13[i])
}
batchNum14 <- c()
dfM_names14 <- row.names(dfM14)
for(i in 1:length(dfM_names14)){
  batchNum14[i] <- sub(".*R", "", dfM_names14[i])
}
#pca based on spectra (without removing any outliers)
pca_dfM13 <- prcomp(dfM13[,3:3327]) 
pca_dfM14 <- prcomp(dfM14[,3:3327]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM13, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum13, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM14, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum14, addEllipses = TRUE, ellipse.level = 0.95) 

######################################################
#based on  svn normalised data
sg_spectra15 <- savitzkyGolay(dfM8[,3:3337], w=11, p=3, m=1)
dfM15 <- CreateObj(dfM8[,1:2], sg_spectra15) #df with transformed specrum values
sg_spectra16 <- savitzkyGolay(dfM8[,3:3337], w=11, p=3, m=0)
dfM16 <- CreateObj(dfM8[,1:2], sg_spectra16) 

waveLength15 <- c(as.numeric(colnames(dfM15[,3:3327])))
spc15 <- new ("hyperSpec", spc = dfM15[,3:3327], wavelength = waveLength15, data = dfM15[,3:3327])
plotspc(spc15)

waveLength16 <- c(as.numeric(colnames(dfM16[,3:3327])))
spc16 <- new ("hyperSpec", spc = dfM16[,3:3327], wavelength = waveLength16, data = dfM16[,3:3327])
plotspc(spc16)

#pca
#split data into batches
batchNum15 <- c()
dfM_names15 <- row.names(dfM15)
for(i in 1:length(dfM_names15)){
  batchNum15[i] <- sub(".*R", "", dfM_names15[i])
}
batchNum16 <- c()
dfM_names16 <- row.names(dfM16)
for(i in 1:length(dfM_names16)){
  batchNum16[i] <- sub(".*R", "", dfM_names16[i])
}
#pca based on spectra (without removing any outliers)
pca_dfM15 <- prcomp(dfM15[,3:3327]) 
pca_dfM16 <- prcomp(dfM16[,3:3327]) 

#grouped PCA scores plot 
# grouped based on batch and condition
# PC1 vs PC2
fviz_pca_ind(pca_dfM15, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum15, addEllipses = TRUE, ellipse.level = 0.95) 
fviz_pca_ind(pca_dfM16, axes = c(1,2), geom = c("point", "text"), 
             col.ind = batchNum16, addEllipses = TRUE, ellipse.level = 0.95) 

#####################################################################################################
################ 3. NORMALISATION (spectra) #########################################################
#####################################################################################################

#range
min_dfM3 <- c() #vector with minimum values
max_dfM3 <- c() #vector with maximum values
for (i in 1:nrow(dfM3)){
  min_dfM3[i] <- min(dfM3[i, 3:3738])
  max_dfM3[i] <- max(dfM3[i, 5:3738])
}
min_dfM3 <- as.matrix(min_dfM3)
max_dfM3 <- as.matrix(max_dfM3)
num_dfM3 <- c(1: length(min_dfM3)) #vector with sample numbers
num_dfM3 <- as.matrix(num_dfM3)

#plot ranges
range_dfM3 <- CreateObj(min_dfM3, max_dfM3) #function to marge to datasets (it's on the github in the separate file)
range_dfM3 <- CreateObj(num_dfM3, range_dfM3)
x11()
ggplot(range_dfM3, aes(range_dfM3[,1])) + geom_line(aes(y=range_dfM3[,2]), colour="red")+
  labs( title="SPECTRUM RANGE")+
  geom_line(aes(y=range_dfM3[,3]), colour="blue") +
  xlab("sample") + ylab('spectra') 

###############################################
# remove columnns with 5.0 result
dfM4 <- dfM3
dfM4 <- dfM4[,-c(3:400)]
#range
min_dfM4 <- c() #vector with minimum values
max_dfM4 <- c() #vector with maximum values
for (i in 1:nrow(dfM4)){
  min_dfM4[i] <- min(dfM4[i, 3:3340])
  max_dfM4[i] <- max(dfM4[i, 3:3340])
}
min_dfM4 <- as.matrix(min_dfM4)
max_dfM4 <- as.matrix(max_dfM4)
num_dfM4 <- c(1: length(min_dfM4)) #vector with sample numbers
num_dfM4 <- as.matrix(num_dfM4)

#plot ranges
range_dfM4 <- CreateObj(min_dfM4, max_dfM4) #function to marge to datasets (it's on the github in the separate file)
range_dfM4 <- CreateObj(num_dfM4, range_dfM4)
x11()
ggplot(range_dfM4, aes(range_dfM4[,1])) + geom_line(aes(y=range_dfM4[,2]), colour="red")+
  labs( title="SPECTRUM RANGE")+
  geom_line(aes(y=range_dfM4[,3]), colour="blue") +
  xlab("sample") + ylab('spectra') 


##############################################
# transform spectra
spectraT_dfM4 <- log2(dfM4[,3:3340] + 1)
dfM5 <- CreateObj(dfM4[,1:2], spectraT_dfM4) #df with transformed specrum values
# range
min_dfM5 <- c()
max_dfM5 <- c()
for (i in 1:nrow(dfM5)){
  min_dfM5[i] <- min(dfM5[i, 3:3340])
  max_dfM5[i] <- max(dfM5[i, 3:3340])
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

#############################################
# normalisation
spectraT_dfM6<-t(dfM5[,3:3340])
spectraT_dfM6_snv<-scale(spectraT_dfM6,center=TRUE,scale=TRUE)
wavelengths<-as.numeric(colnames(dfM5[, 3:3340]))
# spectrum plot based on dfM5
matplot(wavelengths,(spectraT_dfM6),lty=1,pch=21,
        xlab="data_points(nm)" ,ylab="log(1/R)")
# spectrum plot based on dfM6 (SVN, scaled and centered)
matplot(wavelengths,(spectraT_dfM6_snv),lty=1,pch=21,
        xlab="data_points(nm)" ,ylab="log(1/R)")
# create dataframe with normalised data
spectraT_dfM6_snv_t <- t(spectraT_dfM6_snv)
dfM6 <- CreateObj(dfM5[,1:2], spectraT_dfM6_snv_t)
# range
min_dfM6 <- c()
max_dfM6 <- c()
for (i in 1:nrow(dfM6)){
  min_dfM6[i] <- min(dfM6[i, 3:3340])
  max_dfM6[i] <- max(dfM6[i, 3:3340])
}
min_dfM6 <- as.matrix(min_dfM6)
max_dfM6 <- as.matrix(max_dfM6)
num_dfM6 <- c(1: length(min_dfM6))
num_dfM6 <- as.matrix(num_dfM6)

#plot ranges
range_dfM6 <- CreateObj(min_dfM6, max_dfM6)
range_dfM6 <- CreateObj(num_dfM6, range_dfM6)
x11()
ggplot(range_dfM6, aes(range_dfM6[,1])) + geom_line(aes(y=range_dfM6[,2]), colour="red", show.legend = T)+ 
  labs( title="SPECTRUM RANGE \n(after normalisation)")+
  geom_line(aes(y=range_dfM6[,3]), colour="blue", show.legend = T) +
  xlab("sample") + ylab('spectra')

###########################################################################################
############ Savitzky-Golay
sg_spectra7 <- savitzkyGolay(dfM6[,3:3340], w=11, p=3, m=1)
dfM7 <- CreateObj(dfM6[,1:2], sg_spectra7) #df with transformed specrum values

waveLength7 <- c(as.numeric(colnames(dfM7[,3:3330])))
spc7 <- new ("hyperSpec", spc = dfM7[,3:3330], wavelength = waveLength7, data = dfM7[,3:3330])
plotspc(spc7)
waveLength6 <- c(as.numeric(colnames(dfM6[,3:3340])))
spc6 <- new ("hyperSpec", spc = dfM6[,3:3340], wavelength = waveLength6, data = dfM6[,3:3340])
plotspc(spc6)

#the next try w=11, 

sg_spectra8 <- savitzkyGolay(dfM6[,3:3340], p=3,  w= 11,m=0)
dfM8 <- CreateObj(dfM6[,1:2], sg_spectra8) #df with transformed specrum values
waveLength8 <- c(as.numeric(colnames(dfM8[,3:3330])))
spc8 <- new ("hyperSpec", spc = dfM8[,3:3330], wavelength = waveLength8, data = dfM8[,3:3330])
plotspc(spc8)
###############################################
# BASIC STATISTICS ############################
###############################################

# calulate summmary stats for all bacterial counts
summary(object = bacterialCounts)

################################
# statistics of raw spectra data 
#the smallest and the biggest values
range(spectra)
max(spectra)
# to check how many samples have maximusm value
num <- 0
for (i in 1:ncol(spectra)){
  for(j in 1: nrow(spectra)){
    if(as.numeric(spectra[j,i]) == 5.00000){
      num = num + 1
    } else {
      num = num
    }
  }
}

# sum
sum(spectra)
# mean
mean(as.matrix(spectra))
# median
median(as.matrix(spectra))
# sd
sd(as.matrix(spectra))
# 1st and 3rd quantile 
quantile(as.matrix(spectra), probs = c(.25, .75))


###############################################
###### SPECTRA STANDARISATION #################
###############################################

#apply log transform to magnify/scale spectra data
spectraLogTrans <- log2(spectra + 5)

#statistics of tranformed spectra data
#the smallest and the biggest values
range(spectraLogTrans)
# sum
sum(spectraLogTrans)
# mean
mean(as.matrix(spectraLogTrans))
# median
median(as.matrix(spectraLogTrans))
# sd
sd(as.matrix(spectraLogTrans))
# 1st and 3rd quantile 
quantile(as.matrix(spectraLogTrans), probs = c(.25, .75))

#####################################################
###### Prepare the data by creating objectives ######
#####################################################

# funcrion to create combined data sets (bacterial counts + spectrum)
obj <- CreateObj(bacterialCounts, spectra) # Objective raw data (bacterial counts + spectrum)
objTrans <- CreateObj(bacterialCounts, spectraLogTrans) # Objective raw data (bacterial counts + log transformed spectrum)

###############################################
############### PCA ###########################
###############################################

## Visualise FTIR data set by using PCA
# create groups based on batches
batch <- c(rep("1st batch", 100), rep("2nd batch", 98))

##########################################################
# pca based on raw spectra 
pca <- prcomp(obj[,3:3738], center = TRUE, scale = TRUE) # for raw spectra

# PC1 vs PC2
ggbiplot(pca, choices = c(1,2), obs.scale = 1, var.scale = 1, 
         groups = batch, ellipse = TRUE, 
         circle = TRUE) + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', 
                                                                  legend.position = 'top') 

# PC2 vs PC3
ggbiplot(pca, choices = c(2,3), obs.scale = 1, var.scale = 1, 
         groups = batch, ellipse = TRUE, 
         circle = TRUE) + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', 
                                                                  legend.position = 'top')

# PC1 vs PC3
ggbiplot(pca, choices = c(1,3), obs.scale = 1, var.scale = 1, 
         groups = batch, ellipse = TRUE, 
         circle = TRUE) + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', 
                                                                  legend.position = 'top')
#grouped PCA scores plot (PC1 vs PC2)
# grouped based on batch
fviz_pca_ind(pca, axes = c(1,2), geom = c("point", "text"), col.ind = batch, addEllipses = TRUE)

#scree plot
fviz_eig(pca, addlabels = TRUE)


###########################################################################################################
####################### write a csv with pre-processed data ###############################################
###########################################################################################################

#write down raw data with removed 3 outliers
batchNum2df <- data.frame(batchNum)
row.names(batchNum2df) <- dfM2_names
dfM2b <- CreateObj(batchNum2df, dfM2)
dfM2_samplesNames <- row.names(dfM2) 
dfM2_samplesNames <- as.matrix(dfM2_samplesNames )
row.names(dfM2_samplesNames) <- dfM2_samplesNames[,1]
dfM2_withSampleNames <- CreateObj(dfM2_samplesNames, dfM2b)
write_xlsx(dfM2_withSampleNames, path = "chickenThighFillet_FTIR_raw_v3.xlsx", col_names = TRUE)
#write down log2 transformed data with removed 3 outliers
batchNum5df <- data.frame(batchNum5)
row.names(batchNum5df) <- dfM5_names
dfM5b <- CreateObj(batchNum5df, dfM5)
dfM5_samplesNames <- row.names(dfM5) 
dfM5_samplesNames <- as.matrix(dfM5_samplesNames )
row.names(dfM5_samplesNames) <- dfM5_samplesNames[,1]
dfM5_withSampleNames <- CreateObj(dfM5_samplesNames, dfM5b)
write_xlsx(dfM5_withSampleNames, path = "chickenThighFillet_FTIR_log2_v3.xlsx", col_names = TRUE)

######################### write batches ##################################################################
### raw data
# batch A
dfA_samplesNames <- row.names(dfA) 
dfA_samplesNames <- as.matrix(dfA_samplesNames )
row.names(dfA_samplesNames) <- dfA_samplesNames[,1]
dfA_withSampleNames <- CreateObj(dfA_samplesNames, dfA)
write_xlsx(dfA_withSampleNames, path = "chickenThighFillet_FTIR_raw_batchA.xlsx", col_names = TRUE)
# batch B
dfB_samplesNames <- row.names(dfB) 
dfB_samplesNames <- as.matrix(dfB_samplesNames )
row.names(dfB_samplesNames) <- dfB_samplesNames[,1]
dfB_withSampleNames <- CreateObj(dfB_samplesNames, dfB)
write_xlsx(dfB_withSampleNames, path = "chickenThighFillet_FTIR_raw_batchB.xlsx", col_names = TRUE)

### log2 data
# batch A
dfA_log2_samplesNames <- row.names(dfA_log2) 
dfA_log2_samplesNames <- as.matrix(dfA_log2_samplesNames )
row.names(dfA_log2_samplesNames) <- dfA_log2_samplesNames[,1]
dfA_log2_withSampleNames <- CreateObj(dfA_log2_samplesNames, dfA_log2)
write_xlsx(dfA_log2_withSampleNames, path = "chickenThighFillet_FTIR_log2_batchA.xlsx", col_names = TRUE)
# batch B
dfB_log2_samplesNames <- row.names(dfB_log2) 
dfB_log2_samplesNames <- as.matrix(dfB_log2_samplesNames )
row.names(dfB_log2_samplesNames) <- dfB_log2_samplesNames[,1]
dfB_log2_withSampleNames <- CreateObj(dfB_log2_samplesNames, dfB_log2)
write_xlsx(dfB_log2_withSampleNames, path = "chickenThighFillet_FTIR_log2_batchB.xlsx", col_names = TRUE)



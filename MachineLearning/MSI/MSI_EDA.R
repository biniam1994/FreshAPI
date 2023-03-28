##############################
#                            #
#   Group_project MSI Data   #
#   Machine Learning         #
##############################
# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()


library(mixOmics)
library(openxlsx)
#install.packages("mixomics")

# Read datasets
MSI   <- read.xlsx("MSIdata_CT_R1_2_updated.xlsx", sheet=1, rowName=T)
bacteria<- MSI[,37:38]##Extracting bacterial column
##MSI Spectra 
MSI <- MSI[,1:18]
##boxplots_Based on Bacterial counts####################################################
boxplot(data.frame(bacteria), main = "Bacterial count Boxplot", ylab= "Bacterial counts", col = c("green","red"))####Bacteria boxplot.
######################################################PCA_raw#############################################

# create groups based on batches
batch <- c(rep("1st batch", 100), rep("2nd batch", 96))
x11()
pcaForEnose<-pca(MSI,ncomp = 8, scale = T, center=T)
plotIndiv(pcaForEnose,ind.names=T, comp = c(1, 2), ellipse = TRUE, ellipse.level = 0.95,legend = T, group=batch, style="ggplot2",title = "PCA of MSI Raw data 95% cl ellipse")
##Data standardization 
MSI_Nor<- preProcess(MSI, scale = T, center=T)
MSI_NP<- predict(MSI_Nor, MSI)
x11()
pcaForEnose<-pca(MSI_NP,ncomp = 8, scale = T, center=T)
plotIndiv(pcaForEnose,ind.names=T, comp = c(1, 2), ellipse = TRUE, ellipse.level = 0.95,legend = T, group=batch, style="ggplot2",title = "PCA of MSI PreProcessed data 95% cl ellipse")
##################################

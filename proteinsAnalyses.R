#packages used
library(devtools)
library(Biobase)

library(dplyr)



library(mice)
library(ggplot2)
library(corrplot)
library("Hmisc")
library("PerformanceAnalytics")
#library(ggbiplot)
library(reshape2)



#reading-------------------------------------------
source("dataLoader.R")

# data treatment --------------------------------------


#remove data with half of the set missing
filterData <- function() {
  comb_mat <<- filter(comb_mat, rowSums(is.na(comb_mat)) < 14)
}

#imputes data via MICE package
imputeData <- function() {
  comb_mat_mice <- mice(comb_mat, m = 5, method = 'pmm')
  comb_mat <<- complete(comb_mat_mice)
}


filterData()
imputeData()

#plotting distributions -----------------------------------------------------------------

boxplot(comb_mat, col=2, range = 0, ylab="Abundance")

#plots histograms per each gene and protein
plotHistograms <- function() {
  for (i in 1:24) {
    hist(comb_mat[,i], col=2, main = paste("Histogram of", colnames(comb_mat)[i]))
  }
}


#pca-----------------------------------------------------------------
#from https://www.statology.org/principal-components-analysis-in-r/




#linear discriminant analysis --------------------------------

#https://little-book-of-r-for-multivariate-analysis.readthedocs.io/en/latest/src/multivariateanalysis.html
#https://richardlent.github.io/post/multivariate-analysis-with-r/



#linear correlations-------------------------------------------

getCorrelationMatrices <- function() {
  pear_corr <<- cor(as.matrix(comb_mat), method = c("pearson"))  
  spear_corr <<- cor(as.matrix(comb_mat), method = c("spearman"))
  ken_corr <<- cor(as.matrix(comb_mat), method = c("kendall"))
  
}

getPearsonList <- function() {
  pear_corr_list <<- rep(0, 12)
  for (i in 1:12) {
    pear_corr_list[i] <<- c(pear_corr[i,i+12])
  }  
}


getSpearmanList <- function() {
  spear_corr_list <<- rep(0, 12)
  for (i in 1:12) {
    spear_corr_list[i] <<- c(spear_corr[i,i+12])
  }  
}

getKendallList <- function() {
  ken_corr_list <<- rep(0, 12)
  for (i in 1:12) {
    ken_corr_list[i] <<- c(ken_corr[i,i+12])
  }  
}


plotCorrelationPlots <- function() {
  corrplot(pear_corr, type = "upper", title = "Pearson Correlation", is.corr = TRUE)
  corrplot(spear_corr, type = "upper", title = "Spearman Correlation", is.corr = TRUE)
  corrplot(ken_corr, type = "upper", title = "Kendall Correlation", is.corr = TRUE)
}


plotCorrelationList <- function(cl) {
  xLabels <- c("uterus", "kidney", "testis", "panceas", "stomach", "prostate", "ovary", "thyroid gl.", "adrenal gl.", "salivary gl.", "spleen", "esophagus")
  #bar chart
  barplot(cl, names.arg = xLabels, xlab = "Tissue", ylab = "Correlation Coefficient", col = "blue", main = "Correlation Chart")
}

#regression-------------------------------------------

lr_data <- lm(comb_mat)

#linear
plotLR <- function() {
  for (i in 1:12) {
    plot(comb_mat[,i+12]~comb_mat[,i],main = paste("RNA vs Protein Scatter of ", colnames(comb_mat)[i+12]))
    lines(comb_mat[,i], fitted(lm(comb_mat[,i+12]~comb_mat[,i])), col="red")
  }
}


#plot residuals
#https://www.statology.org/residual-plot-r/#:~:text=How%20to%20Create%20a%20Residual%20Plot%20in%20R,normal%20distribution.%20...%204%20Produce%20a%20density%20plot.

plotResidualsLR <- function() {
  for (i in 1:12) {
    lr <- lm(comb_mat[,i+12]~comb_mat[,i], data = comb_mat)
    residuals <- residuals(lr)
    predicted <- predict(lr)
    
    
    plot(fitted(lr), residuals, main = paste("RNA vs Protein Fitted Scatter of ", colnames(comb_mat)[i+12]))
    abline(0,0)
    qqnorm(residuals, main = paste("QQ Plot of RNA vs Protein of ", colnames(comb_mat)[i+12]))
    qqline(residuals)
    
    
    #ggplot(comb_mat, aes(x = comb_mat[,1], y = comb_mat[,1+12])) +  # Set up canvas with outcome variable on y-axis
    #  geom_segment(aes(xend = comb_mat[,1], yend = predicted)) +
    #  geom_point() + 
    #  geom_point(aes(y = predicted), shape = 1)
  }  
}


#logistic
#https://www.r-bloggers.com/2016/08/visualising-residuals/








#k-means clustering ---------------------------------------
 




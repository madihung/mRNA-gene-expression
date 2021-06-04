# Packages ----
library(shiny)  # Required to run any Shiny app
library(ggplot2)  # For creating pretty plots
library(dplyr)  # For filtering and manipulating data
library(tidyr)
library(tidyverse)
library(data.table)
library(mice)
library(corrplot)
library("FactoMineR")
library("factoextra")

source("dataLoader.R")
source("proteinsAnalyses.R")
source("geneAnalysis.R")

# Preparing data ----------
geneUN <- read.csv("https://raw.githubusercontent.com/gcohenfr/Neonatal-Vaccination/main/data/geneUN.csv?token=AMFFFE76ZJVXMSCF2JVQTD3AW3JV2")
proteinUN <- read.csv("https://raw.githubusercontent.com/gcohenfr/Neonatal-Vaccination/main/data/proteinUN.csv?token=AMFFFE77VWMLYQG5QGQQZYLAX3TF4")
#Replacing NA values using MICE package
geneUN.mice <- mice(geneUN, method = 'pmm')
proteinUN.mice <- mice(proteinUN, method = 'pmm')
geneUN <- complete(geneUN.mice)
proteinUN <- complete(proteinUN.mice)
#Logarithmic Scale
geneUN.log <- geneUN
geneUN.log[,2:13] <- log(geneUN[,2:13]+ 0.0000000001)
proteinUN.log <- proteinUN
proteinUN.log[,2:13] <- log(proteinUN[,2:13]+ 0.0000000001)
#Merging Data
geneUN.log <- setNames(geneUN.log,paste(colnames(geneUN.log), sep="", "_ge"))
colnames(geneUN.log)[colnames(geneUN.log) == "X_ge"] <- "X"
mergedData <- full_join(proteinUN.log, geneUN.log, by = intersect(colnames(geneUN.log), colnames(proteinUN.log)))



# Functions -------------

#For plotting tissue-specific protein vs. gene expression 
tissue_plot <- function(tissue) {
  if(tissue == "uterus") {
    ggplot(data=mergedData,aes(uterus_ge, uterus)) +
      xlab("log2(Gene Expression Level)") +
      ylab("log2(Protein Expression Level)") +
      coord_cartesian(xlim = c(0, 0.002)) + 
      coord_cartesian(ylim = c(-20, 0)) +
      geom_point(size=0.5, fill="white") +
      geom_smooth(method='lm')
  } else if (tissue == "kidney") {
    ggplot(data=mergedData,aes(kidney_ge, kidney)) +
      xlab("log2(Gene Expression Level)") +
      ylab("log2(Protein Expression Level)") +
      coord_cartesian(xlim = c(0, 0.005)) + 
      coord_cartesian(ylim = c(-20, 0)) + 
      geom_point(size=0.7, fill="white") + 
      geom_smooth(method='lm')
  } else if (tissue == "testis") {
    ggplot(data=mergedData,aes(testis_ge, testis)) +
      xlab("log2(Gene Expression Level)") +
      ylab("log2(Protein Expression Level)") +
      coord_cartesian(xlim = c(0, 0.005)) + 
      coord_cartesian(ylim = c(-20, 0)) + 
      geom_point(size=0.7, fill="white") + 
      geom_smooth(method='lm')  
  } else if (tissue == "pancreas") {
    ggplot(data=mergedData,aes(pancreas_ge, pancreas)) +
      xlab("log2(Gene Expression Level)") +
      ylab("log2(Protein Expression Level)") +
      coord_cartesian(xlim = c(0, 0.005)) + 
      coord_cartesian(ylim = c(-20, 0)) + 
      geom_point(size=0.7, fill="white") + 
      geom_smooth(method='lm')  
  } else if (tissue == "stomach") {
    ggplot(data=mergedData,aes(stomach_ge, stomach)) +
      xlab("log2(Gene Expression Level)") +
      ylab("log2(Protein Expression Level)") +
      coord_cartesian(xlim = c(0, 0.005)) + 
      coord_cartesian(ylim = c(-20, 0)) + 
      geom_point(size=0.7, fill="white") + 
      geom_smooth(method='lm')
  }
}
  
#For Correlations
pear_corr <<- cor(mergedData[sapply(mergedData,is.numeric)], method = c("pearson"))
spear_corr <<- cor(mergedData[sapply(mergedData,is.numeric)], method = c("spearman"))
ken_corr <<- cor(mergedData[sapply(mergedData,is.numeric)], method = c("kendall"))

md.pca <- PCA(mergedData[,c(2:13)], scale.unit = TRUE, ncp = 5, graph = FALSE) 

scree_plot <- function() {
  eig.val <- get_eigenvalue(md.pca)
  fviz_eig(md.pca, addlabels = TRUE, ylim = c(0, 50))  #Scree Plot - plot of eigenvalues ordered from largest to the smallest
}

cos_plot <- function() {
  fviz_cos2(md.pca, choice = "var", axes = 1:2) 
}

corr_circle <- function() {
  fviz_pca_var(md.pca, col.var = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
               repel = TRUE) # Avoid text overlapping
               
}





# unused code ----------
tissue_plot1 <- function(ge, pe) {
  ggplot(data=mergedData,aes(ge, pe)) +
    coord_cartesian(xlim = c(-25, -8)) + 
    coord_cartesian(ylim = c(-25, 0)) +
    geom_point() +
    geom_smooth(method='lm')  #Add a regression line
}

tissue_plot2 <- function(ge, pe) {
  lm <- lm(formula = ge ~ pe, data = mergedData)
  plot(x=mergedData$ge, 
       y=mergedData$pe,  
       xlab = "log2(gene expression levels)",
       ylab = "log2(protein expression levels)",
       main = "Regression fits of protein and gene expression for tissue type")
  abline(lm) #Success
}
#----------








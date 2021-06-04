#Literatures referenced:
#https://www.nature.com/articles/nature23293
#https://www.nature.com/articles/nature13319

##Packages--------
install.packages("dplyr")
install.packages("ggpubr")
install.packages("ggplot2")
install.packages("ggcorrplot")
install.packages("mice")
library(dplyr)
library(tidyr)
library(tidyverse)
library("data.table")
library(mice)

##Preparing data------
#Load Data
geneUN <- read.csv("https://raw.githubusercontent.com/gcohenfr/Neonatal-Vaccination/main/data/geneUN.csv?token=AMFFFE76ZJVXMSCF2JVQTD3AW3JV2")
proteinUN <- read.csv("https://raw.githubusercontent.com/gcohenfr/Neonatal-Vaccination/main/data/proteinUN.csv?token=AMFFFE5A5EVFXMS4V56YUIDAW3JYI")
#1) Replace NA values with median
#https://statisticsglobe.com/replace-missing-values-by-column-mean-in-r
for(i in 1:ncol(geneUN)) {   
  geneUN[ , i][is.na(geneUN[ , i])] <- median(geneUN[ , i], na.rm = TRUE)
}
for(i in 1:ncol(proteinUN)) {   
  proteinUN[ , i][is.na(proteinUN[ , i])] <- median(proteinUN[ , i], na.rm = TRUE)
}
#2) Replace NA values using MICE package
#Replace NA values - replaces data via MICE package  
  geneUN.mice <- mice(geneUN, m=5, method = 'pmm')
  proteinUN.mice <- mice(proteinUN, m=5, method = 'pmm')
  geneUN <- complete(geneUN.mice)
  proteinUN <- complete(proteinUN.mice)

#Logarithmic Scale
#1) longer 
proteinUN.log <- proteinUN
for(i in 2:ncol(proteinUN.log)) {  
 proteinUN.log[ ,i] <- log2(proteinUN.log[ ,i] + 0.0000000001)
}
geneUN.log <- geneUN
for(i in 2:ncol(geneUN.log)) {  
  geneUN.log[ ,i] <- log2(geneUN.log[ ,i] + 0.0000000001)
}
#2) shorter
geneUN.log <- geneUN
geneUN.log[,2:13] <- log(geneUN[,2:13]+ 0.0000000001)
proteinUN.log <- proteinUN
proteinUN.log[,2:13] <- log(proteinUN[,2:13]+ 0.0000000001)

#Viewing expression levels on plots
plot(geneUN[,i],ylab="gene expression level",type="l",xlab="tissue type index",ylim=c(0,0.02))
plot(geneUN.log[,i],ylab="gene expression level",type="l",xlab="tissue type index",ylim=c(-30,0.08))
plot(proteinUN[,i],ylab="protein expression level",type="l",xlab="tissue type index",ylim=c(0,0.08))
plot(proteinUN.log[,i],ylab="log2(protein expression level)",type="l",xlab="tissue type index",ylim=c(-30,0.08))


##Merge Datasets---------
#1) renaming column and row names
#attempt1: iterate through column names and rename
new_names1 <- paste(colnames(geneUN),"_ge")
for(i in geneUN.log[i]){
  rename(geneUN.log[i] = new_names[i])
}
#attempt2: using data.tables package
#https://stackoverflow.com/questions/44064983/whats-the-best-way-to-add-a-specific-string-to-all-column-names-in-a-dataframe/44065175
setnames(geneUN.log,paste(colnames(geneUN.log), sep="", "_ge"))
#need to rename just column X_ge -> X
#https://statisticsglobe.com/rename-column-name-in-r-data-frame/
colnames(geneUN.log)[colnames(geneUN) == "X_ge"] <- "X"


##2) joining dataframes
#attempt 1: rbind
rbind_test <- rbind(geneUN, proteinUN)  #does not workUp
# Error in match.names(clabs, names(xi)) : names do not match previous names
#attempt2: cbind
cbind_test <- cbind(geneUN, proteinUN) #only combines data set does not merge column X
#attempt3: dplyr full join
mergedData <- full_join(proteinUN.log, geneUN.log, by = intersect(colnames(geneUN.log), colnames(proteinUN.log))) #successfully merges datasets 


##Using Linear Regression Model---------
#attempt 1:
#https://www.datacamp.com/community/tutorials/linear-regression-R#what
lmExpression = lm(uterus_pe~uterus_ge, data = mergedData)   #create linear regression
summary(lmExpression) #review results
#Scatter plot ge vs. pe of one tissue
#uterus tissue 
ggplot(data=mergedData,aes(uterus_ge, uterus)) +
  coord_cartesian(xlim = c(0, 0.002)) + 
  coord_cartesian(ylim = c(-25, 0)) +
  geom_point(size=0.5) +
  geom_smooth(method='lm')  #Add a regression line
#kidney tissue
ggplot(data=mergedData,aes(kidney_ge, kidney)) +
  coord_cartesian(xlim = c(0, 0.005)) + 
  coord_cartesian(ylim = c(-25, 0)) + 
  geom_point(size=0.5) + 
  geom_smooth(method='lm')  #Add a regression line - success
#attempt 2:
#https://bookdown.org/ndphillips/YaRrr/linear-regression-with-lm.html
# Regression line formula
uterus.lm <- lm(formula = uterus_ge ~ uterus_pe,
                  data = mergedData)
plot(x=mergedData$uterus_ge, 
     y=mergedData$uterus_pe,  
     xlab = "log2(gene expression levels)",
     ylab = "log2(protein expression levels)",
     main = "Regression fits of protein and gene expression for uterus tissue")
abline(uterus.lm) #Success
#attempt 3:
ggplot(mergedData, aes(x = uterus_ge, y = uterus_pe)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") #Success

#Using Multiple Linear Regression
#https://stackoverflow.com/questions/15633714/adding-a-regression-line-on-a-ggplot

##Correlation Coefficient----------
#https://machinelearningmastery.com/how-to-use-correlation-to-understand-the-relationship-between-variables/#:~:text=The%20Pearson's%20correlation%20coefficient%20is,to%20give%20an%20interpretable%20score.
#Test for Guassian-like distribution
library("dplyr")
library("ggpubr")
ggqqplot(geneUN$uterus_ge)      #Can only test one column or whole dataset
shapiro.test(geneUN$uterus_ge[0:5000])
#data:  geneUN$uterus_ge[0:5000]
#W = 0.3239, p-value < 2.2e-16
#not normally distributed
##Rotate table
#attempt 1: transposing data frame -> turns it into a matrix which loses some properties
#https://www.dataanalytics.org.uk/rotating-or-transposing-r-objects/#:~:text=Rotating%20or%20transposing%20R%20objects&text=frame%20so%20that%20the%20rows,use%20the%20t()%20command.&text=The%20result%20of%20the%20t()%20command%20is%20always%20a%20matrix%20object.
tgeneUN <- t(geneUN)
#Make gene names header
colnames(tgeneUN) <- tgeneUN[1, ] 
#Delete first row 
tgeneUN <- tgeneUN[-1,]
#attempt 2:  -----

##Spearman Correlation
scorr_uterus <-cor.test(mergedData$uterus_ge, mergedData$uterus_pe,  method = "spearman", exact = FALSE) 
scorr_kidney <-cor.test(mergedData$kidney_ge, mergedData$kidney_pe,  method = "spearman", exact = FALSE) 


##Creating a Corrplot-------------
#1) converting data frame to a matrix
data.matrix(mergedData, rownames.force = NA) #unsure if this is required
#2) Calculated correlation coefficients
library(corrplot)
mergedData.scor <- cor(mergedData[sapply(mergedData,is.numeric)], method = c("spearman"))
mergedData.pcor <- cor(mergedData[sapply(mergedData,is.numeric)], method = c("pearson"))
#3) plotted correlation coefficients
corrplot(mergedData.pcor, method = "square")
corrplot(mergedData.scor, method = "square")

##Creating a Histogram of Correlations
#Calculating correlation coefficients
pear_corr <<- cor(mergedData[sapply(mergedData,is.numeric)], method = c("pearson"))
spear_corr <<- cor(mergedData[sapply(mergedData,is.numeric)], method = c("spearman"))
ken_corr <<- cor(mergedData[sapply(mergedData,is.numeric)], method = c("kendall"))
#replicating coeff into a vector?
pear_corr_list <<- rep(pear_corr)
spear_corr_list <<- rep(spear_corr)
ken_corr_list <<- rep(ken_corr)
#plotting values



# Network Plot -----------
install.packages("corrr")
library(corrr)
spear_corr %>% correlate() %>% network_plot(min_cor=0.4)


#library
install.packages("igraph")
library(igraph)
# build the graph object
network <- graph_from_adjacency_matrix(spear_corr)
# plot it
plot(network)

# Principle Component Analysis
pca <- prcomp(mergedData[,c(2:13)], center = TRUE,scale. = TRUE)
install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")
 
fviz_cos2(md.pca, choice = "var", axes = 1:2)
fviz_pca_var(md.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
fviz_contrib(md.pca, choice = "var", axes = 1:2)
fviz_pca_ind(md.pca, pointsize = 0.5, 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)




# Correlation Box Plots ----------
subsetData <- mergedData  #subsetData <- comb_mat_t 
subsetDataRNA <- subsetData[1:12,] #subsetDataRNA <- subsetData[1:12,]
subsetDataPro <- subsetData[13:24,]  #subsetDataPro <- subsetData[13:24,]
rm(subsetData)

correlateAcross <- function(pro_data,rna_data) {
  across_corrs <- rep(0,nrow(pro_data))
  for (i in 1:nrow(pro_data)) {
    #across_corrs[i] <- cor(rna_data[i,], pro_data[i,], method = c("spearman"))
    cor(rna_data[sapply(rna_data,is.numeric)], pro_data[sapply(pro_data,is.numeric)], method = c("spearman"))
  }
  return (across_corrs)
}

estimateOriginalMethod <- function(pro_data, rna_data) {
  
  transl_rates <- rep(2,ncol(pro_data))
  for (j in 2: ncol(pro_data)) {
    rsf <- rep(0,12)
    for(i in 1:12) { 
      rsf[i] <- (pro_data[i,j]/rna_data[i,j])
    }
    #dotchart(rsf, main = paste("protein_gt/mRNA_gt ", colnames(subsetData)[j]))
    
    transl_rates[j] <- median(rsf)
    #plot(rsf, xaxt = "n",col= "orange", ylab = "Ratio Protein/mRNA", xlab = "Tissue", main = paste("Median Protein/mRNA rate of ", rownames(subsetData)[j]))
    #axis(1, at = 1:12, labels = rownames(subsetDataPro))
    #abline(a = median(rsf), b = 0)
  }
  return(transl_rates)
}

#in: translation rates list and rna data matrix. Genes are cols
#INVARIANT: len(transl_rates) == len(colnames(rna_data))
#out: matrix of predicted values
predictOriginalMethod <- function(transl_rates, rna_data) {
  predict_mat <- matrix(nrow = 12,ncol = ncol(rna_data))
  for (j in 2:ncol(rna_data)) {
    for (i in 2:12) {
      predict_mat[i, j] <- transl_rates[j] * rna_data[i,j]
    }
  }
  return (predict_mat)
}



orig_data_across_corrs <- correlateAcross(subsetDataPro,subsetDataRNA)  #for FigA

orig_transl_rates <- estimateOriginalMethod(subsetDataPro, subsetDataRNA) 
orig_predictions <- predictOriginalMethod(orig_transl_rates, subsetDataRNA)
orig_predict_across_corrs <- correlateAcross(subsetDataPro,orig_predictions)  #for FigA

orig_predict_within_corrs <- correlateWithin(subsetDataPro,orig_predictions)

free_predictions <- estimateRNAFree(subsetDataPro)
#free_predict_within_corrs <- correlateWithin(subsetDataPro,free_predictions) #sdev is 0
free_predict_across_corrs <- correlateAcross(subsetDataPro, free_predictions) #for FigA

rand_predictions <- estimateRandomMethod(subsetDataPro,subsetDataRNA)
rand_predict_within_corrs <- correlateWithin(subsetDataPro,rand_predictions)
rand_predict_across_corrs <- correlateAcross(subsetDataPro, rand_predictions) #for FigA

plotFigureA <- function(data_corrs, orig_across, free_across, rand_across) {
  correlations_df <- as.data.frame(cbind(data_corrs,
                                         orig_across,
                                         free_across,
                                         rand_across))
  colnames(correlations_df) <- c("Wilhelm et al. (mRNA)", "Wilhelm et al. (predictions)", "mRNA-free (predictions)", "Random genes (predictions)")
  correlations_m <- melt(correlations_df, measure.vars = colnames(correlations_df))
  
  myplot <- ggplot(data = correlations_m,aes(x = correlations_m[,1], y = correlations_m[,2])) +
    geom_boxplot(color="darkgrey", fill="grey") +
    geom_point() +
    labs(x = "", y = "Correlation across genes") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  return(myplot)
}

figA_across <- plotFigureA(orig_data_across_corrs,
                           orig_predict_across_corrs,
                           free_predict_across_corrs,
                           rand_predict_across_corrs)

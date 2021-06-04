#packages used
library(devtools)
library(Biobase)


library(mice)
library(dplyr)
library(ggplot2)
library(reshape2)


#reading-------------------------------------------
source("dataLoader.R")

comb_mat <- dplyr::filter(comb_mat, rowSums(is.na(comb_mat)) <= 0)
comb_mat <- log2(comb_mat+0.000000000001) #transformation
comb_mat_t = t(comb_mat)
rm(comb_mat)
# data treatment --------------------------------------

#FN3K <- comb_mat_t[,"ENSG00000167363"] #strong pos correlation rate = 0.397617016
#HNRNPH2 <- comb_mat_t[,"ENSG00000126945"] #strong neg correlation r = -0.082980364
#UFC1 <- comb_mat_t[,"ENSG00000143222"] #-0.026910727
#PPIA <- comb_mat_t[,"ENSG00000196262"] #wrong 1.661287122
#PPP1CC <- comb_mat_t[,"ENSG00000186298"] #-0.430659209
#DIABLO <- comb_mat_t[,"ENSG00000184047"]#wrong 0.475429122
#SAMHD1 <- comb_mat_t[,"ENSG00000101347"] #0.667574619
#ALDH2 <- comb_mat_t[,"ENSG00000111275"] #0.88018996
#CORO1B <- comb_mat_t[,"ENSG00000172725"]#wrong 0.465624524
#CAPN1 <- comb_mat_t[,"ENSG00000014216"] #0.817671031
#LAP3 <- comb_mat_t[,"ENSG00000002549"] #low correlation within, high across ;10
#ITGA3 <- comb_mat_t[,"ENSG00000005884"] #low correlation within, high across ;14

#subsetData <- cbind(FN3K,HNRNPH2,UFC1,PPIA,PPP1CC,DIABLO,SAMHD1,ALDH2,CORO1B,CAPN1, LAP3,ITGA3) #sample of 11
subsetData <- comb_mat_t 
subsetDataRNA <- subsetData[1:12,]
subsetDataPro <- subsetData[13:24,]
rm(subsetData)
#original method***************************************************


#in: protein and rna matrix data. Genes are x axis
#out: list of translation rates per gene
estimateOriginalMethod <- function(pro_data, rna_data) {
  
  #for (i in 1:12) { rsf[i] <-  (roughVecRNA[i]/(roughVec[i]))}
  transl_rates <- rep(0,ncol(pro_data))
  for (j in 1: ncol(pro_data)) {
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
  for (j in 1:ncol(rna_data)) {
    for (i in 1:12) {
      predict_mat[i, j] <- transl_rates[j] * rna_data[i,j]
    }
  }
  return (predict_mat)
}


#linear
plotOriginalMethod <- function(pro_data,rna_data,pro_predict) {
  for (i in 1:ncol(pro_data)) {
    plot(pro_data[,i]~rna_data[,i],main = paste ("protein_gt = r_g * miRNA_gt ", colnames(subsetData)[i]))
    lines(x = rna_data[,i], y = pro_predict[,i], col="red")
  }
}
#plotOriginalMethod(subsetDataPro,subsetDataRNA,orig_predictions)


#first control experiment ************************************

#in: protein data matrix
#out: matrix of median protein abundance in each gene
estimateRNAFree <- function(pro_data) {
  mRNA_free_list <- rep(0,ncol(pro_data))
  for (i in 1:ncol(pro_data)) {
    mRNA_free_list[i] <- median(pro_data[,i])
  }
  mRNA_free_mat <- rbind(mRNA_free_list)
  for (i in 2:nrow(pro_data)) {
    mRNA_free_mat <- rbind(mRNA_free_mat, mRNA_free_list)
  }
  return (mRNA_free_mat)
}



#second control experiment **********************************


estimateRandomMethod <- function(pro_data,rna_data) {
  
  rand_transl_rates <- rep(0,ncol(pro_data))
  rand_transl_cols <- rep(0,ncol(pro_data))
  for (j in 1: ncol(pro_data)) {
    rsf <- rep(0,12)
    randCol <- floor(runif(1,min = 1, max = ncol(pro_data)))
    for(i in 1:12) { 
      rsf[i] <- (pro_data[i,j]/rna_data[i,randCol])
    }
    rand_transl_rates[j] <- median(rsf)
    rand_transl_cols[j] <- randCol
    
    #plot(rsf, xaxt = "n",col= "orange", ylab = "Ratio Protein/mRNA", xlab = "Tissue", main = paste("Median Protein/mRNA rate of ", rownames(subsetData)[j]))
    #axis(1, at = 1:12, labels = rownames(subsetDataPro))
    #abline(a = median(rsf), b = 0)
  }
  
  
  rand_predict_mat <- matrix(nrow = 12,ncol = ncol(pro_data))
  for (j in 1:ncol(pro_data)) {
    for (i in 1:12) {
      rand_predict_mat[i, j] <- rand_transl_rates[j] * rna_data[i,rand_transl_cols[j]]
    }
  }
  
  
  #plot
  #for (i in 1:ncol(pro_data)) {
  #  plot(pro_data[,i]~rna_data[,rand_transl_cols[i]],main = paste ("protein_gt = r_g * miRNA_gt ", colnames(subsetData)[i]))
  #  lines(x = rna_data[,rand_transl_cols[i]], y = rand_predict_mat[,i], col="red")
  #}
  return (rand_predict_mat)
}




#shared---------------------------------------------------------------------------


#what is shown on gene plots of shiny
correlateWithin <- function(pro_data,pro_predict) {
  within_corrs <- rep(0,ncol(pro_data))
  for (i in 1:ncol(pro_data)) {
    within_corrs[i] <- cor(pro_predict[,i], pro_data[,i], method = c("spearman"))
  }
  return (within_corrs)
}



correlateAcross <- function(pro_data,rna_data) {
  across_corrs <- rep(0,nrow(pro_data))
  for (i in 1:nrow(pro_data)) {
    across_corrs[i] <- cor(rna_data[i,], pro_data[i,], method = c("spearman"))
  }
  return (across_corrs)
}

orig_transl_rates <- estimateOriginalMethod(subsetDataPro, subsetDataRNA)
orig_predictions <- predictOriginalMethod(orig_transl_rates, subsetDataRNA)
orig_data_across_corrs <- correlateAcross(subsetDataPro,subsetDataRNA)
orig_predict_within_corrs <- correlateWithin(subsetDataPro,orig_predictions)
orig_predict_across_corrs <- correlateAcross(subsetDataPro,orig_predictions)


free_predictions <- estimateRNAFree(subsetDataPro)
#free_predict_within_corrs <- correlateWithin(subsetDataPro,free_predictions) #sdev is 0
free_predict_across_corrs <- correlateAcross(subsetDataPro, free_predictions)



rand_predictions <- estimateRandomMethod(subsetDataPro,subsetDataRNA)
rand_predict_within_corrs <- correlateWithin(subsetDataPro,rand_predictions)
rand_predict_across_corrs <- correlateAcross(subsetDataPro, rand_predictions)


#figureA-----------------------------------------------------------------------

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



plotFigureB <- function(pro_data,pro_predict) {
  df <- as.data.frame(cbind(pro_data,pro_predict))
  #11,4
  firstGroup <- 12
  secondGroup <- 11
  thirdGroup <- 5
  fourthGroup <- 4
  
  myplot <- ggplot() +
    geom_point(data = df, aes(x = df[,firstGroup], y = df[,ncol(pro_data)+firstGroup]), colour = 'orange', size = 2) +
    stat_ellipse(aes(x = df[,firstGroup], y = df[,ncol(pro_data)+firstGroup]),type = "norm", colour = 'orange') +
    
    geom_point(data = df, aes(x = df[,secondGroup], y = df[,ncol(pro_data)+secondGroup]), colour = 'red', size = 2) +
    stat_ellipse(aes(x = df[,secondGroup], y = df[,ncol(pro_data)+secondGroup]),type = "norm", colour = 'red') +
    
    geom_point(data = df, aes(x = df[,thirdGroup], y = df[,ncol(pro_data)+thirdGroup]), colour = '#4DB6D0', size = 2) +
    stat_ellipse(aes(x = df[,thirdGroup], y = df[,ncol(pro_data)+thirdGroup]),type = "norm", colour = '#4DB6D0') +
    
    geom_point(data = df, aes(x = df[,fourthGroup], y = df[,ncol(pro_data)+fourthGroup]), colour = 'mediumorchid1', size = 2) +
    stat_ellipse(aes(x = df[,fourthGroup], y = df[,ncol(pro_data)+fourthGroup]),type = "norm", colour = 'mediumorchid1') +
    
    labs(x = "Log2 Predicted Protein", y = "Log2 Observed Protein") +
    geom_abline(intercept = 0,slope = 1,size = 2, colour = "grey") + 
    
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
  return (myplot)
}
figB_orig_predict <- plotFigureB(subsetDataPro,orig_predictions)





#plot FigureC----------------------------------------------------------------------------------------------

#in: protein and rna matrix no n/a
#out: list of translation rates
predictRatesOriginal <- function(pro_mat, rna_mat) {
  ratesList <- rep(0,nrow(pro_mat))
  for (j in 1: nrow(pro_mat)) {
    rsf <- rep(0,12)
    for(i in 1:12) { 
      rsf[i] <- (pro_mat[j,i]/rna_mat[j,i])
    }
    ratesList[j] <- median(rsf)
  }
  return (ratesList)
}

#in: list of transl. rates, rna matrix no n/a
#out: matrix of predicted values
predictDataOriginal <- function (ratesList, rna_mat) {
  predict_mat <- matrix(nrow = nrow(rna_mat),ncol = 12)
  for (j in 1:nrow(rna_mat)) {
    for (i in 1:12) {
      predict_mat[j,i] <- ratesList[j] * rna_mat[j,i]
    }
  }
  return (predict_mat)
}


#in: protein matrix, rna matrix no n/a
#out: list of correlations per gene
correlateWithinFigureC <- function(pro_mat, rna_mat) {
  pro_mat_t <- t(pro_mat)
  predictedData <- t(predictDataOriginal(predictRatesOriginal(pro_mat,rna_mat),rna_mat))
  withinCorrelations <- rep(0,nrow(pro_mat))
  for (i in 1:nrow(pro_mat)) {
    withinCorrelations[i] <- cor(predictedData[,i], pro_mat_t[,i], method = c("spearman"))
  }
  return (withinCorrelations)
}


#in: availTissues in [1,12]
#out: plots figure C from the article
plotFigureC <- function(availTissues) {
  
  pro_df_f <- dplyr::filter(pro_df, rowSums(is.na(pro_df)) <= 12 - availTissues)
  rna_df_f <- rna_df[(rna_df$X %in% pro_df_f$X),]
  
  pro_mat <- pro_df_f[,-1]
  rownames(pro_mat) <- pro_df_f[,1]
  pro_mat <- as.matrix(log(pro_mat + 0.00000000001))
  rna_mat <- rna_df_f[,-1]
  rownames(rna_mat) <- rna_df_f[,1]
  rna_mat <- as.matrix(log(rna_mat + 0.00000000001))
  
  pro_mat_mice <- mice(pro_mat, m = 5, method = 'pmm')
  pro_mat <- complete(pro_mat_mice)
  rna_mat_mice <- mice(rna_mat, m = 5, method = 'pmm')
  rna_mat <- complete(rna_mat_mice)
  
  corrList <- correlateWithinFigureC(pro_mat, rna_mat)
  corrList <- as.data.frame(corrList)
  myplot <- ggplot(data=corrList, aes(corrList)) +
    geom_histogram(color = 'blue', fill = 'lightblue1',bins = 30, alpha = 0.5) +
    labs(x = "Correlation", y = "Count") +
    xlim(c(-1,1)) +
    theme_bw()
  
  return (myplot)
}
figC <- plotFigureC(12)

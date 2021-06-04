#navigate to folder with merged data sets
rna_df <- read.csv("https://raw.githubusercontent.com/gcohenfr/Neonatal-Vaccination/main/data/geneUN.csv?token=AMFFFE76ZJVXMSCF2JVQTD3AW3JV2")
pro_df <- read.csv("https://raw.githubusercontent.com/gcohenfr/Neonatal-Vaccination/main/data/proteinUN.csv?token=AMFFFE77VWMLYQG5QGQQZYLAX3TF4")
comb_df <- merge(rna_df,pro_df,by=c('X'),all.x=T, suffixes = c("_rna","_pro"))


rownames(comb_df) <- comb_df[,1]
comb_mat <- comb_df[,-1]



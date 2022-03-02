### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/');

library('tidyr');
library('dplyr');

### FIND CORRELATION COEFFICIENTS ##############################################

find.correlation.coefficients <- function(linear.models.folder) {
  
  ### LOAD DATA ##################################################################################
  
  #access folder with summary statistics
  linear.model.files <- list.files(path = file.path('outputs', 'statistics', linear.models.folder));
  
  #empty dataframe to hold correlation coefficients
  beta.coefficient.coefficient.dataframe <- NULL; 
  spearmans.rho.coefficient.dataframe <- NULL; 
  
  #linear.models.folder <- '2022-02-10-PRISM-all.features';
  #linear.model <- read.table('outputs/statistics/2022-02-10-PRISM-all.features/2022-02-10-NH.F1-PRISM-linear-model-statistics.txt')
  
  ### FORMAT DATA ################################################################################
  #read in those summary statistics
  for (file in 1:length(linear.model.files)){
    #read file
    linear.model <- read.table(file.path('outputs', 'statistics', linear.models.folder, linear.model.files[file]));
    
    #read summary statistics
    beta.coefficient.value <- linear.model$beta.coefficient;
    spearmans.rho.value <- linear.model$spearmans.rho;
    
    #rbind
    beta.coefficient.coefficient.dataframe <- rbind(beta.coefficient.coefficient.dataframe, beta.coefficient.value);
    spearmans.rho.coefficient.dataframe <- rbind(spearmans.rho.coefficient.dataframe, spearmans.rho.value);
    
    #set column name
    rownames(beta.coefficient.coefficient.dataframe)[file] <- linear.model$nuclear.feature[1];
    rownames(spearmans.rho.coefficient.dataframe)[file] <- linear.model$nuclear.feature[1];
    
  }
  
  beta.coefficient.coefficient.dataframe <- t(beta.coefficient.coefficient.dataframe);
  spearmans.rho.coefficient.dataframe <- t(spearmans.rho.coefficient.dataframe);
  
  #adjust row and column names
  rownames(beta.coefficient.coefficient.dataframe) <- linear.model$drug;
  rownames(spearmans.rho.coefficient.dataframe) <- linear.model$drug;
  
  #set as dataframe
  beta.coefficient.coefficient.dataframe <- as.data.frame(beta.coefficient.coefficient.dataframe);
  spearmans.rho.coefficient.dataframe <- as.data.frame(spearmans.rho.coefficient.dataframe);
  
  ### SAVE DATA ##################################################################################
  write.table(beta.coefficient.coefficient.dataframe, file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'beta.coefficients.txt', sep = '-')), sep = '\t');
  write.table(spearmans.rho.coefficient.dataframe, file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'spearmans.rho.txt', sep = '-')), sep = '\t');
  
}

### DATA ANALYSIS ################################################################################

find.correlation.coefficients('2022-02-27-PRISM-all.features');
find.correlation.coefficients('2022-02-26-AUC.CTD2-all.features');
find.correlation.coefficients('2022-02-26-AUC.GDSC1-all.features');
find.correlation.coefficients('2022-02-26-AUC.GDSC2-all.features');
find.correlation.coefficients('2022-02-26-AUC.PRISM-all.features');
find.correlation.coefficients('2022-02-26-dose.level.PRISM-all.features');
find.correlation.coefficients('2022-02-26-IC50-all.features');
find.correlation.coefficients('2022-02-26-replicate.CTD2-all.features');
find.correlation.coefficients('2022-02-27-replicate.GDSC1-all.features');
find.correlation.coefficients('2022-02-27-replicate.GDSC2-all.features');
find.correlation.coefficients('2022-02-27-replicate.PRISM-all.features');

### END ##########################################################################################
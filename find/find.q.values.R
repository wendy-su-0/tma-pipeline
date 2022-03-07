### DESCRIPTION ###################################################################################

# Analyze the results from the linear models. Find the models with significance.

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/tma-pipeline/');
library('tidyr');
library('dplyr');
library('BoutrosLab.plotting.general');

### FIND Q VALUES #################################################################################
find.q.values <- function(linear.models.folder) { 
  
  ### LOAD DATA ###################################################################################
  linear.model.files <- list.files(path = file.path('outputs', 'statistics', linear.models.folder));
  
  #create empty dataframes to hold summary statistics
  main.effect.q.values <- NULL;
  spearmans.rho.q.values <- NULL;
  
  ### FORMAT DATA #################################################################################
  #read in those summary statistics
  for (file in 1:length(linear.model.files)){
    
    #read file
    linear.model <- read.table(file.path('outputs', 'statistics', linear.models.folder, linear.model.files[file]));
    
    #read summary statistics
    main.effect.q.value <- linear.model$main.effect.q.value;
    spearmans.rho.q.value <- linear.model$spearmans.rho.q.value;
    
    #rbind
    main.effect.q.values <- rbind(main.effect.q.values, main.effect.q.value);
    spearmans.rho.q.values <- rbind(spearmans.rho.q.values, spearmans.rho.q.value);
    
    #set column name
    rownames(main.effect.q.values)[file] <- linear.model$nuclear.feature[1];
    rownames(spearmans.rho.q.values)[file] <- linear.model$nuclear.feature[1];
  }
  
  main.effect.q.values <- t(main.effect.q.values);
  spearmans.rho.q.values <- t(spearmans.rho.q.values);
  
  
  #adjust row and column names
  rownames(main.effect.q.values) <- linear.model$drug;
  rownames(spearmans.rho.q.values) <-  linear.model$drug;
  
  #set as dataframe
  main.effect.q.values <- as.data.frame(main.effect.q.values);
  spearmans.rho.q.values <- as.data.frame(spearmans.rho.q.values);
  
  ### SAVE DATA ###################################################################################
  write.table(main.effect.q.values, file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'main.effect.q.values.txt', sep = '-')), sep = '\t');
  write.table(spearmans.rho.q.values, file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'spearmans.rho.q.values.txt', sep = '-')), sep = '\t');
  
  
};

### DATA ANALYSIS #################################################################################

#find.q.values('2022-02-26-AUC.CTD2-all.features');

### END ###########################################################################################
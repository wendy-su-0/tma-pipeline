### DESCRIPTION ###################################################################################

# Analyze the results from the linear models. Find the models with significance.

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/');
library('tidyr');
library('dplyr');
library('BoutrosLab.plotting.general');

### FIND SIGNIFICANT MODELS #######################################################################

find.significant.models <- function(linear.models.folder) {
  
  ### LOAD DATA ###################################################################################
  main.effect.q.values <- read.table(file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'main.effect.q.values.txt', sep = '-')));
  spearmans.rho.q.values <- read.table(file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'spearmans.rho.q.values.txt', sep = '-')));
  
  
  ### FORMAT DATA #################################################################################
  #filter main.effect.q.values
  filtered.main.effect.q.values <- main.effect.q.values %>% filter_all(any_vars(. < 0.1));
  filtered.main.effect.q.values[filtered.main.effect.q.values > 0.1] <- NA;
  filtered.main.effect.q.values <- filtered.main.effect.q.values[,colSums(is.na(filtered.main.effect.q.values))<nrow(filtered.main.effect.q.values)];
  
  #filter spearmans.rho.q.values
  filtered.spearmans.rho.q.values <- spearmans.rho.q.values %>% filter_all(any_vars(. < 0.1));
  filtered.spearmans.rho.q.values[filtered.spearmans.rho.q.values > 0.1] <- NA;
  filtered.spearmans.rho.q.values <- filtered.spearmans.rho.q.values[,colSums(is.na(filtered.spearmans.rho.q.values))<nrow(filtered.spearmans.rho.q.values)];
  
  ### SAVE DATA ###################################################################################
  dir.create(file.path('outputs', 'analysis', linear.models.folder, 'significant'));
  
  write.table(filtered.main.effect.q.values, file.path('outputs', 'analysis', linear.models.folder, 'significant', paste(linear.models.folder, 'significant.main.effect.q.values.txt', sep = '-')), sep = '\t');
  write.table(filtered.spearmans.rho.q.values, file.path('outputs', 'analysis', linear.models.folder, 'significant', paste(linear.models.folder, 'significant.spearmans.rho.q.values.txt', sep = '-')), sep = '\t');
  
};

### DATA ANALYSIS #################################################################################

find.significant.models('2022-02-27-PRISM-all.features');
find.significant.models('2022-02-26-AUC.CTD2-all.features');
find.significant.models('2022-02-26-AUC.GDSC1-all.features');
find.significant.models('2022-02-26-AUC.GDSC2-all.features');
find.significant.models('2022-02-26-AUC.PRISM-all.features');
find.significant.models('2022-02-26-dose.level.PRISM-all.features');
find.significant.models('2022-02-26-IC50-all.features');
find.significant.models('2022-02-26-replicate.CTD2-all.features');
find.significant.models('2022-02-27-replicate.GDSC1-all.features');
find.significant.models('2022-02-27-replicate.GDSC2-all.features');
find.significant.models('2022-02-27-replicate.PRISM-all.features');

### END ###########################################################################################
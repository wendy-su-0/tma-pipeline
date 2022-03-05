### DESCRIPTOIN #####################################################

# iterates through the nuclear features of interest to make lin models of

### PREAMBLE ########################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/tma-pipeline/');

library(tidyr);
library(dplyr);
library(lme4);
library(sjstats);
library(effsize);
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);
source('generate/generate.linear.models.R');
#check necessary libraries

### WHICH FEATURES ##################################################

generate.selected.linear.models <- function(tma.medians.file.name, depmap.data.file.name, nuclear.feature.type, depmap.data.type){
  
  ### LOAD DATA #####################################################
  tma.medians.data <- read.table(file.path('outputs','data', tma.medians.file.name));
  
  depmap.data <- read.table(file.path('outputs', 'data', depmap.data.file.name));
  
  ### MAKE FOLDER FOR DATA ##########################################
  dir.create(file.path('outputs', 'statistics', paste(Sys.Date(), depmap.data.type, nuclear.feature.type, sep = '-')));
  dir.create(file.path('outputs', 'plots', paste(Sys.Date(), depmap.data.type, nuclear.feature.type, sep = '-')));
  dir.create(file.path('outputs', 'analysis', paste(Sys.Date(), depmap.data.type, nuclear.feature.type, sep = '-')));
  
  ### ITERATE THROUGH FEATURES ######################################
  if (nuclear.feature.type == 'all.features') {
    for(tma.col in 1:ncol(tma.medians.data)) {
      generate.linear.model(tma.col, tma.medians.data = tma.medians.data, depmap.data, nuclear.feature.type, depmap.data.type);
      #generate.linear.model(tma.col, depmap.data.type);
      cat(paste0('Done running model...', tma.col, ' out of ', ncol(tma.medians.data), '. ', '\n'));
    }
  }
  
};

### DATA ANALYSIS ##############################################################

select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-AUC.CTD2.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'AUC.CTD2');

### END ########################################################################
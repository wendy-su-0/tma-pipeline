### DESCRIPTION ###################################################################################

# Formats data from TMA. Removes rows w/o associated DepMap data
# Renames the TMA rows so it lines up with DepMap names. Saves the formatted data

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/tma-pipeline/');
library('tidyr');
library('dplyr');

#### FORMAT DATA #########################################################################################

format.TMA.data <- function(TMA.file.name, TMA.data.type) {
  ### LOAD DATA ###################################################################################
  tma.data <- read.csv(file.path('raw-data', 'nuclear-features', TMA.file.name), fileEncoding = 'UTF-8-BOM');
  
  ### REMOVE ROWS W/O DEPMAPDATA ##########################################################################
  tma.data <- subset(tma.data, tma.data$RegionID == 'Cell_lines_TMA_Core_CaoV3-R11' | 
                       tma.data$RegionID == 'Cell_lines_TMA_Core_IGROV-R21' |
                       tma.data$RegionID == 'Cell_lines_TMA_Core_OV90-R9' | 
                       tma.data$RegionID == 'Cell_lines_TMA_Core_OVCAR8-R8' |
                       tma.data$RegionID == 'Cell_lines_TMA_Core_TOV112D-R17' | 
                       tma.data$RegionID == 'Cell_lines_TMA_Core_TOV21G-R20');
  
  ### REMOVE COLUMN ##########################################################################
  tma.data = subset(tma.data, select = -c(CaseID));
  
  ### RENAME ROWS ###################################################################################
  tma.data[tma.data == 'Cell_lines_TMA_Core_CaoV3-R11'] <- 'CAOV3';
  tma.data[tma.data == 'Cell_lines_TMA_Core_IGROV-R21'] <- 'IGROV1';
  tma.data[tma.data == 'Cell_lines_TMA_Core_OV90-R9'] <- 'OV90';
  tma.data[tma.data == 'Cell_lines_TMA_Core_OVCAR8-R8'] <- 'OVCAR8';
  tma.data[tma.data == 'Cell_lines_TMA_Core_TOV112D-R17'] <- 'TOV112D';
  tma.data[tma.data == 'Cell_lines_TMA_Core_TOV21G-R20'] <- 'TOV21G';
  
  ### SAVE DATA ####################################################################################
  write.table(tma.data, file.path('outputs', 'data', paste(Sys.Date(), TMA.data.type, 'formatted.txt', sep = '-')), sep = '\t');
  
};

### DATA ANALYSIS ##################################################################################

format.TMA.data('Cell_lines_TMA_HE.csv','tma.data.six.ovarian');

### END ############################################################################################
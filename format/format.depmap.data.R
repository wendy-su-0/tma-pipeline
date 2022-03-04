### DESCRIPTION ###################################################################################

# Formats data downloaded from DepMap. Removes unnecessary columns: depmapid, cell 
# lineage information past cancer type. Saves the formatted data

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/tma-pipeline/');
library('tidyr');
library('dplyr');

### TITLE #########################################################################################

format.CCLE.data <- function(depmap.file.name, depmap.data.type) {
  
  ### LOAD DATA #################################################################################
  depmap.outputs <- read.csv(file.path('raw-data', 'depmap', depmap.file.name));
  
  depmap.sample.info <- read.csv(file.path('raw-data', 'depmap', 'sample_info.csv'));
  
  ### FORMAT DATA ###############################################################################
  depmap.data <- merge(
    x = depmap.sample.info[, c(1,22)],
    y = depmap.outputs,
    by.x = 'DepMap_ID',
    by.y = 'depmap_id',
    all.x = FALSE,
    all.y = TRUE
  );
  
  #sort alphabetically
  depmap.data <- depmap.data[order(depmap.data$cell_line_display_name), ];
  
  ### SAVE DATA #################################################################################
  write.table(depmap.data, file.path('outputs', 'data', paste(Sys.Date(), depmap.data.type, 'formatted.txt', sep = '-')), sep = '\t');

};

### DATA ANALYSIS #################################################################################

format.CCLE.data('Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4_subsetted.csv', 'PRISM.drug.sensitivity.six.ovarian', return.data = FALSE);

### END ###########################################################################################
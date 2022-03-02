### DESCRIPTION ###################################################################################

# Formats data downloaded from DepMap. Removes unnecessary columns: depmapid, cell 
# lineage information past cancer type. Saves the formatted data

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/');
library('tidyr');
library('dplyr');

### TITLE #########################################################################################

format.CCLE.data <- function(depmap.file.name, depmap.data.type, return.data = TRUE) {
  
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
  #date, depmap.data.type. number of cell lines?
  
  if(return.data){
    return(depmap.data);
  }
};

### DATA ANALYSIS #################################################################################

format.CCLE.data('Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4_subsetted.csv', 'PRISM.drug.sensitivity.six.ovarian', return.data = FALSE);

format.CCLE.data('Drug_sensitivity_AUC_(CTD^2)_subsetted_NAsdropped.csv', 'AUC.CTD2.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_19Q4_subsetted_NAsdropped.csv', 'AUC.PRISM.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_AUC_(Sanger_GDSC1)_subsetted_NAsdropped.csv', 'AUC.GDSC1.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted_NAsdropped.csv', 'AUC.GDSC2.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_dose-level_(PRISM_Repurposing_Secondary_Screen)_19Q4_subsetted.csv', 'dose.level.PRISM.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_IC50_(Sanger_GDSC1)_subsetted_NAsdropped.csv', 'IC50.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_replicate-level_dose_(CTD^2)_subsetted_NAsdropped.csv', 'replicate.CTD2.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_replicate-level_dose_(PRISM_Repurposing_Secondary_Screen)_19Q4_subsetted_NAsdropped.csv', 'replicate.PRISM.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_replicate-level_dose_(Sanger_GDSC1)_subsetted_NAsdropped.csv', 'replicate.GDSC1.drug.sensitivity.six.ovarian', return.data = FALSE);
format.CCLE.data('Drug_sensitivity_replicate-level_dose_(Sanger_GDSC2)_subsetted_NAsdropped.csv', 'replicate.GDSC2.drug.sensitivity.six.ovarian', return.data = FALSE);

### END ###########################################################################################

temp <- read.table('outputs/data/2022-01-25-AUC.CTD2.drug.sensitivity.six.ovarian-formatted.txt')

### DESCRIPTOIN #####################################################

# iterates through the nuclear features of interest to make lin models of

### PREAMBLE ########################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/');

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

select.nuclear.feature.models <- function(tma.data.file, depmap.data.file, nuclear.feature.type, depmap.data.type){
  
  ### LOAD DATA #####################################################
  tma.medians.data <- read.table(file.path('outputs','data', tma.data.file));
  
  depmap.data <- read.table(file.path('outputs', 'data', depmap.data.file));
  
  #plement read text file for run specified nuclear features
  
  #spread out from a few categories
  ##add txt file
  #val = vars
  #top 20 perc - nuclear feature IQR
  # if (nuclear.feature.type == 'all') {
  #   #do nothing
  # } else if (nuclear.feature.type = 'GLME'){
  #   #segment out those nuclear features
  #   #not the right numbers
  #   tma.data <- tma.data[, c(20:30, 70:90)]
  # }
  
  #indices of the tma data columns
  
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
  # 
  # else if (nuclear.feature.type = '#filepath to txt file') {
  #   #read the text file
  #   nuclear.feature.type = 'hematoxylin-features.txt'
  #   tomatch <- scan(read.table(file.path('raw-data/nuclear-features', nuclear.feature.type)), character(), quote = '')
  #   #this one
  #   tomatch <- read.table(file.path('raw-data/nuclear-features', nuclear.feature.type))
  #   tomatch
  #   #match the indices
  #   #nh.36 --> tma.medians column 5, save 5
  #   #ncol(index)
  #   #frerf in to lapply and parse
  #   
  #   #feed in feature name
  #   
  #   #helpr function to search for inxef
  #   
  #   #for(tma.col in 1:ncol(tma.medians.data)) {
  #     #generate.linear.model(tma.col, tma.medians.data = tma.medians.data, depmap.data.type);
  #     #cat(paste0('Done running model...', tma.col, ' out of ', ncol(tma.medians.data), '. ', '\n'));
  #   #}
  # }
  # 
  
}

#figure out cut dataset name out. date|dataset|drug.sens.
#could just delete drug sensitvy
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-AUC.CTD2.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'AUC.CTD2');

#out <- read.table('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/outputs/statistics/01-19 stats/01-19-22-NE.F8-drug-linear-model-statistics.txt');


select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-AUC.GDSC1.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'AUC.GDSC1');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-AUC.GDSC2.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'AUC.GDSC2');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-AUC.PRISM.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'AUC.PRISM');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-dose.level.PRISM.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'dose.level.PRISM');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-IC50.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'IC50');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-replicate.CTD2.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'replicate.CTD2');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-replicate.GDSC1.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'replicate.GDSC1');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-PRISM.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'PRISM');

#not all 6 cell lines
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-replicate.GDSC2.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'replicate.GDSC2');
select.nuclear.feature.models('2022-02-26-six.cell.lines-tma.medians.txt', '2022-02-26-replicate.PRISM.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'replicate.PRISM');
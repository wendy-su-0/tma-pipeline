### DESCRIPTION ###################################################################################

# Finds the median for all nuclear features
# Still need to edit

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/tma-pipeline/');
library('tidyr');
library('dplyr');

### TITLE #########################################################################################

format.tma.medians <- function(tma.data.filename, which.cell.lines) {
  
  ### LOAD DATA #################################################################################
  tma.data <- read.table(file.path('outputs', 'data', tma.data.filename));
  
  cell.lines <- unique(tma.data$RegionID);
  
  tma.medians <- NULL;
  
  ### FORMAT DATA ###############################################################################
  
  for (cell.line.index in 1:length(cell.lines)) {
    
    temp.tma.data <- filter(tma.data, tma.data$RegionID == cell.lines[cell.line.index]);
    temp.median.row <- NULL;
    
    for(col in 2:ncol(tma.data)) {
      
      median.nuclear.feature <- median(as.vector(temp.tma.data[ , col]));
      temp.median.row <-c(temp.median.row, median.nuclear.feature)
      
    }
    
    tma.medians <- rbind(tma.medians, temp.median.row);
    
  }
  
  tma.medians <- as.data.frame(tma.medians);
  rownames(tma.medians) <- cell.lines;
  colnames(tma.medians) <- colnames(tma.data)[-1];
  
  #sort alphabetically
  tma.medians <- tma.medians[ order(row.names(tma.medians)), ]
  
  ### SAVE DATA #################################################################################
  write.table(tma.medians, file.path('outputs', 'data', paste(Sys.Date(), which.cell.lines, 'tma.medians.txt', sep = '-')), sep = '\t');
  
};

### DATA ANALYSIS #################################################################################

generate.tma.medians('2022-02-26-tma.data.six.ovarian-formatted.txt', 'six.cell.lines', return.data = FALSE);

### END ###########################################################################################

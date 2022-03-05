### DESCRIPTION ###################################################################################

# Assesses the medians of the formatted TMA data. Assesses other summary statistics of the TMA data.
# Saves those medians. Saves summary statistics.

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/')
library(tidyr);
library(dplyr);
library(Hmisc);

### GENERATE HISTOGRAMS ############################################################################

generate.histograms <- function(tma.file.name, cell.lines) {
  
  ### READ DATA ####################################################################################
  tma.data.formatted <- read.table(file.path('outputs', 'data', tma.file.name));
  
  ### MAKE PDF #####################################################################################
  pdf(file = file.path('outputs', 'plots', paste(Sys.Date(), cell.lines, 'nuclear.feature.histograms.pdf', sep = '-')), onefile = TRUE);
  
  ### MAKE HISTOGRAMS ##############################################################################
  hist.data.frame(tma.data.formatted[c(-1)]);
  
  ### SAVE PDF #####################################################################################
  dev.off(); 
  
}

### GENERATE SUMMARY STATISTICS ####################################################################

generate.summary.statistics <- function(tma.file.name, depmap.file.name, cell.lines) {
  
  ### LOAD DATA ####################################################################################
  tma.data.formatted <- read.table(file.path('outputs', 'data', tma.file.name));
  depmap.data.formatted <- read.table(file.path('outputs', 'data', depmap.file.name));
  
  ### PREPARE MEDIAN DATAFRAME #####################################################################
  column.names <- colnames(tma.data.formatted)[2:length(tma.data.formatted)];
  cell.lines <- depmap.data.formatted$cell_line_display_name;
  tma.medians <- tma.data.formatted[cell.lines, column.names];
  rownames(tma.mediansns) <- cell.lines;
  
  ### PREPARE SUMMARY DATAFRAME ####################################################################
  tma.summary.statsitics <- NULL;
  
  ### ASSESS SUMMARY STATISTICS #####################################################################
  for (col in 2:ncol(tma.data.formatted)) {
    
    for (row in 1:length(cell.lines)) {
      
      cell.line <- filter(tma.data.formatted, RegionID == cell.lines[row]);
      cell.line.median <- median(as.vector(cell.line[, col]));
      tma.medians[row, col-1] <- cell.line.median;
      
    }
    
    tma.mean <- mean(as.vector(tma.data.formatted[, col]));
    tma.median <- median(as.vector(tma.data.formatted[, col]));
    tma.iqr <- IQR(as.vector(tma.data.formatted[, col]));
    
    temp.df <- c(tma.mean, tma.median, tma.iqr);
    
    tma.summary.statsitics <- rbind(tma.summary.statsitics, temp.df);
    
  }
  
  ### RENAME ROWS OF SUMMARY STATISTICS ############################################################
  colnames(tma.summary.statsitics) <-  c('mean', 'median', 'IQR');
  rownames(tma.summary.statsitics) <- colnames(tma.data.formatted[ , 2:length(tma.data.formatted)]);
  
  ### SAVE DATA ####################################################################################
  write.table(tma.summary.statsitics, file.path('outputs', 'statistics', paste(Sys.Date(), cell.lines, 'tma.summary.statsitics.txt', sep = '-')), sep = '\t');
  
};

### DATA ANALYSIS ##################################################################################



### END ############################################################################################
### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/tma-pipeline/');

library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);


### PRINT DOT MAP #################################################################################
print.corr.dotmap <- function(linear.models.folder, correlation.type) {
  
  ### LOAD DATA ###################################################################################
  
  q.values <- NULL;
  correlation.coefficients <- NULL;
  
  if (correlation.type == 'beta.coefficient'){
    
    q.values <- read.table( file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'main.effect.q.values.txt', sep = '-') ) );
    correlation.coefficients <- read.table( file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'beta.coefficients.txt', sep = '-') ) );
    
    
  } else if (correlation.type == 'spearmans.rho') {
    
    q.values <- read.table( file.path('outputs', 'analysis', linear.models.folder,  paste(linear.models.folder, 'spearmans.rho.q.values.txt', sep = '-') ) );
    correlation.coefficients <- read.table( file.path('outputs', 'analysis', linear.models.folder, paste(linear.models.folder, 'spearmans.rho.txt', sep = '-') ) );
    
  } else {
    #return an error
  }
  
  ### FILTER DATA ##################################################################################
  
  # filter the q value data so only significant models
  q.values.filtered <- q.values %>% 
    filter_all(any_vars(. < 0.1)); 
  
  q.values.filtered <- t(q.values.filtered);
  q.values.filtered <- as.data.frame(q.values.filtered);
  
  q.values.filtered <- q.values.filtered %>% 
    filter_all(any_vars(. < 0.1));
  
  q.values.filtered <- t(q.values.filtered)
  q.values.filtered <- as.data.frame(q.values.filtered)
  
  log.q.values.filtered <- -log10(q.values.filtered);
  
  #match the beta coefficient data rows and columns 
  correlation.coefficients.filtered <- filter(correlation.coefficients, rownames(correlation.coefficients) %in% rownames(q.values.filtered))
  correlation.coefficients.filtered <- select(correlation.coefficients.filtered, colnames(q.values.filtered))
  
  ### SET UP VARIABLES #############################################################################
  
  #zscore transform beta coefficients
  correlation.coefficients.filtered <- apply(correlation.coefficients.filtered, 1, zscore);
  correlation.coefficients.filtered <- as.data.frame(t(correlation.coefficients.filtered));
  
  beta.key <- seq(-2, 2, 0.5);
  
  spot.col <- function(x) { # specifies the colours for the dots +ve vs -ve vs 0
    colours <- rep('white', length(x));
    colours[x == 0]      <- 'transparent';
    colours[x > 0]       <- 'firebrick1';
    colours[x < 0]       <- 'dodgerblue';
    return(colours);
  }
  
  spot.size.function   <- function(x) { 0.1 + 1.5 * abs(x); } # Change the 1.5 number to make the dots visible
  
  ### PLOT DATA ####################################################################################
  
  significant.models.dotmap <- create.dotmap(
    x = correlation.coefficients.filtered,
    #xlab.label = 'Nuclear Features',
    ylab.label = 'Drugs',
    xlab.cex = 0.5,
    ylab.cex = 1.5,
    yaxis.cex = 0.75,
    xaxis.cex = 0.70,
    xaxis.tck = 0.5,
    yaxis.tck = 0.5,
    col.lwd = 0.25,
    row.lwd = 0.25,
    na.spot.size = 1,
    spot.colour.function = spot.col,
    spot.size.function = spot.size.function,
    bg.data = log.q.values.filtered,
    colour.scheme = c('white' , 'gray75', 'gray50', 'black'),
    at = c(0, 1,  2 , 3),
    
    colourkey = TRUE,
    colourkey.labels = c(1, 0.1,  0.01, 0.001),
    colourkey.labels.at = c(0, 1,  2 , 3),
    
    top.padding = 5,
    bottom.padding = 5,
    left.padding = 5,
    right.padding = 5,
    xaxis.rot = 45,
    use.legacy.settings = TRUE,
    
    key = list(
      space = 'right',
      points = list(
        cex = spot.size.function(beta.key),
        col = spot.col(beta.key),
        pch = 19
      ),
      text = list(
        lab = as.character(beta.key),
        cex = 1,
        adj = 1
      ),
      padding.text = 3,
      background = 'white'
    )
    
    
  );
  
  ### SAVE PLOT ##################################################################################  
  png(file = file.path('outputs', 'plots', linear.models.folder, paste(linear.models.folder, correlation.type, 'dotmap.png', sep = '-')), width = 1080); 
  #png(file = file.path('outputs', 'plots', paste(linear.models.folder, correlation.type, 'dotmap.png', sep = '-')), width = 1080, height = 1080); 
  print(significant.models.dotmap);
  dev.off();
  
}

### DATA ANALYSIS ################################################################################

#dir.create(file.path('outputs', 'plots', '2022-02-09-AUC.GDSC1-all.features'));
print.corr.dotmap('2022-02-26-AUC.CTD2-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-AUC.CTD2-all.features', 'spearmans.rho');

### END ##########################################################################################
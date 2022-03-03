# BIND ALL THE ADJUSTED P-VALUES OF ALL NUCLEAR FEATURES INTO A DATA FRAME
# SORT THE THING SO THAT ONLY ROWS AND COLUMNS WITH SIGNIFICANTS ARE KEPT
# PLOT IT

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/');

library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);


### PRINT DOT MAP #################################################################################
print.corr.dotmap <- function(linear.models.folder, correlation.type) {
  
  #linear.models.folder <- '2022-02-10-PRISM-all.features';
  #correlation.type <- 'spearmans.rho';
  
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
  
  #notes
  # 1 data frame.
  #nef12 beta
  #nef12 qvalue
  #grab all the column names
  # 2. x is all the column names
  # 3. 1 dataframe. x = the beta column names, bg.data = the qvalue column names
  
  # USE LOG P VALUES
  
  # strip the drug name
  
  # qualiy control with the chi square dotmap that is idential 
  
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
  #png(file = file.path('outputs', 'plots', linear.models.folder, paste(linear.models.folder, correlation.type, 'dotmap.png', sep = '-')), width = 1080); 
  png(file = file.path('outputs', 'plots', paste(linear.models.folder, correlation.type, 'dotmap.png', sep = '-')), width = 1080, height = 1080); 
  print(significant.models.dotmap);
  dev.off();
  
}

### DATA ANALYSIS ################################################################################

#dir.create(file.path('outputs', 'plots', '2022-02-10-PRISM-all.features'));
print.corr.dotmap('2022-02-27-PRISM-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-27-PRISM-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-09-AUC.GDSC1-all.features'));
print.corr.dotmap('2022-02-26-AUC.CTD2-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-AUC.CTD2-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-09-AUC.GDSC1-all.features'));
print.corr.dotmap('2022-02-26-AUC.GDSC1-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-AUC.GDSC1-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-09-AUC.GDSC2-all.features'));
print.corr.dotmap('2022-02-26-AUC.GDSC2-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-AUC.GDSC2-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-09-AUC.PRISM-all.features'));
print.corr.dotmap('2022-02-26-AUC.PRISM-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-AUC.PRISM-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-09-dose.level.PRISM-all.features'));
print.corr.dotmap('2022-02-26-dose.level.PRISM-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-dose.level.PRISM-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-09-IC50-all.features'));
print.corr.dotmap('2022-02-26-IC50-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-IC50-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-09-replicate.CTD2-all.features'));
print.corr.dotmap('2022-02-26-replicate.CTD2-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-26-replicate.CTD2-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-10-replicate.GDSC1-all.features'));
print.corr.dotmap('2022-02-27-replicate.GDSC1-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-27-replicate.GDSC1-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-10-replicate.GDSC2-all.features'));
print.corr.dotmap('2022-02-27-replicate.GDSC2-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-27-replicate.GDSC2-all.features', 'spearmans.rho');

#dir.create(file.path('outputs', 'plots', '2022-02-10-replicate.PRISM-all.features'));
print.corr.dotmap('2022-02-27-replicate.PRISM-all.features', 'beta.coefficient');
print.corr.dotmap('2022-02-27-replicate.PRISM-all.features', 'spearmans.rho');

### END ##########################################################################################

p.values <- read.table('2021-08-12_NM_p.values.txt'); #there's a code error going from p values to the adjusted
beta.coefficients <- read.table('2021-8-12_NM_beta.coefficients.txt')

create.dotmap(
  x = beta.coefficients
);

NF.corr.data <- data.frame(
  spearman.rho = beta.coefficients$NE.F16,
  pval = p.values$NE.F16
);

NF.beta.data <- data.frame(
  'NH.F12' = beta.coefficients$NH.F12,
  'NH.F23' = beta.coefficients$NH.F23,
  'NH.F26' = beta.coefficients$NH.F26,
  'NH.F41' = beta.coefficients$NH.F41,
  'NH.F43' = beta.coefficients$NH.F43,
  'NH.F48' = beta.coefficients$NH.F48,
  'NH.F49' = beta.coefficients$NH.F49,
  'NH.F51' = beta.coefficients$NH.F51,
  'NH.F52' = beta.coefficients$NH.F52,
  'NH.F54' = beta.coefficients$NH.F54,
  'NE.F10' = beta.coefficients$NE.F10,
  'NE.F14' = beta.coefficients$NE.F14,
  'NE.F16' = beta.coefficients$NE.F16,
  'NE.F21' = beta.coefficients$NE.F21,
  'NE.F23' = beta.coefficients$NE.F23,
  'NE.F26' = beta.coefficients$NE.F26,
  'NE.F28' = beta.coefficients$NE.F28,
  'NE.F29' = beta.coefficients$NE.F29,
  'NE.F41' = beta.coefficients$NE.F41,
  'NE.F48' = beta.coefficients$NE.F48,
  'NE.F49' = beta.coefficients$NE.F49,
  'NE.F52' = beta.coefficients$NE.F52,
  'NE.F53' = beta.coefficients$NE.F53,
  'NE.F54' = beta.coefficients$NE.F54,
  'NE.F57' = beta.coefficients$NE.F57
);

NF.pval.data <- data.frame(
  'NH.F12' = adjusted.p.values$NH.F12,
  'NH.F23' = adjusted.p.values$NH.F23,
  'NH.F26' = adjusted.p.values$NH.F26,
  'NH.F41' = adjusted.p.values$NH.F41,
  'NH.F43' = adjusted.p.values$NH.F43,
  'NH.F48' = adjusted.p.values$NH.F48,
  'NH.F49' = adjusted.p.values$NH.F49,
  'NH.F51' = adjusted.p.values$NH.F51,
  'NH.F52' = adjusted.p.values$NH.F52,
  'NH.F54' = adjusted.p.values$NH.F54,
  'NE.F10' = adjusted.p.values$NE.F10,
  'NE.F14' = adjusted.p.values$NE.F14,
  'NE.F16' = adjusted.p.values$NE.F16,
  'NE.F21' = adjusted.p.values$NE.F21,
  'NE.F23' = adjusted.p.values$NE.F23,
  'NE.F26' = adjusted.p.values$NE.F26,
  'NE.F28' = adjusted.p.values$NE.F28,
  'NE.F29' = adjusted.p.values$NE.F29,
  'NE.F41' = adjusted.p.values$NE.F41,
  'NE.F48' = adjusted.p.values$NE.F48,
  'NE.F49' = adjusted.p.values$NE.F49,
  'NE.F52' = adjusted.p.values$NE.F52,
  'NE.F53' = adjusted.p.values$NE.F53,
  'NE.F54' = adjusted.p.values$NE.F54,
  'NE.F57' = adjusted.p.values$NE.F57
);
?row.names

rm(tma.data.zscore)

row.names(NF.pval.data) <- row.names(adjusted.p.values)
row.names(NF.beta.data) <- row.names(beta.coefficients)

zscore.transform <- function() {
  
  ### LOAD DATA ####################################################################################
  tma.data <- read.table('2021-08-10_NM_tma.medians.txt');
  
  ### TRANSFORM DATA ###############################################################################
  NF.beta.data <- apply(NF.beta.data, 1, zscore);
  NF.beta.data <- as.data.frame(t(NF.beta.data));
  
  ### SAVE DATA ####################################################################################
  write.table(tma.data.zscore, '2021-08-10_NM_tma.medians.zscore.txt', sep = '\t');
  
}

NF.beta.data <- subset(NF.beta.data, row.names(NF.beta.data) %in% rownames(filtered.p.adjust))

NF.pval.data <- subset(NF.pval.data, row.names(NF.pval.data) %in% rownames(filtered.p.adjust))

NF.pval.log.data <- -log10(NF.pval.data);

spot.size.function   <- function(x) { 0.1 + 1.25*abs(x); }

drug.names <- c('Betamethasone', 'Enilconazole', 'Methocarbamol', 'Heptaminol', 'Colfosceril Palmitate',
                'Foretinib', 'D Delta Tocopherol', 'Ibuprofen', 'Etacrynic Acid', 'Asenapine')
drug.names <- t(drug.names)

feature.names <- c('H - Energy Intensity', 'H - Difference Variance', 'H - Information Correlation', 
                   'H - Low Gray-Level Zone Emphasis', 'H - Small Zone / Low Gray Emphasis',
                   'H - Zone Size Non-Uniformity', 'H - Zone Size Percentage', 'H - Long Run Emphasis',
                   'H - Gray Level Non-Uniformity', 'H - Run Length Non-Uniformity',
                   'E - Skewness', 'E - Angular Second Moment', 'E - Correlation', 'E - Sum Entropy',
                   'E- Difference Variance', 'E - Information Correlation', 'E - Autocorrelation', 
                   'E - Dissimilarity', 'E - Low Gray-Level Zone Emphasis', 'E - Zone Size Non-Uniformity',
                   'E - Zone Size Percentage', 'E - Gray Level Non-Uniformity', 'E - Run Percentage',
                   'E - Run Length Non-Uniformity', 'E - Chromatin Heterogenity' 
)

row.names(NF.pval.data) <- (drug.names)
row.names(NF.beta.data) <- (drug.names)
colnames(NF.pval.data) <- (feature.names)
colnames(NF.beta.data) <- (feature.names)

beta.key <- seq(-2, 2, 0.5);

spot.col <- function(x) { # specifies the colours for the dots +ve vs -ve vs 0
  colours <- rep('white', length(x));
  colours[x == 0]      <- 'transparent';
  colours[x > 0]       <- 'firebrick1';
  colours[x < 0]       <- 'dodgerblue';
  return(colours);
}

create.dotmap(
  x = NF.beta.data,
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
  bg.data = NF.pval.log.data,
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

?create.dotmap

##only sig models is done here


### LOAD DATA #################################################################################
NF.corr.data <- data.frame(
  spearman.rho = runif(n = 100, min = -1, max = 1),
  pval = runif(n = 100, min = 0.000000000000001, max = 0.2)
);

### ADJUST FOR MULTIPLE TESTING ###############################################################
NF.corr.data$FDR        <- p.adjust(NF.corr.data$pval, method = 'fdr');
NF.corr.data$log.FDR    <- -log10(NF.corr.data$FDR);

min(NF.corr.data$FDR)

### FORMAT PLOT DATA ##########################################################################
# Subset to show top 25
plot.data <- subset(NF.corr.data, NF.corr.data$log.FDR > 0.5); # Play around with this number til you get something that looks good
# Cluster correlation rho to get row order
row.order <- diana(plot.data$spearman.rho)$order; # Very important to double check that you have this correct for the labels

### CREATE DOT MAP #############################################################################
rho.key              <- seq(-1, 1, 0.2);

spot.col <- function(x) { # specifies the colours for the dots +ve vs -ve vs 0
  colours <- rep('white', length(x));
  colours[x == 0]      <- 'transparent';
  colours[x > 0]       <- 'firebrick1';
  colours[x < 0]       <- 'dodgerblue';
  return(colours);
}

border.colours.for.dotmap <- function(x) {
  colours <- rep('transparent', length(x));
  colours[abs(x) > 0]  <- 'black';
  colours[x == 0] <- 'transparent';
  return(colours);
}

spot.size.function   <- function(x) { 0.1 + 1.5 * abs(x); } # Change the 1.5 number to make the dots visible

create.dotmap(
  plot.data[row.order, 'spearman.rho'],
  bg.data = plot.data[row.order, 'log.FDR'],
  spot.size.function = spot.size.function,
  spot.colour.function = spot.col,
  yaxis.cex = 1,
  yaxis.lab = row.names(plot.data)[row.order], #If you have the NF labels available
  ylab.label = NULL,
  xaxis.lab = NULL,
  xaxis.rot = 90,
  xaxis.cex = 1,
  xaxis.tck = 0,
  ylab.cex = 1,
  xlab.cex = 1,
  top.padding = 5,
  bottom.padding = 5,
  # Changing background color scheme
  colour.scheme = c('white' , 'gray75', 'gray50', 'black'),
  at = c(0, 1,  2 , 3),
  colourkey.labels = c(1, 0.1,  0.01, 0.001),
  colourkey.labels.at = c(0, 1,  2 , 3),
  colourkey = FALSE,
  bg.alpha = 1,
  colourkey.cex = 1,
  na.pch = 4,
  na.spot.size = 1,
  resolution = 1500,
  yaxis.tck = 0,
  key = list(
    space = 'right',
    points = list(
      cex = spot.size.function(rho.key),
      col = spot.col(rho.key),
      pch = 19
    ),
    text = list(
      lab = as.character(rho.key),
      cex = 1,
      adj = 1
    ),
    padding.text = 3,
    background = 'white'
  )
);
}

print.corr.dotmap()
### END ###########################################################################################


print.dotmap <- function() {
  
  ###load data ###
  p.values <- read.table('2021-08-12_NM_adjusted.p.values.txt'); #there's a code error going from p values to the adjusted
  beta.coefficients <- read.table('2021-8-12_NM_beta.coefficients.txt');
  
  #adjust the pvalues
  for (i in 1:ncol(p.values)) {
    
    p.values[,i] <- p.adjust(p.values[,i], method = "fdr",  n = (nrow(p.values)));
    
  };
  
  #grab the rows you want
  NF.beta.data <- beta.coefficients %>% select(36:53, 89:106);
  NF.pval.data <- p.values %>% select(36:53, 89:106);
  
  #grab the significant models
  NF.sig.beta.data <- subset(NF.beta.data, row.names(NF.beta.data) %in% rownames(filtered.p.adjust))
  NF.sig.pval.data <- subset(NF.pval.data, row.names(NF.pval.data) %in% rownames(filtered.p.adjust))
  
  #grab the platin models
  NF.platin.beta.data <- NF.beta.data[c(2724, 3648, 3793, 3975),];
  NF.platin.pval.data <- NF.pval.data[c(2724, 3648, 3793, 3975),];
  
  NF.beta.data <- rbind(NF.sig.beta.data, NF.platin.beta.data);
  NF.beta.data <- apply(NF.beta.data, 1, zscore);
  NF.beta.data <- as.data.frame(t(NF.beta.data));
  
  NF.pval.data <- rbind(NF.sig.pval.data, NF.platin.pval.data);
  
  NF.pval.log.data <- -log10(NF.pval.data);
  
  spot.size.function   <- function(x) { 0.1 + 1.25*abs(x); }
  
  drug.names <- c('Betamethasone', 'Enilconazole', 'Methocarbamol', 'Heptaminol', 'Colfosceril Palmitate',
                  'Foretinib', 'D Delta Tocopherol', 'Ibuprofen', 'Etacrynic Acid', 'Asenapine', 'Satraplatin', "Cisplatin", 'Nedaplatin', 'Oxaliplatin')
  drug.names <- t(drug.names)
  row.names(NF.pval.data) <- (drug.names)
  row.names(NF.beta.data) <- (drug.names)
  
  beta.key <- seq(-2, 2, 0.5);
  
  spot.col <- function(x) { # specifies the colours for the dots +ve vs -ve vs 0
    colours <- rep('white', length(x));
    colours[x == 0]      <- 'transparent';
    colours[x > 0]       <- 'firebrick1';
    colours[x < 0]       <- 'dodgerblue';
    return(colours);
  }
  
  create.dotmap(
    x = NF.beta.data,
    xlab.label = 'Nuclear Features',
    ylab.label = 'Drugs',
    xlab.cex = 1.5,
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
    bg.data = NF.pval.log.data,
    colour.scheme = c('white' , 'gray75', 'gray50', 'black'),
    at = c(0, 1,  2 , 3),
    colourkey.labels = c(1, 0.1,  0.01, 0.001),
    colourkey.labels.at = c(0, 1,  2 , 3),
    top.padding = 5,
    bottom.padding = 5,
    left.padding = 5,
    right.padding = 5,
    xaxis.rot = 90,
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
}
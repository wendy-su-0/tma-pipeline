### DESCRIPTION ###################################################################################

# Print a heatmap of all the nuclear features

### PREAMBLE ######################################################################################

setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/linear-models/');

library('tidyr');
library('dplyr');
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);

### LIN MODELS ####################################################################################

print.heatmap.all <- function(date = '11-21-21', heatmap.covariates = c('Hematoxylin', 'Eosin')){
  
  ### LOAD DATA ###################################################################################
  
  tma.medians <- read.table('outputs/data/2021-08-10_NM_tma.medians.txt');
  
  all.heatmap<- cor(tma.medians)
  #text file going to the excel sheet
  #argument is subsection
  #filter tma medians by this text in the text file
  
  #I am manually dividing the features up into H/E using the colors rn
  subtype.cov <- rep(c('orange', 'chartreuse4'),times= c(60, 53))
  
  subtype.covariate <- list(
    rect = list(
      col = "black",
      fill = subtype.cov,
      lwd = 1.5
    )
  );
  
  #how do i automatically coerce the two colors to the two different features
  #default.colours
  subtype.cov.legend <- list(
    legend = list(
      colours =  default.colours(12),
      labels = c('Hematoxylin', 'Eosin'),
      title = 'Nuclear Feature Stain '
    )
  );
  
  tma.heatmap <- create.heatmap(
    all.heatmap,
    clustering.method = 'none',
    cluster.dimensions = 'both',
    grid.col = TRUE,
    grid.row = TRUE,
    plot.dendrograms = TRUE,
    print.colour.key = TRUE,
    use.legacy.settings = TRUE,
    
    colourkey.cex = 1,
    xaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.cex = 1,
    yaxis.fontface = 1,
    
    covariates = subtype.covariate,
    covariate.legends = subtype.cov.legend,
    legend.side = "right",
    legend.cex = 1,
    covariates.grid.border = list(col = 'black', lwd = 1)
  );
  
  ### SAVE HEATMAP #################################################################
  png(file= 'outputs/plots/11-21-21-heatmap.all.png');
  tma.heatmap;
  dev.off();
  
}

### DATA ANALYSIS ############################################################

print.heatmap.all();

# using boutrous.plotting.general, how do i functionize the colors of this

### END ###################################################################

### SET THE FOLDER YOUR FILES ARE ##############################################
setwd('/Users/wsu31/OneDrive/Desktop/R Stuff/tma-pipeline')

### R PACKAGES AND SCRIPTS #####################################################

  ### INSTALL PACKAGES: ONLY NEED TO DO ONCE ###################################
  install.packages('tidyr');
  install.packages('dplyr');
  install.packages('Hmisc');
  install.packages('lme4');
  install.packages('lsr');
  install.packages('insight');
  install.packages('sjstats');
  install.packages('effsize');
  install.packages(dependencies/BoutrosLab.utilities_1.9.10.tar.gz);
  install.packages('BoutrosLab.plotting.general');
  install.packages(dependencies/BoutrosLab.statistics.general_2.1.3.tar.gz);

  ### LOAD PACKAGES: EVERY R SESSION ###########################################
  library(tidyr);
  library(dplyr);
  library(Hmisc);
  library(lme4);
  library(lsr);
  library(insight);
  library(sjstats);
  library(effsize);
  library(BoutrosLab.utilities);
  library(BoutrosLab.plotting.general);
  library(BoutrosLab.statistics.general);

  ### LOAD SCRIPTS: EVERY R SESSION ############################################
  source('format/format.depmap.data.R');
  source('format/format.tma.data.R');
  source('format/format.tma.medians.R');
  source('generate/generate.linear.model.R');
  source('generate/generate.distributions.R');
  source('generate/generate.selected.linear.models.R');
  source('find/find.q.values.R');
  source('find/find.correlation.coefficients.R');
  source('find/find.significant.models.R');
  source('print/print.dotmap.R');

### PREPROCESS THE DATA ########################################################

  ### FORMAT NUCLEAR FEATURE DATA ##############################################
  format.depmap.data('Drug_sensitivity_AUC_(CTD^2)_subsetted_NAsdropped.csv', 'AUC.CTD2.drug.sensitivity.six.ovarian');
  #format.depmap.data(depmap.file.name, depmap.data.type)  

  ### FORMAT DEPMAP DATA #######################################################
  format.tma.data('Cell_lines_TMA_HE.csv','tma.data.six.ovarian');
  #format.tma.data(tma.file.name, tma.data.type)

  ### GENERATE NUCLEAR FEATURE MEDIANS #########################################
  format.tma.medians('2022-03-07-tma.data.six.ovarian-formatted.txt', 'six.cell.lines');
  #format.tma.median(tma.data.filename, which.cell.lines)


### GENERATE DISTRIBUTIONS #####################################################

  ### GENERATE HISTOGRAMS ######################################################
  generate.histograms('2022-03-07-tma.data.six.ovarian-formatted.txt', 'six.cell.lines');
  #generate.histograms(tma.file.name, cell.lines)


  ### GENERATE SUMMARY STATISTICS ##############################################
  generate.summary.statistics('2022-03-07-tma.data.six.ovarian-formatted.txt', '2022-03-07-AUC.CTD2.drug.sensitivity.six.ovarian-formatted.txt', 'six.cell.lines');
  #generate.summary.statistics(tma.file.name, depmap.file.name, cell.lines);

### GENERATE LINEAR MODELS #####################################################
generate.selected.linear.models('2022-03-07-six.cell.lines-tma.medians.txt', '2022-03-07-AUC.CTD2.drug.sensitivity.six.ovarian-formatted.txt', 'all.features', 'AUC.CTD2');
#generate.selected.linear.models(tma.medians.file.name, depmap.data.file.name, nuclear.feature.type, depmap.data.type)


### FIND KEY STATISTICS ########################################################

  ### FIND CORRELATION COEFFICIENTS ############################################
  find.correlation.coefficients('2022-03-07-AUC.CTD2-all.features');
  #find.correlation.coefficients(linear.models.folder);
  
  ### FIND Q VALUES ############################################################
  find.q.values('2022-03-07-AUC.CTD2-all.features');
  #find.q.values(linear.models.folder);
                          
  ### FIND SIGNIFICANT MODELS ##################################################
  find.significant.models('2022-03-07-AUC.CTD2-all.features');
  #find.significant.models(linear.models.folder);

### PRINT DOTMAP ###############################################################
print.corr.dotmap('2022-03-07-AUC.CTD2-all.features', 'beta.coefficient');
#print.corr.dotmap(linear.models.folder, correlation.type);

### END ########################################################################

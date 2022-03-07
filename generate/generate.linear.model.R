### DESCRIPTOIN #####################################################

# generate all the linear models for depmap data ~ 1 nucler feature
# full/null models with age as a fixed effect
# separate script iterates through the nuclear features of interest

### PREAMBLE ########################################################

library(tidyr);
library(dplyr);
library(lme4);
library(lsr);
library(insight);
library(sjstats);
library(effsize);
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);

### LIN MODELS ######################################################
generate.linear.model <- function(nuclear.feature.index, tma.medians.data, depmap.data, nuclear.feature.type, depmap.data.type) {
  
  nuclear.feature.name <- colnames(tma.medians.data[nuclear.feature.index]);
  
  if (nrow(tma.medians.data) != nrow(depmap.data)) {
    tma.medians.data <- subset(tma.medians.data, (rownames(tma.medians.data) %in% depmap.data$cell_line_display_name))
    
  }
  
  tma.data <- tma.medians.data[ ,nuclear.feature.index];
  
  depmap.data <- depmap.data;
  ### FORMAT DATA ###################################################
  
  result.summary.stats <- NULL;
  
  ### CREATE LIN MODELS #############################################
  
  #for loop. loop through the drug sensitivity. apply to one singular nuclear feature
  for(depmapCol in 8:ncol(depmap.data)) {
    
    # If there is NULL results for drug of interest, mark result summary data frame with NA
    #if (sum(depmap.data[, depmapCol], na.rm = TRUE) == 0) 
    if (any(is.na(depmap.data[, depmapCol]))) {
      
      temp.summary.stats <- data.frame(
        nuclear.feature = nuclear.feature.name,
        drug = colnames(depmap.data)[depmapCol],
        sample.size = ncol(depmap.data),
        main.effect.p.value = NA,
        main.effect.q.value = NA,
        main.effect.log.10.fdr = NA,
        chi.squared.p.value = NA,
        chi.squared.q.value = NA,
        chi.squared.log.10.fdr = NA,
        beta.coefficient = NA,
        spearmans.rho = NA,
        spearmans.rho.p.value = NA,
        spearmans.rho.q.value = NA,
        spearmans.rho.log.10.fdr = NA,
        eta.coefficient = NA,
        cohens.d = NA,
        variance = NA,
        confidence.intervals = NA
      );
      
      result.summary.stats <- rbind(result.summary.stats, temp.summary.stats);
      next
    }
    
    # If there are data for the gene of interest, run the LMEM model
    #else if (sum(depmap.data[, depmapCol], na.rm = TRUE) != 0) 
    else if (any(!is.na(depmap.data[, depmapCol]))){
      
      #create full and null models
      # dependent ~ independent
      
      # see the as fomula difference 
      AB.LMEM <-  lm(depmap.data[ ,depmapCol] ~ tma.data + depmap.data$age, REML = FALSE);
      A.LMEM <-  lm(depmap.data[ ,depmapCol] ~ depmap.data$age, REML = FALSE);
      # AB.LMEM <-  lm(depmap.data[ ,depmapCol] ~ tma.data + depmap.data$age);
      # A.LMEM <-  lm(depmap.data[ ,depmapCol] ~ depmap.data$age);
      
      #for chisquare, we want rmel to be false
      #for p-values, we awnt rmel to be true
      
      ### STATISTICS ####################################################
      
      # Assess p-value
      #is it just age? that is the driver
      LM.chi.squared.p.value <- as.data.frame(anova(AB.LMEM, A.LMEM))$`Pr(>F)`[2];
      
      ###maybe uncomment
      AB.LMEM <-  lm(depmap.data[ ,depmapCol] ~ tma.data + depmap.data$age, REML = TRUE);
      
      # Assess AB p value 
      #does the model fit
      LM.main.effect.p.value <- summary(AB.LMEM)$coefficients[2,4];
      
      # Assess beta coefficient
      LM.beta.coefficient <- summary(AB.LMEM)$coefficients[2,1];
      
      # Assess spearman's rho
      LM.spearmans.rho.test <- cor.test(x = tma.data, y = depmap.data[, depmapCol], method = 'spearman');
      LM.spearmans.rho <- LM.spearmans.rho.test[["estimate"]][["rho"]];
      
      # Assess spearman's rho p value
      LM.spearmans.rho.p.value <- LM.spearmans.rho.test$p.value;
      
      # Assess eta.coefficient
      LM.eta.coefficient <- etaSquared(AB.LMEM)[1,1];
      
      # Assess cohen's d
      LM.cohens.d <- cohen.d(as.numeric(tma.data), as.numeric(depmap.data[, depmapCol]))$sd;
      
      # Assess variance
      #between.cluster.var <- get_variance_intercept(AB.LMEM);
      #within.cluster.var  <- get_variance_residual(AB.LMEM);
      #LM.variance <- between.cluster.var / (between.cluster.var + within.cluster.var);
      LM.variance <- (summary(AB.LMEM)$sigma)**2;
      
      # Assess confidence intervals - double check. beta coef estimates. summary(ab.lmer)
      # use the stnardard error
      #LM.confidence.intervals <- summary(AB.LMEM)[['sigma']];
      LM.confidence.intervals <-(summary(AB.LMEM)$coefficients[2,2])**2;
      
      #add p value of ab model, 
      #bind the model's statistics to the results df
      temp.summary.stats <- data.frame(
        nuclear.feature = nuclear.feature.name,
        drug = colnames(depmap.data)[depmapCol],
        sample.size = ncol(depmap.data),
        main.effect.p.value = LM.main.effect.p.value,
        main.effect.q.value = NA,
        main.effect.log.10.fdr = NA,
        chi.squared.p.value = LM.chi.squared.p.value,
        chi.squared.q.value = NA,
        chi.squared.log.10.fdr = NA,
        beta.coefficient = LM.beta.coefficient,
        spearmans.rho = LM.spearmans.rho,
        spearmans.rho.p.value = LM.spearmans.rho.p.value,
        spearmans.rho.q.value = NA,
        spearmans.rho.log.10.fdr = NA,
        eta.coefficient = LM.eta.coefficient,
        cohens.d = LM.cohens.d,
        variance = LM.variance,
        confidence.intervals = LM.confidence.intervals
      );
      
      result.summary.stats <- rbind(result.summary.stats, temp.summary.stats);
      
    };
  };
  
  ###ADJUST P VAL ###################################################
  #main effect
  result.summary.stats$main.effect.q.value <- p.adjust(result.summary.stats$main.effect.p.value, method = 'fdr')
  result.summary.stats$main.effect.log.10.fdr <- -log(result.summary.stats$main.effect.q.value, base = 10)
  
  
  #chi square
  result.summary.stats$chi.squared.q.value <- p.adjust(result.summary.stats$chi.squared.p.value, method = 'fdr')
  result.summary.stats$chi.squared.log.10.fdr <- -log(result.summary.stats$chi.squared.q.value, base = 10)
  
  
  #adjust spearman's rho p values
  result.summary.stats$spearmans.rho.q.value <- p.adjust(result.summary.stats$spearmans.rho.p.value, method = 'fdr')
  result.summary.stats$spearmans.rho.log.10.fdr <- -log(result.summary.stats$spearmans.rho.q.value, base = 10)
  
  ### SAVE DATA #####################################################
  write.table(result.summary.stats, file.path('outputs', 'statistics', paste(Sys.Date(), depmap.data.type, nuclear.feature.type, sep = '-'), paste(Sys.Date(), nuclear.feature.name, depmap.data.type, 'linear-model-statistics.txt', sep = '-')), sep = '\t');
  
  #file.path('outputs/statistics', paste(Sys.Date(), depmap.data.type, nuclear.feature.type, sep = '-'))
};


### DATA ANALYSIS ###################################################

### END #############################################################
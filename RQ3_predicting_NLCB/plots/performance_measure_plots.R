###########################
### prepare environment ### 
###########################

# clear environment
rm(list = ls())

# load packages
library(ggplot2) # create plots
library(RColorBrewer) # use nice colors
library(reshape2) # manipulate dataframes
library(dplyr)  # manipulate dataframes
library(stringr)  # manipulate dataframes

# set global parameters
t0 <- 4
minAUC <- 0.4
maxAUC <- 0.9
minBrier <- 0.04
maxBrier <- 0.3

# function to summarize the C-indices
CI.func <- function(CI.list){
  ## calculate summarizing quantities
  CI_mean <- as.data.frame(simplify2array(lapply(CI.list, function(x) mean(as.matrix(x[,6]), na.rm = T))))
  CI_SD <- as.data.frame(simplify2array(lapply(CI.list, function(x) mean(as.matrix(x[,7]), na.rm = T))))
  CI_NA <- as.data.frame(simplify2array(lapply(CI.list, function(x) sum(is.na(as.matrix(x[,1:5]))))))
  ## clean up dataframes
  names(CI_mean) <- 'CI_mean'
  names(CI_SD) <- 'CI_SD'
  names(CI_NA) <- 'CI_nrmissing'
  CI_mean$model <- CI_SD$model <- CI_NA$model <- str_sub(rownames(CI_mean), 8, nchar(rownames(CI_mean))-5)
  rownames(CI_mean) <- rownames(CI_SD) <- rownames(CI_NA) <- 1:nrow(CI_mean)
  ## combine dataframes
  CI_sum <- left_join(CI_mean, CI_SD, by = 'model') %>%
    left_join(CI_NA, by = 'model') %>%
    mutate(subtype = str_sub(gsub('^.*?_', '_', model), 2, 100),
           model.int = model,
           model = NA,
           model = ifelse(grepl('pencal', model.int), 'pencal', model),
           model = ifelse(grepl('SPJM', model.int), 'JM', model),
           model = ifelse(grepl('strictlandmarking', model.int), 'landmarking', model),
           CI_mean = round(CI_mean, 2),
           CI_SD = round(CI_SD, 2),
           CI = paste0(CI_mean, ' (', CI_SD, ')')) %>%
    select(model, subtype, CI, CI_mean, CI_SD, CI_nrmissing)  %>%
    filter(!(subtype == 'LMM')) %>%
    filter(!(grepl('lmark', subtype))) %>%
    group_by(model) %>%
    mutate(best_CI = ifelse(CI_mean == max(CI_mean), 'yes', 'no'),
           best_CI = ifelse(model == 'pencal' & subtype %in% c('enet', 'lasso'), 'no', best_CI)) %>%
    ungroup()
  
  return(CI_sum)
}

#####################################################
### load data: only PSA as longitudinal predictor ### 
#####################################################

# only PSA as longitudinal predictor
## load dataframes
derived_data <- list.files(path = "Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/onevariable/",
                           pattern = ".csv$", full.names = T)
list2env(
  lapply(setNames(derived_data, 
                  make.names(gsub("\\.csv$", "", 
                                  gsub("\\Z..Documents.Scripts.RQ3_prediction.saved_dataframes.onevariable.", "", derived_data)))), 
         read.csv), envir = .GlobalEnv)
rm(derived_data)

## create dataframe with all the predictive performance values
### AUC and Brier score
pp_onepred <- left_join(pp_pencal_elasticnet_20x5, pp_pencal_lasso_20x5, by = 'times') %>%
  left_join(pp_pencal_ridge_20x5, by = 'times') %>%
  left_join(pp_SPJM_lmm_slope_20x5, by = 'times') %>%
  left_join(pp_SPJM_lmm_tv_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm1_slope_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm1_tv_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm2_slope_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm2_tv_20x5, by = 'times') %>%
  left_join(pp_strictlandmarking_LOCF_20x5, by = 'times')
names(pp_onepred) <- c('times', 'auc_elasticnet', 'brier_elasticnet',
                       'auc_lasso', 'brier_lasso',
                       'auc_ridge', 'brier_ridge',
                       'auc_lmm_slope', 'brier_lmm_slope',
                       'auc_lmm_tv', 'brier_lmm_tv',
                       'auc_mm1_slope', 'brier_mm1_slope', 
                       'auc_mm1_tv', 'brier_mm1_tv',
                       'auc_mm2_slope', 'brier_mm2_slope',
                       'auc_mm2_tv', 'brier_mm2_tv',
                       'auc_landmarking_LOCF', 'brier_landmarking_LOCF')
### C-index
CI.list <- ls()[grepl('cindex', ls())]
CI.list <- mget(CI.list)
CI_onepred <- CI.func(CI.list)

# clean environment 
rm(list = setdiff(ls(), c('pp_onepred', 'CI_onepred',
                          't0', 'minAUC', 'maxAUC', 'minBrier', 'maxBrier', 'CI.func')))

# multiple longitudinal predictors ---------------------------------------------
## load dataframes
derived_data <- list.files(path = "Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/mult_longvariable/",
                           pattern = ".csv$", full.names = T)
list2env(
  lapply(setNames(derived_data, 
                  make.names(gsub("\\.csv$", "", 
                                  gsub("\\Z..Documents.Scripts.RQ3_prediction.saved_dataframes.mult_longvariable.", "", derived_data)))), 
         read.csv), envir = .GlobalEnv)
rm(derived_data)

## create dataframe with all the predictive performance values
### AUC and Brier score
pp_multpred <- left_join(pp_pencal_enet_20x5, pp_pencal_lasso_20x5, by = 'times') %>%
  left_join(pp_pencal_ridge_20x5, by = 'times') %>%
  left_join(pp_SPJM_lmm_slope_20x5, by = 'times') %>%
  left_join(pp_SPJM_lmm_tv_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm1_slope_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm1_tv_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm2_slope_20x5, by = 'times') %>%
  left_join(pp_SPJM_mm2_tv_20x5, by = 'times') %>%
  left_join(pp_strictlandmarking_LOCF_20x5, by = 'times')
names(pp_multpred) <- c('times', 'auc_elasticnet', 'brier_elasticnet',
                        'auc_lasso', 'brier_lasso',
                        'auc_ridge', 'brier_ridge',
                        'auc_lmm_slope', 'brier_lmm_slope',
                        'auc_lmm_tv', 'brier_lmm_tv',
                        'auc_mm1_slope', 'brier_mm1_slope', 
                        'auc_mm1_tv', 'brier_mm1_tv',
                        'auc_mm2_slope', 'brier_mm2_slope',
                        'auc_mm2_tv', 'brier_mm2_tv',
                        'auc_landmarking_LOCF', 'brier_landmarking_LOCF')
### C-index
CI.list <- ls()[grepl('cindex', ls())]
CI.list <- mget(CI.list)
CI_multpred <- CI.func(CI.list)

# clean environment 
rm(list = setdiff(ls(), c('pp_onepred', 'CI_onepred', 'pp_multpred', 'CI_multpred',
                          't0', 'minAUC', 'maxAUC', 'minBrier', 'maxBrier', 'CI.func')))

# multiple baseline predictors -------------------------------------------------
## load dataframes
derived_data <- list.files(path = "Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/mult_baseline/",
                           pattern = ".csv$", full.names = T)
list2env(
  lapply(setNames(derived_data, 
                  make.names(gsub("\\.csv$", "", 
                                  gsub("\\Z..Documents.Scripts.RQ3_prediction.saved_dataframes.mult_baseline.", "", derived_data)))), 
         read.csv), envir = .GlobalEnv)
rm(derived_data)

## create dataframe with all the predictive performance values
### AUC and Brier score
pp_multbase <- pp_pencal_ridge_20x5 %>% 
  left_join(pp_SPJM_mm1_slope_20x5, by = 'times') %>%
  left_join(pp_strictlandmarking_LOCF_20x5, by = 'times')
names(pp_multbase) <- c('times', 
                        'auc_ridge', 'brier_ridge',
                        'auc_mm1_slope', 'brier_mm1_slope',
                        'auc_landmarking_LOCF', 'brier_landmarking_LOCF')
### C-index
CI.list <- ls()[grepl('cindex', ls())]
CI.list <- mget(CI.list)
CI_multbase <- CI.func(CI.list)

## clean environment 
rm(list = setdiff(ls(), c('pp_onepred', 'CI_onepred', 'pp_multpred', 'CI_multpred', 'pp_multbase', 'CI_multbase',
                          't0', 'minAUC', 'maxAUC', 'minBrier', 'maxBrier')))

##########################################
### seperate AUC and Brier information ### 
##########################################

# functions to clean up the AUC and Brier score dataframes
## AUC
AUC.func <- function(df, CI_df){
  melt(df[-1,], id.vars = "times") %>%
    mutate(variable = gsub("auc_", "", variable)) %>%
    mutate(model = NA,
           model = ifelse(variable %in% c('lasso', 'elasticnet', 'ridge'), 'pencal', model),
           model = ifelse(variable %in% c('lmm_slope', 'lmm_tv', 'mm1_slope', 'mm1_tv', 'mm2_slope', 'mm2_tv'), 'JM', model),
           model = ifelse(grepl('landmarking', variable), 'landmarking', model),
           variable = ifelse(variable == 'elasticnet', 'enet', variable),
           variable = ifelse(grepl('landmarking', variable), str_sub(variable, 13, 100), variable),
           model_subtype = paste0(model, '_', variable)) %>%
    rename(subtype = variable,
           auc = value) %>%
    left_join(select(CI_df, model, subtype, best_CI), by = c('model', 'subtype')) %>%
    select(times, auc, model, subtype, model_subtype, best_CI)
}
## Brier 
brier.func <- function(df, CI_df){
  melt(df[-1,], id.vars = "times")%>%
    mutate(variable = gsub("brier_", "", variable)) %>%
    mutate(model = NA,
           model = ifelse(variable %in% c('lasso', 'elasticnet', 'ridge'), 'pencal', model),
           model = ifelse(variable %in% c('lmm_slope', 'lmm_tv', 'mm1_slope', 'mm1_tv', 'mm2_slope', 'mm2_tv'), 'JM', model),
           model = ifelse(grepl('landmarking', variable), 'landmarking', model),
           variable = ifelse(variable == 'elasticnet', 'enet', variable),
           variable = ifelse(grepl('landmarking', variable), str_sub(variable, 13, 100), variable),
           model_subtype = paste0(model, '_', variable)) %>%
    rename(subtype = variable,
           brier = value) %>%
    left_join(select(CI_df, model, subtype, best_CI), by = c('model', 'subtype')) %>%
    select(times, brier, model, subtype, model_subtype, best_CI)
}

# create dataframes for the AUC and Brier score separately
## wide format
auc_onepred <- pp_onepred[,c(1,seq(2,20,2))] # one longitudinal predictor
auc_multpred <- pp_multpred[,c(1,seq(2,20,2))] # multiple longitudinal predictors
auc_multbase <- pp_multbase[,c(1,seq(2,6,2))] # multiple baseline predictors
brier_onepred <- pp_onepred[,c(1,seq(3,21,2))] # one longitudinal predictor
brier_multpred <- pp_multpred[,c(1,seq(3,21,2))] # multiple longitudinal predictors
brier_multbase <- pp_multbase[,c(1,seq(3,7,2))]# multiple baseline predictors
## long format
auc_onepred_long <- AUC.func(auc_onepred, CI_onepred) # one longitudinal predictor
auc_multpred_long <- AUC.func(auc_multpred, CI_multpred) # multiple longitudinal predictors
auc_multbase_long <- AUC.func(auc_multbase, CI_multbase)# multiple baseline predictors
brier_onepred_long <-  brier.func(brier_onepred, CI_onepred) # one longitudinal predictor
brier_multpred_long <-  brier.func(brier_multpred, CI_multpred) # multiple longitudinal predictors
brier_multbase_long <-  brier.func(brier_multbase, CI_multbase)# multiple baseline predictors

# combine data 
## AUC
auc_onepred_long$pred <- 'uni'
auc_multpred_long$pred <- 'mult_long'
auc_multbase_long$pred <- 'mult_base'
auc_long <- rbind(auc_onepred_long, auc_multpred_long, auc_multbase_long) 
## Brier
brier_onepred_long$pred <- 'uni'
brier_multpred_long$pred <- 'mult_long'
brier_multbase_long$pred <- 'mult_base'
brier_long <- rbind(brier_onepred_long, brier_multpred_long, brier_multbase_long) 
# C-index
CI_onepred$pred <- 'uni'
CI_multpred$pred <- 'mult_long'
CI_multbase$pred <- 'mult_basic'
CI <- rbind(CI_onepred, CI_multpred, CI_multbase)

# clean environment 
rm(list = setdiff(ls(), c('t0', 'minAUC', 'maxAUC', 'minBrier', 'maxBrier',
                          'auc_long', 'brier_long', 'CI')))

#######################
### plot perfomance ### 
#######################

# time-varying AUC's
#for(i in c('uni', 'mult_long')){
#  for(j in unique(auc_long$model)){
#    print(ggplot(auc_long[auc_long$pred == i & auc_long$model == j,], 
#                 aes(x = times, y = auc, color = subtype)) +
#            geom_line() +
#            labs(title = "AUC", x = "Months since baseline", y = "AUC") +
#            theme_minimal())
#  }
#}
#for(i in c('uni', 'mult_long')){
#  for(j in unique(brier_long$model)){
#    print(ggplot(brier_long[brier_long$pred == i &  brier_long$model == j,], 
#                 aes(x = times, y = brier, color = subtype)) +
#            geom_line() +
#            labs(title = "Brier score", x = "Months since baseline", y = "Brier score") +
#            theme_minimal())
#  }
#}

nr.preds <- unique(auc_long$pred)

# AUC
ggs.auc <- vector(mode = 'list', length(nr.preds))
for(i in 1:length(nr.preds)){
  ggs.auc[[i]] <- ggplot(auc_long[auc_long$pred == nr.preds[i],], 
                         aes(x = times, y = auc, group = model_subtype,
                             color = model, alpha = best_CI, linetype = model)) +
    geom_line(aes(size = model)) +
    scale_alpha_manual(labels = c('yes' = 'Optimal', 'no' = 'Non-optimal'),
                       values = c('yes' = 1, 'no' = 0.3)) +
    labs(x = "Months since baseline", y = "AUC") +
    #ggtitle('AUC') +
    ylim(c(minAUC, maxAUC)) +
    geom_vline(xintercept = t0, linetype = "dashed", color = "grey75") +
    theme_minimal() +
    theme(legend.position = 'none',
          axis.text=element_text(size=20),
          axis.title=element_text(size=18)) +
    scale_color_manual(labels=c('Landmarking', 'PRC', 'JM'),
                       values = brewer.pal(name ="Paired", n = 9)[c(8,2,4)]) +
    scale_linetype_manual(labels=c('Landmarking', 'PRC', 'JM'),
                          values = c('landmarking' = 'dotdash', 'pencal' = 'solid', 'JM' = 'longdash')) +
    scale_size_manual(labels=c('Landmarking', 'PRC', 'JM'),
                      values = c('landmarking' = 1.5, 'pencal' = 0.8, 'JM' = 0.8)) +
    guides(color = guide_legend(
      override.aes = list(linetype = c('dotdash', 'solid', 'longdash'),
                          color = brewer.pal(name ="Paired", n = 9)[c(2,4,8)],
                          linewidth = c(1, 0.7, 0.7),
                          size = c(3,3,3))),
      linetype = 'none',
      size = 'none',
      alpha = 'none') +
    labs(color = 'Model type')
  
  print(ggs.auc[[i]])
}

# Brier
ggs.brier <- vector(mode = 'list', length(nr.preds))
for(i in 1:length(nr.preds)){
  ggs.brier[[i]] <- ggplot(brier_long[brier_long$pred == nr.preds[i],], 
                           aes(x = times, y = brier, group = model_subtype, 
                               color = model, alpha = best_CI, linetype = model)) +
    geom_line(aes(size = model)) +
    scale_alpha_manual(labels = c('yes' = 'Optimal', 'no' = 'Non-optimal'),
                       values = c('yes' = 1, 'no' = 0.25)) +
    labs(x = "Months since baseline", y = "Brier score") +
    # ggtitle('Brier score') +
    ylim(c(minBrier, maxBrier)) +
    geom_vline(xintercept = t0, linetype = "dashed", color = "grey75") +
    theme_minimal() +
    theme(legend.position = 'none',
          axis.text=element_text(size=20),
          axis.title=element_text(size=18)) +
    scale_color_manual(labels=c('Landmarking', 'PRC', 'JM'),
                       values = brewer.pal(name ="Paired", n = 9)[c(8,2,4)]) +
    scale_linetype_manual(labels=c('Landmarking', 'PRC', 'JM'),
                          values = c('landmarking' = 'dotdash', 'pencal' = 'solid', 'JM' = 'longdash')) +
    scale_size_manual(labels=c('Landmarking', 'PRC', 'JM'),
                      values = c('landmarking' = 1.5, 'pencal' = 0.8, 'JM' = 0.8)) +
    guides(color = guide_legend(
      override.aes = list(
        linetype = c('dotdash', 'solid', 'longdash'),
        color = brewer.pal(name ="Paired", n = 9)[c(2,4,8)],
        linewidth = c(1, 0.7, 0.7),
        size = c(3,3,3))),
      linetype = 'none',
      size = 'none',
      alpha = 'none') +
    labs(color = 'Model type')
  
  print(ggs.brier[[i]])
}

###################
### export data ### 
###################

# plots
plots <- list(ggs.auc[[1]], ggs.auc[[2]], ggs.auc[[3]],
              ggs.brier[[1]], ggs.brier[[2]], ggs.brier[[3]]) 

filenames <- c(paste0('auc_', nr.preds),
               paste0('brier_', nr.preds))

for (i in 1:length(plots)){  
  file_name = paste('Z:/Documents/Scripts/RQ3_prediction/saved_plots/original/', filenames[i], '.pdf', sep='')
  pdf(file_name, height = 5, width = 8)
  print(plots[[i]])
  dev.off()
}

# dataframe CI
write.csv(CI, 'Z:/Documents/Scripts/RQ3_prediction/saved_plots/CI_original.csv', row.names = FALSE)

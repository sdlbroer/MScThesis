###########################
### prepare environment ###
###########################

# load pre-processed data
source('Z:/Documents/Scripts/RQ3_data_management.R')

# load libraries 
library(nlme) # fit the mixed model
library(survival) # fit the survival model
library(JMbayes2) # fit the joint model
library(splines) # use splines
library(riskRegression) # calculate prediction performance
library(survcomp) # calculate C-index

#################################
### prepare for model fitting ### 
#################################

# set hyperparameters
t0 <- 4 # moment of prediction 
tn <- 10 # maximum prediction time
times <- seq(t0, tn, length.out = (tn-t0)*2+1) # number of predictions 
n.folds <- 5 # number of folds in the cross-validation
n.RCV <- 20 # number of repeats of the cross-validation
## specify model-specific hyperparameters
time.model <- 'lmm' # lmm/mm1
assoc.param <- 'slope' # slope/tv
if(time.model == 'lmm' & assoc.param == 'slope') {seed <- floor(1803158/116) ; seeds <- seed:(seed + n.RCV)} # set seeds
if(time.model == 'lmm' & assoc.param == 'tv') {seed <- floor(1803158/117) ; seeds <- seed:(seed + n.RCV)} # set seeds
if(time.model == 'mm1' & assoc.param == 'slope') {seed <- floor(1803158/118) ; seeds <- seed:(seed + n.RCV)} # set seeds
if(time.model == 'mm1' & assoc.param == 'tv') {seed <- floor(1803158/119) ; seeds <- seed:(seed + n.RCV)} # set seeds

# add baseline date to longitudinal measurements 
longdata <- long_meas_train %>%
  filter(!is.na(PSA)) %>%
  filter(!is.na(LDH)) %>%
  filter(!is.na(ALP)) %>%
  left_join(select(long_meas_train[long_meas_train$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId')  %>%
  select(id_num, therapy_received, fup_time = time, time_obs, NLCB, 
         log2PSA, LDH, ALP) %>%
  mutate(LDH = log(LDH + 0.01, 2), ALP = log(ALP + 0.01, 2),
         therapy_received = as.factor(as.character(therapy_received))) %>%
  rename(log2LDH = LDH, log2ALP = ALP) 

# create dataframe with only the survival data 
survdata <- baseline_data_table1_train[baseline_data_table1_train$id_num %in% longdata$id_num,] %>%
  distinct(id_num, therapy_received, NLCB_overall_num, time_obs, age,
           Gleason, LocationMetastases, ecog, treatmentline, RT, ST) %>%
  # filter(treatmentline != 3) %>%
  mutate(treatmentline = as.character(treatmentline),
         therapy_received = as.factor(as.character(therapy_received)),
         Gleason = as.factor(as.character(Gleason)),
         LocationMetastases = as.factor(as.character(LocationMetastases)),
         ecog = as.factor(as.character(ecog)),
         treatmentline = as.factor(as.character(treatmentline))) 
#longdata <- longdata[longdata$id_num %in% survdata$id_num,]

# perform strict landmarking 
long_strict_landmarked <- longdata[longdata$fup_time <= t0 & longdata$time_obs > t0,]
surv_strict_landmarked <- survdata[survdata$id_num %in% long_strict_landmarked$id_num,]

# time forms 
if(time.model == 'lmm'){ # linear
  fixed.formula.PSA <- as.formula('log2PSA ~ fup_time*therapy_received')
  fixed.formula.LDH <- as.formula('log2LDH ~ fup_time*therapy_received')
  fixed.formula.ALP <- as.formula('log2ALP ~ fup_time*therapy_received')
  random.formula <- as.formula('~ fup_time | id_num')
}
if(time.model == 'mm1'){ # splines with 1 internal knot
  fixed.formula.PSA <- as.formula('log2PSA ~ ns(fup_time, k = c(2), B = c(0, 3.8))*therapy_received')
  fixed.formula.LDH <- as.formula('log2LDH ~ ns(fup_time, k = c(2), B = c(0, 3.8))*therapy_received')
  fixed.formula.ALP <- as.formula('log2ALP ~ ns(fup_time, k = c(2), B = c(0, 3.8))*therapy_received')
  random.formula <- as.formula('~ ns(fup_time, k = c(2), B = c(0, 3.8)) | id_num')
}

# functional forms
if(assoc.param == 'slope') def.form <- ~ value(log2PSA) + slope(log2PSA, eps = 1, direction = 'back') +
  value(log2LDH) + slope(log2LDH, eps = 1, direction = 'back') +
  value(log2ALP) + slope(log2ALP, eps = 1, direction = 'back') # slope
if(assoc.param == 'tv') def.form <- ~ value(log2PSA) * ns(fup_time, k = c(2), B = c(0, 3.8)) +
  value(log2LDH) * ns(fup_time, k = c(2), B = c(0, 3.8)) +
  value(log2ALP) * ns(fup_time, k = c(2), B = c(0, 3.8)) # time-varying 

##########################
### dynamic prediction ###
##########################

# m-times k-fold cross-validation
#all_auc <- all_brier <- vector('list', n.RCV)
auc_reps <- brier_reps <- as.data.frame(times)
all_cindex <- data.frame(matrix(NA, nrow = n.RCV, ncol = n.folds))
names(all_cindex) <- paste0('fold', 1:n.folds)

for(m in 1:n.RCV){
  # print progress
  print(paste('Starting repition', m))
  
  # set seed 
  set.seed(seeds[m])
  
  # create k random folds
  ## repeated measurements  
  folds_long <- create_folds(longdata, V = n.folds, id_var = 'id_num', seed = seeds[m])
  ## survival data
  folds_surv <- list('training' = vector('list', n.folds), 
                     'testing' = vector('list', n.folds))
  for(i in 1:n.folds){
    ids <- unique(folds_long$training[[i]]$id_num)
    folds_surv$training[[i]] <- survdata[survdata$id_num %in% ids,]
    folds_surv$testing[[i]] <- survdata[!(survdata$id_num %in% ids),]
  }
  
  # perform prediction k times
  auc_df <- brier_df <- as.data.frame(times)
  
  for(i in 1:n.folds){
    # print progress
    print(paste0('Repetition ', m, ', fold ', i))
    
    # train and test set
    ## train (including imputation of missing values)
    df_train_long <- folds_long$training[[i]]
    df_train_surv <- folds_surv$training[[i]]
    df_train <- left_join(df_train_long, distinct(df_train_surv, id_num, NLCB_overall_num), by = 'id_num')
    ## test
    df_test_long <- folds_long$testing[[i]]
    df_test_surv <- folds_surv$testing[[i]] 
    df_test <- left_join(df_test_long, distinct(df_test_surv, id_num, NLCB_overall_num), by = 'id_num') %>%
      filter(fup_time <= t0) %>% # keep only observations up until t0 in the test set
      mutate(NLCB_overall_num = 0, time_obs = t0)
    df_test_surv <- df_test_surv[df_test_surv$id_num %in% df_test$id_num,]
    
    # fit models 
    ## survival submodel
    train_surv <- coxph(Surv(time_obs, NLCB_overall_num) ~ therapy_received, 
                        data = df_train_surv)
    ## longitudinal submodel: time modelled linearly
    train_longit_PSA <- lme(fixed = fixed.formula.PSA, 
                            random = random.formula,
                            data = df_train_long,
                            control = lmeControl(opt = 'optim'))
    train_longit_LDH <- try(lme(fixed = fixed.formula.LDH, 
                                random = random.formula,
                                data = df_train_long,
                                control = lmeControl(opt = 'optim')), silent = TRUE)
    if(inherits(train_longit_LDH, "try-error")) print(paste0('For repetition ', m, ', fold ', i, ': the MEM for LDH did not converge'))
    train_longit_ALP <- try(lme(fixed = fixed.formula.ALP, 
                                random = random.formula,
                                data = df_train_long,
                                control = lmeControl(opt = 'optim')) , silent = TRUE)
    if(inherits(train_longit_ALP, "try-error")) print(paste0('For repetition ', m, ', fold ', i, ': the MEM for ALP did not converge'))
    
    if(!inherits(train_longit_LDH, "try-error") & !inherits(train_longit_ALP, "try-error")){
      ## joint model, time-varying
      train_joint <- jm(train_surv, list(train_longit_PSA, train_longit_LDH, train_longit_ALP), 
                        time_var = "fup_time", functional_forms = assoc.param)
      
      # validate model on test set
      ## predict survival outcome
      fail_prob <- predict(train_joint, newdata = df_test, 
                           process = "event", times = times,
                           return_newdata = TRUE) %>%
        select(id_num, fup_time, pred_CIF) %>%
        reshape(idvar = 'id_num', timevar = 'fup_time', direction = "wide")
      fail_prob <- fail_prob[,2:ncol(fail_prob)]
      colnames(fail_prob) <- times
      
      # predictive performance
      acc_measures <- Score(as.list(fail_prob), 
                            formula = Surv(time_obs, NLCB_overall_num) ~ 1, 
                            data = df_test_surv,
                            times = times, 
                            cens.model = 'km',
                            metrics = c('auc','brier'), 
                            conf.int = FALSE, 
                            exact = FALSE, 
                            split.method	= 'none', 
                            B = 0)
      ## extract AUC
      auc <- acc_measures$AUC$score
      auc <- auc[auc$model == auc$times,-1]
      ## extract Brier
      brier <- acc_measures$Brier$score
      brier <- brier[brier$model == brier$times,-1]
      ## safe accuracy measures
      auc_df <- left_join(auc_df, auc, by = 'times')
      brier_df <- left_join(brier_df, brier, by = 'times')
      
      # C-index
      c_index <- concordance.index(x = fail_prob[,2], method = 'noether',
                                   surv.time = df_test_surv$time_obs, 
                                   surv.event = df_test_surv$NLCB_overall_num)$c.index
      ## save all C-indexes
      all_cindex[m,i] <- c_index
    }
  }
  # summarize predictive performance
  ## calculate mean performance
  auc_df$mean_auc <- apply(auc_df[,-1], 1, mean, na.rm = T)
  brier_df$mean_brier <- apply(brier_df[,-1], 1, mean, na.rm = T)
  ## safe predictive performance of each fold
  #all_auc[[m]] <- auc_df[,-1]
  #all_brier[[m]] <- brier_df[,-1]
  ## safe mean predictive performance
  auc_reps <- left_join(auc_reps, auc_df[,c('times', 'mean_auc')], by = 'times')
  brier_reps <- left_join(brier_reps, brier_df[,c('times', 'mean_brier')], by = 'times')
  
  # clean up environment
  rm(auc_df, brier_df, auc, brier, acc_measures, fail_prob,
     df_train_long, df_train_surv, df_test_long, df_test_surv, df_test, 
     train_surv, train_longit_PSA, train_longit_LDH, train_longit_ALP, train_joint)
}

# summarize predictive performance
## change row/column names
names(auc_reps) <- c('times', paste0(rep('repetition', n.RCV), 1:n.RCV))
names(brier_reps) <- c('times', paste0(rep('repetition', n.RCV), 1:n.RCV))
## calculate mean predictive performance of all RCV
auc_reps$mean_auc <- apply(auc_reps[,-1], 1, mean, na.rm = T)
brier_reps$mean_brier <- apply(brier_reps[,-1], 1, mean, na.rm = T)
## create dataframe with the predictive performance
pp <- data.frame(times = times,
                 auc = auc_reps$mean_auc,
                 brier = brier_reps$mean_brier)
## mean cindex + sd
all_cindex$mean_cindex <- apply(all_cindex, 1, mean, na.rm = T)
all_cindex$sd_cindex <- apply(all_cindex[,-ncol(all_cindex)], 1, sd, na.rm = T)
## export dataframe
write.csv(pp, paste0('Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/multlong_pp_SPJM_lmark_', time.model,'_', assoc.param, '_20x5.csv'), row.names = FALSE)
write.csv(all_cindex, paste0('Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/multlong_cindex_SPJM_lmark_', time.model,'_', assoc.param, '_20x5.csv'), row.names = FALSE)

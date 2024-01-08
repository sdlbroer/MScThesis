###########################
### prepare environment ###
###########################

# load pre-processed data
source('Z:/Documents/Scripts/RQ3_data_management.R')

# load libraries 
library(survival) # fit the survival model
library(JMbayes2) # create_folds function
library(riskRegression) # calculate prediction performance
library(pec) # predict survival probability based on Cox PH model
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
seed <- floor(1803158/102) # set seed 
seeds <- seed:(seed + n.RCV) # create new seed for each repetition

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
  print(paste('Starting repetition', m))
  
  # set seed 
  set.seed(seeds[m])
  
  # create k random folds
  folds_surv <- create_folds(surv_strict_landmarked, V = n.folds, id_var = 'id_num', seed = seeds[m])
  ## survival data
  folds_long <- list('training' = vector('list', n.folds), 
                     'testing' = vector('list', n.folds))
  for(i in 1:n.folds){
    ids <- unique(folds_surv$training[[i]]$id_num)
    folds_long$training[[i]] <- long_strict_landmarked[long_strict_landmarked$id_num %in% ids,]
    folds_long$testing[[i]] <- long_strict_landmarked[!(long_strict_landmarked$id_num %in% ids),]
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
    ## test (including removal of missing values)
    df_test_long <- folds_long$testing[[i]]
    df_test_surv <- folds_surv$testing[[i]]
    
    # fit landmark model (Cox PH)
    ## fit LME's
    model_lme_PSA <- lme(fixed = log2PSA ~ fup_time*therapy_received, 
                         random = ~ fup_time | id_num,
                         data = df_train_long,
                         na.action = na.omit,
                         control = lmeControl(opt = 'optim'))
    model_lme_ALP <- try(lme(fixed = log2ALP ~ fup_time*therapy_received, 
                             random = ~ fup_time | id_num,
                             data = df_train_long,
                             na.action = na.omit,
                             control = lmeControl(opt = 'optim')), silent = TRUE)
    model_lme_LDH <- try(lme(fixed = log2LDH ~ fup_time*therapy_received, 
                             random = ~ fup_time | id_num,
                             data = df_train_long,
                             na.action = na.omit,
                             control = lmeControl(opt = 'optim')), silent = TRUE)
    if(inherits(model_lme_ALP, "try-error")) print(paste0('For repetition ', m, ', fold ', i, ': the MEM for LDH did not converge'))
    if(inherits(model_lme_LDH, "try-error")) print(paste0('For repetition ', m, ', fold ', i, ': the MEM for ALP did not converge'))
    if(!inherits(model_lme_ALP, "try-error") & !inherits(model_lme_LDH, "try-error")){
      ## extract predicted value at landmark time
      df_train_surv <- cbind(df_train_surv,
                             predict(model_lme_PSA,
                                     newdata = data.frame(df_train_surv[,c('id_num', 'therapy_received')], fup_time = t0))) 
      names(df_train_surv)[12] <- 'predPSA'
      df_train_surv <- cbind(df_train_surv,
                             predict(model_lme_ALP,
                                     newdata = data.frame(df_train_surv[,c('id_num', 'therapy_received')], fup_time = t0))) 
      names(df_train_surv)[13] <- 'predALP'
      df_train_surv <- cbind(df_train_surv,
                             predict(model_lme_LDH,
                                     newdata = data.frame(df_train_surv[,c('id_num', 'therapy_received')], fup_time = t0))) 
      names(df_train_surv)[14] <- 'predLDH'
      ## fit Cox model 
      model_cox <- coxph(Surv(time_obs, NLCB_overall_num) ~ therapy_received + predPSA + predALP + predLDH, 
                         data = df_train_surv, x = TRUE)
      
      # validate model on test set 
      ## predict longitudinal outcome (new patients don't have random effects)
      df_test_surv <- cbind(df_test_surv,
                            predict(model_lme_PSA,
                                    newdata = data.frame(df_test_surv[,c('id_num', 'therapy_received')], fup_time = t0),
                                    level = 0)) 
      names(df_test_surv)[12] <- 'predPSA'
      df_test_surv <- cbind(df_test_surv,
                            predict(model_lme_ALP,
                                    newdata = data.frame(df_test_surv[,c('id_num', 'therapy_received')], fup_time = t0),
                                    level = 0)) 
      names(df_test_surv)[13] <- 'predALP'
      df_test_surv <- cbind(df_test_surv,
                            predict(model_lme_LDH,
                                    newdata = data.frame(df_test_surv[,c('id_num', 'therapy_received')], fup_time = t0),
                                    level = 0)) 
      names(df_test_surv)[14] <- 'predLDH'
      ## predict survival outcome
      fail_prob <- as.data.frame(1 - predictSurvProb(model_cox, newdata = df_test_surv, times = times))
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
     df_train_long, df_train_surv, df_test_long, df_test_surv)
}

# summarize predictive performance
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
write.csv(pp, 'Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/mult_longvariable/pp_strictlandmarking_LMM_20x5.csv', row.names = FALSE)
write.csv(all_cindex, 'Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/mult_longvariable/cindex_strictlandmarking_LMM_20x5.csv', row.names = FALSE)

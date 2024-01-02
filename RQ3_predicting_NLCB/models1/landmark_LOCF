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
seed <- floor(1803158) # set seed 
seeds <- seed:(seed + n.RCV) # create new seed for each repetition

# add baseline date to psa_long dataframe
longdata <- psa_long_train %>%
  filter(!is.na(PSA)) %>% 
  left_join(select(psa_long[psa_long$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId')  %>%
  select(id_num, therapy_received, fup_time = time, log2PSA, NLCB, time_obs)

# create dataframe with only the survival data 
survdata <- baseline_data_table1[baseline_data_table1$id_num %in% longdata$id_num,] %>%
  distinct(id_num, therapy_received, NLCB_overall_num, time_obs)

# perform strict landmarking 
long_strict_landmarked <- longdata[longdata$time_obs > t0,]
long_strict_landmarked <- long_strict_landmarked[long_strict_landmarked$fup_time <= t0,] 

surv_strict_landmarked <- survdata[survdata$id_num %in% long_strict_landmarked$id_num,] 

# make time-varying information LOCF
long_strict_landmarked <- long_strict_landmarked %>%
  arrange(id_num, desc(fup_time)) %>%
  filter(!duplicated(id_num))

surv_strict_landmarked <- left_join(surv_strict_landmarked,
                                    select(long_strict_landmarked, id_num, log2PSA),
                                    by = 'id_num')

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
  
  # perform prediction k times
  auc_df <- brier_df <- as.data.frame(times)
  
  for(i in 1:n.folds){
    # print progress
    print(paste0('Repetition ', m, ', fold ', i))
    
    # train and test set
    df_train_surv <- folds_surv$training[[i]]
    df_test_surv <- folds_surv$testing[[i]]
    
    # fit landmark model (Cox PH)
    model_landmark_LOCF <- coxph(Surv(time_obs, NLCB_overall_num) ~ therapy_received + log2PSA, 
                                 data = df_train_surv, x = TRUE)
    
    # validate model on test set 
    ## predict survival outcome
    fail_prob <- as.data.frame(1 - predictSurvProb(model_landmark_LOCF, newdata = df_test_surv[,c('therapy_received', 'log2PSA')], times = times))
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
  # summarize predictive performance
  ## change row/column names
  names(auc_df) <- c('times', paste0(rep('fold', n.folds), 1:n.folds))
  rownames(auc_df) <- auc_df[,1]
  names(brier_df) <- c('times', paste0(rep('fold', n.folds), 1:n.folds))
  rownames(brier_df) <- brier_df[,1]
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
     df_train_surv, df_test_surv)
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
write.csv(pp, 'Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/pp_strictlandmarking_LOCF_20x5.csv', row.names = FALSE)
write.csv(all_cindex, 'Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/cindex_strictlandmarking_LOCF_20x5.csv', row.names = FALSE)

###########################
### prepare environment ###
###########################

# load pre-processed data
source('Z:/Documents/Scripts/RQ3_data_management.R')

# load libraries 
library(survival) # fit the survival model
library(JMbayes2) # fit the joint model
library(pencal) # fit PRC model
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
seed <- floor(1803158/4) # set seed 
seeds <- seed:(seed + n.RCV) # create new seed for each repetition

pen <- 'ridge' # penalty to be used in the penalized regression, options: ridge, lasso, enet

# add baseline date to psa_long dataframe
longdata <- psa_long_train %>%
  filter(!is.na(PSA)) %>% 
  left_join(select(psa_long[psa_long$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId') %>%
  distinct(id = id_num, fup_time = time, log2PSA, therapy_received)

# create dataframe with only the survival data 
survdata <- baseline_data_table1[baseline_data_table1$id_num %in% long_meas_train$id_num,] %>%
  distinct(id = id_num, time = time_obs, event = NLCB_overall_num, therapy_received)

# perform strict landmarking 
survdata <- survdata[survdata$time > t0,] %>% 
  mutate(therapy_received = as.factor(as.character(therapy_received)))
longdata <- longdata[longdata$id %in% survdata$id,]
longdata <- longdata[longdata$fup_time <= t0,] %>% 
  mutate(therapy_received = as.factor(as.character(therapy_received)))

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
  folds_long <- create_folds(longdata, V = n.folds, id_var = 'id', seed = seeds[m])
  ## survival data
  folds_surv <- list('training' = vector('list', n.folds), 
                     'testing' = vector('list', n.folds))
  for(i in 1:n.folds){
    ids <- unique(folds_long$training[[i]]$id)
    folds_surv$training[[i]] <- survdata[survdata$id %in% ids,]
    folds_surv$testing[[i]] <- survdata[!(survdata$id %in% ids),]
  }
  
  # perform prediction k times
  auc_df <- brier_df <- as.data.frame(times)
  
  for(i in 1:n.folds){
    # print progress
    print(paste0('Repetition ', m, ', fold ', i))
    
    # train and test set
    df_train_long <- folds_long$training[[i]]
    df_train_surv <- folds_surv$training[[i]]
    df_test_long <- folds_long$testing[[i]]
    df_test_surv <- folds_surv$testing[[i]]

    # fit models 
    ## step 1: estimate mixed-effects model
    model_lme <- fit_lmms(y.names = c('log2PSA'),
                          fixefs = ~ therapy_received*fup_time, 
                          ranefs = ~ fup_time | id,
                          long.data = df_train_long, surv.data = df_train_surv,
                          t.from.base = fup_time,
                          n.boots = 0, n.cores = 8, verbose = F)
    ## step 2: compute predicted random effects
    sum_lme <- summarize_lmms(object = model_lme, n.cores = 8, verbose = F)
    ## step 3: estimate penalized Cox model
    model_pcr <- fit_prclmm(object = sum_lme, surv.data = df_train_surv,
                            baseline.covs = ~ therapy_received,
                            penalty = pen, standardize = F,
                            n.cores = 8, verbose = F)
    
    # validate model on test set
    ## predict survival outcome
    surv_prob <- survpred_prclmm(step1 = model_lme, step2 = sum_lme, 
                                 step3 = model_pcr, times = times,
                                 new.longdata = df_test_long,
                                 new.basecovs = df_test_surv)
    fail_prob <- 1 - surv_prob$predicted_survival[,2:ncol(surv_prob$predicted_survival)]
    colnames(fail_prob) <- times
    
    # predictive performance
    acc_measures <- Score(as.list(fail_prob), 
                          formula = Surv(time, event) ~ 1, 
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
                                 surv.time = df_test_surv$time, 
                                 surv.event = df_test_surv$event)$c.index
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
  rm(auc_df, brier_df, auc, brier, acc_measures, fail_prob, surv_prob,
     df_train_long, df_train_surv, df_test_long, df_test_surv, 
     model_lme, sum_lme, model_pcr)
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
write.csv(pp, 'Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/pp_pencal_ridge_20x5.csv', row.names = FALSE)
write.csv(all_cindex, 'Z:/Documents/Scripts/RQ3_prediction/saved_dataframes/cindex_pencal_ridge_20x5.csv', row.names = FALSE)

###########################
### prepare environment ###
###########################

# load pre-processed data
source('Z:/Documents/Scripts/data_management_RQ3.R')

# load libraries 
library(survival) # fit the survival model
library(pencal) # fit PRC model
library(riskRegression) # calculate prediction performance
library(survcomp) # calculate C-index
library(RColorBrewer) # use nice colors

#################################
### prepare for model fitting ### 
#################################

# set hyperparameters
t0 <- 4 # moment of prediction 
tn <- 10 # maximum prediction time
times <- seq(t0, tn, length.out = (tn-t0)*2+1) # number of predictions 
seed <- floor(1803158/1000) # set seed

# add baseline date to longitudinal measurements 
longdata <- psa_long_train %>%
  filter(!is.na(PSA)) %>% 
  left_join(select(psa_long[psa_long$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId') %>%
  distinct(id = id_num, fup_time = time, time = time_obs, log2PSA, therapy_received) %>%
  mutate(therapy_received = as.factor(as.character(therapy_received))) %>%
  filter(fup_time <= t0 & time > t0) # landmark data

# create dataframe with only the survival data 
survdata <- baseline_data_table1[baseline_data_table1$id_num %in% longdata$id,] %>%
  distinct(id = id_num, time = time_obs, event = NLCB_overall_num, therapy_received) %>%
  mutate(therapy_received = as.factor(as.character(therapy_received)))

# prepare validation set 
## obtain patients that are not in the training set 
longdata_validation <- psa_long[!(psa_long$id_num %in% psa_long_train$id_num),] %>%
  filter(!is.na(PSA)) %>% 
  left_join(select(psa_long[psa_long$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId') %>%
  mutate(time = as.numeric(date_lab - baseline_date)/(365.25/12),
         therapy_received = as.factor(as.character(therapy_received))) %>%
  distinct(id = id_num, fup_time = time, time = time_obs, log2PSA, therapy_received) %>%
  filter(fup_time <= t0 & time > t0) # landmark data

survdata_validation <- baseline_data_table1[baseline_data_table1$id_num %in% longdata_validation$id,] %>%
  distinct(id = id_num, time = time_obs, event = NLCB_overall_num, therapy_received) %>% 
  mutate(therapy_received = as.factor(as.character(therapy_received)))

#####################
### model fitting ###
#####################

# step 1: estimate mixed-effects model
model_lme <- fit_lmms(y.names = c('log2PSA'),
                      fixefs = ~ therapy_received*fup_time, 
                      ranefs = ~ fup_time | id,
                      long.data = as.data.frame(longdata), surv.data = survdata,
                      t.from.base = fup_time,
                      n.boots = 0, n.cores = 8, verbose = F)
# step 2: compute predicted random effects
sum_lme <- summarize_lmms(object = model_lme, n.cores = 8, verbose = F)
# step 3: estimate penalized Cox model
model_pcr <- fit_prclmm(object = sum_lme, surv.data = survdata,
                        baseline.covs = ~ therapy_received,
                        penalty = 'ridge', standardize = F,
                        n.cores = 8, verbose = F)

########################
### model validation ###
########################

# predict survival outcome
surv_prob <- survpred_prclmm(step1 = model_lme, step2 = sum_lme, 
                             step3 = model_pcr, times = times,
                             new.longdata = as.data.frame(longdata_validation),
                             new.basecovs = survdata_validation)
fail_prob <- 1 - surv_prob$predicted_survival[,2:ncol(surv_prob$predicted_survival)]
colnames(fail_prob) <- times

# predictive performance
acc_measures <- Score(as.list(fail_prob), 
                      formula = Surv(time, event) ~ 1, 
                      data = survdata_validation,
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

# C-index
c_index <- concordance.index(x = fail_prob[,2], method = 'noether',
                             surv.time = survdata_validation$time, 
                             surv.event = survdata_validation$event)$c.index

##########################
### data visualization ### 
##########################

# dynamic prediction
## obtain sample of patients from the validation set 
#pts <- c(9, 134, 262, 133,58,90,264,132,255) 
pts <- c(9, 58, 134)

## calculate estimates for new patients
new_est <- survpred_prclmm(step1 = model_lme, step2 = sum_lme, 
                                  step3 = model_pcr, times = times, keep.ranef = TRUE,
                                  new.longdata = as.data.frame(longdata_validation[longdata_validation$id %in% pts,]),
                                  new.basecovs = survdata_validation[survdata_validation$id %in% pts,])

### survival probabilities
plot.surv_prob <- new_est$predicted_survival
names(plot.surv_prob) <- c('id', seq(4, 10, 0.5))
plot.surv_prob <- reshape(plot.surv_prob, 
                          direction = 'long',
                          varying = list(names(plot.surv_prob)[2:14]),
                          v.names = 'surv_prob',
                          idvar = 'id',
                          timevar = 'times',
                          times = seq(4, 10, 0.5)) %>%
  mutate(id = factor(id, levels = pts))

### longitudinal trajectories
#### random effects 
random_effects <- new_est$predicted_ranefs 
names(random_effects) <- c('intercept', 'fup_time')
#### fixed effects
fixed_effects <- coef(model_lme$lmm.fits.orig$log2PSA)[1,]
names(fixed_effects) <- c('intercept.fix', 'PARPi', 'platinum', 'taxane', 'fup_time.fix',
                          'PARPi.fup_time', 'platinum.fup_time', 'taxane.fup_time')
fixed_effects$intercept.fix <- 4.205480274
fixed_effects$fup_time.fix <- -1.174017740
#### coefficients per patient
coefs <- cbind(random_effects, fixed_effects) %>%
  mutate(intercept = intercept + intercept.fix,
         fup_time = fup_time + fup_time.fix) 
coefs$intercept.fix <- coefs$fup_time.fix <- NULL
#### design matrix
eval.time <- seq(0, 4, 0.1)
X.prep <- survdata_validation[survdata_validation$id %in% pts, c('id', 'therapy_received')]
X <- data.frame(matrix(NA, nrow = 0, ncol = length(pts) + 1))
for(j in 1:length(eval.time)){
  X.int <- data.frame(matrix(NA, nrow = 0, ncol = ncol(coefs) + 1))
  for(i in 1:nrow(X.prep)){
    if(X.prep[i,]$therapy_received == 'ARSi') X.int <- rbind(X.int, c(pts[i], 1,eval.time[j],0,0,0,0,0,0))
    if(X.prep[i,]$therapy_received == 'PARPi') X.int <- rbind(X.int, c(pts[i], 1,eval.time[j],1,0,0,eval.time[j],0,0))
    if(X.prep[i,]$therapy_received == 'Platinum') X.int <- rbind(X.int, c(pts[i], 1,eval.time[j],0,1,0,0,eval.time[j],0))
    if(X.prep[i,]$therapy_received == 'Taxane') X.int <- rbind(X.int, c(pts[i], 1,eval.time[j],0,0,1,0,0,eval.time[j]))
  } 
  X <- rbind(X, c(diag(as.matrix(X.int[,-1]) %*% as.matrix(t(coefs))), eval.time[j]))
}
names(X) <- c(pts, 'time')
#### longitudinal trajectories
Y <- data.frame(matrix(NA, nrow = 0, ncol = 3))
for(i in 1:length(pts)) Y <- rbind(Y, cbind('log2PSA' = X[,i,ncol(X)],
                                            'time' = X[,ncol(X)], 
                                            'id' = rep(pts[i], length(eval.time))))
rm(eval.time, X.prep, X.int, coefs, random_effects, fixed_effects)

### save event time
vline.event <- survdata_validation[survdata_validation$id %in% pts,] %>%
  filter(event != 0)
vline.cens <- survdata_validation[survdata_validation$id %in% pts,] %>%
  filter(event != 1)
coef <- 1/max(longdata_validation$log2PSA[longdata_validation$id %in% pts])

## plot longitudinal and survival data 
gg.dynpred <- ggplot(plot.surv_prob, aes(x = times, y = surv_prob/(coef), group = id)) +
  geom_line() +
  geom_point(data = longdata_validation[longdata_validation$id %in% pts,],
             mapping = aes(x = fup_time, y = log2PSA, group = id),
             colour = brewer.pal(name = 'Paired', n = 12)[2]) +
  geom_line(data = Y, mapping = aes(x = time, y = log2PSA, group = id),
            colour = brewer.pal(name = 'Paired', n = 12)[2]) +
  scale_y_continuous(name = expression(log[2](PSA)),
                     limits = c(0, 5.7),
                     breaks = seq(0, 5.7, 1),
                     sec.axis = sec_axis(~ .*coef, 
                                         name = 'Survival probability', 
                                         breaks = seq(0, 1, 0.2))) +
  labs(x = 'Months since baseline', y = 'Brier score') +
  # ggtitle('') +
  xlim(c(0,10)) +
  geom_vline(xintercept = 4, linetype = 'dashed', color = 'grey60') +
  geom_vline(aes(xintercept = time), linetype = 'dotted', vline.event, 
             colour = brewer.pal(name = 'Paired', n = 12)[6], linewidth = 0.75) +
  geom_vline(aes(xintercept = time), linetype = 'dotted', vline.cens,
             colour = brewer.pal(name = 'Paired', n = 12)[4], linewidth = 0.75) +
  theme_bw() +
  facet_wrap(~ id) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        strip.text = element_text(size=14))

###################
### export data ### 
###################

# plots
plots <- list(gg.dynpred) 

filenames <- c('dyn_pred')

for (i in 1:length(plots)){  
  file_name = paste('Z:/Documents/Scripts/RQ3_prediction/saved_plots/validation/', filenames[i], '.pdf', sep='')
  pdf(file_name, height = 5, width = 11)
  print(plots[[i]])
  dev.off()
}

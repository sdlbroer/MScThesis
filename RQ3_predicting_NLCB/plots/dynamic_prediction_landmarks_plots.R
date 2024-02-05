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
tn <- 10 # maximum prediction time
seed <- floor(1803158/1000) # set seed

# add baseline date to longitudinal measurements 
longdata <- psa_long_train %>%
  filter(!is.na(PSA)) %>% 
  left_join(select(psa_long[psa_long$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId') %>%
  distinct(id = id_num, fup_time = time, time = time_obs, log2PSA, therapy_received) %>%
  mutate(therapy_received = as.factor(as.character(therapy_received))) 

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
  distinct(id = id_num, fup_time = time, time = time_obs, log2PSA, therapy_received)

survdata_validation <- baseline_data_table1[baseline_data_table1$id_num %in% longdata_validation$id,] %>%
  distinct(id = id_num, time = time_obs, event = NLCB_overall_num, therapy_received) %>% 
  mutate(therapy_received = as.factor(as.character(therapy_received)))

#####################
### model fitting ###
#####################

# t0 = 4 -----------------------------------------------------------------------
longdata.4 <- longdata %>% filter(fup_time <= 4 & time > 4)
survdata.4 <- survdata %>% filter(time > 4)
# step 1: estimate mixed-effects model
model_lme4 <- fit_lmms(y.names = c('log2PSA'),
                      fixefs = ~ therapy_received*fup_time, 
                      ranefs = ~ fup_time | id,
                      long.data = as.data.frame(longdata.4), surv.data = survdata.4,
                      t.from.base = fup_time,
                      n.boots = 0, n.cores = 8, verbose = F)
# step 2: compute predicted random effects
sum_lme4 <- summarize_lmms(object = model_lme4, n.cores = 8, verbose = F)
# step 3: estimate penalized Cox model
model_pcr4 <- fit_prclmm(object = sum_lme4, surv.data = survdata.4,
                        baseline.covs = ~ therapy_received,
                        penalty = 'ridge', standardize = F,
                        n.cores = 8, verbose = F)

# t0 = 5 -----------------------------------------------------------------------
longdata.5 <- longdata %>% filter(fup_time <= 5 & time > 5)
survdata.5 <- survdata %>% filter(time > 5)
# step 1: estimate mixed-effects model
model_lme5 <- fit_lmms(y.names = c('log2PSA'),
                      fixefs = ~ therapy_received*fup_time, 
                      ranefs = ~ fup_time | id,
                      long.data = as.data.frame(longdata.5), surv.data = survdata.5,
                      t.from.base = fup_time,
                      n.boots = 0, n.cores = 8, verbose = F)
# step 2: compute predicted random effects
sum_lme5 <- summarize_lmms(object = model_lme5, n.cores = 8, verbose = F)
# step 3: estimate penalized Cox model
model_pcr5 <- fit_prclmm(object = sum_lme5, surv.data = survdata.5,
                        baseline.covs = ~ therapy_received,
                        penalty = 'ridge', standardize = F,
                        n.cores = 8, verbose = F)

# t0 = 6 -----------------------------------------------------------------------
longdata.6 <- longdata %>% filter(fup_time <= 6 & time > 6)
survdata.6 <- survdata %>% filter(time > 6)
# step 1: estimate mixed-effects model
model_lme6 <- fit_lmms(y.names = c('log2PSA'),
                      fixefs = ~ therapy_received*fup_time, 
                      ranefs = ~ fup_time | id,
                      long.data = as.data.frame(longdata.6), surv.data = survdata.6,
                      t.from.base = fup_time,
                      n.boots = 0, n.cores = 8, verbose = F)
# step 2: compute predicted random effects
sum_lme6 <- summarize_lmms(object = model_lme6, n.cores = 8, verbose = F)
# step 3: estimate penalized Cox model
model_pcr6 <- fit_prclmm(object = sum_lme6, surv.data = survdata.6,
                        baseline.covs = ~ therapy_received,
                        penalty = 'ridge', standardize = F,
                        n.cores = 8, verbose = F)

########################
### model validation ###
########################

# t0 = 4 -----------------------------------------------------------------------
longdata_validation.4 <- longdata_validation %>% filter(fup_time <= 4 & time > 4)
survdata_validation.4 <- survdata_validation %>% filter(time > 4)
times4 <- seq(4, tn, length.out = (tn-4)*2+1)
# predict survival outcome
surv_prob <- survpred_prclmm(step1 = model_lme4, step2 = sum_lme4, 
                             step3 = model_pcr4, times = times4,
                             new.longdata = as.data.frame(longdata_validation.4),
                             new.basecovs = survdata_validation.4)
fail_prob <- 1 - surv_prob$predicted_survival[,2:ncol(surv_prob$predicted_survival)]
colnames(fail_prob) <- times4

# predictive performance
acc_measures <- Score(as.list(fail_prob), 
                      formula = Surv(time, event) ~ 1, 
                      data = survdata_validation.4,
                      times = times4, 
                      cens.model = 'km',
                      metrics = c('auc','brier'), 
                      conf.int = FALSE, 
                      exact = FALSE, 
                      split.method	= 'none', 
                      B = 0)
## extract AUC
auc4 <- acc_measures$AUC$score
auc4 <- auc4[auc4$model == auc4$times,-1]
## extract Brier
brier4 <- acc_measures$Brier$score
brier4 <- brier4[brier4$model == brier4$times,-1]

# C-index
c_index4 <- concordance.index(x = fail_prob[,2], method = 'noether',
                             surv.time = survdata_validation.4$time, 
                             surv.event = survdata_validation.4$event)$c.index

# t0 = 5 -----------------------------------------------------------------------
longdata_validation.5 <- longdata_validation %>% filter(fup_time <= 5 & time > 5)
survdata_validation.5 <- survdata_validation %>% filter(time > 5)
times5 <- seq(5, tn, length.out = (tn-5)*2+1)
# predict survival outcome
surv_prob <- survpred_prclmm(step1 = model_lme5, step2 = sum_lme5, 
                             step3 = model_pcr5, times = times5,
                             new.longdata = as.data.frame(longdata_validation.5),
                             new.basecovs = survdata_validation.5)
fail_prob <- 1 - surv_prob$predicted_survival[,2:ncol(surv_prob$predicted_survival)]
colnames(fail_prob) <- times5

# predictive performance
acc_measures <- Score(as.list(fail_prob), 
                      formula = Surv(time, event) ~ 1, 
                      data = survdata_validation.5,
                      times = times5, 
                      cens.model = 'km',
                      metrics = c('auc','brier'), 
                      conf.int = FALSE, 
                      exact = FALSE, 
                      split.method	= 'none', 
                      B = 0)
## extract AUC
auc5 <- acc_measures$AUC$score
auc5 <- auc5[auc5$model == auc5$times,-1]
## extract Brier
brier5 <- acc_measures$Brier$score
brier5 <- brier5[brier5$model == brier5$times,-1]

# C-index
c_index5 <- concordance.index(x = fail_prob[,2], method = 'noether',
                             surv.time = survdata_validation.5$time, 
                             surv.event = survdata_validation.5$event)$c.index

# t0 = 6 -----------------------------------------------------------------------
longdata_validation.6 <- longdata_validation %>% filter(fup_time <= 6 & time > 6)
survdata_validation.6 <- survdata_validation %>% filter(time > 6)
times6 <- seq(6, tn, length.out = (tn-6)*2+1)
# predict survival outcome
surv_prob <- survpred_prclmm(step1 = model_lme6, step2 = sum_lme6, 
                             step3 = model_pcr6, times = times6,
                             new.longdata = as.data.frame(longdata_validation.6),
                             new.basecovs = survdata_validation.6)
fail_prob <- 1 - surv_prob$predicted_survival[,2:ncol(surv_prob$predicted_survival)]
colnames(fail_prob) <- times6

# predictive performance
acc_measures <- Score(as.list(fail_prob), 
                      formula = Surv(time, event) ~ 1, 
                      data = survdata_validation.6,
                      times = times6, 
                      cens.model = 'km',
                      metrics = c('auc','brier'), 
                      conf.int = FALSE, 
                      exact = FALSE, 
                      split.method	= 'none', 
                      B = 0)
## extract AUC
auc6 <- acc_measures$AUC$score
auc6 <- auc6[auc6$model == auc6$times,-1]
## extract Brier
brier6 <- acc_measures$Brier$score
brier6 <- brier6[brier6$model == brier6$times,-1]

# C-index
c_index6 <- concordance.index(x = fail_prob[,2], method = 'noether',
                             surv.time = survdata_validation.6$time, 
                             surv.event = survdata_validation.6$event)$c.index

##########################
### combine dataframes ### 
##########################

pt <- 9 # 9, 58

# t0 = 4 -----------------------------------------------------------------------
## calculate estimates for new patients
new_est4 <- survpred_prclmm(step1 = model_lme4, step2 = sum_lme4, 
                            step3 = model_pcr4, times = times4, keep.ranef = TRUE,
                            new.longdata = as.data.frame(longdata_validation[longdata_validation$id == pt,]),
                            new.basecovs = survdata_validation[survdata_validation$id == pt,])

### survival probabilities
plot.surv_prob4 <- new_est4$predicted_survival
names(plot.surv_prob4) <- c('id', seq(4, 10, 0.5))
plot.surv_prob4 <- reshape(plot.surv_prob4, 
                           direction = 'long',
                           varying = list(names(plot.surv_prob4)[2:14]),
                           v.names = 'surv_prob',
                           idvar = 'id',
                           timevar = 'times',
                           times = seq(4, 10, 0.5))

### longitudinal trajectories
#### random effects 
random_effects4 <- new_est4$predicted_ranefs 
names(random_effects4) <- c('intercept', 'fup_time')
#### fixed effects
fixed_effects4 <- coef(model_lme4$lmm.fits.orig$log2PSA)[1,]
names(fixed_effects4) <- c('intercept.fix', 'PARPi', 'platinum', 'taxane', 'fup_time.fix',
                           'PARPi.fup_time', 'platinum.fup_time', 'taxane.fup_time')
fixed_effects4$intercept.fix <- 4.205480274
fixed_effects4$fup_time.fix <- -1.174017740
#### coefficients per patient
coefs4 <- cbind(random_effects4, fixed_effects4) %>%
  mutate(intercept = intercept + intercept.fix,
         fup_time = fup_time + fup_time.fix) 
coefs4$intercept.fix <- coefs4$fup_time.fix <- NULL
#### design matrix
eval.time4 <- seq(0, 4, 0.1)
X4.prep <- survdata_validation[survdata_validation$id == pt, c('id', 'therapy_received')]
X4 <- data.frame(matrix(NA, nrow = 0, ncol = length(pt) + 1))
for(j in 1:length(eval.time4)){
  X4.int <- data.frame(matrix(NA, nrow = 0, ncol = ncol(coefs4) + 1))
  for(i in 1:nrow(X4.prep)){
    if(X4.prep[i,]$therapy_received == 'ARSi') X4.int <- rbind(X4.int, c(pt[i], 1,eval.time4[j],0,0,0,0,0,0))
    if(X4.prep[i,]$therapy_received == 'PARPi') X4.int <- rbind(X4.int, c(pt[i], 1,eval.time4[j],1,0,0,eval.time4[j],0,0))
    if(X4.prep[i,]$therapy_received == 'Platinum') X4.int <- rbind(X4.int, c(pt[i], 1,eval.time4[j],0,1,0,0,eval.time4[j],0))
    if(X4.prep[i,]$therapy_received == 'Taxane') X4.int <- rbind(X4.int, c(pt[i], 1,eval.time4[j],0,0,1,0,0,eval.time4[j]))
  } 
  X4 <- rbind(X4, c(diag(as.matrix(X4.int[,-1]) %*% as.matrix(t(coefs4))), eval.time4[j]))
}
names(X4) <- c(pt, 'time')
#### longitudinal trajectories
Y4 <- data.frame(matrix(NA, nrow = 0, ncol = 3))
for(i in 1:length(pt)) Y4 <- rbind(Y4, cbind('log2PSA' = X4[,i,ncol(X4)],
                                             'time' = X4[,ncol(X4)], 
                                             'lmark' = rep(4, length(eval.time4))))
rm(eval.time4, X4.prep, X4.int, coefs4, random_effects4, fixed_effects4)
# t0 = 5 -----------------------------------------------------------------------
## calculate estimates for new patients
new_est5 <- survpred_prclmm(step1 = model_lme5, step2 = sum_lme5, 
                            step3 = model_pcr5, times = times5, keep.ranef = TRUE,
                            new.longdata = as.data.frame(longdata_validation[longdata_validation$id == pt,]),
                            new.basecovs = survdata_validation[survdata_validation$id == pt,])

### survival probabilities
plot.surv_prob5 <- new_est5$predicted_survival
names(plot.surv_prob5) <- c('id', seq(5, 10, 0.5))
plot.surv_prob5 <- reshape(plot.surv_prob5, 
                           direction = 'long',
                           varying = list(names(plot.surv_prob5)[2:12]),
                           v.names = 'surv_prob',
                           idvar = 'id',
                           timevar = 'times',
                           times = seq(5, 10, 0.5))

### longitudinal trajectories
#### random effects 
random_effects5 <- new_est5$predicted_ranefs 
names(random_effects5) <- c('intercept', 'fup_time')
#### fixed effects
fixed_effects5 <- coef(model_lme5$lmm.fits.orig$log2PSA)[1,]
names(fixed_effects5) <- c('intercept.fix', 'PARPi', 'platinum', 'taxane', 'fup_time.fix',
                           'PARPi.fup_time', 'platinum.fup_time', 'taxane.fup_time')
fixed_effects5$intercept.fix <- 3.8640221
fixed_effects5$fup_time.fix <- -0.9207029
#### coefficients per patient
coefs5 <- cbind(random_effects5, fixed_effects5) %>%
  mutate(intercept = intercept + intercept.fix,
         fup_time = fup_time + fup_time.fix) 
coefs5$intercept.fix <- coefs5$fup_time.fix <- NULL
#### design matrix
eval.time5 <- seq(0, 5, 0.1)
X5.prep <- survdata_validation[survdata_validation$id == pt, c('id', 'therapy_received')]
X5 <- data.frame(matrix(NA, nrow = 0, ncol = length(pt) + 1))
for(j in 1:length(eval.time5)){
  X5.int <- data.frame(matrix(NA, nrow = 0, ncol = ncol(coefs5) + 1))
  for(i in 1:nrow(X5.prep)){
    if(X5.prep[i,]$therapy_received == 'ARSi') X5.int <- rbind(X5.int, c(pt[i], 1,eval.time5[j],0,0,0,0,0,0))
    if(X5.prep[i,]$therapy_received == 'PARPi') X5.int <- rbind(X5.int, c(pt[i], 1,eval.time5[j],1,0,0,eval.time5[j],0,0))
    if(X5.prep[i,]$therapy_received == 'Platinum') X5.int <- rbind(X5.int, c(pt[i], 1,eval.time5[j],0,1,0,0,eval.time5[j],0))
    if(X5.prep[i,]$therapy_received == 'Taxane') X5.int <- rbind(X5.int, c(pt[i], 1,eval.time5[j],0,0,1,0,0,eval.time5[j]))
  } 
  X5 <- rbind(X5, c(diag(as.matrix(X5.int[,-1]) %*% as.matrix(t(coefs5))), eval.time5[j]))
}
names(X5) <- c(pt, 'time')
#### longitudinal trajectories
Y5 <- data.frame(matrix(NA, nrow = 0, ncol = 3))
for(i in 1:length(pt)) Y5 <- rbind(Y5, cbind('log2PSA' = X5[,i,ncol(X5)],
                                             'time' = X5[,ncol(X5)], 
                                             'lmark' = rep(5, length(eval.time5))))
rm(eval.time5, X5.prep, X5.int, coefs5, random_effects5, fixed_effects5)
# t0 = 6 -----------------------------------------------------------------------
## calculate estimates for new patients
new_est6 <- survpred_prclmm(step1 = model_lme6, step2 = sum_lme6, 
                            step3 = model_pcr6, times = times6, keep.ranef = TRUE,
                            new.longdata = as.data.frame(longdata_validation[longdata_validation$id == pt,]),
                            new.basecovs = survdata_validation[survdata_validation$id == pt,])

### survival probabilities
plot.surv_prob6 <- new_est6$predicted_survival
names(plot.surv_prob6) <- c('id', seq(6, 10, 0.5))
plot.surv_prob6 <- reshape(plot.surv_prob6, 
                           direction = 'long',
                           varying = list(names(plot.surv_prob6)[2:10]),
                           v.names = 'surv_prob',
                           idvar = 'id',
                           timevar = 'times',
                           times = seq(6, 10, 0.5))

### longitudinal trajectories
#### random effects 
random_effects6 <- new_est6$predicted_ranefs 
names(random_effects6) <- c('intercept', 'fup_time')
#### fixed effects
fixed_effects6 <- coef(model_lme6$lmm.fits.orig$log2PSA)[1,]
names(fixed_effects6) <- c('intercept.fix', 'PARPi', 'platinum', 'taxane', 'fup_time.fix',
                           'PARPi.fup_time', 'platinum.fup_time', 'taxane.fup_time')
fixed_effects6$intercept.fix <- 3.8640221
fixed_effects6$fup_time.fix <- -0.9207029
#### coefficients per patient
coefs6 <- cbind(random_effects6, fixed_effects6) %>%
  mutate(intercept = intercept + intercept.fix,
         fup_time = fup_time + fup_time.fix) 
coefs6$intercept.fix <- coefs6$fup_time.fix <- NULL
#### design matrix
eval.time6 <- seq(0, 6, 0.1)
X6.prep <- survdata_validation[survdata_validation$id == pt, c('id', 'therapy_received')]
X6 <- data.frame(matrix(NA, nrow = 0, ncol = length(pt) + 1))
for(j in 1:length(eval.time6)){
  X6.int <- data.frame(matrix(NA, nrow = 0, ncol = ncol(coefs6) + 1))
  for(i in 1:nrow(X6.prep)){
    if(X6.prep[i,]$therapy_received == 'ARSi') X6.int <- rbind(X6.int, c(pt[i], 1,eval.time6[j],0,0,0,0,0,0))
    if(X6.prep[i,]$therapy_received == 'PARPi') X6.int <- rbind(X6.int, c(pt[i], 1,eval.time6[j],1,0,0,eval.time6[j],0,0))
    if(X6.prep[i,]$therapy_received == 'Platinum') X6.int <- rbind(X6.int, c(pt[i], 1,eval.time6[j],0,1,0,0,eval.time6[j],0))
    if(X6.prep[i,]$therapy_received == 'Taxane') X6.int <- rbind(X6.int, c(pt[i], 1,eval.time6[j],0,0,1,0,0,eval.time6[j]))
  } 
  X6 <- rbind(X6, c(diag(as.matrix(X6.int[,-1]) %*% as.matrix(t(coefs6))), eval.time6[j]))
}
names(X6) <- c(pt, 'time')
#### longitudinal trajectories
Y6 <- data.frame(matrix(NA, nrow = 0, ncol = 3))
for(i in 1:length(pt)) Y6 <- rbind(Y6, cbind('log2PSA' = X6[,i,ncol(X6)],
                                             'time' = X6[,ncol(X6)], 
                                             'lmark' = rep(6, length(eval.time6))))
rm(eval.time6, X6.prep, X6.int, coefs6, random_effects6, fixed_effects6)

# combine data -----------------------------------------------------------------
plot.surv_prob <- rbind(cbind(plot.surv_prob4, 'lmark' = 4), 
                        cbind(plot.surv_prob5, 'lmark' = 5), 
                        cbind(plot.surv_prob6, 'lmark' = 6))
Y <- rbind(Y4, Y5, Y6)

##########################
### data visualization ### 
##########################

lmark.times <- data.frame('id' = pt, lmark = c(4,5,6), time = c(4,5,6))
vline.event <- survdata_validation[survdata_validation$id == pt,] %>%
  filter(event != 0)
vline.cens <- survdata_validation[survdata_validation$id == pt,] %>%
  filter(event != 1)
coef <- 1/max(longdata_validation$log2PSA[longdata_validation$id == pt])

longdata_validation.4$lmark <- 4
longdata_validation.5$lmark <- 5
longdata_validation.6$lmark <- 6

gg.dynpred <- ggplot(plot.surv_prob, aes(x = times, y = surv_prob/(coef), group = lmark)) +
  geom_line() +
  geom_point(data = longdata_validation.4[longdata_validation.4$id == pt,], 
            mapping = aes(x = fup_time, y = log2PSA, group = lmark),
            colour = brewer.pal(name = 'Paired', n = 12)[2]) +
  geom_point(data = longdata_validation.5[longdata_validation.5$id == pt,], 
            mapping = aes(x = fup_time, y = log2PSA, group = lmark),
            colour = brewer.pal(name = 'Paired', n = 12)[2]) +
  geom_point(data = longdata_validation.6[longdata_validation.6$id == pt,], 
            mapping = aes(x = fup_time, y = log2PSA, group = lmark),
            colour = brewer.pal(name = 'Paired', n = 12)[2]) +
  geom_line(data = Y, mapping = aes(x = time, y = log2PSA, group = lmark),
            colour = brewer.pal(name = 'Paired', n = 12)[2]) + 
  scale_y_continuous(name = expression(log[2](PSA)),
                     limits = c(0, 5.8),
                     breaks = seq(0, 5.8, 1),
                     sec.axis = sec_axis(~ .*coef, 
                                         name = "Survival probability", 
                                         breaks = seq(0, 1, 0.2))) +
  labs(x = 'Months since baseline', y = 'Brier score') +
  # ggtitle('') +
  xlim(c(0,10)) +
  geom_vline(xintercept = vline.event$time, linetype = 'dotted', linewidth = 0.75,
             colour = brewer.pal(name = 'Paired', n = 12)[6]) +
  geom_vline(xintercept = vline.cens$time, linetype = 'dotted', linewidth = 0.75,
             colour = brewer.pal(name = 'Paired', n = 12)[4]) +
  geom_vline(aes(xintercept = time), linetype = 'dashed', lmark.times, 
             colour = 'grey60') +
  theme_bw() +
  facet_wrap(~ lmark) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        strip.text = element_text(size=14))


###################
### export data ### 
###################

# plots
plots <- list(gg.dynpred) 

filenames <- c(paste0('dyn_pred_', pt))

for (i in 1:length(plots)){  
  file_name = paste('Z:/Documents/Scripts/RQ3_prediction/saved_plots/validation/', filenames[i], '.pdf', sep='')
  pdf(file_name, height = 5, width = 11)
  print(plots[[i]])
  dev.off()
}

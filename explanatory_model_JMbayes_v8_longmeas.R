# first run the data_management file written by me

###########################
### prepare environment ###
###########################

# load libraries 
library(nlme) # to fit the mixed model
library(survival) # to fit the survival model
library(JMbayes2) # to fit the joint model
library(car) # to use the qqPlot function for model diagnostics
library(splines) # to use splines

# set seed
set.seed(1803158)

######################################
### prepare data for model fitting ### 
######################################

# create dataframe with only the survival data 
survdata <- baseline_data_table1[baseline_data_table1$id_num %in% long_meas_train$id_num,] %>%
  distinct(id_num, therapy_received, NLCB_overall, NLCB_overall_num, time_obs)

################################################
### model 1: lmm + current value association ### 
################################################

# fit the longitudinal sub-model
model1_longit <- lme(fixed = log2PSA ~ months_in_followup + therapy_received:months_in_followup, 
                         random = ~ months_in_followup | id_num,
                         data = long_meas_train,
                         na.action = na.omit)

# fit the survival sub-model 
model1_surv <- coxph(Surv(time_obs, NLCB_overall_num) ~ therapy_received, 
                         data = survdata)

# fit joint model 
model1_joint <- jm(model1_surv, model1_longit, 
                 time_var = "months_in_followup")

#############################################################
### model 2: lmm with splines + current value association ### 
#############################################################

# fit the longitudinal sub-model
model2_longit <- lme(fixed = log2PSA ~ ns(months_in_followup, k = c(2, 6), B = c(0, 35.5))*therapy_received, 
                     random = ~ ns(months_in_followup, k = c(2, 6), B = c(0, 35.5)) | patientId,
                     data = long_meas_train,
                     control = lmeControl(opt = 'optim'))

# fit the survival sub-model 
model2_surv <- coxph(Surv(time_obs, NLCB_overall_num) ~ therapy_received, 
                     data = survdata)

# fit joint model 
model2_joint <- jm(model2_surv, model2_longit, time_var = "months_in_followup")

##########################################################################
### model 3: lmm with splines + time-varying current value association ### 
##########################################################################

form_splines3 <- ~ value(log2PSA * ns(months_in_followup, k = c(2, 6), B = c(0, 35.5)))
model3_joint <- update(model2_joint, functional_forms = form_splines3)

#######################################################################
### model 4: lmm with splines + current value and slope association ### 
#######################################################################

form_splines4 <- ~ value(log2PSA) + slope(log2PSA)
model4_joint <- update(model2_joint, functional_forms = form_splines4)

form_splines5 <- ~ value(log2PSA) + slope(log2PSA, eps = 1, direction = 'back')
model5_joint <- update(model2_joint, functional_forms = form_splines5)

################################
### comparing model outcomes ### 
################################

summary(model1_joint)$Survival
summary(model2_joint)$Survival
summary(model3_joint)$Survival
summary(model4_joint)$Survival
summary(model5_joint)$Survival

compare_jm(model1_joint, model2_joint, model3_joint, model4_joint, model5_joint)

#########################
### model diagnostics ###
#########################















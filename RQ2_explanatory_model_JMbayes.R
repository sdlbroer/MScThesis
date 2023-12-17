###########################
### prepare environment ###
###########################

# load data
source('C:/Users/lanbro/OneDrive - Karolinska Institutet/Dokument/Scripts/data_management_RQ2.R')

# load libraries 
library(nlme) # fit the mixed model
library(survival) # fit the survival model
library(JMbayes2) # fit the joint model
library(car) # use the qqPlot function for model diagnostics
library(splines) # use splines
library(HLMdiag) # find outlying and influential points
library(survminer) # diagnostics of the survival submodel
library(RColorBrewer) # plot colours

# set seed
set.seed(1803158)

######################################
### prepare data for model fitting ### 
######################################

# add baseline date to psa_long dataframe
psa_long_train <- psa_long_train %>%
  filter(!is.na(PSA)) %>% 
  left_join(select(psa_long_train[psa_long_train$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId')

# create dataframe with only the survival data 
survdata <- baseline_data_table1_train[baseline_data_table1_train$id_num %in% psa_long_train$id_num,] %>%
  distinct(id_num, therapy_received, NLCB_overall, NLCB_overall_num, time_obs)

##########################
### survival sub-model ### 
##########################

model_surv <- coxph(Surv(time_obs, NLCB_overall_num) ~ therapy_received, 
                    data = survdata)

#######################
### mixed sub-model ### 
#######################

# lmm
model1_longit <- lme(fixed = log2PSA ~ therapy_received*time, 
                     random = ~ time | id_num,
                     data = psa_long_train)

# mm with splines and 1 knot
model2_longit <- lme(fixed = log2PSA ~ ns(time, k = c(2), B = c(0, 18.5))*therapy_received, 
                     random = ~ ns(time, k = c(2), B = c(0, 18.5)) | id_num,
                     data = psa_long_train,
                     control = lmeControl(opt = 'optim'))

# mm with splines and 2 knots
model3_longit <- lme(fixed = log2PSA ~ ns(time, k = c(2, 4), B = c(0, 18.5))*therapy_received, 
                     random = ~ ns(time, k = c(2, 4), B = c(0, 18.5)) | id_num,
                     data = psa_long_train,
                     control = lmeControl(opt = 'optim'))

# mm with splines and 3 knots
model4_longit <- lme(fixed = log2PSA ~ ns(time, k = c(2, 4, 6), B = c(0, 18.5))*therapy_received, 
                     random = ~ ns(time, k = c(2, 4, 6), B = c(0, 18.5)) | id_num,
                     data = psa_long_train,
                     control = lmeControl(opt = 'optim'))

# mm with quadratic time  
model5_longit <- lme(fixed = log2PSA ~ poly(time, 2)*therapy_received, 
                     random = ~ poly(time, 2) | id_num,
                     data = psa_long_train,
                     control = lmeControl(opt = 'optim'))

####################
### joint models ### 
####################

# functional forms 
## time-varying association 1
form_splines1 <- ~ value(log2PSA) * ns(time, k = c(4), B = c(0, 18.5))
## time-varying association 2
form_splines2 <- ~ value(log2PSA) * ns(time, k = c(2, 4), B = c(0, 18.5))
## time-varying association 3
form_splines3 <- ~ value(log2PSA) * ns(time, k = c(2, 4, 6), B = c(0, 18.5))
## current-value + slope
form_slope <- ~ value(log2PSA) + slope(log2PSA, eps = 1, direction = 'back')

# lmm
## current-value association
model1_lmm_joint <- jm(model_surv, model1_longit, time_var = "time")
## current-value + slope association 
model2_lmm_joint <- update(model1_lmm_joint, functional_forms = form_slope)
## time-varying association 1
model3_lmm_joint <- update(model1_lmm_joint, functional_forms = form_splines1)
## time-varying association 2
model4_lmm_joint <- update(model1_lmm_joint, functional_forms = form_splines2)
## time-varying association 3
model5_lmm_joint <- update(model1_lmm_joint, functional_forms = form_splines3)

# mm with splines and 1 internal knot 
## current-value association
model1_mm1int_joint <- jm(model_surv, model2_longit, time_var = "time")
## current-value + slope association 
model2_mm1int_joint <- update(model1_mm1int_joint, functional_forms = form_slope)
## time-varying association 1
model3_mm1int_joint <- update(model1_mm1int_joint, functional_forms = form_splines1)
## time-varying association 2
model4_mm1int_joint <- update(model1_mm1int_joint, functional_forms = form_splines2)
## time-varying association 3
model5_mm1int_joint <- update(model1_mm1int_joint, functional_forms = form_splines3)

# mm with splines and 2 internal knots
## current-value association
model1_mm2int_joint <- jm(model_surv, model3_longit, time_var = "time")
## current-value + slope association 
model2_mm2int_joint <- update(model1_mm2int_joint, functional_forms = form_slope)
## time-varying association 1
model3_mm2int_joint <- update(model1_mm2int_joint, functional_forms = form_splines1)
## time-varying association 2
model4_mm2int_joint <- update(model1_mm2int_joint, functional_forms = form_splines2)
## time-varying association 3
model5_mm2int_joint <- update(model1_mm2int_joint, functional_forms = form_splines3)

# mm with splines and 3 internal knots
## current-value association
model1_mm3int_joint <- jm(model_surv, model4_longit, time_var = "time")
## current-value + slope association 
model2_mm3int_joint <- update(model1_mm3int_joint, functional_forms = form_slope)
## time-varying association 1
model3_mm3int_joint <- update(model1_mm3int_joint, functional_forms = form_splines1)
## time-varying association 2
model4_mm3int_joint <- update(model1_mm3int_joint, functional_forms = form_splines2)
## time-varying association 3
model5_mm3int_joint <- update(model1_mm3int_joint, functional_forms = form_splines3)

# mm with time modelled quadratic
## current-value association
model1_quad_joint <- jm(model_surv, model5_longit, time_var = "time")
## current-value + slope association 
model2_quad_joint <- update(model1_quad_joint, functional_forms = form_slope)
## time-varying association 1
model3_quad_joint <- update(model1_quad_joint, functional_forms = form_splines1)
## time-varying association 2
model4_quad_joint <- update(model1_quad_joint, functional_forms = form_splines2)
## time-varying association 3
model5_quad_joint <- update(model1_quad_joint, functional_forms = form_splines3)

################################
### comparing model outcomes ### 
################################

compare_jm(model1_lmm_joint, model2_lmm_joint, model3_lmm_joint,
           model4_lmm_joint, model5_lmm_joint,
           model1_mm1int_joint, model2_mm1int_joint, model3_mm1int_joint,
           model4_mm1int_joint, model5_mm1int_joint,
           model1_mm2int_joint, model2_mm2int_joint, model3_mm2int_joint,
           model4_mm2int_joint, model5_mm2int_joint,
           model1_mm3int_joint, model2_mm3int_joint, model3_mm3int_joint,
           model4_mm3int_joint, model5_mm3int_joint,
           model1_quad_joint, model2_quad_joint, model3_quad_joint,
           model4_quad_joint, model5_quad_joint)

compare_jm(model1_lmm_joint, model2_lmm_joint, model3_lmm_joint,
           model4_lmm_joint, model5_lmm_joint)
compare_jm(model1_mm1int_joint, model2_mm1int_joint, model3_mm1int_joint,
           model4_mm1int_joint, model5_mm1int_joint)
compare_jm(model1_mm2int_joint, model2_mm2int_joint, model3_mm2int_joint,
           model4_mm2int_joint, model5_mm2int_joint)
compare_jm(model1_mm3int_joint, model2_mm3int_joint, model3_mm3int_joint,
           model4_mm3int_joint, model5_mm3int_joint)
compare_jm(model1_quad_joint, model2_quad_joint, model3_quad_joint,
           model4_quad_joint, model5_quad_joint)

compare_jm(model5_lmm_joint, # time-dependent association with splines (3 knots)
           model5_mm1int_joint, # time-dependent association with splines (2 knots)
           model5_mm2int_joint, # time-dependent association with splines (2 knots)
           model4_mm3int_joint, # time-dependent association with splines (2 knots)
           model5_quad_joint) # time-dependent association with splines (2 knots)

round(summary(model5_lmm_joint)$Survival, 2) # linear mixed, alpha(t), 3 knots
round(summary(model5_mm1int_joint)$Survival, 2) # splines mixed (1 knot), alpha(t), 2 knots
round(summary(model5_quad_joint)$Survival, 2) # quadratic mixed, alpha(t), 2 knots
round(summary(model4_mm3int_joint)$Survival, 2) # splines mixed (3 knots), alpha(t), 2 knots
round(summary(model5_mm2int_joint)$Survival, 2) # splines mixed (2 knots), alpha(t), 2 knots

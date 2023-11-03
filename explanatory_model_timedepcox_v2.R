###########################
### prepare environment ###
###########################

# load libraries 
library(survival) # to fit the Cox model

#######################################
### fit time-dependent Cox PH model ###
#######################################

# dataframe with only the necessary survival data 
survdata <- left_join(
  distinct(long_meas_train, id, patientId, PSA, NLCB, months_in_followup),
  distinct(baseline_data_table1_train, patientId, therapy_received, NLCB_overall, 
           NLCB_overall_num, time_obs),
  by = 'patientId'
) %>%
  mutate(therapy_received = as.character(therapy_received)) %>%
  arrange(patientId, months_in_followup)

# create the tstart variable 
for(i in 1:nrow(survdata)){
  if(survdata[i, 'months_in_followup'] == 0) survdata[i, 'tstart'] = 0
  else if(survdata[i-1, 'patientId'] == survdata[i, 'patientId']){
    survdata[i, 'tstart'] = survdata[i-1, 'months_in_followup']
  } 
}

# time-dependent Cox
coxmodel <- coxph(Surv(tstart, months_in_followup, NLCB) ~ therapy_received + 
                    log(PSA + 0.01, 2), 
      data = survdata,
      id = patientId)

summary(coxmodel)

###########################
### prepare environment ###
###########################

# load data
source('C:/Users/lanbro/Documents/Scripts/data_management_RQ2.R')

# load libraries 
library(survival) # to fit the Cox model

#######################################
### fit time-dependent Cox PH model ###
#######################################

# dataframe with only the necessary survival data 
survdata <- left_join(
  distinct(long_meas_train, id, patientId, PSA, NLCB, time),
  distinct(baseline_data_table1_train, patientId, therapy_received, NLCB_overall, 
           NLCB_overall_num, time_obs),
  by = 'patientId'
) %>%
  mutate(therapy_received = as.character(therapy_received)) %>%
  arrange(patientId, time)

# create the tstart variable 
for(i in 1:nrow(survdata)){
  if(survdata[i, 'time'] == 0) survdata[i, 'tstart'] = 0
  else if(survdata[i-1, 'patientId'] == survdata[i, 'patientId']){
    survdata[i, 'tstart'] = survdata[i-1, 'time']
  } 
}

# time-dependent Cox
coxmodel <- coxph(Surv(tstart, time, NLCB) ~ therapy_received + 
                    log(PSA + 0.01, 2), 
      data = survdata,
      id = patientId)

summary(coxmodel)

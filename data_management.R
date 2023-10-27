###########################
### prepare environment ###
###########################

rm(list = ls())

# load libraries 
library(tidyverse)

# load all data
derived_data <- list.files(path = "C:/Users/lanbro/Documents/Data/2023-10-17",
                           pattern = ".RData$", full.names = T)
invisible(lapply(derived_data, load, .GlobalEnv))
rm(derived_data)

###################################################################
### new dataframe with months of follow-up + experiencing event ###
###################################################################

# longitudinal PSA measurements, with randomization as baseline 
# the code for the psa_long dataframe was written by Alessio Crippa, 20231023

subjects_cr$id_num <- consecutive_id(subjects_cr$id)

psa_long <- bind_rows(
  # treatment start as new baseline
  select(visit0_cr, id, patientId, date_visit = E3_F17_dateOfVisitMonth0,
         date_lab = E3_F7_dateOfRoutineLab, PSA = E3_F7_PSA) %>% 
    mutate(visit = "Treatment start", psa_prog = 0, radio_prog = 0, clinic_prog = 0, NLCB = 0),
  # at first visit, check treatment initiation
  select(visit1_cr, id, patientId, date_visit = E4_F21_dateOfVisitMonth1,
         date_lab = E4_F7_dateOfRoutineLab, PSA = E4_F7_PSA) %>% 
    mutate(visit = "Treatment check", psa_prog = 0, radio_prog = 0, clinic_prog = 0, NLCB = 0),
  # regular or unscheduled visits (follow-up with treatment evaluation)
  # all planned visits are included, exclude those not happened yet (missing date_lab)
  filter(visits_cr, status != "Not Active") %>% 
    select(id, patientId, visit, date_visit = F23_dateOfVisit,
           date_lab = F7_dateOfRoutineLab, PSA = F7_PSA,
           psa_prog, radio_prog, clinic_prog, NLCB = F18_noLongerClinBenefit)
) %>% 
  # add info from the main object (subjects_cr)
  left_join(select(subjects_cr, id, id_num, status_probio, therapy_class,
                   date_rand = E1_F16_dateRand, therapy_received, locked, 
                   prog, time_obs, last_dte), by = "id") %>% 
  # add info from from baseline evaluation
  left_join(select(baseline_cr, id, baseline_PSA = E1_F7_PSA), by = "id") %>% 
  # considering only patients in subjects_cr
  filter(id %in% subjects_cr$id) %>% 
  mutate(
    # replacing date of treatment start with date of randomization, to be consistent with survival model
    date_lab = replace(date_lab, visit == "Treatment start", date_rand[visit == "Treatment start"]),
    # replacing missing dates of lab with visits date (not actually needed since PSA was not measured then)
    date_lab = replace(date_lab, is.na(date_lab), date_visit[is.na(date_lab)]),
    # some dates of lab are done after discontinuation, replacing those with date of visit instead 
    date_lab = replace(date_lab, !is.na(date_lab) & date_lab > last_dte & date_visit <= last_dte,
                       date_visit[!is.na(date_lab) & date_lab > last_dte & date_visit <= last_dte]),
    # replacing missing PSA at treatment start with PSA at study baseline (7 subjects)
    PSA = replace(PSA, visit == "Treatment start" & is.na(PSA), baseline_PSA[visit == "Treatment start" & is.na(PSA)]),
    # log transformations
    log2PSA = log2(PSA + 0.001), 
    logPSA = log(PSA + 0.001),
    log10PSA = log10(PSA + 0.001),
    therapy_class = droplevels(therapy_class),
    therapy_received = droplevels(therapy_received)
  ) %>% 
  # keep only those where a PSA measure is available
  filter(!is.na(date_lab)) %>%  # none of these visits has a PSA measure
  #filter(!is.na(PSA)) %>%       # excluding about 37 PSA measurements
  # NB: remove after checks
  filter(date_lab >= date_rand & date_lab <= last_dte) %>% 
  arrange(patientId, date_lab) %>% 
  group_by(patientId) %>% 
  mutate(
    date_last = max(date_lab, na.rm = T),
    num_psa_value = sum(!is.na(PSA)),
    seq_psa_value = seq(num_psa_value)
  ) %>% 
  ungroup()

# adding info on number of PSA measurements to main data
subjects_cr <- left_join(
  subjects_cr,
  select(filter(psa_long, seq_psa_value == 1), id, num_psa_value),
  by = "id"
) %>% 
  mutate(
    therapy_class = droplevels(therapy_class),
    therapy_received = droplevels(therapy_received)
  )

# baseline measurements 
baseline_data_table1 <- left_join(
  select(visit0_cr, id, patientId, date_lab = E3_F7_dateOfRoutineLab, PSA = E3_F7_PSA,
         Hb = E3_F7_Hb, RBC = E3_F7_RBC, leukocytes = E3_F7_leukocytes, 
         neutrophils = E3_F7_neutrophils, lymphocytes = E3_F7_lymphocytes, 
         monocytes = E3_F7_monocytes, platelets = E3_F7_platelets, 
         testosterone = E3_F7_testosterone, LDH = E3_F7_LDH, ALP = E3_F7_ALP, 
         albumin = E3_F7_albumin, crea = E3_F7_crea, AST = E3_F7_AST, ALT = E3_F7_ALT, 
         bilirubin = E3_F7_bilirubin, glucose = E3_F7_glucose, potassium = E3_F7_potassium,
         sbp = E3_F20_sbp, dbp = E3_F20_dbp, hr = E3_F20_hr, weight = E3_F20_weight),
  select(baseline_cr, patientId, country, dateOfBirth, height = E1_F4_height, 
         PSADT = E1_F3_PSADT, dateOfmetastaticDisease = E2_F5_dateOfmetastaticDisease, 
         dateOfPCdiagnosis = E2_F5_dateOfPCdiagnosis, PSAatDiagnosis = E2_F5_PSAatDiagnosis, 
         TNMstagingT = E2_F5_TNMstagingT, TNMstagingN = E2_F5_TNMstagingN, 
         TNMstagingM = E2_F5_TNMstagingM, Gleason = E2_F5_ISUP, 
         LocationMetastases = E2_F26_location_metastases)) %>%
  left_join(select(subjects_cr, id_num, patientId, therapy_class, therapy_received), by = 'patientId') %>%
  left_join(distinct(psa_long[psa_long$NLCB == 1,], patientId, prog), 
            by = 'patientId') %>%
  mutate(dateOfBirth = as.numeric(date_lab - dateOfBirth)/365.25,
         dateOfmetastaticDisease = as.numeric(date_lab - dateOfmetastaticDisease)/365.25,
         dateOfPCdiagnosis = as.numeric(date_lab - dateOfPCdiagnosis)/365.25,
         therapy_class = as.character(therapy_class),
         therapy_received = as.character(therapy_received),
         NLCB_overall_num = ifelse(is.na(prog), 0, prog),
         NLCB_overall = ifelse(NLCB_overall_num == 0, 'No event', 'Event')) %>% 
  rename(age = dateOfBirth,
         TimeSinceMetastaticDisease = dateOfmetastaticDisease,
         TimeSincePCdiagnosis = dateOfPCdiagnosis)

# longitudinal measurements

long_meas <- bind_rows(
  # at baseline (defined as moment of randomization)
  select(visit0_cr, id, patientId, date_lab = E3_F7_dateOfRoutineLab, PSA = E3_F7_PSA,
         date_visit = E3_F17_dateOfVisitMonth0, Hb = E3_F7_Hb, RBC = E3_F7_RBC, 
         leukocytes = E3_F7_leukocytes, neutrophils = E3_F7_neutrophils, 
         lymphocytes = E3_F7_lymphocytes, monocytes = E3_F7_monocytes, 
         platelets = E3_F7_platelets, testosterone = E3_F7_testosterone, LDH = E3_F7_LDH, 
         ALP = E3_F7_ALP, albumin = E3_F7_albumin, crea = E3_F7_crea, AST = E3_F7_AST, 
         ALT = E3_F7_ALT, bilirubin = E3_F7_bilirubin, glucose = E3_F7_glucose, 
         potassium = E3_F7_potassium, sbp = E3_F20_sbp, dbp = E3_F20_dbp, 
         hr = E3_F20_hr, weight = E3_F20_weight) %>% 
    mutate(visit = 'Treatment start', psa_prog = 0, radio_prog = 0, clinic_prog = 0, NLCB = 0),
  # at first visit, check treatment initiation
  select(visit1_cr, id, patientId, date_lab = E4_F7_dateOfRoutineLab, PSA = E4_F7_PSA,
         Hb = E4_F7_Hb, RBC = E4_F7_RBC, leukocytes = E4_F7_leukocytes, 
         date_visit = E4_F21_dateOfVisitMonth1, neutrophils = E4_F7_neutrophils, 
         lymphocytes = E4_F7_lymphocytes, monocytes = E4_F7_monocytes, 
         platelets = E4_F7_platelets, testosterone = E4_F7_testosterone, 
         LDH = E4_F7_LDH, ALP = E4_F7_ALP, albumin = E4_F7_albumin, crea = E4_F7_crea, 
         AST = E4_F7_AST, ALT = E4_F7_ALT, bilirubin = E4_F7_bilirubin, 
         glucose = E4_F7_glucose, potassium = E4_F7_potassium) %>% 
    mutate(visit = 'Treatment check', psa_prog = 0, radio_prog = 0, clinic_prog = 0, NLCB = 0),
  # regular or unscheduled visits (follow-up with treatment evaluation)
  select(visits_cr, id, patientId, visit, date_lab = F7_dateOfRoutineLab, PSA = F7_PSA,
         date_visit = F23_dateOfVisit, psa_prog, radio_prog, clinic_prog, 
         NLCB = F18_noLongerClinBenefit, Hb = F7_Hb, RBC = F7_RBC, 
         leukocytes = F7_leukocytes, neutrophils = F7_neutrophils, 
         lymphocytes = F7_lymphocytes, monocytes = F7_monocytes, platelets = F7_platelets, 
         testosterone = F7_testosterone, LDH = F7_LDH, ALP = F7_ALP, 
         albumin = F7_albumin, crea = F7_crea, AST = F7_AST, ALT = F7_ALT, 
         bilirubin = F7_bilirubin, glucose = F7_glucose, potassium = F7_potassium,
         sbp = F20_sbp, dbp = F20_dbp, hr = F20_hr, weight = F20_weight)
) %>% 
  # add info from from baseline evaluation
  left_join(select(baseline_cr, id, baseline_PSA = E1_F7_PSA), by = "id") %>% 
  # add info from the main object (subjects_cr)
  left_join(select(subjects_cr, id, id_num, therapy_class, date_rand = E1_F16_dateRand,
                   therapy_received, locked, prog, time_obs, country, last_dte), by = 'id') %>% 
  left_join(select(psa_long[psa_long$visit == 'Treatment start',], patientId, 
                   baseline_date = date_lab), by = 'patientId') %>%
  left_join(distinct(psa_long[psa_long$NLCB == 1,], patientId, NLCB_overall = NLCB), 
            by = 'patientId') %>%
  mutate(
    # replacing date of treatment start with date of randomization, to be consistent with survival model
    date_lab = replace(date_lab, visit == 'Treatment start', 
                       date_rand[visit == 'Treatment start']),
    # replacing missing dates of lab with visits date (not actually needed since PSA was not measured then)
    date_lab = replace(date_lab, is.na(date_lab), date_visit[is.na(date_lab)]),
    # some dates of lab are done after discontinuation, replacing those with date of visit instead 
    date_lab = replace(date_lab, !is.na(date_lab) & date_lab > last_dte & date_visit <= last_dte,
                       date_visit[!is.na(date_lab) & date_lab > last_dte & date_visit <= last_dte]),
    # replacing missing PSA at treatment start with PSA at study baseline (7 subjects)
    PSA = replace(PSA, visit == "Treatment start" & is.na(PSA), baseline_PSA[visit == "Treatment start" & is.na(PSA)]),
    NLCB_overall_num = ifelse(is.na(prog), 0, prog),
    NLCB_overall = ifelse(NLCB_overall_num == 0, 'No event', 'Event')
  ) %>% 
  # keep only those where a PSA measure is available
  filter(!is.na(date_lab)) %>%  # none of these visits has a PSA measure
  #filter(!is.na(PSA)) %>%       # excluding about 37 PSA measurements
  # NB: remove after checks
  filter(date_lab >= date_rand & date_lab <= last_dte) %>% 
  group_by(patientId) %>% 
  mutate(date_last = max(date_lab, na.rm = T),
         num_psa_value = sum(!is.na(PSA)),
         time_obs = as.numeric(max(date_lab, na.rm = T) - min(date_lab))/(365.25/12),
         seq_psa_value = seq(num_psa_value)
  ) %>% 
  ungroup() %>% 
  # remove after checks
  mutate(therapy_class = as.character(therapy_class),
         therapy_received = factor(therapy_received, levels = c('ARSi', 'Taxane', 'PARPi', 'Platinum')),
         months_in_followup = as.numeric(date_lab - baseline_date)/(365.25/12),
         months_in_followup_rounded = round(months_in_followup, 0),
         reverse_time = months_in_followup - time_obs,
         reverse_time_rounded = months_in_followup_rounded - round(time_obs, 0),
         log2PSA = log(PSA + 0.01, 2)
  ) %>%
  arrange(patientId, date_lab)

# add observed time to baseline for easy use during modeling 
baseline_data_table1 <- baseline_data_table1 %>%
  left_join(distinct(long_meas, patientId, time_obs), by = 'patientId')

# remove patients that received either PI3Kinhibitor or other treatment
removeId <- baseline_data_table1$patientId[baseline_data_table1$therapy_received %in% 
                                             c('PI3Kinhibitor', 'Other')]

baseline_data_table1 <- baseline_data_table1[!(baseline_data_table1$patientId %in% removeId),]
long_meas <- long_meas[!(long_meas$patientId %in% removeId),]

# remove patients with one measurement 
removeId <- table(long_meas$patientId[!is.na(long_meas$PSA)])
removeId <- names(removeId[removeId < 2])
long_meas <- long_meas[!(long_meas$patientId %in% removeId),]

# remove patients with no baseline measurement
removeId <- long_meas$patientId[is.na(long_meas$PSA) & long_meas$date_lab == long_meas$baseline_date]
long_meas <- long_meas[!(long_meas$patientId %in% removeId),]

rm(removeId)

# 
psa_long <- psa_long[psa_long$patientId %in% long_meas$patientId,]

## training datasets ----
psa_long_train <- psa_long %>% 
  filter(locked == "ARPI-all") %>% 
  mutate(
    prog = replace(prog, is.na(prog), 0),
    time = time_length(difftime(date_lab, date_rand), "months"),
    # for time-reversed plot
    time_reversed = time_length(difftime(date_lab, date_last), "months")
  ) %>% 
  filter(num_psa_value > 1) %>% 
  # # temporary filters
  # filter(therapy_received %in% c("ARSi", "Taxane")) %>% 
  # filter(patientId %in% sub_id) %>% 
  # safely ordering data for joint models
  arrange(id_num, date_lab)

# corresponding data for survival submodel
subjects_train <- filter(subjects_cr, id %in% psa_long_train$id) %>% 
  mutate(prog = replace(prog, is.na(prog), 0)) %>% 
  arrange(id_num)

# corresponding data for longitudinal dataframe
long_meas_train <- filter(long_meas, id %in% psa_long_train$id) %>%
  arrange(id_num, date_lab)

# corresponding data for baseline dataframe
baseline_data_table1_train <- filter(baseline_data_table1, id %in% psa_long_train$id)

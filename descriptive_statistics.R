###########################
### prepare environment ###
###########################

# load data
source('C:/Users/lanbro/OneDrive - Karolinska Institutet/Dokument/Scripts/data_management_RQ2.R')
# source('C:/Users/lanbro/OneDrive - Karolinska Institutet/Dokument/Scripts/data_management_RQ3.R')

# load libraries 
library(tableone) # to create table 1
library(xtable) # create LateX from table 
library(RColorBrewer) # for the plot colours

############################
### numeric descriptives ###
############################

# nr. of patients
length(unique(long_meas_train$patientId)) 

# distribution of therapy given
plot(as.factor(baseline_data_table1_train$therapy_received)) 
table(baseline_data_table1_train$therapy_received)
table(baseline_data_table1_train$therapy_received, 
      baseline_data_table1_train$NLCB_overall)

# nr. of events
table((long_meas_train %>% distinct(patientId, NLCB_overall))$NLCB_overall) 

# counts per visit type
table(long_meas_train$visit) 

# nr. of locked observations
table(long_meas_train$locked) 

# summary of the follow-up times 
table(round((arrange(long_meas_train, patientId, desc(time)) %>%
               filter(!duplicated(patientId)))$time))
summary((arrange(long_meas_train, patientId, desc(time)) %>%
           filter(!duplicated(patientId)))$time)
hist((arrange(long_meas_train, patientId, desc(time)) %>%
        filter(!duplicated(patientId)))$time, 25)

##################################
### Metastases differentiation ### 
##################################

# create dataframe with metastasis information
metastasis_info <- baseline_cr[baseline_cr$patientId %in% unique(long_meas_train$patientId),] %>%
  select(patientId,
         bone_disease = 'E1_F3_boneDisease0', 
         lymph_disease = 'E1_F3_LNdisease', 
         lymph_pelvic = 'E1_F3_whichLNsite_0', 
         lymph_retroperitoneal = 'E1_F3_whichLNsite_1', 
         lymph_mediastinal = 'E1_F3_whichLNsite_2',
         lymph_thoracic = 'E1_F3_whichLNsite_3', 
         lymph_other = 'E1_F3_whichLNsite_4', 
         visc_disease = 'E1_F3_VisceralDisease', 
         visc_liver = 'E1_F3_VisceralSpec_0',
         visc_lung = 'E1_F3_VisceralSpec_1',  
         visc_adrenalgland = 'E1_F3_VisceralSpec_2', 
         visc_CNS = 'E1_F3_VisceralSpec_3', 
         visc_other = 'E1_F3_VisceralSpec_4', 
         visc_other_name = 'E1_F3_VisceralOther',
         locationmet = 'E2_F26_location_metastases', 
         amount_bone_met = 'E2_F26_TC99boneMetsAtBA', 
         amount_LN_met = 'E2_F26_measurableLN', 
         amount_visc_met = 'E2_F26_measurableViscera', 
         other_sites = 'E2_F26_otherSites') %>%
  mutate(bone_disease = ifelse(bone_disease == 1, 'bone', NA),
         lymph_disease = ifelse(lymph_disease == 1, 'lymph', NA),
         visc_disease = ifelse(visc_disease == 1, 'visc', NA))  %>%
  unite(col = 'met_place', bone_disease, lymph_disease, visc_disease, na.rm = TRUE) %>%
  mutate(met_place = ifelse(met_place == '', 'Undetectable', met_place),
         met_place = ifelse(met_place == 'bone', 'Only bone', met_place),
         met_place = ifelse(met_place == 'lymph', 'Only lymph', met_place),
         met_place = ifelse(met_place == 'visc', 'Only visceral', met_place),
         met_place = ifelse(met_place == 'bone_lymph', 'Bone + lymph', met_place),
         met_place = ifelse(met_place == 'bone_visc', 'Bone + visceral', met_place),
         met_place = ifelse(met_place == 'bone_lymph_visc', 'Bone + lymph + visceral', met_place),
         lymph_pelvic = ifelse(is.na(lymph_pelvic), 'Not measured', lymph_pelvic), 
         lymph_retroperitoneal = ifelse(is.na(lymph_retroperitoneal), 'Not measured', lymph_retroperitoneal), 
         lymph_mediastinal = ifelse(is.na(lymph_mediastinal), 'Not measured', lymph_mediastinal), 
         lymph_thoracic = ifelse(is.na(lymph_thoracic), 'Not measured', lymph_thoracic), 
         lymph_other = ifelse(is.na(lymph_other), 'Not measured', lymph_other),
         lymph_pelvic = ifelse(lymph_pelvic == 0, NA, lymph_pelvic), 
         lymph_retroperitoneal = ifelse(lymph_retroperitoneal == 0, NA, lymph_retroperitoneal), 
         lymph_mediastinal = ifelse(lymph_mediastinal == 0, NA, lymph_mediastinal), 
         lymph_thoracic = ifelse(lymph_thoracic == 0, NA, lymph_thoracic), 
         lymph_other = ifelse(lymph_other == 0, NA, lymph_other),
         lymph_pelvic = ifelse(lymph_pelvic == 1, 'pelvic', lymph_pelvic), 
         lymph_retroperitoneal = ifelse(lymph_retroperitoneal == 1, 'retroperitoneal', lymph_retroperitoneal), 
         lymph_mediastinal = ifelse(lymph_mediastinal == 1, 'mediastinal', lymph_mediastinal), 
         lymph_thoracic = ifelse(lymph_thoracic == 1, 'thoracic', lymph_thoracic), 
         lymph_other = ifelse(lymph_other == 1, 'other', lymph_other),
         visc_liver = ifelse(is.na(visc_liver), 'Not measured', visc_liver), 
         visc_lung = ifelse(is.na(visc_lung), 'Not measured', visc_lung), 
         visc_adrenalgland = ifelse(is.na(visc_adrenalgland), 'Not measured', visc_adrenalgland), 
         visc_CNS = ifelse(is.na(visc_CNS), 'Not measured', visc_CNS), 
         visc_other = ifelse(is.na(visc_other), 'Not measured', visc_other),
         visc_liver = ifelse(visc_liver == 0, NA, visc_liver), 
         visc_lung = ifelse(visc_lung == 0, NA, visc_lung), 
         visc_adrenalgland = ifelse(visc_adrenalgland == 0, NA, visc_adrenalgland), 
         visc_CNS = ifelse(visc_CNS == 0, NA, visc_CNS), 
         visc_other = ifelse(visc_other == 0, NA, visc_other),
         visc_liver = ifelse(visc_liver == 1, 'liver', visc_liver), 
         visc_lung = ifelse(visc_lung == 1, 'lung', visc_lung), 
         visc_adrenalgland = ifelse(visc_adrenalgland == 1, 'adrenalgland', visc_adrenalgland), 
         visc_CNS = ifelse(visc_CNS == 1, 'CNS', visc_CNS), 
         visc_other = ifelse(visc_other == 1, 'other', visc_other),
         amount_bone_met = as.character(amount_bone_met),
         amount_bone_met = ifelse(is.na(amount_bone_met), 'Not evaluated', amount_bone_met),
         amount_bone_met = as.factor(amount_bone_met)) %>%
  left_join(distinct(baseline_data_table1_train, patientId, therapy_received,
                     event = NLCB_overall)) 
metastasis_info$visc_other_name <- metastasis_info$locationmet <- metastasis_info$patientId <- NULL

# distribution of metastasis types
table(metastasis_info$met_place)
table(metastasis_info$met_place, metastasis_info$therapy_received)

# plot with number of bone metastases
## dataframe with summary of number of bone metastases
bone_count <- metastasis_info[grepl('bone', metastasis_info$met_place, ignore.case = T),] %>%
  select(amount_bone_met, met_place) %>%
  group_by(amount_bone_met, met_place) %>%
  mutate(Freq = length(amount_bone_met),
         Freq = ifelse(is.na(Freq), 0, Freq),
         amount_bone_met = factor(amount_bone_met, levels = c('Undetected', '0-1', '2-4', '5-9', '10-19','>=20', 'Diffuse', 'Not evaluated'))) %>%
  rename(Var1 = amount_bone_met, 
         Var2 = met_place)
bone_count <- unique(bone_count)

## dataframe with summary of number of bone metastases, split by treatment
bone_count_treat <- metastasis_info[grepl('bone', metastasis_info$met_place, ignore.case = T),] %>%
  select(amount_bone_met, met_place, therapy_received) %>%
  group_by(amount_bone_met, met_place, therapy_received) %>%
  mutate(Freq = length(amount_bone_met),
         Freq = ifelse(is.na(Freq), 0, Freq),
         amount_bone_met = factor(amount_bone_met, levels = c('Undetected', '0-1', '2-4', '5-9', '10-19','>=20', 'Diffuse', 'Not evaluated'))) %>%
  rename(Var1 = amount_bone_met, 
         Var2 = met_place)
bone_count_treat <- unique(bone_count_treat)

## dataframe with summary of number of bone metastases, split by event
bone_count_event <- metastasis_info[grepl('bone', metastasis_info$met_place, ignore.case = T),] %>%
  select(amount_bone_met, met_place, event) %>%
  group_by(amount_bone_met, met_place, event) %>%
  mutate(Freq = length(amount_bone_met),
         Freq = ifelse(is.na(Freq), 0, Freq),
         amount_bone_met = factor(amount_bone_met, levels = c('Undetected', '0-1', '2-4', '5-9', '10-19','>=20', 'Diffuse', 'Not evaluated'))) %>%
  rename(Var1 = amount_bone_met, 
         Var2 = met_place)
bone_count_event <- unique(bone_count_event)

## plot
bone_count_plot <- ggplot(bone_count, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(angle=45, hjust = 1, size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab(NULL) +
  ylab('Frequency') +
  scale_fill_manual(values = brewer.pal(name = 'Paired', n = 12)[c(1,3,5,7)]) + 
  scale_x_discrete(labels = c('Undetected', '0-1', '2-4', '5-9', '10-19', 
                              expression(phantom(x) >=20), 'Diffuse', 'Not evaluated')) +
  guides(fill = guide_legend(title = 'Metastasis location', nrow = 2))

bone_count_plot_treat <- ggplot(bone_count_treat, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") + # position="fill"
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(angle=45, hjust = 1, size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab(NULL) +
  ylab('Frequency') +
  scale_fill_manual(values = brewer.pal(name = 'Paired', n = 12)[c(1,3,5,7)]) + 
  scale_x_discrete(labels = c('Undetected', '0-1', '2-4', '5-9', '10-19', 
                              expression(phantom(x) >=20), 'Diffuse', 'Not evaluated')) +
  guides(fill = guide_legend(title = 'Metastasis location', nrow = 2)) +
  facet_wrap(~ therapy_received)

bone_count_plot_event <- ggplot(bone_count_event, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") + # position="fill"
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(angle=45, hjust = 1, size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=18)) +
  xlab(NULL) +
  ylab('Frequency') +
  scale_fill_manual(values = brewer.pal(name = 'Paired', n = 12)[c(1,3,5,7)]) + 
  scale_x_discrete(labels = c('Undetected', '0-1', '2-4', '5-9', '10-19', 
                              expression(phantom(x) >=20), 'Diffuse', 'Not evaluated')) +
  guides(fill = guide_legend(title = 'Metastasis location', nrow = 2)) +
  facet_wrap(~ event)

# export plots 
pdf('C:/Users/lanbro/OneDrive - Karolinska Institutet/Dokument/Figures/Descriptive/bone_count.pdf', height = 5, width = 8)
print(bone_count_plot)
dev.off()

pdf('C:/Users/lanbro/OneDrive - Karolinska Institutet/Dokument/Figures/Descriptive/bone_count_treat.pdf', height = 5, width = 8)
print(bone_count_plot_treat)
dev.off()

pdf('C:/Users/lanbro/OneDrive - Karolinska Institutet/Dokument/Figures/Descriptive/bone_count_event.pdf', height = 5, width = 8)
print(bone_count_plot_event)
dev.off()

###############
### Table 1 ###
###############

# select variables for the table 
listVars <- setdiff(names(baseline_data_table1_train[,]), 
                    c('id', 'patientId', 'date_lab', 
                      'TimeSinceMetastaticDisease', 'TimeSincePCdiagnosis',
                      'therapy_class'))
catVars <- c('country', 'TNMstagingT', 'TNMstagingN', 'TNMstagingM', 'Gleason', 
             'LocationMetastases', 'therapy_received', 'ecog')

# create table 
table1 <- CreateTableOne(data = baseline_data_table1_train, vars = listVars, 
                         factorVars = catVars, includeNA = T)
table1

table1.strat <- CreateTableOne(data = baseline_data_table1_train, vars = listVars, 
                               factorVars = setdiff(catVars, 'therapy_received'), 
                               includeNA = T, strata = c('therapy_received')) 
table1.strat

# make LaTeX file from table 1 
tabAsStringMatrix <- print(table1.strat, nonnormal = setdiff(listVars, catVars),
                           printToggle = FALSE, noSpaces = TRUE)
xtable(tabAsStringMatrix)
# print(xtable(xtable(tabAsStringMatrix), type = 'latex'), file = 'Table1.tex')

# clean environment
rm(table1, listVars, catVars, tabAsStringMatrix)

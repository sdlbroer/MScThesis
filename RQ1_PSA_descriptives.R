# first run the data_management_v6 file

# load libraries 
library(ggplot2) # create plots
library(GGally) # create pairs matrix with ggplot
library(RColorBrewer) # use nice colors

# dataframe without missing PSA values 
long_meas_train_noNA <- long_meas_train %>%
  filter(!is.na(PSA))

############################
### numeric descriptives ###
############################

# nr. of patients
length(unique(long_meas_train$patientId[!is.na(long_meas_train$PSA)])) 

# nr. of repeated measurements
table(distinct(long_meas_train, patientId, num_psa_value)$num_psa_value)
summary(distinct(long_meas_train, patientId, num_psa_value)$num_psa_value)

# summary of the PSA values (incl. nr. missing)
summary(long_meas_train$PSA)
tapply(long_meas_train$PSA, long_meas_train$therapy_received, summary)

# summary of follow-up times
summary((arrange(long_meas_train, patientId, desc(months_in_followup)) %>% 
  filter(!duplicated(patientId)))$months_in_followup)

#######################################
### plot(s) 1: histograms + boxplot ###
#######################################

# histograms
hist1 <- ggplot(long_meas_train_noNA, aes(x=PSA)) + 
  geom_histogram(color=brewer.pal(name="Paired", n = 12)[2], 
                 fill=brewer.pal(name="Paired", n = 12)[1], 
                 bins = 100) + 
  theme_bw() + 
  xlab('PSA (ng/ml)') +
  ylab('Count')

hist2 <- ggplot(long_meas_train_noNA, aes(x=log2PSA)) + 
  geom_histogram(color=brewer.pal(name="Paired", n = 12)[2], 
                 fill=brewer.pal(name="Paired", n = 12)[1], 
                 bins = 100) + 
  theme_bw() + 
  xlab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) + 
  ylab('Count')

# box plot 
box <- ggplot(data = long_meas_train_noNA, aes(x = months_in_followup, y = log2PSA, 
                                   fill = therapy_received)) +
  geom_boxplot() +
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  scale_fill_manual(values = brewer.pal(name="Paired", n = 12)[c(1,3,5,7)]) +
  labs(fill = 'Treatment') +
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) # split on treatment

##################################
### plot(s) 2: spaghetti plots ### 
##################################

# with LOESS smoothing
spag.tot <- ggplot(data = long_meas_train_noNA, aes(x = months_in_followup, y = PSA)) + 
  geom_line(aes(group = patientId, colour = patientId)) +
  geom_smooth(colour = 'black') +
  theme_bw() + 
  theme(legend.position = 'none') +
  scale_colour_manual(values = brewer.pal(name="PuBu", n = 9)[rep(4:9, 31)]) + 
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) # all observations

spag.tot.log <- ggplot(data = long_meas_train_noNA, aes(x = months_in_followup, y = log2PSA)) + 
  geom_line(aes(group = patientId, colour = patientId)) +
  geom_smooth(colour = 'black') +
  scale_colour_manual(values = brewer.pal(name="PuBu", n = 9)[rep(4:9, 31)]) + 
  theme_bw() + 
  theme(legend.position = 'none') +
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) # all observations

spag.treat <- ggplot(data = long_meas_train_noNA, aes(x = months_in_followup, y = log2PSA)) + 
  geom_line(aes(group = patientId, colour = therapy_received)) +
  geom_smooth(colour = 'black') +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values = brewer.pal(name="Paired", n = 12)[c(1,3,5,7)]) +
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) +
  facet_wrap(~ therapy_received) # split on treatment

# reverse time scale 

reverse.event <- ggplot(data = long_meas_train_noNA, aes(x = reverse_time, y = log2PSA)) + 
  geom_line(aes(group = patientId, colour = NLCB_overall)) +
  geom_smooth(colour = 'black') + 
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values = brewer.pal(name="PRGn", n = 11)[c(4,2)]) +
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) +
  facet_wrap(~ NLCB_overall) # split on outcome

#####################################
### plot(s) 3: mean profile plots ### 
#####################################

# dataframe with the mean information per treatment
statsMeanTreat <- long_meas_train_noNA %>%
  group_by(therapy_received, months_in_followup_rounded) %>%
  summarise(
    count = n(),
    meanPSA = mean(log2PSA,na.rm=TRUE),
    sdPSA = sd(log2PSA, na.rm=TRUE),
    sePSA = sdPSA/sqrt(count),
    ci95lower = meanPSA - sePSA*qnorm(0.975),
    ci95upper = meanPSA + sePSA*qnorm(0.975)
  ) %>%
  rename(therapy_received = therapy_received)

long_meas_train_noNA_int <- long_meas_train_noNA %>%
  mutate(treat_num = as.character(as.numeric(therapy_received)))

# mean profile plot
mean.treat_withleg <- ggplot(data = long_meas_train_noNA_int, aes(x = months_in_followup_rounded, y = log2PSA)) +
  geom_line(aes(group = patientId, color = treat_num), alpha = 0.4) + 
  geom_line(data = statsMeanTreat, 
            mapping = aes(x = months_in_followup_rounded, y = meanPSA, 
                          color = therapy_received), linewidth = 0.8) +
  geom_point(data = statsMeanTreat, 
             mapping = aes(x = months_in_followup_rounded, y = meanPSA, 
                           color = therapy_received, shape = therapy_received), size = 2) +
  geom_errorbar(data = statsMeanTreat,
                mapping = aes(x = months_in_followup_rounded, y = meanPSA,
                              ymin = ci95lower, ymax = ci95upper, color = therapy_received), 
                width = 0.5, alpha = 0.5, linewidth = 0.7) +
  scale_color_manual(values = c('grey80', 'grey60', 'grey40', 'grey30',
                                brewer.pal(name="Paired", n = 12)[c(1,3,5,7)]),
                     labels = c('ARSi', 'Taxane', 'PARPi', 'Platinum', 
                                'ARSi', 'Taxane', 'PARPi', 'Platinum')) +
  scale_shape_manual(values = c(15, 16, 17, 18)) + 
  guides(colour = guide_legend(title = 'Treatment',
                               override.aes = list(shape = c(rep(NA, 4), 15, 16, 17, 18),
                                                   color = c('grey80', 'grey60', 'grey40', 'grey30',
                                                             brewer.pal(name="Paired", n = 12)[c(1,3,5,7)]))),
         shape = 'none') + 
  theme_bw() + 
  theme(legend.position = 'bottom') +
  labs(color = 'Treatment') +
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) 

mean.treat_withoutleg <- mean.treat_withleg + theme(legend.position = 'none')

mean.treat.facet <- ggplot(statsMeanTreat, aes(x = months_in_followup_rounded, y = meanPSA, 
                           color = therapy_received, shape = therapy_received)) +
  geom_line(aes(x = months_in_followup_rounded, y = log2PSA, group = patientId, 
                color = ''), 
            data = long_meas_train_noNA) +
  geom_line() + 
  geom_errorbar(aes(ymin = ci95lower, ymax = ci95upper), width = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = c('grey87', brewer.pal(name="Paired", n = 12)[c(1,3,5,7)])) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  guides(shape = 'none') + 
  theme_bw() + 
  theme(legend.position = 'none') +
  labs(color = 'Treatment') +
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) + 
  facet_wrap(~ therapy_received)

rm(statsMeanTreat, long_meas_train_noNA_int)

########################################
### plot(s) 4: correlation structure ###
########################################

# dataframe with rounded measurement times 
meas <- long_meas_train_noNA %>%
  select(id_num, log2PSA, months_in_followup) %>% 
  mutate(months_in_followup = round(months_in_followup, 0)) %>%
  filter(months_in_followup %in% c(0,1,2,4,6,9,12,15))

full.df <- data.frame('id_num' = sort(rep(unique(meas$id_num), 8)),
                      'months_in_followup' = rep(c(0,1,2,4,6,9,12,15), 184))

psa_wide_train <- unique(left_join(full.df, meas)) %>%
  mutate(unique_id = paste(id_num, months_in_followup)) %>%
  filter(!duplicated(unique_id)) %>%
  select(id_num, months_in_followup, log2PSA)

# wide format 
psa_wide_train <- reshape(psa_wide_train, idvar = 'id_num', 
                          timevar = 'months_in_followup', direction = "wide")

paircor <- ggpairs(psa_wide_train, columns = 2:9, 
        upper = list(continuous = wrap("cor", size=3, stars = F)),
        columnLabels = c('Baseline', paste('Month', c(1,2,4,6,9,12,15)))) + 
  #xlab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) +
  #ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) +
  theme_bw() +
  theme(panel.grid = element_blank())

rm(meas,full.df,psa_wide_train)

#####################################################
### plot(s) 5: categorization of PSA trajectories ###
#####################################################

# category 1: flat trajectories
## category 1a: flat throughout
flat <- c('UZG3006', 'UZG3007', 'UZG3025', 'UZG3027', 'OLVA4005', 
          'OLVA4015', 'ZOLG3401', 'CHUL3709', 'ZOLG3404', 'GZA3902', 
          'AHUS5001', 'STAV5108','SG1001', 'SG1004', 
          'SG1005', 'SG1017',  'AK1208', 'AK1212', 'NU1310', 'NU1314', 
          'SU2009', 'SU2016', 'SU2018', 'RY1720', 'RY1727', 'RY1729', 'KAR2402', 
          'KAR2404','KAR2405', 'VX2211', 'VX2212', 'KAR2407', 'KAR2414',
          'FALU1502', 'RY1728', 'NU1317', 'SG1045', 'SG1009',
          'KS1108', 'AZSJ3152')
## category 1b: litte decrease, then flat
decr_to_flat <- c('UZG3004', 'STAV5107', 'SG1031','SG1068', 'KS1114', 'NU1304',
                  'NU1316', 'SU2014', 'RY1734', 'KAR2406','KAR2416', 'SG1035', 
                  'CHUL3712', 'AZGR3501', 'AZSJ3155', 'SG1014', 'SU2030')
## category 1c: flat cup
flatincr <- c('UZG3012', 'CHUL3703', 'JESS3301', 'AHUS5006', 'SG1049', 'KS1101',
              'KS1116', 'KS1118', 'AK1205', 'AK1213', 'AK1221', 'NU1311', 'SU2007', 
              'SU2013', 'RY1710', 'RY1718', 'RY1726', 'VX2204', 'SG1027', 'SG1006',
              'SG1008')
## category 1d: flat into weak increase
flatcup <- c('ZOLG3403', 'STAV5112', 'KS1120', 'KS1121', 'NU1309', 'SG1018', 
             'AK1223', 'KAR2412', 'RY1712', 'RY1707', 'SG1043', 'KS1113', 'VX2207', 
             'RY1717', 'AHUS5010', 'CHUL3713', 'UZG3017', 
             'UZG3021', 'AHUS5009')

# category 2: increasing trajectories
## category 2a: flat, then strong increase
flat_to_incr <- c('UZG3013', 'ALES5302', 'KS1112', 'AK1209', 'SU2008', 'SU2011', 
                  'RY1704', 'RY1711', 'RY1719', 'RY1721', 'NU1302', 'SG1028', 
                  'AK1219', 'OLVA4010')
## category 2b: weak decrease into strong increase
decr_to_INCR <- c('AHUS5008', 'AHUS5012', 'AK1218', 'KAR2411', 'RY1724', 
                  'KS1106', 'NU1312', 'KS1103', 'AZSJ3160')
## category 2c: strong decrease into strong increase
DECR_to_INCR <- c('KS1117', 'SG1025', 'AK1214', 'RY1715', 'SU2015', 'NU1320')
## category 2d: strictly monotonic increase
monotincr <- c('STAV5104', 'SG1019', 'NU1301', 'SU2004', 'RY1708', 
               'SG1007', 'AZSJ3157', 'VX2203')

# category 3: decreasing trajectories
## category 3a: strong decrease into flat
DECR_to_flat <- c('KAR2413', 'UZG3018', 'UZG3028', 'KS1104')
## category 3b: increase into decrease (cap)
incr_to_decr <- c('CHUL3710', 'SG1073', 'NU1307', 'SG1047')

# category 4
## category 4a: s-shaped (decrease, increase, decrease, increase)
sshape <- c('AZGR3504', 'SG1012', 'AZSJ3101', 'UZG3030', 'AK1202', 'AK1220')

# category 5: not enough information 
## category 5a: only two measurements
two_meas <- c('UZG3005', 'OLVA4012', 'AZSLB3603', 'STAV5113', 'SG1015', 'SG1026', 
              'AK1204')
## category 5b: only three measurements
three_meas <- c('UZG3022', 'SG1002', 'KS1123', 'AK1226', 'NU1303', 'NU1306', 
                'SU2029', 'VX2205', 'AK1222',  # flat
                'UZG3023', # decr_to_flat
                'NU1321', 'FALU1505', # flatincr
                'ZOLG3402', 'STAV5103', 'STAV5106', 'AK1224', 'NU1315', 'RY1705', 
                'RY1733', 'VX2209', 'KS1109', 'KS1111', 'STAV5114', 'SG1016', 
                'SG1011', 'SG1010', # monotincr
                'NIKO3201', 'SU2010', # incr_to_decr
                'RY1730') # noncategorizable

# check that all patients are categorized
types <- list(flat, decr_to_flat, flatincr, flatcup, # flat
              monotincr, flat_to_incr, decr_to_INCR, DECR_to_INCR, # increasing
              DECR_to_flat, incr_to_decr, # decreasing
              sshape) # odd-shaped
nontypes <- list(two_meas, three_meas)

setdiff(long_meas_train_noNA$patientId, c(unlist(types), unlist(nontypes)))
setdiff(c(unlist(types), unlist(nontypes)), long_meas_train_noNA$patientId)

# plot 
pts <- c('UZG3006', 'UZG3007', 'UZG3025', 'UZG3027', 'STAV5107', 'SG1068', 'SG1031', 'UZG3004',  
         'ZOLG3403', 'STAV5112', 'KS1120', 'KS1121', 'UZG3012', 'CHUL3703', 'JESS3301', 'AHUS5006',
         'AZSJ3157', 'KS1111', 'RY1708','SG1010', 'UZG3013', 'ALES5302', 'KS1112', 'AK1209',
         'AHUS5008', 'AHUS5012', 'AK1218', 'RY1724', 'KS1117', 'SG1025', 'AK1214', 'RY1715',
         'CHUL3710', 'SG1073', 'NU1307', 'SG1047', 'KAR2413', 'UZG3018', 'UZG3028',
         'AZGR3504', 'SG1012', 'AZSJ3101', 'UZG3030')

trajectories.df <- left_join(
  select(long_meas_train_noNA, patientId, PSA, log2PSA, months_in_followup, date_lab),
  rbind(data.frame('patientId' = flat, 'trajectory' = 'Flat: type a', 'category' = 'Flat'),
        data.frame('patientId' = decr_to_flat, 'trajectory' = 'Flat: type b', 'category' = 'Flat'),
        data.frame('patientId' = flatincr, 'trajectory' = 'Flat: type c', 'category' = 'Flat'),
        data.frame('patientId' = flatcup, 'trajectory' = 'Flat: type d', 'category' = 'Flat'),
        data.frame('patientId' = monotincr, 'trajectory' = 'Increasing: type d', 'category' = 'increasing'),
        data.frame('patientId' = flat_to_incr, 'trajectory' = 'Increasing: type a', 'category' = 'increasing'),
        data.frame('patientId' = decr_to_INCR, 'trajectory' = 'Increasing: type b', 'category' = 'increasing'),
        data.frame('patientId' = DECR_to_INCR, 'trajectory' = 'Increasing: type c', 'category' = 'increasing'),
        data.frame('patientId' = DECR_to_flat, 'trajectory' = 'Decreasing: type a', 'category' = 'decreasing'),
        data.frame('patientId' = incr_to_decr, 'trajectory' = 'Decreasing: type b', 'category' = 'decreasing'),
        data.frame('patientId' = sshape, 'trajectory' = 'S-shape', 'category' = 'odd'),
        data.frame('patientId' = two_meas, 'trajectory' = 'cat5a_2meas', 'category' = 'too_little_info'),
        data.frame('patientId' = three_meas, 'trajectory' = 'cat5b_3meas', 'category' = 'too_little_info')),
  by = 'patientId'
) %>% 
  mutate(trajectory = ifelse(patientId %in% c('KS1111', 'SG1010'), 'Increasing: type d', trajectory),
         category = ifelse(patientId %in% c('KS1111', 'SG1010'), 'increasing', category),
         trajectory = factor(trajectory, levels = c(
           'Flat: type a', 'Increasing: type a', 'Decreasing: type a',
           'Flat: type b', 'Increasing: type b', 'Decreasing: type b',
           'Flat: type c', 'Increasing: type c', 'S-shape',
           'Flat: type d', 'Increasing: type d'))) %>%
  filter(!is.na(log2PSA)) %>%
  filter(patientId %in% pts) %>%
  arrange(trajectory, patientId, date_lab) %>% 
  mutate(patientId = factor(patientId, levels = pts))

PSA_traj_types <- ggplot(trajectories.df[trajectories.df$patientId %in% pts,], 
       aes(x = months_in_followup, y = PSA, group = patientId, colour = patientId)) +
  facet_wrap(~ trajectory, nrow = 4, ncol = 3) +
  geom_point() +
  geom_line() +
  coord_cartesian(ylim = c(0, 800),
                  xlim = c(0, 36)) +
  xlab('Time in follow-up (months)') +
  ylab('PSA (ng/ml)') +
  theme_bw() +
  theme(legend.position="none") +
  scale_color_manual(values = brewer.pal(name="Paired", n = 12)[rep(1:4, 11)])

rm(flat, decr_to_flat, flatincr, flatcup, monotincr, flat_to_incr, decr_to_INCR, 
   DECR_to_INCR, DECR_to_flat, incr_to_decr, sshape, two_meas, three_meas,
   types, nontypes, trajectories.df, pts, df.int)

####################
### export plots ### 
#################### 

plots <- list(hist1, hist2, box, spag.tot, spag.tot.log, spag.treat, reverse.event, 
           mean.treat_withoutleg, mean.treat_withleg, mean.treat.facet, paircor, 
           PSA_traj_types) 

filenames <- c('hist_PSA', 'hist_PSA_log', 'box_PSA_treat', 'spaghetti_PSA',
               'spaghetti_PSA_log', 'spaghetti_PSA_log_treat', 'reverse_PSA_event',
               'mean_PSA_treat', 'mean_PSA_treat_legend', 'mean_PSA_treat_facet', 'pairs',
               'trajectories_PSA')

for (i in 1:length(plots)){  
  file_name = paste("C:/Users/lanbro/Documents/Figures/Q1/", filenames[i], ".pdf", sep="")
  pdf(file_name, height=5,width=8)
  print(plots[[i]])
  dev.off()
}

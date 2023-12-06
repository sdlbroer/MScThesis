###########################
### prepare environment ###
###########################

# load data
source('C:/Users/lanbro/Documents/Scripts/data_management_RQ2.R')

# load libraries
library('ggpubr')

################################
### disregard missing values ### 
################################

trajectories.df <- long_meas_train %>% filter(!is.na(PSA))

############################################
### plot all individual PSA trajectories ### 
############################################

#ptids <- unique(trajectories.df$patientId)
#for(i in 1:27){
#  print(arrange(trajectories.df[trajectories.df$patientId %in% ptids[(1 + (i-1)*9):((i-1)*9 + 9)],], patientId, date_lab) %>%
#          ggplot(aes(x = time, y = PSA)) +
#          geom_point() +
#          geom_line() +
#          # coord_cartesian(ylim = c(0, 500)) +
#          theme_bw() +
#          facet_wrap(~ patientId))
#}
#rm(ptids)

###############################################################
### categorize all patients based on their PSA trajectories ### 
###############################################################

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
## category 1c: flat into weak increase
flatcup <- c('ZOLG3403', 'STAV5112', 'KS1120', 'KS1121', 'NU1309', 'SG1018', 
             'AK1223', 'KAR2412', 'RY1712', 'RY1707', 'SG1043', 'KS1113', 'VX2207', 
             'RY1717', 'AHUS5010', 'CHUL3713', 'UZG3017', 
             'UZG3021', 'AHUS5009')
## category 1d: flat cup
flatincr <- c('UZG3012', 'CHUL3703', 'JESS3301', 'AHUS5006', 'SG1049', 'KS1101',
              'KS1116', 'KS1118', 'AK1205', 'AK1213', 'AK1221', 'NU1311', 'SU2007', 
              'SU2013', 'RY1710', 'RY1718', 'RY1726', 'VX2204', 'SG1027', 'SG1006',
              'SG1008')

# category 2: increasing trajectories
## category 2a: strictly monotonic increase
monotincr <- c('STAV5104', 'SG1019', 'NU1301', 'SU2004', 'RY1708', 
               'SG1007', 'AZSJ3157', 'VX2203')
## category 2b: flat, then strong increase
flat_to_incr <- c('UZG3013', 'ALES5302', 'KS1112', 'AK1209', 'SU2008', 'SU2011', 
                  'RY1704', 'RY1711', 'RY1719', 'RY1721', 'NU1302', 'SG1028', 
                  'AK1219', 'OLVA4010')
## category 2c: weak decrease into strong increase
decr_to_INCR <- c('AHUS5008', 'AHUS5012', 'AK1218', 'KAR2411', 'RY1724', 
                  'KS1106', 'NU1312', 'KS1103', 'AZSJ3160')
## category 2d: strong decrease into strong increase
DECR_to_INCR <- c('KS1117', 'SG1025', 'AK1214', 'RY1715', 'SU2015', 'NU1320')

# category 3: decreasing trajectories
## category 3a: strong decrease into flat
DECR_to_flat <- c('KAR2413', 'UZG3018', 'UZG3028', 'KS1104')
## category 3b: increase into decrease (cap)
incr_to_decr <- c('CHUL3710', 'SG1073', 'NU1307', 'SG1047')

# category 4: no clear direction
## category 4a: s-shaped (decrease, increase, decrease, increase)
sshape <- c('AZGR3504', 'SG1012', 'AZSJ3101', 'UZG3030', 'AK1202', 'AK1220')

# category 5: not enough information 
## category 5a: only two measurements
two_meas <- c('UZG3005', 'OLVA4012', 'AZSLB3603', 'STAV5113', 'SG1015', 'SG1026', 
              'AK1204')
## category 5b: only three measurements
three_meas <- c('UZG3022', 'SG1002', 'KS1123', 'AK1226', 'NU1303', 'NU1306', 
                'SU2029', 'VX2205', 'AK1222',  # flat,
                'UZG3023', # decr_to_flat
                'NU1321', 'FALU1505', # flatincr
                'ZOLG3402', 'STAV5103', 'STAV5106', 'AK1224', 'NU1315', 'RY1705', 
                'RY1733', 'VX2209', 'KS1109', 'KS1111', 'STAV5114', 'SG1016', 
                'SG1011', 'SG1010', # monotincr 
                'RY1730', # DECR_to_INCR
                'NIKO3201', 'SU2010') # incr_to_decr

###############################################
### check that all patients are categorized ### 
###############################################

types <- list(flat, decr_to_flat, flatincr, flatcup, # flat
              monotincr, flat_to_incr, decr_to_INCR, DECR_to_INCR, # increasing
              DECR_to_flat, incr_to_decr, # decreasing
              sshape) # odd-shaped

nontypes <- list(two_meas, three_meas)

setdiff(long_meas_train$patientId, c(unlist(types), unlist(nontypes)))
setdiff(c(unlist(types), unlist(nontypes)), long_meas_train$patientId)

#######################################################
### plot three individual trajectories per category ### 
#######################################################

plottitles <- c('Flat', 'Decrease into flat', 'Flat increase', 'Cup-shape: flat', # flat
                'Strictly monotic increase', 'Flat into strong increase', 
                'Cup-shape: weak decrease', 'Cup-shape: pronounced', # increasing
                'Strong decrease into flat', 'Cap-shape', # decreasing
                'S-shape')

# without lowess smoothing
ggs <- vector(mode = 'list', length = length(types))

for(i in 1:length(types)){
  ggs[[i]] <- arrange(trajectories.df[trajectories.df$patientId %in% types[[i]],], 
                      patientId, date_lab) %>%
    ggplot(aes(x = time, y = PSA, group = patientId, colour = patientId)) +
    geom_point() +
    geom_line() +
    coord_cartesian(ylim = c(0, 700),
                    xlim = c(0, 36)) +
    xlab('Time in follow-up (months)') +
    ylab('PSA') +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle(plottitles[i])
}

ggarrange(ggs[[1]], ggs[[2]], ggs[[3]], ggs[[4]], legend = F) # flat
ggarrange(ggs[[5]], ggs[[6]], ggs[[7]], ggs[[8]], legend = F) # increasing
ggarrange(ggs[[9]], ggs[[10]], legend = F) # decreasing
ggarrange(ggs[[11]], legend = F) # odd-shaped

# with lowess smoothing
ggs.lowess <- vector(mode = 'list', length = length(types))

for(i in 1:length(types)){
  ggs.lowess[[i]] <- arrange(trajectories.df[trajectories.df$patientId %in% types[[i]],], 
                      patientId, date_lab) %>%
    ggplot(aes(x = time, y = PSA)) +
    geom_line(aes(colour = patientId), alpha = .3) +
    geom_point(aes(colour = patientId), alpha = .3) +
    geom_smooth(colour = 'black') +
    coord_cartesian(ylim = c(0, 700),
                    xlim = c(0, 36)) +
    xlab('Time in follow-up (months)') +
    ylab('PSA') +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle(plottitles[i])
}

ggarrange(ggs.lowess[[1]], ggs.lowess[[2]], ggs.lowess[[3]], ggs.lowess[[4]], legend = F) # flat
ggarrange(ggs.lowess[[5]], ggs.lowess[[6]], ggs.lowess[[7]], ggs.lowess[[8]], legend = F) # increasing
ggarrange(ggs.lowess[[9]], ggs.lowess[[10]], legend = F) # decreasing
ggarrange(ggs.lowess[[11]], legend = F) # odd-shaped

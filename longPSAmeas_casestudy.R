library('ggpubr')

################################
### disregard missing values ### 
################################

trajectories.df <- long_meas[!is.na(long_meas$PSA),]

############################################
### plot all individual PSA trajectories ### 
############################################

#ptids <- unique(trajectories.df$patientId)
#for(i in 1:27){
#  print(arrange(trajectories.df[trajectories.df$patientId %in% ptids[(1 + (i-1)*9):((i-1)*9 + 9)],], patientId, date_lab) %>%
#          ggplot(aes(x = months_in_followup, y = PSA)) +
#          geom_point() +
#          geom_line() +
#          # coord_cartesian(ylim = c(0, 500)) +
#          theme(panel.grid = element_blank()) +
#          theme_bw() +
#          facet_wrap(~ patientId))
#}
#rm(ptids)

###############################################################
### categorize all patients based on their PSA trajectories ### 
###############################################################

# category 1: flat trajectories

## category 1a: flat throughout
flat <- c('UZG3006', 'UZG3007', 'UZG3025', 'UZG3027', 'UZG3034','OLVA4005', 
          'OLVA4015', 'ZOLG3401', 'CHUL3709', 'ZOLG3404', 'CHUL3716', 'GZA3902', 
          'GZA3910', 'AHUS5001', 'STAV5108', 'ALES5304', 'SG1001', 'SG1004', 
          'SG1005', 'SG1017', 'SG1080', 'AK1208', 'AK1212', 'NU1310', 'NU1314', 
          'SU2009', 'SU2016', 'SU2018', 'RY1720', 'RY1727', 'RY1729', 'KAR2402', 
          'KAR2404','KAR2405', 'VX2211', 'VX2212', 'KAR2407', 'KAR2414', 'BOR1801',
          'FALU1502', 'KAL1609', 'KAL1612', 'RY1728', 'NU1317', 'SG1045', 'SG1009',
          'KS1108', 'AZSJ3152', 'STAV5119')
#flat3 <- c('UZG3022', 'ZOLG3410', 'CHUL3721', 'GZA3908', 'SG1002', 'KS1123', 
#           'AK1226', 'NU1303', 'NU1306', 'SU2029', 'VX2205', 'BOR1802', 'KAL1608', 
#           'AK1222')

## category 1b: litte decrease, then flat
decr_to_flat <- c('UZG3004', 'STAV5107', 'STAV5121', 'SG1031','SG1068', 'SG1078', 
                  'SG1083', 'KS1114', 'NU1304','NU1316', 'SU2014', 'SU2057', 
                  'RY1734', 'KAR2406','KAR2416', 'AZSJ3164', 'SG1035', 'CHUL3712',
                  'AZGR3501', 'AZSJ3155', 'SG1014', 'SU2030', 'KAR2418')
#decr_to_flat3 <- c('UZG3023', 'UHBA7012', 'GZA3919')

## category 1c: flat into weak increase
flatcup <- c('ZOLG3403', 'STAV5112', 'KS1120', 'KS1121', 'NU1309', 'SG1018', 
             'AK1223', 'KAR2412', 'RY1712', 'RY1707', 'SG1043', 'KS1113', 'VX2207', 
             'STAV5116', 'RY1735', 'RY1717', 'AHUS5010', 'CHUL3713', 'UZG3017', 
             'UZG3021', 'AHUS5009', 'SG1091')

## category 1d: flat cup
flatincr <- c('UZG3012', 'CHUL3703', 'JESS3301', 'AHUS5006', 'SG1049', 'KS1101',
              'SG1090', 'KS1116', 'KS1118', 'AK1205', 'AK1213', 'AK1221', 'NU1311',
              'SU2007', 'SU2013', 'RY1710', 'RY1718', 'RY1726', 'VX2204', 'SG1027')

# category 2: increasing trajectories

## category 2a: strictly monotonic increase

monotincr <- c('STAV5104', 'SG1019', 'AK1206', 'NU1301', 'SU2004', 'RY1708', 
               'SG1007', 'AZSJ3157')
#monotincr3 <- c('ZOLG3402', 'ZOLG3409', 'STAV5103', 'STAV5106', 'STAV5123', 
#                'ALES5307', 'AK1207', 'AK1224', 'NU1315', 'RY1701', 'RY1705', 
#                'RY1733', 'VX2209', 'KAL1607')

## category 2b: flat, then strong increase
flat_to_incr <- c('UZG3013', 'ALES5302', 'KS1112', 'AK1209', 'SU2008', 'SU2011', 
                  'RY1704', 'RY1711', 'RY1719', 'RY1721', 'NU1302', 'SG1028', 
                  'AK1219', 'OLVA4010', 'KS1119')
#flat_to_incr3 <- c('NU1321', 'FALU1505')

## category 2c: weak decrease into strong increase
decr_to_INCR <- c('NU1302', 'AHUS5008', 'AHUS5012', 'AK1218', 'KAR2411', 
                  'RY1724', 'KS1106', 'NU1312', 'KS1103', 'AZSJ3160')
#decr_to_INCR3 <- c('GZA3911')

## category 2d: strong decrease into strong increase
DECR_to_INCR <- c('KS1117', 'SG1025', 'UHBA7006', 'AK1214', 'RY1715', 'SU2015', 
                  'FALU1506', 'NU1320')

# category 3: decreasing trajectories

## category 3a: strong decrease into flat
DECR_to_flat <- c('KAR2413', 'UZG3018', 'UZG3028', 'KS1104', 'TROM5202')

## category 3b: increase into decrease (cap)
incr_to_decr <- c('CHUL3710', 'SG1073', 'KAR2417', 'NU1307')
#incr_to_decr3 <- c('NIKO3201', 'SU2010')

# category 4

## category 4a: s-shaped (decrease, increase, decrease, increase)
sshape <- c('AZGR3504', 'SG1012', 'AZSJ3101', 'UZG3030', 'AK1225')

## category 4b: not categorizable
noncategorizable <- c('AK1202', 'AK1220', 'VX2203', 'SG1047')
#noncategorizable3 <- c('RY1730')

# category 5: not enough information 

## category 5a: only two measurements
two_meas <- c('AZSJ3105', 'UZG3005', 'DAMO3804', 'OLVA4012', 'AZSLB3603', 
              'STAV5113', 'STAV5118', 'STAV5126', 'SG1015', 'SG1026', 'AK1204')

## category 5b: only three measurements
three_meas <- c('UZG3023', 'UHBA7012', 'GZA3919', # decr_to_flat
                
                'UZG3022', 'ZOLG3410', 'CHUL3721', 'GZA3908', 'SG1002', 'KS1123', 
                'AK1226', 'NU1303', 'NU1306', 'SU2029', 'VX2205', 'BOR1802', 
                'KAL1608', 'AK1222',  # flat
                
                'ZOLG3402', 'ZOLG3409', 'STAV5103', 'STAV5106', 'STAV5123', 
                'ALES5307', 'AK1207', 'AK1224', 'NU1315', 'RY1701', 'RY1705', 
                'RY1733', 'VX2209', 'KAL1607', 'KS1109', 'KS1111', 'STAV5114', 
                'SG1016', 'SG1011', # monotincr
                
                'NIKO3201', 'SU2010', # incr_to_decr
                'NU1321', 'FALU1505', # flatincr
                'GZA3911', # decr_to_INCR
                'RY1730') # noncategorizable

###############################################
### check that all patients are categorized ### 
###############################################

types <- list(flat, decr_to_flat, flatincr, flatcup, # flat
              monotincr, flat_to_incr, decr_to_INCR, DECR_to_INCR, # increasing
              DECR_to_flat, incr_to_decr, # decreasing
              sshape, noncategorizable) # odd-shaped

nontypes <- list(two_meas, three_meas)

setdiff(trajectories.df$patientId, c(unlist(types), unlist(nontypes)))
setdiff(c(unlist(types), unlist(nontypes)), trajectories.df$patientId)

#######################################################
### plot three individual trajectories per category ### 
#######################################################

plottitles <- c('Flat', 'Decrease into flat', 'Flat increase', 'Cup-shape: flat', # flat
                'Strictly monotic increase', 'Flat into strong increase', 
                'Cup-shape: weak decrease', 'Cup-shape: pronounced', # increasing
                'Strong decrease into flat', 'Cap-shape', # decreasing
                'S-shape', 'Non-categorized')

# without lowess smoothing
ggs <- vector(mode = 'list', length = length(types))

for(i in 1:length(types)){
  ggs[[i]] <- arrange(trajectories.df[trajectories.df$patientId %in% types[[i]][1:4],], 
                      patientId, date_lab) %>%
    ggplot(aes(x = months_in_followup, y = PSA, group = patientId, colour = patientId)) +
    geom_point() +
    geom_line() +
    coord_cartesian(ylim = c(0, 700),
                    xlim = c(0, 36)) +
    xlab('Months') +
    ylab('PSA') +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle(plottitles[i])
}

ggarrange(ggs[[1]], ggs[[2]], ggs[[3]], ggs[[4]], legend = F) # flat
ggarrange(ggs[[5]], ggs[[6]], ggs[[7]], ggs[[8]], legend = F) # increasing
ggarrange(ggs[[9]], ggs[[10]], legend = F) # decreasing
ggarrange(ggs[[11]], ggs[[12]], legend = F) # odd-shaped

# with lowess smoothing
ggs.lowess <- vector(mode = 'list', length = length(types))

for(i in 1:length(types)){
  ggs.lowess[[i]] <- arrange(trajectories.df[trajectories.df$patientId %in% types[[i]],], 
                      patientId, date_lab) %>%
    ggplot(aes(x = months_in_followup, y = PSA)) +
    geom_line(aes(colour = patientId), alpha = .3) +
    geom_point(aes(colour = patientId), alpha = .3) +
    geom_smooth(colour = 'black') +
    coord_cartesian(ylim = c(0, 700),
                    xlim = c(0, 36)) +
    xlab('Months') +
    ylab('PSA') +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle(plottitles[i])
}

ggarrange(ggs.lowess[[1]], ggs.lowess[[2]], ggs.lowess[[3]], ggs.lowess[[4]], legend = F) # flat
ggarrange(ggs.lowess[[5]], ggs.lowess[[6]], ggs.lowess[[7]], ggs.lowess[[8]], legend = F) # increasing
ggarrange(ggs.lowess[[9]], ggs.lowess[[10]], legend = F) # decreasing
ggarrange(ggs.lowess[[11]], ggs.lowess[[12]], legend = F) # odd-shaped

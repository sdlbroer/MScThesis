###########################
### prepare environment ###
###########################

# load data
source('C:/Users/lanbro/Documents/Scripts/data_management_RQ2.R')

# load models 
source('C:/Users/lanbro/Documents/Scripts/RQ2_explanatory/explanatory_model_JMbayes_v22.R')
models <- list(model1_longit, model2_longit, model3_longit, model4_longit,
               model5_longit)

#################
### Cox model ###
#################

# proportional hazards assumptions 
cox.zph(model_surv)

# influential diagnostics
outliers.Cox <- ggcoxdiagnostics(model_surv, type = 'deviance',
                                 linear.predictions = FALSE, ggtheme = theme_bw(),
                                 xlab = 'Patient ID', ylab = 'Residuals (deviance)')

###################
### mixed model ###
###################

# linear relationship between dependent and independent variable
linrel <- vector(mode = 'list', length = 5)

for(i in 1:5){
  linrel[[i]] <- ggplot(data.frame(months_in_followup = models[[i]]$data$months_in_followup, 
                                   pearson = residuals(models[[i]],type = 'pearson')),
                        aes(x = months_in_followup,y = pearson)) +
    geom_point() +
    theme_bw() +
    xlab('Time in follow-up (months)') + 
    ylab('Residuals (pearson)')
}

# constant variance (homoskedasticity)
homosked <- vector(mode = 'list', length = 5)

for(i in 1:5){
  homosked[[i]] <- plot(models[[i]], col = 'black')
}

# normally distributed errors 
for(i in 1:5){
  pdfname <- paste0('C:/Users/lanbro/Documents/Figures/Q2/mm/model', i, '_qqplot.pdf')
  pdf(pdfname, height = 5, width = 8)
  qqPlot(residuals(models[[i]]),
         col.lines = brewer.pal(name = 'Paired', n = 12)[2],
         xlab = 'Normal quantiles', ylab = 'Standardized residuals')
  dev.off()
}
rm(pdfname)

# outliers and influential points
cooksdist <- vector(mode = 'list', length = 5)

for(i in 1:5){
  infl <- hlm_influence(models[[i]])
  cooksdist[[i]] <- dotplot_diag(infl$cooksd, name = 'cooks.distance', cutoff = 'internal') +
    theme_bw() + 
    theme(legend.position = 'none') + 
    ylab('Cooks distance') + 
    xlab('Count')
}
rm(infl)

# observed vs. fitted
obs_vs_fitted <- vector(mode = 'list', length = 5)

for(i in 1:5){
  obs_vs_fitted[[i]] <- ggplot(data = data.frame(observed = psa_long_train$log2PSA,
                                                 fittedval = fitted(models[[i]])),
                               aes(x = fittedval, y = observed)) + 
    geom_point() + 
    geom_abline(colour = brewer.pal(name = 'Paired', n = 12)[2]) + 
    theme_bw() + 
    xlab('Fitted values') + 
    ylab('Observed values')
}

# plot of random sample of patients with all five models
set.seed(31102023)
random.pts <- sample(unique(psa_long_train$id_num), 4)

plots_longit.df <- psa_long_train %>%
  select(id_num, therapy_received, months_in_followup, log2PSA) %>%
  cbind(fitted(model1_longit)) %>%
  cbind(fitted(model2_longit)) %>%
  cbind(fitted(model3_longit)) %>%
  cbind(fitted(model4_longit)) %>%
  cbind(fitted(model5_longit)) %>%
  filter(id_num %in% random.pts) %>%
  rename(log2PSA.pred1 = 'fitted(model1_longit)',
         log2PSA.pred2 = 'fitted(model2_longit)',
         log2PSA.pred3 = 'fitted(model3_longit)',
         log2PSA.pred4 = 'fitted(model4_longit)',
         log2PSA.pred5 = 'fitted(model5_longit)')

pt_sample_model_fits <- ggplot(plots_longit.df, aes(x = months_in_followup, y = log2PSA)) +
  facet_wrap(~ id_num) +
  geom_point() +
  geom_line() +
  geom_line(data = plots_longit.df, 
            mapping = aes(x = months_in_followup, y = log2PSA.pred1, 
                          color = 'linear', linetype = 'linear2')) +
  geom_line(data = plots_longit.df, 
            mapping = aes(x = months_in_followup, y = log2PSA.pred5, 
                          color = 'quadratic', linetype = 'quadratic2')) +
  geom_line(data = plots_longit.df, 
            mapping = aes(x = months_in_followup, y = log2PSA.pred2, 
                          color = 'ns with 1 knot', linetype = '1knot2')) +
  geom_line(data = plots_longit.df, 
            mapping = aes(x = months_in_followup, y = log2PSA.pred3, 
                          color = 'ns with 2 knots', linetype = '2knot2')) +
  geom_line(data = plots_longit.df, 
            mapping = aes(x = months_in_followup, y = log2PSA.pred4, 
                          color = 'ns with 3 knots', linetype = '3knot2')) +
  theme_bw() +
  xlab('Time in follow-up (months)') +
  ylab(expression(paste(log[2](PSA), phantom(x),'(ng/ml)'))) + 
  theme(legend.position = 'bottom') +
  scale_color_manual(labels=c('linear', 'quadratic', 'ns with 1 knot', 'ns with 2 knots',
                              'ns with 3 knots'),
                     values = c('linear' = brewer.pal(name='Paired', n = 12)[c(4)], 
                                'quadratic' = brewer.pal(name='Paired', n = 12)[c(8)],
                                'ns with 1 knot' = brewer.pal(name='Paired', n = 12)[c(2)],
                                'ns with 2 knots' = brewer.pal(name='Paired', n = 12)[c(6)],
                                'ns with 3 knots' = brewer.pal(name='Paired', n = 12)[c(10)])) +
  scale_linetype_manual(labels=c('linear', 'quadratic', 'ns with 1 knot', 'ns with 2 knots',
                                 'ns with 3 knots'), 
                        values = c('linear2' = 'solid', 
                                   'quadratic2' = 'twodash',
                                   '1knot2' = 'dashed',
                                   '2knot2' = 'dotdash',
                                   '3knot2' = 'longdash')) +
  guides(color = guide_legend(override.aes = list(
    linetype = c('solid', 'twodash', 'dashed', 'dotdash', 'longdash'),
    color = brewer.pal(name = 'Paired', n = 12)[c(4,8,2,6,10)])),
    linetype = 'none') +
  labs(color = 'Model type')
rm(random.pts, plots_longit.df)

# mean trajectories per treatment, split on mixed model
fitted_per_treat_ARSi <- data.frame(
  months_in_followup = seq(0, 35.5, 0.01),
  ARSi1 = predict(model1_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 35.5, 0.01),
                                                                    therapy_received = 'ARSi'), level = 0),
  ARSi2 = predict(model2_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 35.5, 0.01),
                                                                    therapy_received = 'ARSi'), level = 0),
  ARSi3 = predict(model3_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 35.5, 0.01),
                                                                    therapy_received = 'ARSi'), level = 0),
  ARSi4 = predict(model4_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 35.5, 0.01),
                                                                    therapy_received = 'ARSi'), level = 0),
  ARSi5 = predict(model5_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 35.5, 0.01),
                                                                    therapy_received = 'ARSi'), level = 0),
  therapy_received = 'ARSi') 
fitted_per_treat_Taxane <- data.frame(
  months_in_followup = seq(0, 21.5, 0.01),
  Taxane1 = predict(model1_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 21.5, 0.01),
                                                                      therapy_received = 'Taxane'), level = 0),
  Taxane2 = predict(model2_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 21.5, 0.01),
                                                                      therapy_received = 'Taxane'), level = 0),
  Taxane3 = predict(model3_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 21.5, 0.01),
                                                                      therapy_received = 'Taxane'), level = 0),
  Taxane4 = predict(model4_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 21.5, 0.01),
                                                                      therapy_received = 'Taxane'), level = 0),
  Taxane5 = predict(model5_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 21.5, 0.01),
                                                                      therapy_received = 'Taxane'), level = 0),
  therapy_received = 'Taxane') 
fitted_per_treat_PARPi <- data.frame(
  months_in_followup = seq(0, 16, 0.01),
  PARPi1 = predict(model1_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 16, 0.01),
                                                                     therapy_received = 'PARPi'), level = 0),
  PARPi2 = predict(model2_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 16, 0.01),
                                                                     therapy_received = 'PARPi'), level = 0),
  PARPi3 = predict(model3_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 16, 0.01),
                                                                     therapy_received = 'PARPi'), level = 0),
  PARPi4 = predict(model4_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 16, 0.01),
                                                                     therapy_received = 'PARPi'), level = 0),
  PARPi5 = predict(model5_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 16, 0.01),
                                                                     therapy_received = 'PARPi'), level = 0),
  therapy_received = 'PARPi') 
fitted_per_treat_Platinum <- data.frame(
  months_in_followup = seq(0, 12, 0.01),
  Platinum1 = predict(model1_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 12, 0.01),
                                                                        therapy_received = 'Platinum'), level = 0),
  Platinum2 = predict(model2_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 12, 0.01),
                                                                        therapy_received = 'Platinum'), level = 0),
  Platinum3 = predict(model3_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 12, 0.01),
                                                                        therapy_received = 'Platinum'), level = 0),
  Platinum4 = predict(model4_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 12, 0.01),
                                                                        therapy_received = 'Platinum'), level = 0),
  Platinum5 = predict(model5_longit, re.form = NA, newdata = data.frame(months_in_followup = seq(0, 12, 0.01),
                                                                        therapy_received = 'Platinum'), level = 0),
  therapy_received = 'Platinum') 

meanmm_treat <- ggplot(data = psa_long_train, aes(x = months_in_followup, y = log2PSA)) + 
  facet_wrap(~ factor(therapy_received, levels = c('ARSi', 'Taxane', 'PARPi', 'Platinum'))) +
  geom_line(aes(group = id_num, colour = therapy_received)) +
  geom_smooth(aes(linetype = 'dotted'), se = F, colour = 'black') +
  geom_line(aes(y = ARSi1, linetype = 'linear2'), 
            data = fitted_per_treat_ARSi,linewidth = 0.8, colour = 'grey90') + 
  geom_line(aes(y = ARSi2, linetype = 'quadratic2'), 
            data = fitted_per_treat_ARSi,linewidth = 0.8, colour = 'grey70') + 
  geom_line(aes(y = ARSi3, linetype = '1knot2'), 
            data = fitted_per_treat_ARSi,linewidth = 0.8, colour = 'grey50') + 
  geom_line(aes(y = ARSi4, linetype = '2knot2'), 
            data = fitted_per_treat_ARSi,linewidth = 0.8, colour = 'grey30') + 
  geom_line(aes(y = ARSi5, linetype = '3knot2'), 
            data = fitted_per_treat_ARSi,linewidth = 0.8, colour = 'grey10') +
  geom_line(aes(y = Taxane1, linetype = 'linear2'), 
            data = fitted_per_treat_Taxane,linewidth = 0.8, colour = 'grey90') + 
  geom_line(aes(y = Taxane2, linetype = 'quadratic2'), 
            data = fitted_per_treat_Taxane,linewidth = 0.8, colour = 'grey70') + 
  geom_line(aes(y = Taxane3, linetype = '1knot2'), 
            data = fitted_per_treat_Taxane,linewidth = 0.8, colour = 'grey50') + 
  geom_line(aes(y = Taxane4, linetype = '2knot2'), 
            data = fitted_per_treat_Taxane,linewidth = 0.8, colour = 'grey30') + 
  geom_line(aes(y = Taxane5, linetype = '3knot2'), 
            data = fitted_per_treat_Taxane,linewidth = 0.8, colour = 'grey10') +  
  geom_line(aes(y = PARPi1, linetype = 'linear2'), 
            data = fitted_per_treat_PARPi,linewidth = 0.8, colour = 'grey90') + 
  geom_line(aes(y = PARPi2, linetype = 'quadratic2'), 
            data = fitted_per_treat_PARPi,linewidth = 0.8, colour = 'grey70') + 
  geom_line(aes(y = PARPi3, linetype = '1knot2'), 
            data = fitted_per_treat_PARPi,linewidth = 0.8, colour = 'grey50') + 
  geom_line(aes(y = PARPi4, linetype = '2knot2'), 
            data = fitted_per_treat_PARPi,linewidth = 0.8, colour = 'grey30') + 
  geom_line(aes(y = PARPi5, linetype = '3knot2'), 
            data = fitted_per_treat_PARPi,linewidth = 0.8, colour = 'grey10') + 
  geom_line(aes(y = Platinum1, linetype = 'linear2'), 
            data = fitted_per_treat_Platinum,linewidth = 0.8, colour = 'grey90') + 
  geom_line(aes(y = Platinum2, linetype = 'quadratic2'), 
            data = fitted_per_treat_Platinum,linewidth = 0.8, colour = 'grey70') + 
  geom_line(aes(y = Platinum3, linetype = '1knot2'), 
            data = fitted_per_treat_Platinum,linewidth = 0.8, colour = 'grey50') + 
  geom_line(aes(y = Platinum4, linetype = '2knot2'), 
            data = fitted_per_treat_Platinum,linewidth = 0.8, colour = 'grey30') + 
  geom_line(aes(y = Platinum5, linetype = '3knot2'), 
            data = fitted_per_treat_Platinum,linewidth = 0.8, colour = 'grey10') + 
  theme_bw() +
  theme(legend.position = 'bottom') + 
  scale_linetype_manual(labels=c('linear', 'quadratic', 'ns with 1 knot', 'ns with 2 knots',
                                 'ns with 3 knots'), 
                        values = c('linear2' = 'solid', 
                                   'quadratic2' = 'dotted',
                                   '1knot2' = 'dashed',
                                   '2knot2' = 'dotdash',
                                   '3knot2' = 'longdash')) +
  scale_color_manual(values = brewer.pal(name = 'Paired', n = 12)[c(1,3,5,7)]) +
  guides(linetype = guide_legend(override.aes = list(color = c('grey90', 'grey70', 'grey50', 'grey30', 'grey10'),
                                                     linetype = c('solid', 'dotted', 'dashed', 
                                                                  'dotdash', 'longdash'))),
         color = 'none') +
  labs(linetype = 'Model type')

rm(fitted_per_treat_ARSi, fitted_per_treat_Taxane, fitted_per_treat_PARPi, fitted_per_treat_Platinum)

###################
### joint model ###
###################

models_joint <- list(model3_lmm_joint, model3_mm3int_joint, model3_mm2int_joint,
                     model3_mm1int_joint, model3_quad_joint,
                     model4_lmm_joint, model4_mm3int_joint, model4_mm2int_joint,
                     model4_mm1int_joint, model4_quad_joint, 
                     model5_lmm_joint, model5_mm3int_joint, model5_mm2int_joint,
                     model5_mm1int_joint, model5_quad_joint)
models_joint_compare <- list(model1_lmm_joint, model1_mm3int_joint, model1_mm2int_joint,
                             model1_mm1int_joint, model1_quad_joint)

# time-varying proportional hazards
time_var_plot <- vector(mode = 'list', length = 15)

for(i in 1:length(models_joint)){
  x_times <- seq(0.001, 18, length = 500)
  
  if(i <= 5) X <- cbind(1, ns(x_times, k = c(4), B = c(0, 35.5)))
  if(i > 5 & i <= 10) X <- cbind(1, ns(x_times, k = c(2, 4), B = c(0, 35.5)))
  if(i > 10) X <- cbind(1, ns(x_times, k = c(2, 4, 6), B = c(0, 35.5)))
  
  mcmc_alphas <- do.call('rbind', models_joint[[i]]$mcmc$alphas)
  log_hr <- X %*% t(mcmc_alphas)
  log_hr_mean <- rowMeans(log_hr)
  log_hr_low <- apply(log_hr, 1, quantile, probs = 0.025)
  log_hr_upp <- apply(log_hr, 1, quantile, probs = 0.975)
  
  time_var_plot.df <- data.frame('months_in_followup' = x_times, 
                                 'HR_mean' = exp(log_hr_mean), 
                                 'HR_low' = exp(log_hr_low), 
                                 'HR_up' = exp(log_hr_upp))
  
  time_var_plot[[i]] <- ggplot(time_var_plot.df) + 
    geom_line(aes(x = months_in_followup, y = HR_mean), colour = 'red', linewidth = 0.8) + 
    geom_line(aes(x = months_in_followup, y = HR_low), 
              colour = 'black', linetype = 'dashed', linewidth = 0.8) + 
    geom_line(aes(x = months_in_followup, y = HR_up), 
              colour = 'black', linetype = 'dashed', linewidth = 0.8) + 
    geom_hline(yintercept = exp(coef(models_joint_compare[[(i-1)%%5 + 1]])$association), colour = 'red', 
               linetype = 'dashed', linewidth = 0.5) + 
    geom_hline(yintercept = 1, colour = 'grey', linetype = 'dashed', linewidth = 0.5) + 
    ylim(c(0,2)) + 
    theme_bw() + 
    xlab('Time in follow-up (months)') +
    ylab(expression(paste('Hazard ratio of', phantom(x), log[2](PSA)))) 
}

####################
### export plots ### 
####################

# survival model 
plots <- list(outliers.Cox) 

filenames <- c('outliers_Cox')

for (i in 1:length(plots)){  
  file_name = paste('C:/Users/lanbro/Documents/Figures/Q2/surv/', filenames[i], '.pdf', sep='')
  pdf(file_name, height = 5, width = 8)
  print(plots[[i]])
  dev.off()
}

# mixed model
plots <- list(linrel[[1]], linrel[[2]], linrel[[3]], linrel[[4]], linrel[[5]],
              homosked[[1]], homosked[[2]], homosked[[3]], homosked[[4]], homosked[[5]],
              cooksdist[[1]], cooksdist[[2]], cooksdist[[3]], cooksdist[[4]], cooksdist[[5]],
              obs_vs_fitted[[1]], obs_vs_fitted[[2]], obs_vs_fitted[[3]], obs_vs_fitted[[4]], 
              obs_vs_fitted[[5]], pt_sample_model_fits, meanmm_treat) 

filenames <- c('model1_linear_relationship', 'model2_linear_relationship',
               'model3_linear_relationship', 'model4_linear_relationship',
               'model5_linear_relationship',
               'model1_homoskedacity', 'model2_homoskedacity', 'model3_homoskedacity',
               'model4_homoskedacity', 'model1_homoskedacity',
               'model1_Cooks_distance', 'model2_Cooks_distance', 'model3_Cooks_distance', 
               'model4_Cooks_distance', 'model1_Cooks_distance', 
               'model1_obs_vs_fitted', 'model2_obs_vs_fitted', 'model3_obs_vs_fitted', 
               'model4_obs_vs_fitted', 'model1_obs_vs_fitted', 
               'fitted_trajectories_random_pts', 'fitted_trajectories_treat')

for (i in 1:length(plots)){  
  file_name = paste('C:/Users/lanbro/Documents/Figures/Q2/mm/', filenames[i], '.pdf', sep='')
  pdf(file_name, height = 5, width = 8)
  print(plots[[i]])
  dev.off()
}

# joint model 
plots <- list(time_var_plot[[1]], time_var_plot[[2]], time_var_plot[[3]], 
              time_var_plot[[4]], time_var_plot[[5]], time_var_plot[[6]],
              time_var_plot[[7]], time_var_plot[[8]], time_var_plot[[9]], 
              time_var_plot[[10]], time_var_plot[[11]], time_var_plot[[12]],
              time_var_plot[[13]], time_var_plot[[14]], time_var_plot[[15]]) 

filenames <- c('lmm_1knot_timevarying', 'mm1knot_1knot_timevarying',
               'mm2knots_1knot_timevarying', 'mm3knots_1knot_timevarying',
               'mmquad_1knot_timevarying',
               
               'extra_time_varying_prop_hazards', 'model4_time_varying_prop_hazards',
               'model3_time_varying_prop_hazards', 'model2_time_varying_prop_hazards',
               'model5_time_varying_prop_hazards',
               
               'lmm_3knots_timevarying', 'mm1knot_3knots_timevarying',
               'mm2knots_3knots_timevarying', 'mm3knots_time_varying_prop_hazards',
               'mmquad_3knots_timevarying')

for (i in 1:length(plots)){  
  file_name = paste('C:/Users/lanbro/Documents/Figures/Q2/joint/', filenames[i], '.pdf', sep='')
  pdf(file_name, height = 5, width = 8)
  print(plots[[i]])
  dev.off()
}

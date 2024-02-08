# MScThesis

This GitHub repository contains all R-code for the MSc Thesis in Statistics and Data Science of Lana Broer, completed in February 2024, at Leiden University; supervised by dr. Alessio Crippa (Karolinska Insitutet) and dr. Mirko Signorelli (Leiden University).

The data used in this project is part of the ProBio trial (Karolinska Insititutet, Sweden), an international, multi-center, randomized controlled phase 3 platform trial, which aims to assess treatment-biomarker combinations in metastatic prostate cancer.

The repository (categorized by folder) contains:
* **main**
  * **data_management**: loads all relevant dataframes and performs initial preparation of the most used data;
  * **descriptive_statistics**: general descriptive statistics on the patient population;
    
* **RQ1_PSA_dynamics**
  * **PSA_descriptives**: produces all descriptive statistics (including plots) to answer the first research question, which focuses on exploring the Prostate-Specific Antigen (PSA) trajectories;
  * **exploration_trajectories**: explores the different types of PSA trajectories for patients in the ProBio trial;
    
* **RQ2_explanatory_model**
  * **explanatory_model_JMbayes**: fits a number of joint models to explore the relationship between repeated measurements of PSA and the survival outcome No-Longer Clinical Benefit (NLCB);
  * **explanatory_model_timedependentCox**: fits a time-dependent Cox model to explore the relationship between PSA and NLCB;
  * **model_diagnostics**: performs model diagnostics for the survival and longitudinal sub-models, as well as the joint models fitted in the *RQ2_explanatory_model_JMbayes* file;
  * **survival_plots**: performs model diagnostics for the survival and longitudinal sub-models, as well as the joint models fitted in the *RQ2_explanatory_model_JMbayes* file;
    
* **RQ3_predicting_NLCB**
  * **create_folds_prop_function**: adaptation of the create_folds-function of the JMbayes2 package which allows for stratified train/test-splitting of the data on the landmark time;
  * **data_management_RQ3**: updated version of the *data_management* file that includes 2 more months of data;
  * **basic_models** (folder): models that predict NLCB using PSA as a longitudinal biomarker and treatment as a baseline covariate. It contains
    * **SPJM**: fits the predictive shared-parameter joint models on the non-landmarked data;
    * **landmark_LOCF**: fits a landmark model using the last-observation carrief forward approach;
    * **pencal**: fits a two-stage landmark model which uses penalized regression calibration;
  * **multiple_longitudinal_cov** (folder): models that predict NLCB using PSA, LDH and ALP as longitudinal biomarkers and treatment as a baseline covariate. It contains files of the same name as in **basic_models**, in which the predictors were suitably updated;
  * **multiple_baseline_cov** (folder): models that predict NLCB using PSA as a longitudinal biomarker and treatment, previous treatmentline, Gleason score, ECOG and location of metastases as baseline covariates. It contains files of the same name as in **basic_models**, in which the predictors were suitably updated;
  * **plots** (folder): plots 
    * **ALP_LDH_plots**: descriptive statistics of ALP and LDH;
    * **dynamic_prediction_landmark_plots**: create dynamic prediction plots for the landmark times 4, 5 and 6 months;
    * **dynamic_prediction_plots**: create dynamic prediction plots at landmark time 4 months;
    * **performance_measure_plots**: creates plots of the predicted AUC and Brier scores.

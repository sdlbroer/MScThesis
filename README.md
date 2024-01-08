# MScThesis

This GitHub repository contains all R-code for the MSc Thesis in Statistics and Data Science of Lana Broer, completed ..., at Leiden University; supervised by dr. Alessio Crippa (Karolinska Insitutet) and dr. Mirko Signorelli (Leiden University).

The data used in this project is part of the ProBio trial (Karolinska Insititutet, Sweden), an international, multi-center, randomized controlled phase 3 platform trial, which aims to investigate novel treatments for metastatic prostate cancer.

The repository contains 3 folders and ... files:
* **main**
  * **data_management**: loads all relevant dataframes and performs initial preparation of the most used data;
  * **descriptive_statistics**: general descriptive statistics on the patient population;
* **RQ1_PSA_dynamics**
  * **RQ1_PSA_descriptives**: produces all descriptive statistics (including plots) to answer the first research question, which focuses on exploring the Prostate-Specific Antigen (PSA) trajectories;
  * **RQ1_exploration_trajectories**: explores the different types of PSA trajectories for patients in the ProBio trial;
* **RQ2_explanatory_model**
  * **RQ2_explanatory_model_JMbayes**: fits a number of joint models to explore the relationship between repeated measurements of PSA and the survival outcome No-Longer Clinical Benefit (NLCB);
  * **RQ2_explanatory_model_timedependentCox**: fits a time-dependent Cox model to explore the relationship between PSA and NLCB;
  * **RQ2_model_diagnostics**: performs model diagnostics for the survival and longitudinal sub-models, as well as the joint models fitted in the *RQ2_explanatory_model_JMbayes* file;
* **RQ3_predicting_NLCB**
  * **RQ3_data_management**: updated version of the *data_management* file that includes 2 more months of data;
  * **create_folds_prop_function**: adaptation of the create_folds-function of the JMbayes2 package which allows for stratified train/test-splitting of the data on the landmark time;
  * **TBD1**: folder that contains models that predict NLCB using PSA as a longitudinal biomarker and treatment as a baseline covariate. It contains
    * **SPJM**: fits the predictive shared-parameter joint models on the *non-landmarked* data;
    * **SPJM_landmarked**: fits the predictive shared-parameter joint models on the *landmarked* data;
    * **landmark_LMM**: fits a two-stage landmark model in which the longitudinal covariate was captured using a linear mixed-effects model;
    * **landmark_LOCF**: fits a landmark model using the last-observation carrief forward approach;
    * **pencal**: fits a two-stage landmark model which uses penalized regression calibration. 
  * **TBD2**: folder that contains models that predict NLCB using PSA, LDH and ALP as longitudinal biomarkers and treatment as a baseline covariate. It contains files of the same name as in **TBD1**, in which the predictors were suitably updated;
  * **TBD3**: folder that contains models that predict NLCB using PSA as a longitudinal biomarker and treatment, previous treatmentline, Gleason score, ECOG and location of metastases as baseline covariates. It contains files of the same name as in **TBD1**, in which the predictors were suitably updated.


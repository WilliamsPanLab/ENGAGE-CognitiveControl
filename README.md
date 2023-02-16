# ENGAGE Cognitive Control Circuit paper

This repo contains steps and code to replicate analysis included in manuscript entitled "Adaptive Changes in the Cognitive Control Brain Circuit Underlie and Predict Behavioral Outcomes for Depression over Two Years". In thi project, we focused on the cognitive control  circuit as a putative neural mechanism of action for a novel behavioral intervention with five repeat measures over two years and explored the possibility of using early changes in this circuit to predict future treatment outcomes. The whole-brain voxel-wise Linear Mixed Model analysis was conducted in Matlab while the cross-validated prediction was conducted in R.


## High-level steps




## Pre-requisite
### Hardware requirements
All stpes could be done on a standard research computer with reasonable CPUs and RAM.

### Software requirements

#### OS requirements

The analysis was conducted and only tested for running on macOS Mojave (10.14.1) and Monterey (12.2.1).

#### Software
- [Matlab_R2020b](https://www.mathworks.com/products/new_products/release2020b.html) for neuroimaging analysis
  -  Matlab dependencies: [SPM8](https://www.fil.ion.ucl.ac.uk/spm/software/spm8/), [DPABI V6.0_210501](http://rfmri.org/content/dpabi-v60-and-dpabinet-v10-were-released).
- [R version 4.0.5](https://www.r-project.org/) for non-neuroimaging analysis and mediation analysis
  - R dependencies: rio, ggplot2, lme4, tidyverse, sjPlot, coefplot2, performance, see, broom.mixed, kableExtra, janitor, ggeffects, dplyr, gridExtra, qqplotr, emmeans, pbkrtest, knitr,ggpubr, here, table1, psych, broom,lsr, rstatix, formatR, RVAideMemoire, labelled, cowplot, readr, svglite, rmcorr, cowplot, grid, gtable, RColorBrewer, extrafont, corrplot, grDevices,icesTAF, gganimate

## Installation guide
  ```
  cd ${where_you_would_like_to_save_the_code}
  git clone https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl
  ```
  Simulated Demo data has been provided in `ENGAGE-CognitiveControl/demo_data` for running the analysis in R.

## Run code with demo data in R

  1. We recommend to use [Rstudio](https://support--rstudio-com.netlify.app/products/rstudio/download/) to open `Ketamine-FEET-Mediation/LME_mediation/FEET_CADSS_5DASC_Analysis.Rproj`   
  
  2. Open `Ketamine-FEET-Mediation/LME_mediation/RMD/P50_LME_visualization.Rmd` and press `Knit` in Rstudio to generate .html results report. The results should be comparable to `LME_mediation/Results/P50_LME_visualization.html`
  
  3. Open `Ketamine-FEET-Mediation/LME_mediation/RMD/P50_correlation_mediation.Rmd` and press `Knit` in Rstudio to generate .html results report. The results should be comparable to `LME_mediation/Results/P50_correlation_mediation.html`

  If preferred, scripts in Step 2 and 3 could also be `Run` in Rstudio instead of `Knit`. Please note the data was simulated so the results generated from these scripts won't match the manuscript.












## Documentation for methods and corresponding scripts  

#### Demographic Table 1

Demographic and symptom data at baseline was summarized in Table 1 using the below code.

[Generate_Baseline_Demographics.Rmd]()


#### Examining effects of I-CARE on behavioral outcomes over U-CARE and delineating individuals’ behavioral outcome trajectories

[ENGAGE_LME_behavioral_outcome.Rmd](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/fd145b7a7260efa2d74b24a0a4728ddf6321dd23/LME_mediation/R/read_redcap_data.R)

This R script examined the effect of I-CARE over U-CARE on improvingSPSI and SCL-20 including six, 12-, and 24-months data using an LMM. Expanding the parent RAINBOW study (8), we examined the effect of I-CARE over U-CARE on improving SPSI and SCL-20 including six, 12-, and 24-months data using an LMM. This LMM included SPSI or SCL-20 change relative to baseline as the dependent variable, intervention (I-CARE or U-CARE), time (six, 12, or 24 months), and their interaction as the fixed effects and subjects as the random effect. Baseline SPSI or SCL-20 was included as a regressor.


It also include delineating individual's behavioral outcome trajectories shown in Fig. S2. To delineate how an individual’s behavioral outcomes (SCL-20 and SPSI) changed across intervention phases and to decide whether the brain-behavior association should be modeled differently across intervention phases, we plotted individuals’ trajectories and assessed the stability across different time points via correlation coefficients and intra-class coefficient (ICC(2, 1)).

#### whole-brain LMMs in identifying mechanism and predictive markers within the cognitive control circuit

The code for generating whole-brain results in the paper is here: [STEP1_ENGAGE_CCC_lme_mechanism_predict.m](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/f125b3186f867c03f4dc09b90eae07e72044b225/fmri/Task_repeated_ANOVA_spm.m), this generates the brain results in Fig. 2A, 3A, 4A, 4C as well as Fig. S5A, S6-S8.

We utilized the fitlme function of the Statistics and Machine Learning Toolbox in Matlab (https://www.mathworks.com/products/matlab.html) and customized scripts to conduct the voxel-wise whole-brain LMM with clinical measures (here ∆SCL-20 or ∆SPSI) as the dependent variable, as detailed below. The significance of fixed effects was tested under the type III hypotheses using the analysis of variance (ANOVA) function in Matlab. The degrees of freedom for ANOVA were estimated via the Satterthwaite method. All results were corrected multiple comparisons with a voxel threshold of p < 0.001 and a Gaussian random field theory (GRF) familywise error cluster-level correction at p < 0.05 using DPABI V5.1 (http://rfmri.org/dpabi). Data smoothness for GRF correction was estimated on the statistical image following a similar procedure to FSL easythresh.

The script for generating scatter plots for each regions is here: [STEP1_ENGAGE_CCC_lme_mechanism_predict.m](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/f125b3186f867c03f4dc09b90eae07e72044b225/fmri/Task_repeated_ANOVA_spm.m). This script also generated the comparison of models using cognitive control circuit as the independent variables and alternative control models using only baseline characteristics in Table S1.

The script for estimating beta coefficient per each intervention group or timepoint is [ENGAGE_LME_Task_Activation_SCL20.Rmd](...) for SCL-20 depression severity and [ENGAGE_LME_Task_Activation_SPSI.Rmd](....) for SPSI. This script produced Fig. S3 and Table S2.

The script that runs cross-validation for the predicting of 6,12, and 24 months behavioral changes using 2 months cognitive control circuit activation is [ENGAGE_LME_Task_Activation.Rmd]
()
##### Identifying cognitive control circuit as a neural mechanisam engaged by the I-CARE intervention

We constructed two linear mixed models (LMMs) with the change of each behavioral outcome – ∆SPSI and ∆SCL-20 – at six, 12, and 24 months relative to baseline as the dependent variable. Fixed-effect terms included the cognitive control circuit activity change at the same time point relative to baseline  (∆Circuit, quantified as activation change of the NoGo > Go contrast), the interaction of the circuit activity change and intervention groups (I-CARE or U-CARE; Intervention x ∆Circuit), the interaction of circuit change and the time (six, 12, and 24 months; Time x ∆Circuit), and the interaction of all three factors (Intervention x Time x ∆Circuit). The interaction of circuit change and intervention groups (Intervention x ∆Circuit) measured how I-CARE modulated the neural mechanism that underlies the improvement of intervention outcomes. The main effect of cognitive control circuit engagement (∆Circuit) assessed the shared neural mechanism that underlies behavior improvement, regardless of I-CARE or U-CARE. Additionally, these two models examined how intervention-dependent and intervention-independent effects were different between intervention phases via the interaction effects of Time.

##### Identifying cognitive control circuit engagement as a predictor of later behavioral improvement

Similar LMMs were conducted with the same SPSI and SCL-20 outcomes at six, 12, and 24 months relative to baseline as the dependent variable in a voxel-wise whole-brain analysis. Fixed effect circuit terms were circuit change from baseline to two months as the behavioral outcome instead of cognitive control circuit activity change at the same timepoint. Beyond cognitive control circuit engagement change at two months relative to baseline (∆Circuit at two months, quantified as activation change of the NoGo > Go contrast at two months), we also modeled the interaction of the circuit change at two months and intervention groups (I-CARE or U-CARE; Intervention x ∆Circuit2MO), the interaction of circuit change at two months and the timepoint (six, 12, and 24 months; Time x ∆Circuit2MO ), and the interaction of all three factors (Intervention x Time x ∆Circuit2MO). The main effect of cognitive control circuit activity (∆Circuit2MO) identified general neural predictors of future behavior improvement, regardless of I-CARE or U-CARE, while the interaction of circuit change at two months and intervention group (Intervention x ∆Circuit2MO) measured how the prediction of future behavioral outcomes using cognitive control circuit change at two months was different in I-CARE versus U-CARE. Time-dependent effects were also examined with time-involved interaction effects.



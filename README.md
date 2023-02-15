# ENGAGE Cognitive Control Circuit paper

This repo contains steps and code to replicate analysis included in manuscript entitled "Adaptive Changes in the Cognitive Control Brain Circuit Underlie and Predict Behavioral Outcomes for Depression over Two Years". In thi project, we focused on the cognitive control  circuit as a putative neural mechanism of action for a novel behavioral intervention with five repeat measures over two years and explored the possibility of using early changes in this circuit to predict future treatment outcomes.


## High-level steps
- [Pre-requisite](#pre-requisite)
- [Installation guide](#Installation-guide)
- [Run code with demo data in R](#Run-code-with-demo-data-in-R)
- [Documentation for methods and corresponding scripts](#Documentation-for-methods-and-corresponding-scripts)
  - [data preprocessing and preparation](#data-preprocessing-and-preparation)
    - [fMRI data](#fmri-data)
    - [CADSS and 5D-ASC](#cadss-and-5d-asc)
  - [data analysis](#data-analysis)
    - [fMRI Analysis of Variance (ANOVA)](#fmri-analysis-of-variance-(anova))
    - [Linear mixed model analysis for CADSS and 5D-ASC](#linear-mixed-model-analysis-for-cadss-and-5d-asc)
    - [Mediation analysis](#mediation-analysis)

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
  Demo data has been provided in `Ketamine-FEET-Mediation/LME_mediation/csv` for running the analysis in R.

## Run code with demo data in R

  1. We recommend to use [Rstudio](https://support--rstudio-com.netlify.app/products/rstudio/download/) to open `Ketamine-FEET-Mediation/LME_mediation/FEET_CADSS_5DASC_Analysis.Rproj`   
  
  2. Open `Ketamine-FEET-Mediation/LME_mediation/RMD/P50_LME_visualization.Rmd` and press `Knit` in Rstudio to generate .html results report. The results should be comparable to `LME_mediation/Results/P50_LME_visualization.html`
  
  3. Open `Ketamine-FEET-Mediation/LME_mediation/RMD/P50_correlation_mediation.Rmd` and press `Knit` in Rstudio to generate .html results report. The results should be comparable to `LME_mediation/Results/P50_correlation_mediation.html`

  If preferred, scripts in Step 2 and 3 could also be `Run` in Rstudio instead of `Knit`. Please note the data was simulated so the results generated from these scripts won't match the manuscript.












## Documentation for methods and corresponding scripts  
### data preprocessing and preparation
#### fMRI data
1. Preprocessing: [fmriprep-20.2.3.job](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/b1ee4f796c71b1707de6bc68edbf99c3c6c7ff38/fmri/Preprocessing/fmriprep-20.2.3.job)

  Results included in this paper come from preprocessing performed using fMRIPrep 20.2.3, details can be found [here](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/5d5a29331a3e20494244e72c544fee06472ac069/neuroimaging_preprocessing.md).

2. Quality control: [summarizing_motion_spikes](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/b1ee4f796c71b1707de6bc68edbf99c3c6c7ff38/fmri/summarizing_motion_spikes.m)
  Quality control diagnostics included visual inspection of the raw fMRI timeseries for artifacts and signal dropout, and a review of the fMRIprep summary report for each participant. Participants’ data were excluded if more than 25% (37/148) of time points were detected as motion spikes. Volumes with frame-wise displacement >0.5 mm or std DVARS >1.5 are defined as motion spikes. One participant’s data for the 0.05 mg/kg was excluded. One participant was not able to complete brain scans due to nausea under the 0.5 mg/kg condition. One participant was unreachable after the completion of the first two scan visits (placebo and 0.5 mg/kg) and is missing the 0.05 mg/kg data. This resulted in n = 13, 11, and 12 for placebo, 0.05mg/kg and 0.5mg/kg conditions, respectively.

3. Definition of region of interest (ROI): As established in our previous work, we defined ROIs with an automated meta-analysis approach using neurosynth.org. Specifically, we used Neurosynth uniformity (previously called forwardinference) maps with a false discovery rate (FDR) threshold of .01 for the search terms of anterior insula and amygdala. We imposed a restriction that the peak of the ROIs should have a minimum z-score of 6. For the anterior insula, we also excluded voxels with a z-score <5 to keep only the most relevant voxels spatially located in the anterior portion of the insula via visual inspection. For the amygdala, neurosynth maps were restricted by anatomically defined boundaries from the Automated Anatomical Labeling atlas. The established ROIs can be found [here](https://github.com/WilliamsPanLab/2021-masks).

4. Generating activation maps: Preprocessed data were entered into a general linear model at the individual level using [SPM8](https://www.fil.ion.ucl.ac.uk/spm/software/spm8/). Each block of emotional expressions was convolved with a canonical hemodynamic response function, and the blocks were used as regressors in the general linear model, as were motion spikes. Activation maps for threat (fear and anger facial expressions) relative to neutral faces, and for happy relative to neutral faces, were estimated to examine ketamine-induced brain activity change in response to negative and positive emotions.


#### Demographic Table 1
[Generate_Baseline_Demographics.Rmd]()



#### Behavioral outcomes

[read_redcap_data.R](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/fd145b7a7260efa2d74b24a0a4728ddf6321dd23/LME_mediation/R/read_redcap_data.R)

Addressing data missingness: for each subject’s questionnaire data under a certain dose condition, if the missing items were fewer than 10% of the total item numbers of the questionnaire, we replaced missing items with the group mean of that dose condition. This brought the sample size of CADSS to be n = 13, n = 12, n = 13 for placebo, 0.05 mg/kg, and 0.5 mg/kg and of 5D-ASC to be n = 13 for all three drug visits.








### data analysis

#### fMRI: Identifying cognitive control circuit as a neural mechanisam engaged by the I-CARE intervention


We constructed two linear mixed models (LMMs) with the change of each behavioral outcome – ∆SPSI and ∆SCL-20 – at six, 12, and 24 months relative to baseline as the dependent variable. Fixed-effect terms included the cognitive control circuit activity change at the same time point relative to baseline  (∆Circuit, quantified as activation change of the NoGo > Go contrast), the interaction of the circuit activity change and intervention groups (I-CARE or U-CARE; Intervention x ∆Circuit), the interaction of circuit change and the time (six, 12, and 24 months; Time x ∆Circuit), and the interaction of all three factors (Intervention x Time x ∆Circuit). The interaction of circuit change and intervention groups (Intervention x ∆Circuit) measured how I-CARE modulated the neural mechanism that underlies the improvement of intervention outcomes. The main effect of cognitive control circuit engagement (∆Circuit) assessed the shared neural mechanism that underlies behavior improvement, regardless of I-CARE or U-CARE. Additionally, these two models examined how intervention-dependent and intervention-independent effects were different between intervention phases via the interaction effects of Time.

#### fMRI: Identifying cognitive control circuit engagement as a predictor of later behavioral improvement

Similar LMMs were conducted with the same SPSI and SCL-20 outcomes at six, 12, and 24 months relative to baseline as the dependent variable in a voxel-wise whole-brain analysis. Fixed effect circuit terms were circuit change from baseline to two months as the behavioral outcome instead of cognitive control circuit activity change at the same timepoint. Beyond cognitive control circuit engagement change at two months relative to baseline (∆Circuit at two months, quantified as activation change of the NoGo > Go contrast at two months), we also modeled the interaction of the circuit change at two months and intervention groups (I-CARE or U-CARE; Intervention x ∆Circuit2MO), the interaction of circuit change at two months and the timepoint (six, 12, and 24 months; Time x ∆Circuit2MO ), and the interaction of all three factors (Intervention x Time x ∆Circuit2MO). The main effect of cognitive control circuit activity (∆Circuit2MO) identified general neural predictors of future behavior improvement, regardless of I-CARE or U-CARE, while the interaction of circuit change at two months and intervention group (Intervention x ∆Circuit2MO) measured how the prediction of future behavioral outcomes using cognitive control circuit change at two months was different in I-CARE versus U-CARE. Time-dependent effects were also examined with time-involved interaction effects.

We utilized the fitlme function of the Statistics and Machine Learning Toolbox in Matlab (https://www.mathworks.com/products/matlab.html) and customized scripts to conduct the voxel-wise whole-brain LMM with clinical measures (here ∆SCL-20 or ∆SPSI) as the dependent variable. The significance of fixed effects was tested under the type III hypotheses using the analysis of variance (ANOVA) function in Matlab. The degrees of freedom for ANOVA were estimated via the Satterthwaite method. All results were corrected multiple comparisons with a voxel threshold of p < 0.001 and a Gaussian random field theory (GRF) familywise error cluster-level correction at p < 0.05 using DPABI V5.1 (http://rfmri.org/dpabi). Data smoothness for GRF correction was estimated on the statistical image following a similar procedure to FSL easythresh.

The code for generating whole-brain results is here: [Task_repeated_ANOVA_spm.m](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/f125b3186f867c03f4dc09b90eae07e72044b225/fmri/Task_repeated_ANOVA_spm.m), this generates the brain results in Fig. 2A, 3A, 4A, 4C as well as Fig. S5A, S6-S8.

The code for generating scatter plots for each regions is here: 

ENGAGE_CCC_lme_mechanism_predict.m






#### Linear mixed model analysis for CADSS and 5D-ASC

[P50_LME_visualization.Rmd](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/fd145b7a7260efa2d74b24a0a4728ddf6321dd23/LME_mediation/RMD/P50_LME_visualization.Rmd)

To examine dose-dependent effects of ketamine on CADSS-assessed subcomponents of dissociation and 5D-ASC-assessed other ASCs, we used linear mixed effects models (LMMs) with dose (placebo, 0.05mg/kg or 0.5mg/kg) as the fixed effect and participant as a random effect using the [lmer package](https://www.rdocumentation.org/packages/lme4/versions/1.1-31/topics/lmer) in R. Time and dose-by-time interaction were added if applicable (Suppl. Methods). Age and biological sex were included as covariates. We implemented an FDR correction to control for the testing of multiple scale sub-components. For significant dose-dependent effects, post-hoc paired t-tests were also run to compare 0.5 mg/kg versus placebo, 0.05 mg/kg versus placebo, and 0.5 mg/kg versus 0.05 mg/kg, to reveal which drug dose condition drove the effect. 

For fMRI clusters that survived multiple comparison corrections, we extracted the peak voxel activation (fMRI beta estimate) for all three conditions and conducted planned contrasts using paired t-tests between each pair of dose conditions, as we did with the ASC data.


#### Mediation analysis

[P50_correlation_mediation.Rmd](https://github.com/WilliamsPanLab/Ketamine-FEET-Mediation/blob/fd145b7a7260efa2d74b24a0a4728ddf6321dd23/LME_mediation/RMD/P50_correlation_mediation.Rmd)

To address our third objective — to test whether specific aspects of ketamine-induced dissociation and other ASCs mediate the effect of dose on acute changes in neural activity during emotional processing — we utilized the Averaged Causal Mediation Effect mediation model because it is powerful for understanding mechanisms of action with dose as the independent variable (0.5mg/kg versus placebo, the X variable), altered states as the mediators (the M variables), and the neural activity of the anterior insula and amygdala in response to emotional faces as the dependent variables (the Y variables). To test our working hypothesis that ketamine will reduce neural activity reflecting relief of negative affective states, mediators included depersonalization and derealization from the CADSS and blissful state from the 5D-ASC. To test our working hypothesis that ketamine will increase neural activity reflecting exacerbation of negative affective states, mediators included dissociative amnesia from the CADSS, as well as anxiety and impaired control and cognition from the 5D-ASC. Mediation models were implemented using the [mediation package](https://cran.r-project.org/web/packages/mediation/index.html) combined with the [lmer package](https://www.rdocumentation.org/packages/lme4/versions/1.1-31/topics/lmer).





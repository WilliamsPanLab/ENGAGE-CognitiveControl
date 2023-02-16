# ENGAGE Cognitive Control Circuit paper

This repo contains steps and code to replicate analysis included in manuscript entitled "Adaptive Changes in the Cognitive Control Brain Circuit Underlie and Predict Behavioral Outcomes for Depression over Two Years". In thi project, we focused on the cognitive control  circuit as a putative neural mechanism of action for a novel behavioral intervention with five repeat measures over two years and explored the possibility of using early changes in this circuit to predict future treatment outcomes. The whole-brain voxel-wise Linear Mixed Model analysis was conducted in Matlab while the cross-validated prediction was conducted in R. Scrpts are listed based on their order in Results section.


## High-level steps
- [Pre-requisite](#pre-requisite)
- [Installation guide](#Installation-guide)
- [Run code with demo data in R](#Run-code-with-demo-data-in-R)
- [Documentation for methods and corresponding scripts](#Documentation-for-methods-and-corresponding-scripts)


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


## Documentation for methods and corresponding scripts  

#### Schematic diagram in Figure 1

Data demonstrating the modulating effect of intervention or time or their interaction were simulated and visualized using [Figure1_Simulation.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/383702627a7a52139e8558ab67e68e19692c61f6/matlab/Figure1_Simulation.m).

#### Demographic Table 1

Demographic and symptom data at baseline was summarized in Table 1 using the below code.

[Generate_Baseline_Demographics.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/13928c1f5a0de9e3a4bebac7b93c48de1789545d/r_script/Rmd/Generate_Baseline_Demographics.Rmd)


#### Examining effects of I-CARE on behavioral outcomes over U-CARE and delineating individuals’ behavioral outcome trajectories

[ENGAGE_LME_behavioral_outcome.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/13928c1f5a0de9e3a4bebac7b93c48de1789545d/r_script/Rmd/ENGAGE_LME_behavioral_outcome.Rmd)

This R script examined the effect of I-CARE over U-CARE on improvingSPSI and SCL-20 including six, 12-, and 24-months data using an LMM. Expanding the parent RAINBOW study, we examined the effect of I-CARE over U-CARE on improving SPSI and SCL-20 including six, 12-, and 24-months data using an LMM. This LMM included SPSI or SCL-20 change relative to baseline as the dependent variable, intervention (I-CARE or U-CARE), time (six, 12, or 24 months), and their interaction as the fixed effects and subjects as the random effect. Baseline SPSI or SCL-20 was included as a regressor.


It also include delineating individual's behavioral outcome trajectories shown in Fig. S2. To delineate how an individual’s behavioral outcomes (SCL-20 and SPSI) changed across intervention phases and to decide whether the brain-behavior association should be modeled differently across intervention phases, we plotted individuals’ trajectories and assessed the stability across different time points via correlation coefficients and intra-class coefficient (ICC(2, 1)).

#### whole-brain LMMs in identifying mechanism and predictive markers within the cognitive control circuit

The code for generating whole-brain results in the paper is here: [STEP1_ENGAGE_CCC_lme_wholebrain_voxelwise.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1b19afde1c0c010ca0b44e19d64cd492a13905d9/matlab/STEP1_ENGAGE_CCC_lme_wholebrain_voxelwise.m), this generates the brain results in Fig. 2A, 3A, 4A, 4C as well as Fig. S5A, S6-S8.

We utilized the fitlme function of the Statistics and Machine Learning Toolbox in [Matlab](https://www.mathworks.com/products/matlab.html) and customized scripts to conduct the voxel-wise whole-brain LMM with clinical measures (here ∆SCL-20 or ∆SPSI) as the dependent variable, as detailed below. The significance of fixed effects was tested under the type III hypotheses using the analysis of variance (ANOVA) function in Matlab. The degrees of freedom for ANOVA were estimated via the Satterthwaite method. All results were corrected multiple comparisons with a voxel threshold of p < 0.001 and a Gaussian random field theory (GRF) familywise error cluster-level correction at p < 0.05 using [DPABI V5.1](http://rfmri.org/dpabi). Data smoothness for GRF correction was estimated on the statistical image following a similar procedure to FSL easythresh.

The script for generating scatter plots for each regions is here: [STEP2_ENGAGE_CCC_lme_ROI_visualization.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1b19afde1c0c010ca0b44e19d64cd492a13905d9/matlab/STEP2_ENGAGE_CCC_lme_ROI_visualization.m). This script also generated the comparison of models using cognitive control circuit as the independent variables and alternative control models using only baseline characteristics in Table S1.

The script for estimating beta coefficient per each intervention group or timepoint is [ENGAGE_LME_Task_Activation_SCL20.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1b19afde1c0c010ca0b44e19d64cd492a13905d9/r_script/Rmd/ENGAGE_LME_Task_Activation_SCL20.Rmd) for SCL-20 depression severity and [ENGAGE_LME_Task_Activation_SPSI.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1b19afde1c0c010ca0b44e19d64cd492a13905d9/r_script/Rmd/ENGAGE_LME_Task_Activation_SPSI.Rmd) for SPSI. This script produced Fig. S3 and Table S2.

The script that runs cross-validation for the continuous and binary prediction of 6,12, and 24 months SCL-20 using 2 months cognitive control circuit activation is [ENGAGE_cross_validation_SCL20.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1b19afde1c0c010ca0b44e19d64cd492a13905d9/r_script/Rmd/ENGAGE_cross_validation_SCL20.Rmd)

Details about each model can be found in the paper.


# ENGAGE Cognitive Control Circuit paper

This repo contains major steps and codes to replicate the analysis included in the manuscript entitled "Adaptive Cognitive Control Circuit Changes Associated with Problem-solving Ability and Depression Symptom Outcomes over 24 months". In this project, we focused on the cognitive control circuit as a putative neural mechanism of action for a novel behavioral intervention with five repeat measures over two years. We explored the possibility of using early changes in this circuit to predict future treatment outcomes. The whole-brain voxel-wise Linear Mixed Model analysis was conducted in Matlab using customized code. 

Scripts are listed based on their order in the Results section.


## High-level steps
- [Pre-requisite](#pre-requisite)
- [Installation guide](#Installation-guide)
- [How to run scripts](#How-to-run-scripts)
- [Documentation for methods and corresponding scripts](#Documentation-for-methods-and-corresponding-scripts)


## Pre-requisite
### Hardware requirements
All steps could be done on a standard research computer with reasonable CPUs and RAM.

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
  
## How to run scripts
  1. When running Matlab scripts, `cd ${where_you_would_like_to_save_the_code}/ENGAGE-CognitiveControl/matlab` and run scripts within this folder.
  2. When running R scripts, open the R project file `ENGAGE-CognitiveControl/r_script/ENGAGE_LME_Task_Activation.Rproj` and run RMD under the `Rmd` subfolder.
  

## Documentation for methods and corresponding scripts  

#### Schematic diagram Figure 1
  
  [Figure1_Simulation.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/383702627a7a52139e8558ab67e68e19692c61f6/matlab/Figure1_Simulation.m) 

- This Matlab script simulted and visualized data demonstrating the modulating effect of intervention or time or their interaction in Fig. 1.

#### Baseline demographic and characteristics Table 1

  [Generate_Baseline_Demographics.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/13928c1f5a0de9e3a4bebac7b93c48de1789545d/r_script/Rmd/Generate_Baseline_Demographics.Rmd)

- This R script summarized baseline demographics and characteristics for each intervention group and the overall group in Table 1.

#### Examining effects of I-CARE on behavioral outcomes over U-CARE and delineating individuals’ behavioral outcome trajectories

  [ENGAGE_LME_behavioral_outcome.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/13928c1f5a0de9e3a4bebac7b93c48de1789545d/r_script/Rmd/ENGAGE_LME_behavioral_outcome.Rmd)

  - This R script examined the effect of I-CARE over U-CARE on improving SPSI and SCL-20 including 6-, 12-, and 24-month data using an LMM. 

  - It also included delineating an individual's behavioral outcome trajectories shown in Fig. S2.

#### whole-brain LMMs in identifying the mechanism and predictive markers within the cognitive control circuit

  [STEP1_ENGAGE_CCC_lme_wholebrain_voxelwise_v4.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1dba4b68dd7b179ed775846183ba001b1fe29ead/matlab/STEP1_ENGAGE_CCC_lme_wholebrain_voxelwise_v4.m)

  - This Matlab script generated whole-brain results in Fig. 2A, 3A, 4A, 4C as well as Fig. S4A, S9-S11.
  
  - We utilized the fitlme function of the Statistics and Machine Learning Toolbox in [Matlab](https://www.mathworks.com/products/matlab.html) and customized scripts to conduct the voxel-wise whole-brain LMM with clinical measures (here ∆SCL-20 or ∆SPSI) as the dependent variable, as detailed below. 
  
  - The significance of fixed effects was tested under the type III hypotheses using the analysis of variance (ANOVA) function in Matlab. 
  
  - The degrees of freedom for ANOVA were estimated via the Satterthwaite method. All results were corrected multiple comparisons with a voxel threshold of p < 0.001 and a Gaussian random field theory (GRF) familywise error cluster-level correction at p < 0.05 using [DPABI V5.1](http://rfmri.org/dpabi).
  - This script could be modified to examine the association between baseline circuit activation and outcomes.

#### ROI LMMs and visualizations

  [STEP2_ENGAGE_CCC_lme_ROI_visualization_mechanistic_v3.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1dba4b68dd7b179ed775846183ba001b1fe29ead/matlab/STEP2_ENGAGE_CCC_lme_ROI_visualization_mechanistic_v3.m) 
  [STEP2_ENGAGE_CCC_lme_ROI_visualization_predictive_v3.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1dba4b68dd7b179ed775846183ba001b1fe29ead/matlab/STEP2_ENGAGE_CCC_lme_ROI_visualization_predictive_v3.m) 

  - These two Matlab scripts generated statistics in Table 2 and scatter plots in Fig. 2B, 3B, 4B, 4D, and Fig. S3 and S4B for significant regions that survived multiple comparison corrections in the above whole-brain analysis, for mechanism markers and predictive markers respectively.
  
  - These scripts also generated the ANOVA comparison of models using cognitive control circuit as the independent variables and alternative control models using only baseline characteristics in Table S1.

#### Post-hoc beta estimation for significant LMMs

  Below scripts produced Table S2.

  - [ENGAGE_LME_Task_Activation_SCL20.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/27a7f97ab58e6f34e76927dd28d81122380142b7/r_script/Rmd/ENGAGE_LME_Task_Activation_SCL20.Rmd) estimated beta coefficient per each intervention group or timepoint for SCL-20 depression severity.
  
  - [ENGAGE_LME_Task_Activation_SPSI.Rmd](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/27a7f97ab58e6f34e76927dd28d81122380142b7/r_script/Rmd/ENGAGE_LME_Task_Activation_SPSI.Rmd) estimated beta coefficient per each intervention group or timepoint for SPSI problem-solving ability.

 Inside these scripts, we tested the confounding effects of go-nogo performance and medication prescription; Additionally, we tested the association between behavioral performance and treatment outcomes.

#### Testing generalizability of two-month neuroimaging predictors using cross-validation

  [STEP3_ENGAGE_CCC_lme_wholebrain_voxelwise_CV_v3.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1dba4b68dd7b179ed775846183ba001b1fe29ead/matlab/STEP3_ENGAGE_CCC_lme_wholebrain_voxelwise_CV_v3.m)

  - The Matlab script ran cross-validation for the continuous prediction of 6,12, and 24 months SCL-20 using 2 months cognitive control circuit activation change and reported its performance against the control model using only baseline characteristics in Fig. 5A, Table 3, and Table S3. It also generated the average feature map in Fig. S6.

  [STEP3_ENGAGE_CCC_lme_wholebrain_voxelwise_CV_binary_v5.m](https://github.com/WilliamsPanLab/ENGAGE-CognitiveControl/blob/1dba4b68dd7b179ed775846183ba001b1fe29ead/matlab/STEP3_ENGAGE_CCC_lme_wholebrain_voxelwise_CV_binary_v5.m)

  - The Matlab script ran cross-validation for the prediction of the 6-month binary response using 2 months cognitive control circuit activation change and reported its performance against the control model using only baseline characteristics in Fig. 5B. It also generated the average feature map in Fig. S6.



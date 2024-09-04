%% Prepare and Visualize Scatterplots for Significant Results
clear; clc; close all

% Load plot parameters
load_plot_parameters

dataDir = fullfile('/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Code/github/matlab');

%% Load non-imaging data and define outcome
load(fullfile(dataDir, 'non_imaging_forgithub.mat'));
outcomeIndex = 1; % Outcome index - 1: SCL-20; 2: SPSI; 12: 'feeling everything an effort'; 7: Go-NoGo reaction time; 8: Go-NoGo false alarms
outcomeVariable = beh_name{outcomeIndex};

%% Define subjects and timepoints to include in the current analysis
numSubjects = 108;
sessionList = {'000', '2MO', '6MO', '12MO', '24MO'};
actualTime = [0, 2, 6, 12, 24];
includedSessions = 1:5; 
includedSessions4Change = setdiff(includedSessions, 1); 

% Construct model timepoints string
modelTimePoints = strjoin(arrayfun(@num2str, actualTime(includedSessions), 'UniformOutput', false), '_');
modelTimePoints = [modelTimePoints, '_MO'];

%% Load imaging data
% Define task and contrast
taskName = '111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO'; % GoNoGo task
contrastName = 'NogovsGo'; % NoGo vs Go

% Load imaging data for ROIs that survived multiple corrections
roiLabel = 'MixedModel_SCL20_mechanism';
load(fullfile(dataDir, taskName, contrastName, ['task_activation_', roiLabel, '.mat']));

%% Define models for hypothesis 1 & 2 (mechanistic markers)
modelNames = {'mechanistic'};
modelNum = length(modelNames);
modelEquations = {
    sprintf('%s ~ %s + %s:%s + %s:%s + %s:%s:%s + %s + %s + %s + %s + (%s|%s)', ...
    outcomeVariable, 'fMRI', 'fMRI', 'Time', 'fMRI', 'Group', 'fMRI', 'Time', 'Group', [beh_name{outcomeIndex}, '_baseline'], 'BMI_baseline', 'age', 'sex', '1', 'Subj')
};
numModelVars = 9;

modelControl = sprintf('%s ~ %s + %s + %s + %s + (%s|%s)', ...
    outcomeVariable, [outcomeVariable, '_baseline'], 'BMI_baseline', 'age', 'sex', '1', 'Subj');

%% Convert data into long format for modeling
timePoints = repmat(includedSessions4Change, numSubjects, 1);
timePoints = timePoints(:);

subjectIds = repmat((1:numSubjects)', 1, length(includedSessions4Change));
subjectIds = subjectIds(:);
subjectIds = nominal(subjectIds);

ageSex = beh_all(:, 4:5, 1);
ageSexLong = repmat(ageSex, length(includedSessions4Change), 1);

% Load intervention and response information
load(fullfile(dataDir, 'Interv_Info.mat')); 
load(fullfile(dataDir, 'Responder_list_6mo.mat')); 
load(fullfile(dataDir, 'Remission_list_6mo.mat'));
load(fullfile(dataDir, 'Responder_list_long.mat')); 
load(fullfile(dataDir, 'Remission_list_long.mat'));

groupInfo = repmat(Interv_Info, 1, length(includedSessions4Change)); 
groupInfo = nominal(groupInfo(:));

responseInfo = 2 - resp_ident_long(:, includedSessions4Change); 
responseInfo = responseInfo(:);

remissionInfo = 2 - remi_ident_long(:, includedSessions4Change); 
remissionInfo = remissionInfo(:);

outcomeChange = squeeze(beh_all(:, outcomeIndex, includedSessions4Change)) - repmat(beh_all(:, outcomeIndex, 1), 1, length(includedSessions4Change));
outcomeChange = outcomeChange(:);

baselineOutcome = repmat(squeeze(beh_all(:, outcomeIndex, 1)), 1, length(includedSessions4Change)); 
baselineOutcome = baselineOutcome(:);

baselineBMI = repmat(squeeze(beh_all(:, 6, 1)), 1, length(includedSessions4Change)); 
baselineBMI = baselineBMI(:);

FDIncludeChange = beh_all(:, 11, includedSessions4Change) - beh_all(:, 11, 1); 
FDIncludeChange = FDIncludeChange(:);

%% Generate statistical results for all linear mixed-effects models (LMMs)
outputDir = fullfile(dataDir, 'SPM', 'Level2', 'MixedModel', taskName, contrastName, outcomeVariable, roiLabel);
mkdir(outputDir);

for iRoi = 1:length(roiList)
    roiName = roiList{iRoi};
    
    % Prepare fMRI data
    fmriChange = squeeze(taskActivation(:, iRoi, includedSessions4Change)) - repmat(squeeze(taskActivation(:, iRoi, 1)), 1, length(includedSessions4Change));
    fmriChange = fmriChange(:);
    
    fmriBaseline = repmat(squeeze(taskActivation(:, iRoi, 1)), 1, length(includedSessions4Change));
    fmriBaseline = fmriBaseline(:);
    
    fmri2Mo = repmat(squeeze(taskActivation(:, iRoi, 2)) - squeeze(taskActivation(:, iRoi, 1)), 1, length(includedSessions4Change));
    fmri2Mo = fmri2Mo(:);
    
    % Create table for LME model
    tblData = table(outcomeChange, baselineOutcome, timePoints, subjectIds, FDIncludeChange, fmriBaseline, fmri2Mo, fmriChange, ...
                    ageSexLong(:, 1), ageSexLong(:, 2), groupInfo, baselineBMI, ...
                    'VariableNames', {outcomeVariable, [outcomeVariable, '_baseline'], 'Time', 'Subj', 'FD', ...
                                      'fMRI_baseline', 'fMRI_2MO', 'fMRI', 'age', 'sex', 'Group', 'BMI_baseline'});
    tblData.sex = categorical(tblData.sex); 
    tblData.Group = nominal(tblData.Group);
    tblData = rmmissing(tblData, 'DataVariables', {'fMRI', outcomeVariable, 'Time', 'Group', ...
        [outcomeVariable, '_baseline'], 'age', 'sex', 'BMI_baseline'});
    
    % Fit LME model
    lmeModel = fitlme(tblData, modelEquations{1}, ...
        'DummyVarCoding', 'effects', 'StartMethod', 'random');
    
    % Perform ANOVA
    anovaTable = anova(lmeModel, 'DFMethod', 'satterthwaite');
    anovaTable = dataset2table(anovaTable);
    
    % Fit control model and compare
    controlModel = fitlme(tblData, modelControl, ...
        'DummyVarCoding', 'effects', 'StartMethod', 'random');
    comparisonTable = dataset2table(compare(controlModel, lmeModel));
    writetable(comparisonTable, fullfile(outputDir, ['likelihood_', roiName, '.xls']), 'Sheet', roiName);
end

%% Creating Scatter Plots for Significant ROIs
% Define group info for plotting
legends = {'U-CARE', 'I-CARE'};
groupInfoLong = repmat(Interv_Info, 1, length(includedSessions4Change)); % for dlpfc replace Interv_Info with timePoints
groupInfoLong = groupInfoLong(:);
groupNames = cell(length(groupInfoLong), 1);
groupNames(groupInfoLong == 0) = {legends{1}};
groupNames(groupInfoLong == 1) = {legends{2}};

% Loop through ROIs to generate scatter plots
for iRoi = 1:length(roiList)
    roiName = roiList{iRoi};
    
    % Prepare fMRI change data
    fmriChange = squeeze(taskActivation(:, iRoi, includedSessions4Change)) - repmat(squeeze(taskActivation(:, iRoi, 1)), 1, length(includedSessions4Change));
    fmriChange = fmriChange(:);

    % Overall correlation
    [r, p] = partialcorr(outcomeChange, fmriChange, ageSexLong, 'Rows', 'pairwise', 'Type', 'Spearman');
    fitModel = fitlm(fmriChange, outcomeChange);
    
    % Group-specific correlations
    [rU, pU] = partialcorr(outcomeChange(groupInfoLong == 0), fmriChange(groupInfoLong == 0), ageSexLong(groupInfoLong == 0), 'Rows', 'pairwise');
    [rI, pI] = partialcorr(outcomeChange(groupInfoLong == 1), fmriChange(groupInfoLong == 1), ageSexLong(groupInfoLong == 1), 'Rows', 'pairwise');
    
    % Plot scatter plot
    xLabel = ['\Delta Activation in ', roiName];
    yLabel = ['\Delta ', beh_name{outcomeIndex}];
    figPosition = [2, 2, 9.5, 6.35];
    scatterPlot = fmri_behavior_plot4(fmriChange, outcomeChange, 1 - groupInfoLong, [183, 90, 101; 166, 166, 172]/255, 1, 0.8, xLabel, yLabel, '', flip(legends), figPosition, '', 0);
    
    if outcomeIndex == 1
        set(gca, 'YDir', 'reverse');
    end
    
    legend off;
    savePath = fullfile(outputDir, [roiLabel,'_', roiName, '_', outcomeVariable, '.tiff']);
    export_fig(savePath, '-r300', '-transparent');
end

clear, clc, close all
parpool
dataDir = fullfile('/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Code/github/matlab');

%% Load non-imaging data and define outcome
load(fullfile(dataDir, 'non_imaging_forgithub.mat'));
outcomeIndex = 1; % Outcome index - 1: SCL-20; 2: SPSI
outcomeVariable = beh_name{outcomeIndex};

%% Define subjects and timepoints to include in the current analysis
numSubjects = 108;
sessionList = {'000'; '2MO'; '6MO'; '12MO'; '24MO'};
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

% Define mask used for extracting imaging data
maskFile = fullfile(dataDir, 'GreyMask_02_91x109x91.img');
[mask3D, ~, ~, header] = y_ReadAll(maskFile); 
maskSize = size(mask3D); 
threshold = 0; 
maskIndices = find(mask3D > threshold);
header.pinfo = [1; 0; 0]; 
header.dt = [16, 0]; 
maskSuffix = ['GM_', num2str(threshold)];

% Make mask for GRF correction
maskGRF = zeros(maskSize); 
maskGRF(maskIndices) = 1; 
[a, b, c] = fileparts(maskFile); 
maskGRFOutName = fullfile(a, [b, '_thresh_', num2str(threshold), c]); 
y_Write(maskGRF, header, maskGRFOutName);

% Load imaging data
load(fullfile(dataDir, [taskName, '_', contrastName, '_task_activation_wholebrain_', maskSuffix, '.mat']));

% Set up GRF correction
voxelPThreshold = 0.001;
isTwoTailed = 0;
clusterPThreshold = 0.05;

% Plot parameters
load_plot_parameters
minNegative = -3; 
minPositive = 3;
maxNegative = -4; 
maxPositive = 4;
colorMap = y_AFNI_ColorMap(12);
flag = 'F';
clusterSize = 0;
connectivityCriterion = 18;
sliceInterval = 6;

% Set up cross-validation
rng(1); % For reproducibility
folds = 5;
cvIndices = crossvalind('Kfold', numSubjects, folds);

% Initialize variables to store final results
ypredAll = []; 
yactualAll = []; 
aicAll = []; 
bicAll = []; 
maeAll = []; 
rAll = []; 
pAll = []; 

% Start cross-validation
for fold = 1:folds
    fprintf('Starting fold %d/%d...\n', fold, folds);
    
    % Split data into training and testing sets
    trainIndices = (cvIndices ~= fold);
    testIndices = ~trainIndices;
    
    % Training data
    fmriTrain = fmri_wb(trainIndices, :, :);
    outcomeChangeTrain = outcomeChange(trainIndices);
    baselineOutcomeTrain = baselineOutcome(trainIndices);
    fdIncludeChangeTrain = FDIncludeChange(trainIndices);
    fdIncludeBaselineTrain = FDIncludeBaseline(trainIndices);
    subjectIdsTrain = subjectIds(trainIndices);
    ageSexLongTrain = ageSexLong(trainIndices, :);
    groupInfoTrain = groupInfo(trainIndices);
    
    % Test data
    fmriTest = fmri_wb(testIndices, :, :);
    outcomeChangeTest = outcomeChange(testIndices);
    baselineOutcomeTest = baselineOutcome(testIndices);
    fdIncludeChangeTest = FDIncludeChange(testIndices);
    fdIncludeBaselineTest = FDIncludeBaseline(testIndices);
    subjectIdsTest = subjectIds(testIndices);
    ageSexLongTest = ageSexLong(testIndices, :);
    groupInfoTest = groupInfo(testIndices);

    % Output directory for the current fold
    outputDir = fullfile(dataDir, 'MixedModel', taskName, contrastName, outcomeVariable, [maskSuffix, '_', modelTimePoints, '_Fold_', num2str(fold)]);
    mkdir(outputDir);

    % Train the model and save F-maps
    modelFStats = nan(length(maskIndices), 9); % 9 variables in the model
    modelDF = nan(length(maskIndices), 2, 9); % 2 degrees of freedom (DF1 and DF2) for each variable

    % Parallel loop for each voxel
    parfor voxelIndex = 1:length(maskIndices)
        % Prepare fMRI data in long format
        fmriDataChange = squeeze(fmriTrain(:, voxelIndex, includedSessions4Change)) - repmat(squeeze(fmriTrain(:, voxelIndex, 1)), 1, length(includedSessions4Change));
        fmriDataChange = fmriDataChange(:);

        fmriBaseline = squeeze(fmriTrain(:, voxelIndex, 1));
        fmriBaseline = repmat(fmriBaseline, 1, length(includedSessions4Change));
        fmriBaseline = fmriBaseline(:);

        fmri2Mo = squeeze(fmriTrain(:, voxelIndex, 2));
        fmri2Mo = repmat(fmri2Mo, 1, length(includedSessions4Change));
        fmri2Mo = fmri2Mo(:) - fmriBaseline;

        % Create table for the model
        tblData = table(outcomeChangeTrain, baselineOutcomeTrain, timePoints(trainIndices), subjectIdsTrain, ...
            fmriBaseline, fmri2Mo, fmriDataChange, ageSexLongTrain(:, 1), ageSexLongTrain(:, 2), groupInfoTrain, ...
            'VariableNames', {outcomeVariable, [outcomeVariable, '_baseline'], 'Time', 'Subj', 'fMRI_baseline', ...
                              'fMRI_2MO', 'fMRI', 'age', 'sex', 'Group'});

        % Fit the predictive model
        mixModel = fitlme(tblData, modelEquations{2}, 'DummyVarCoding', 'effects', 'StartMethod', 'random');
        
        % Extract F-statistics and degrees of freedom
        anovaResults = anova(mixModel);
        modelFStats(voxelIndex, :) = [anovaResults.FStat; nan(9 - length(anovaResults.FStat), 1)]';
        modelDF(voxelIndex, :, :) = [[anovaResults.DF1, anovaResults.DF2]; nan(9 - length(anovaResults.DF1), 2)];
    end

    % Save F-maps for multiple comparison correction
    for varIndex = 1:9
        data = modelFStats(:, varIndex);
        data3D = zeros(maskSize);
        data3D(maskIndices) = data;

        outputFileName = fullfile(outputDir, sprintf('F_model_predict_var%d.nii', varIndex));
        y_Write(data3D, header, outputFileName);

        % Apply GRF correction
        [~, filename, ~] = fileparts(outputFileName);
        outputNameCorr = fullfile(outputDir, [filename, '_GRF_', num2str(voxelPThreshold), '.nii']);
        df1 = max(modelDF(~isnan(modelDF(:, 1, varIndex)), 1, varIndex));
        df2 = max(modelDF(~isnan(modelDF(:, 2, varIndex)), 2, varIndex));
        y_GRF_Threshold(outputFileName, voxelPThreshold, isTwoTailed, clusterPThreshold, outputNameCorr, maskGRFOutName, flag, df1, df2);
    end

    % Extract significant peaks
    significantIndices = [];
    for varIndex = 1:9
        corrFileName = fullfile(outputDir, sprintf('F_model_predict_var%d_GRF_%0.3f.nii', varIndex, voxelPThreshold));
        if exist(corrFileName, 'file')
            peakIndices = extract_peak_ind(corrFileName);
            significantIndices = [significantIndices; peakIndices];
        end
    end
    significantIndices = unique(significantIndices);

    %% Final model with significant voxels
    tblDataTest = table(outcomeChangeTest, baselineOutcomeTest, timePoints(testIndices), subjectIdsTest, ...
        ageSexLongTest(:, 1), ageSexLongTest(:, 2), groupInfoTest, ...
        'VariableNames', {outcomeVariable, [outcomeVariable, '_baseline'], 'Time', 'Subj', 'age', 'sex', 'Group'});

    for voxelIndex = significantIndices'
        % Prepare data from significant voxels for testing
        fmri2MoTest = squeeze(fmriTest(:, voxelIndex, 2));
        fmri2MoTest = fmri2MoTest - fmriBaselineTest;
        fmriBaselineTest = squeeze(fmriTest(:, voxelIndex, 1));
       
        % Add significant voxel data to table
        tblDataTest.(sprintf('fMRI_%d', voxelIndex)) = fmri2MoTest;
        tblDataTest.(sprintf('fMRI_%d_baseline', voxelIndex)) = fmriBaselineTest;
    end

    % Fit the model using significant voxels and predict on testing data
    modelEqTest = sprintf('%s ~ %s + %s + %s + %s + %s', ...
        outcomeVariable, strjoin(tblDataTest.Properties.VariableNames(contains(tblDataTest.Properties.VariableNames, 'fMRI')), ' + '), ...
        [outcomeVariable, '_baseline'], 'age', 'sex', 'Group');

    finalModel = fitlme(tblDataTest, modelEqTest, 'DummyVarCoding', 'effects', 'StartMethod', 'random');
    yPred = predict(finalModel, tblDataTest);

    % Store results
    ypredAll = [ypredAll; yPred]; 
    yactualAll = [yactualAll; outcomeChangeTest];
    aicAll = [aicAll; finalModel.ModelCriterion.AIC]; 
    bicAll = [bicAll; finalModel.ModelCriterion.BIC];
end

% Save final results after all folds
save(fullfile(outputDir, 'final_results.mat'), 'ypredAll', 'yactualAll', 'aicAll', 'bicAll');

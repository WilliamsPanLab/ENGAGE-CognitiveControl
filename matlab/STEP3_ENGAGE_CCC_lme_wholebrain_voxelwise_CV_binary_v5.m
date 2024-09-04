%% Clear Workspace and Set Paths
clear; clc; close all;
% parpool; % Uncomment if parallel processing is needed
% pctRunOnAll javaaddpath '/Users/xuezhang/Downloads/Software/ParforProgMon'

% Define data directory
dataDir = fullfile('/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Code/github/matlab');

%% Load Non-Imaging Data and Define Outcome
load(fullfile(dataDir, 'non_imaging_forgithub.mat'));

% Define outcome variable index and name
outcomeIndex = 1; % 1: SCL-20; 2: SPSI
outcomeVariable = beh_name{outcomeIndex};

%% Define Subjects and Timepoints for Analysis
numSubjects = 108;
sessionList = {'000', '2MO', '6MO', '12MO', '24MO'};
actualTime = [0, 2, 6, 12, 24];

% Specify included sessions to predict response at 6 months
includedSessions = [1, 3]; % Baseline and 6 months
includedSessions4Change = setdiff(includedSessions, 1);

% Construct model timepoints string
modelTimePoints = strjoin(arrayfun(@(x) sprintf('%d', x), actualTime(includedSessions), 'UniformOutput', false), '_');
modelTimePoints = [modelTimePoints, '_MO'];

%% Load Imaging Data
% Define task and contrast names
taskName = '111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO'; % Go/No-Go task
contrastName = 'NogovsGo'; % Contrast between No-Go and Go conditions

% Define mask file for extracting imaging data
maskFile = fullfile(dataDir, 'GreyMask_02_91x109x91.img');
[mask3D, ~, ~, header] = y_ReadAll(maskFile);
maskSize = size(mask3D);
threshold = 0;
maskIndices = find(mask3D > threshold);

% Update header information
header.pinfo = [1; 0; 0];
header.dt = [16, 0];
maskSuffix = ['GM_', num2str(threshold)];

% Create mask for GRF correction
maskGRF = zeros(maskSize);
maskGRF(maskIndices) = 1;
[maskPath, maskName, maskExt] = fileparts(maskFile);
maskGRFOutName = fullfile(maskPath, [maskName, '_thresh_', num2str(threshold), maskExt]);
y_Write(maskGRF, header, maskGRFOutName);

% Load whole-brain fMRI data
load(fullfile(dataDir, [taskName, '_', contrastName, '_task_activation_wholebrain_', maskSuffix, '.mat']));

% Prepare baseline and 6-month change fMRI data
fmriBaseline = squeeze(fmri_wb(:, :, 1)); % Baseline data
fmri2MO = squeeze(fmri_wb(:, :, 2)) - fmriBaseline; % 2-month change data from baseline
fmri2MO = fmri2MO';

%% Define Logistic Regression Model for Predictive Markers
modelNames = {'predictive'};
numModels = length(modelNames);

% Define logistic regression model equation
modelEquation = 'Response ~ fMRI_2MO + SCL_20_baseline + age + sex + BMI_baseline';
numModelVars = 6;

% Control model equation (without fMRI data)
modelEquationControl = 'Response ~ SCL_20_baseline + age + sex + BMI_baseline';

%% Set Up GRF Correction Parameters
voxelPThreshold = 0.01;
isTwoTailed = 1;
clusterPThreshold = 0.05;
flag = 'T';
%% Set Up Cross-Validation Parameters
numFolds = 2;
rng(10001); % Set random seed for reproducibility
allIndices = randperm(numSubjects);
trainPercent = 0.8; % Percentage of data for training

%% Prepare Short Format Data
% Prepare subject IDs and age/sex data
subjectIDs = (1:numSubjects)';
ageSex = beh_all(:, 4:5, 1);
baselineOutcome = beh_all(:, outcomeIndex, 1);
baselineBMI = beh_all(:, 6, 1);

% Load intervention and response information
load(fullfile(dataDir, 'Responder_list_6mo.mat')); 

% Prepare group and response information
responseInfo = 2 - resp_ident; % Binary response variable for 6-month outcome

%% Initialize Result Storage Variables
ypredAll = [];
yactualAll = [];
prAucAll = [];
rocAucAll = [];
sensitivityAll = [];
specificityAll = [];

ypredAllControl = [];
prAucAllControl = [];
rocAucAllControl = [];
sensitivityAllControl = [];
specificityAllControl = [];

%% Cross-Validation Loop
for fold = 1:numFolds
    fprintf('Processing fold %d/%d\n', fold, numFolds);
    
    % Define testing and training indices
    testIndices = floor((1 + (fold - 1) * (1 - trainPercent) * numSubjects) : fold * (1 - trainPercent) * numSubjects);
    testSubjects = allIndices(testIndices);
    trainSubjects = setdiff(allIndices, testSubjects);
    
    %% Prepare Training Data
    responseTrain = responseInfo(trainSubjects);
    baselineOutcomeTrain = baselineOutcome(trainSubjects);
    fmri2MOTrain = fmri2MO(:, trainSubjects);
    ageTrain = ageSex(trainSubjects, 1);
    sexTrain = ageSex(trainSubjects, 2);
    bmiBaselineTrain = baselineBMI(trainSubjects);
    
    %% Prepare Testing Data
    responseTest = responseInfo(testSubjects);
    baselineOutcomeTest = baselineOutcome(testSubjects);
    fmri2MOTest = fmri2MO(:, testSubjects);
    ageTest = ageSex(testSubjects, 1);
    sexTest = ageSex(testSubjects, 2);
    bmiBaselineTest = baselineBMI(testSubjects);
    
    %% Initialize Model Result Storage for Current Fold
    modelTStats = nan(length(maskIndices), numModelVars);
    modelDF = nan(length(maskIndices), numModelVars);
    maskUsed = zeros(length(maskIndices), 1);
%     ppm = ParforProgMon('Progress bar', length(maskIndices));
    
    %% Loop Through All Voxels for Logistic Regression
%     for voxelIdx = 1:length(maskIndices)
%         % Extract fMRI data for current voxel
%         fmri2MOVoxelTrain = fmri2MOTrain(voxelIdx, :);
%         
%         % Prepare data table for modeling
%         tblTrain = table(...
%             responseTrain, ...
%             baselineOutcomeTrain, ...
%             fmri2MOVoxelTrain(:), ...
%             ageTrain, ...
%             sexTrain, ...
%             bmiBaselineTrain, ...
%             'VariableNames', {...
%                 'Response', ...
%                 'SCL_20_baseline', ...
%                 'fMRI_2MO', ...
%                 'age', ...
%                 'sex', ...
%                 'BMI_baseline'});
%         
%         % Fit logistic regression model
%         try
%             logisticModel = fitglm(tblTrain, modelEquation, 'Distribution', 'binomial');
%             
%             % Store t-statistics and degrees of freedom
%             coeffTable = logisticModel.Coefficients;
%             modelTStats(voxelIdx, :) = coeffTable.tStat;
%             modelDF(voxelIdx,:) = repmat(logisticModel.DFE,1, numModelVars);
%             maskUsed(voxelIdx) = 1;
%         catch
%             % Skip voxel if model fitting fails
%             continue;
%         end
% %         ppm.increment();
%     end
%     
%     for varIdx = 1:numModelVars
%         % Prepare data for saving
%         tValueData = modelTStats(:, varIdx);
%         tValue3D = zeros(maskSize);
%         tValue3D(maskIndices) = tValueData;
%         
%         % Define output filenames
%         outFileName = fullfile(dataDir, sprintf('T_model_%s_var%d_fold%d.nii', modelNames{1}, varIdx, fold));
%         y_Write(tValue3D, header, outFileName);
%         
%         % Perform GRF correction
%         df1 = max(modelDF(~isnan(modelDF(:, varIdx)), varIdx));
%         outputNameCorr = fullfile(dataDir, sprintf('T_model_%s_var%d_fold%d_GRF_%.2f.nii', modelNames{1}, varIdx, fold, voxelPThreshold));
%         y_GRF_Threshold(outFileName, voxelPThreshold, isTwoTailed, clusterPThreshold, outputNameCorr, maskGRFOutName, flag, df1);
%     end
    
    %% Save Used Mask
    mask3D = zeros(maskSize);
    mask3D(maskIndices) = maskUsed;
    maskOutName = fullfile(dataDir, 'mask_used.nii');
    y_Write(mask3D, header, maskOutName);
    
    %% Identify Significant Voxels and Finalize Model
    % Extract indices for significant clusters' peaks
    varIdx = 3;
    outname = fullfile(dataDir, sprintf('Z_ClusterThresholded_T_model_%s_var%d_fold%d_GRF_%.2f.nii', modelNames{1}, varIdx, fold, voxelPThreshold));
    try
        [ind] = extract_peak_ind(outname);
        [~, ind_inmask] = ismember(ind, maskIndices);
        peakSigCluster = setdiff(unique(ind_inmask), 0);
        numPeakSigCluster = length(peakSigCluster);
        roiNames = arrayfun(@(x) sprintf('ROI%d', x), 1:numPeakSigCluster, 'UniformOutput', false);
  
        tblTrain = table(...
            responseTrain, ...
            baselineOutcomeTrain, ...
            ageTrain, ...
            sexTrain, ...
            bmiBaselineTrain, ...
            'VariableNames', {...
            'Response', ...
            'SCL_20_baseline', ...
            'age', ...
            'sex', ...
            'BMI_baseline'});
        
        % Construct final logistic regression model equation with significant ROIs
        roiTerms = strjoin(roiNames, ' + ');
        finalModelEquation = sprintf('Response ~ %s + SCL_20_baseline + age + sex + BMI_baseline', roiTerms);
        
        % Prepare training data table with significant voxels
        for roiIdx = 1:numPeakSigCluster
            voxelIdx = peakSigCluster(roiIdx);
            fmri2MOVoxelTrain = fmri2MOTrain(voxelIdx, :);
            tblTrain.(roiNames{roiIdx}) = fmri2MOVoxelTrain(:);
        end
        
        % Prepare testing data table with significant voxels
        tblTest = table(...
            responseTest, ...
            baselineOutcomeTest, ...
            ageTest, ...
            sexTest, ...
            bmiBaselineTest, ...
            'VariableNames', {...
            'Response', ...
            'SCL_20_baseline', ...
            'age', ...
            'sex', ...
            'BMI_baseline'});
        
        for roiIdx = 1:numPeakSigCluster
            voxelIdx = peakSigCluster(roiIdx);
            fmri2MOVoxelTest = fmri2MOTest(voxelIdx, :);
            tblTest.(roiNames{roiIdx}) = fmri2MOVoxelTest(:);
        end
        
        % Fit final logistic regression model and predict on testing data
        try
            finalModel = fitglm(tblTrain, finalModelEquation, 'Distribution', 'binomial');
            yPredProb = predict(finalModel, tblTest);
            
            % Determine optimal threshold to maximize F1 score
            [optimalThreshold, prAuc, rocAuc, sensitivity, specificity] = calculateOptimalThreshold(responseTest, yPredProb);
            
            % Convert probabilities to binary predictions using optimal threshold
            yPredBinary = yPredProb >= optimalThreshold;
            
            % Store prediction results
            ypredAll = [ypredAll; yPredProb];
            yactualAll = [yactualAll; responseTest];
            prAucAll = [prAucAll; prAuc];
            rocAucAll = [rocAucAll; rocAuc];
            sensitivityAll = [sensitivityAll; sensitivity];
            specificityAll = [specificityAll; specificity];
        catch
            % Handle cases where model fitting fails
            ypredAll = [ypredAll; nan(size(responseTest))];
            yactualAll = [yactualAll; responseTest];
            prAucAll = [prAucAll; nan];
            rocAucAll = [rocAucAll; nan];
            sensitivityAll = [sensitivityAll; nan];
            specificityAll = [specificityAll; nan];
        end
        
        %% Fit and Predict with Control Model
        try
            controlModel = fitglm(tblTrain, modelEquationControl, 'Distribution', 'binomial');
            yPredProbControl = predict(controlModel, tblTest);
            
            % Determine optimal threshold to maximize F1 score for control model
            [optimalThresholdControl, prAucControl, rocAucControl, sensitivityControl, specificityControl] = calculateOptimalThreshold(responseTest, yPredProbControl);
            
            % Convert probabilities to binary predictions using optimal threshold for control model
            yPredBinaryControl = yPredProbControl >= optimalThresholdControl;
            
            % Store control model results
            ypredAllControl = [ypredAllControl; yPredProbControl];
            prAucAllControl = [prAucAllControl; prAucControl];
            rocAucAllControl = [rocAucAllControl; rocAucControl];
            sensitivityAllControl = [sensitivityAllControl; sensitivityControl];
            specificityAllControl = [specificityAllControl; specificityControl];
        catch
            % Handle cases where control model fitting fails
            ypredAllControl = [ypredAllControl; nan(size(responseTest))];
            prAucAllControl = [prAucAllControl; nan];
            rocAucAllControl = [rocAucAllControl; nan];
            sensitivityAllControl = [sensitivityAllControl; nan];
            specificityAllControl = [specificityAllControl; nan];
        end
    catch
        % Handle cases where model fitting fails
        ypredAll = [ypredAll; nan(size(responseTest))];
        yactualAll = [yactualAll; responseTest];
        prAucAll = [prAucAll; nan];
        rocAucAll = [rocAucAll; nan];
        sensitivityAll = [sensitivityAll; nan];
        specificityAll = [specificityAll; nan];
        
        % Handle cases where control model fitting fails
        ypredAllControl = [ypredAllControl; nan(size(responseTest))];
        prAucAllControl = [prAucAllControl; nan];
        rocAucAllControl = [rocAucAllControl; nan];
        sensitivityAllControl = [sensitivityAllControl; nan];
        specificityAllControl = [specificityAllControl; nan];
    end
end



%% Recalculate Overall Performance Metrics

% Compute Precision-Recall and ROC curves
[precisionOverall, recallOverall, thresholdsPR] = perfcurve(yactualAll, ypredAll, 1, 'XCrit', 'recall', 'YCrit', 'prec');
[rocXOverall, rocYOverall, thresholdsROC] = perfcurve(yactualAll, ypredAll, 1);

% Compute PR-AUC and ROC-AUC
overallPrAuc = trapz(recallOverall, precisionOverall);
overallRocAuc = trapz(rocXOverall, rocYOverall);

% Compute sensitivity, specificity, PPV (precision), and NPV
[optimalThresholdOverall, ~, ~, sensitivityOverall, specificityOverall] = calculateOptimalThreshold(yactualAll, ypredAll);

tp = sum((ypredAll >= optimalThresholdOverall) & (yactualAll == 1));
fp = sum((ypredAll >= optimalThresholdOverall) & (yactualAll == 0));
tn = sum((ypredAll < optimalThresholdOverall) & (yactualAll == 0));
fn = sum((ypredAll < optimalThresholdOverall) & (yactualAll == 1));

ppvOverall = tp / (tp + fp);
npvOverall = tn / (tn + fn);

% Repeat for Control Model
[precisionOverallControl, recallOverallControl, thresholdsPRControl] = perfcurve(yactualAll, ypredAllControl, 1, 'XCrit', 'recall', 'YCrit', 'prec');
[rocXOverallControl, rocYOverallControl, thresholdsROCControl] = perfcurve(yactualAll, ypredAllControl, 1);

overallPrAucControl = trapz(recallOverallControl, precisionOverallControl);
overallRocAucControl = trapz(rocXOverallControl, rocYOverallControl);

[optimalThresholdOverallControl, ~, ~, sensitivityOverallControl, specificityOverallControl] = calculateOptimalThreshold(yactualAll, ypredAllControl);

tpControl = sum((ypredAllControl >= optimalThresholdOverallControl) & (yactualAll == 1));
fpControl = sum((ypredAllControl >= optimalThresholdOverallControl) & (yactualAll == 0));
tnControl = sum((ypredAllControl < optimalThresholdOverallControl) & (yactualAll == 0));
fnControl = sum((ypredAllControl < optimalThresholdOverallControl) & (yactualAll == 1));

ppvOverallControl = tpControl / (tpControl + fpControl);
npvOverallControl = tnControl / (tnControl + fnControl);

%% Display Overall Results

disp('Overall Predictive Model Results:');
disp(['Optimal Threshold: ', num2str(optimalThresholdOverall)]);
disp(['PR-AUC: ', num2str(overallPrAuc)]);
disp(['ROC-AUC: ', num2str(overallRocAuc)]);
disp(['Sensitivity: ', num2str(sensitivityOverall)]);
disp(['Specificity: ', num2str(specificityOverall)]);
disp(['PPV (Precision): ', num2str(ppvOverall)]);
disp(['NPV: ', num2str(npvOverall)]);

disp('Overall Control Model Results:');
disp(['Optimal Threshold Control: ', num2str(optimalThresholdOverallControl)]);
disp(['PR-AUC Control: ', num2str(overallPrAucControl)]);
disp(['ROC-AUC Control: ', num2str(overallRocAucControl)]);
disp(['Sensitivity Control: ', num2str(sensitivityOverallControl)]);
disp(['Specificity Control: ', num2str(specificityOverallControl)]);
disp(['PPV (Precision) Control: ', num2str(ppvOverallControl)]);
disp(['NPV Control: ', num2str(npvOverallControl)]);


%% Plot PR Curve for Overall Performance
figure;
hold on;
plot(recallOverall, precisionOverall, '-b', 'LineWidth', 2);
plot(recallOverallControl, precisionOverallControl, '--r', 'LineWidth', 2);
xlabel('Recall');
ylabel('Precision');
title('Overall Precision-Recall Curve');
legend('Predictive Model', 'Control Model', 'Location', 'SouthEast');
hold off;

%% Plot ROC Curve for Overall Performance
figure;
hold on;
plot(rocXOverall, rocYOverall, '-b', 'LineWidth', 2);
plot(rocXOverallControl, rocYOverallControl, '--r', 'LineWidth', 2);
plot([0 1], [0 1], '--k'); % Diagonal line for random chance
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('Overall ROC Curve');
legend('Predictive Model', 'Control Model', 'Random Chance', 'Location', 'SouthEast');
hold off;


%% Initialize Variables for Averaging Feature Maps
numFolds = 5;  % Number of folds in cross-validation
varIdx = 3; % Update this to your actual number of model variables if different
averageFeatureMap = zeros(maskSize);  % Initialize an array to store the sum of feature maps
maskIndices = find(mask3D > threshold); % Use maskIndices from your setup

%% Load and Sum GRF-Corrected Feature Maps
for fold = 1:numFolds
    % Define the filename for the GRF-corrected feature map of the current fold and variable index
    grfFileName = fullfile(dataDir, sprintf('Z_ClusterThresholded_T_model_%s_var%d_fold%d_GRF_%.2f.nii', modelNames{1}, varIdx, fold, voxelPThreshold));
    
    % Load the GRF-corrected feature map
    [featureMap, ~] = y_Read(grfFileName);
    
    % Sum the feature maps across all folds and variables within the mask
    averageFeatureMap(maskIndices) = averageFeatureMap(maskIndices) + featureMap(maskIndices);
end

% Compute the average by dividing the sum by the total number of entries (numFolds)
averageFeatureMap(maskIndices) = averageFeatureMap(maskIndices) / (numFolds);

%% Save the Average Feature Map
averageFeatureMapFileName = fullfile(dataDir, sprintf('Average_Z_ClusterThresholded_T_model_%s_GRF_%.2f.nii', modelNames{1}, voxelPThreshold));
y_Write(averageFeatureMap, header, averageFeatureMapFileName);
disp('Average feature map saved successfully.');



% Define a function to calculate the optimal threshold for binary classification
function [optimalThreshold, prAuc, rocAuc, sensitivity, specificity] = calculateOptimalThreshold(trueLabels, predictedScores)
    [precision, recall, thresholdsPR] = perfcurve(trueLabels, predictedScores, 1, 'XCrit', 'recall', 'YCrit', 'prec');
    f1Scores = 2 * (precision .* recall) ./ (precision + recall);
    [~, optimalIdx] = max(f1Scores);
    optimalThreshold = thresholdsPR(optimalIdx);
    prAuc = trapz(recall, precision);
    
    [rocX, rocY, thresholdsROC] = perfcurve(trueLabels, predictedScores, 1);
    rocAuc = trapz(rocX, rocY);
    
    sensitivity = rocY(optimalIdx);
    specificity = 1 - rocX(optimalIdx);
end

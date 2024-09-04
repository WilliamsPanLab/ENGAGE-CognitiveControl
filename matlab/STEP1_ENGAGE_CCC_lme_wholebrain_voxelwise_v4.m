clear, clc, close all
parpool
pctRunOnAll javaaddpath '/Users/xuezhang/Downloads/Software/ParforProgMon'
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

%% Define models for hypothesis 1 & 2 (mechanistic markers) and 3 (predictive markers)
modelNames = {'mechanistic'; 'predictive'};
modelNum = length(modelNames);
modelEquations = {
    % Identifying cognitive control circuit as a neural mechanism engaged by the I-CARE intervention
    sprintf('%s ~ %s + %s:%s + %s:%s + %s:%s:%s + %s + %s + %s + %s + (%s|%s)', ...
    outcomeVariable, 'fMRI', 'fMRI', 'Time', 'fMRI', 'Group', 'fMRI', 'Time', 'Group', [beh_name{outcomeIndex}, '_baseline'], 'BMI_baseline', 'age', 'sex', '1', 'Subj');...
    
    % Identifying cognitive control circuit engagement as a predictor of later outcome improvement
    sprintf('%s ~ %s + %s:%s + %s:%s + %s:%s:%s + %s + %s + %s + %s + (%s|%s)', ...
    outcomeVariable, 'fMRI_2MO', 'fMRI_2MO', 'Time', 'fMRI_2MO', 'Group', 'fMRI_2MO', 'Time', 'Group', [beh_name{outcomeIndex}, '_baseline'], 'BMI_baseline', 'age', 'sex', '1', 'Subj')
    };
numModelVars = [9; 9];

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
load(fullfile(dataDir, 'Responder_list_6mo.mat')); load(fullfile(dataDir, 'Remission_list_6mo.mat'));
load(fullfile(dataDir, 'Responder_list_long.mat')); load(fullfile(dataDir, 'Remission_list_long.mat'));

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

%% Modeling
% Define output folder and prepare data for modeling
outputDir = fullfile(dataDir, 'MixedModel', taskName, contrastName, beh_name{outcomeIndex}, [maskSuffix, '_', modelTimePoints]);
mkdir(outputDir);
modelFStats = nan(length(maskIndices), max(numModelVars), length(modelNames));
modelDF = nan(length(maskIndices), max(numModelVars), 2, length(modelNames));
tic 

% Parallel loop for each region of interest (ROI)
ppm = ParforProgMon( 'Progress bar', length(maskIndices));
parfor roiIndex = 1: length(maskIndices)
    % Prepare fMRI data in long format
    fmriDataChange = squeeze(fmri_wb(:, roiIndex, :));
    fmriDataChange = fmriDataChange(:, includedSessions4Change) - repmat(fmriDataChange(:, 1), 1, length(includedSessions4Change)); 
    fmriDataChange = fmriDataChange(:); 

    fmriBaseline = squeeze(fmri_wb(:, roiIndex, 1)); 
    fmriBaseline = repmat(fmriBaseline, 1, length(includedSessions4Change)); 
    fmriBaseline = fmriBaseline(:); 
    
    fmri2Mo = squeeze(fmri_wb(:, roiIndex, 2)); 
    fmri2Mo = repmat(fmri2Mo, 1, length(includedSessions4Change)); 
    fmri2Mo = fmri2Mo(:) - fmriBaseline; 
    
    %% Mechanism model: △SCL-20 ~ △fMRI + △fMRI * Time + △fMRI * Group + △fMRI * Time * Group + SCL-20_baseline + age + gender + BMI
    % Identify data instances with no missing values for all variables
    validIndices = intersect(intersect(find(~isnan(fmriDataChange)), find(~isnan(outcomeChange))), ...
        intersect(find(~isnan(FDIncludeChange)), intersect(find(~isnan(ageSexLong(:, 1))), find(~isnan(ageSexLong(:, 2))))));
    
    % Only model for voxels where data is in full rank
    if length(validIndices) > 0.3 * length(includedSessions4Change) * numSubjects && ...
       rank([outcomeChange(validIndices), fmriDataChange(validIndices), ageSexLong(validIndices, :), FDIncludeChange(validIndices)]) == 5 && ...
       length(find(fmriBaseline(validIndices) == 0)) < 0.7 * length(includedSessions4Change) * numSubjects && ...
       length(find(fmriDataChange(validIndices) == 0)) < 0.7 * length(includedSessions4Change) * numSubjects
       
        % Form a table that contains all variables needed for modeling
        tblData = table(outcomeChange, baselineOutcome, timePoints, subjectIds, fmriBaseline, fmri2Mo, fmriDataChange, ...
            responseInfo, remissionInfo, ageSexLong(:, 1), ageSexLong(:, 2), groupInfo, baselineBMI, ...
            'VariableNames', {beh_name{outcomeIndex}, [beh_name{outcomeIndex}, '_baseline'], 'Time', 'Subj', ...
                        'fMRI_baseline', 'fMRI_2MO', 'fMRI', 'Responder', 'Remitter', 'age', 'sex', 'Group', 'BMI_baseline'});
        
        % Loop across all models
        for modelIndex = 1: modelNum
            % Fit linear mixed-effects model
            mixModel = fitlme(tblData, modelEquations{modelIndex}, 'DummyVarCoding', 'effects', 'StartMethod', 'random'); 
            
            % Extract ANOVA table for the model
            modelAnova = anova(mixModel); 
            
            % Store degrees of freedom and F-statistics
            modelDF(roiIndex, :, :, modelIndex) = [[modelAnova.DF1, modelAnova.DF2]; nan(max(numModelVars) - length(modelAnova.DF1), 2)]; 
            modelFStats(roiIndex, :, modelIndex) = [modelAnova.FStat; nan(max(numModelVars) - length(modelAnova.DF1), 1)]; 
        end
    end
    ppm.increment();
end
toc % End timing the operation


%% Run multiple comparison correction on F maps and save them as NIfTI files
for modelIndex = 1: length(modelNames)
    for subConIndex = 1: numModelVars(modelIndex)
        data = modelFStats(:, subConIndex, modelIndex); 
        data3D = zeros(maskSize); 
        data3D(maskIndices) = data; % Prepare 3D F-statistic map
        
        outputFileName = fullfile(outputDir, ['F_model', modelNames{modelIndex}, '_var', num2str(subConIndex), '.nii']);
        y_Write(data3D, header, outputFileName); % Write F-statistic map to NIfTI file
        
        [~, filename, ~] = fileparts(outputFileName);
        outputNameCorr = fullfile(outputDir, [filename, '_GRF_', num2str(voxelPThreshold), '.nii']); % Define output file name for GRF correction
        
        % Perform GRF correction
        df1 = max(modelDF(~isnan(modelDF(:, subConIndex, 1, modelIndex)), subConIndex, 1, modelIndex));
        df2 = max(modelDF(~isnan(modelDF(:, subConIndex, 2, modelIndex)), subConIndex, 2, modelIndex));
        y_GRF_Threshold(outputFileName, voxelPThreshold, isTwoTailed, clusterPThreshold, outputNameCorr, maskGRFOutName, flag, df1, df2);
        
        % Save degrees of freedom
        save(fullfile(outputDir, [modelNames{modelIndex}, '_item', num2str(subConIndex), '_DF.mat']), 'df1', 'df2');
        clear df1; clear df2; clear data; clear data3D;
    end
end

%% Visualize results that survived multiple comparison correction in volume space
[DPABIPath, ~, ~] = fileparts(which('DPABI.m')); % Get DPABI path
underlayFileName = fullfile(DPABIPath, 'Templates', 'ch2.nii'); % Define underlay file

for modelIndex = 1: length(modelNames)
    for varIndex = 1:numModelVars(modelIndex)
        imageFile = fullfile(outputDir, ['Z_ClusterThresholded_F_model', modelNames{modelIndex}, '_var', num2str(varIndex), '_GRF_', num2str(voxelPThreshold), '.nii']);
        
        if exist(imageFile, 'file')
            [~, filename, ~] = fileparts(imageFile);
            hView = w_Call_DPABI_VIEW(imageFile, minNegative, minPositive, clusterSize, connectivityCriterion, underlayFileName, colorMap, maxNegative, maxPositive);
            [imageMontage, ~] = w_MontageImage([-46:sliceInterval:-22; -22:sliceInterval:2; 2:sliceInterval:26; 26:sliceInterval:50; 50:sliceInterval:74], 'T', hView);
            imageMontage = flipdim(imageMontage, 1);
            imwrite(imageMontage, fullfile(outputDir, [filename, '.tif']));
            close all
        end
    end
end

%% Render results that survived multiple comparison correction in surface space
maskSurface = [];
space = 'fsaverage'; % Space for surface rendering ('fsaverage' or 'fsaverage5')
surfaceType = 'white'; % Surface type ('white' or 'pial')
DPABISurfPath = fileparts(which('DPABISurf.m')); % Get DPABI Surf path
surfUnderlay = {
    fullfile(DPABISurfPath, 'SurfTemplates', [space, '_lh_', surfaceType, '.surf.gii']);
    fullfile(DPABISurfPath, 'SurfTemplates', [space, '_rh_', surfaceType, '.surf.gii'])
}; % Define underlay files for surface rendering

for modelIndex = 1: length(modelNames)
    for varIndex = 1:numModelVars(modelIndex)
        imageFile = fullfile(outputDir, ['Z_ClusterThresholded_F_model', modelNames{modelIndex}, '_var', num2str(varIndex), '_GRF_', num2str(voxelPThreshold), '.nii']);
        
        if exist(imageFile, 'file')
            [~, filename, ~] = fileparts(imageFile);
            tic % Start timing the rendering
            mkdir(fullfile(outputDir, [space, '_', surfaceType])); % Create output directory for surface images
            jpgFile = fullfile(outputDir, [space, '_', surfaceType, filename, '.tif']); % Define output file name for the surface image
            system(['killall Docker && open /Applications/Docker.app']); % Restart Docker (necessary for surface rendering)
            pause(100); % Pause to ensure Docker starts
            y_Call_DPABISurf_VIEW_FromVolume(imageFile, jpgFile, minNegative, minPositive, maskSurface, '', clusterSize, connectivityCriterion, space, surfUnderlay, colorMap, maxNegative, maxPositive);
            close all
            toc % End timing the rendering
        end
    end
end
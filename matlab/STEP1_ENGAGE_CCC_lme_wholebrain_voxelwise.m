%% if interested in other contrasts than nogovsgo, should modify the code to save other contrasts first
%% !!!!! choose time as catergory or continuous variables
%% can do continous and binary LMM or LM

clear, clc, close all
parpool
% pctRunOnAll javaaddpath '/Users/xuezhang/Downloads/Software/ParforProgMon'

%% load beh_all and beh_name
load(fullfile('..','engage','beh_all_17-Oct-2021.mat'));

%% read fmri performance
load(['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/engage/Behavior/gng/ENGAGE/gng_performance.mat']);
beh_name_itask = cellfun(@(x) ['fmri_', x], beh_name_itask,'UniformOutput',false);
fmri_gonogo_data = task_beh_data;

%% combine with beh_all
beh_all = [beh_all, fmri_gonogo_data];
beh_name = [beh_name; beh_name_itask];

%% need pre-definition please

% behavior - 24: SCL-20; 25: SPSI
i_beh = 24;

% time as catergory or continuous variables
time_catergorical = 1;
actual_time = [0, 2, 6, 12, 24];
include_ses = [1:5]; include_ses_change = setdiff(include_ses, [1]);

% binary prediction or continuous prediction
binary_outcome = 0; 
binary_suffix = '';
y_variable = beh_name{i_beh};

% GRF correction
VoxelPThreshold = 0.001;
IsTwoTailed = 0;
ClusterPThreshold = 0.05;

% plot parameters
NMin = -3; PMin = 3;
NMax = -4; PMax = 4;

ColorMap=y_AFNI_ColorMap(12);%jet(128);
Flag = 'F';
ClusterSize=0;
ConnectivityCriterion=18;
slice_interval = 6;


%% main script
load_plot_parameters

subjects_list = {'LA13012';'LA13272';'LA14016';'MV00878';'MV00962';'MV00992';'MV01113';'MV01438';'MV01836';'MV01950';'MV04661';...
    'MV05158';'MV05953';'MV06084';'MV06904';'MV07296';'MV07303';'MV07572';'MV07647';'MV08032';'MV08112';'MV08176';'MV08645';...
    'MV08712';'MV08866';'MV09122';'MV09305';'MV09434';'MV09441';'MV09560';'MV09586';'MV09876';'MV11065';'MV11133';'MV11135';...
    'MV11150';'MV11202';'PA20147';'PA21728';'PA21991';'PA22014';'PA22518';'PA22544';'PA22561';'PA22568';'PA22594';'PA22725';...
    'PA22728';'PA22772';'PA23284';'PA23955';'PA24195';'PA24326';'PA24603';'PA24859';'PA24876';'PA25084';'PA25119';'PA25306';...
    'PA25642';'PA25692';'PA25870';'PA25894';'PA25960';'PA25994';'PA26039';'PA26203';'PA26376';'PA26623';'PA26650';'PA26904';...
    'PA27040';'PA27394';'PA27432';'PA27434';'PA27493';'PA27541';'PA27578';'PA27784';'PA27793';'PA27962';'PA27995';'PA28033';...
    'PA28219';'PA28336';'PA28460';'PA28464';'PA28564';'PA28985';'PA28989';'PA29385';'PA29661';'PA29685';'PA29689';'PA30071';...
    'PA30104';'PA30563';'PA30677';'PA30861';'PA30862';'PA30895';'PA30973';'SU30700';'SU30734';'SU30816';'SU31067';'SU33550';...
    'SU35282'};
session_list = {'000';'2MO';'6MO';'12MO';'24MO'};
task_list = {'111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO';};
task_contrasts_select = {{'NogovsGo'}}
task_contrasts_select_con = {{'03'}};

%% load FD
load(fullfile('..','engage','FD_alltasks.mat'));% FD_alltasks
FD_alltasks = FD_alltasks(:,:,3);
perc_motion = 0.25;
frame_N = 151; 

%% load demo and treatment related info
age_sex = beh_all(:,19:20, 1);
% load intervention
load(fullfile('..','engage','Interv_Info.mat')); load(fullfile('..','engage','Responder_list_6mo.mat')); load(fullfile('..','engage','Remission_list_6mo.mat'));
load(fullfile('..','engage','Responder_list_long.mat')); load(fullfile('..','engage','Remission_list_long.mat'));

% more parameters for modeling
% included timepoints
timepoints = repmat(include_ses_change, length(subjects_list), 1); timepoints = timepoints(:);
model_timepoints = [];
for i_include_ses = include_ses
    model_timepoints = [model_timepoints, num2str(actual_time(i_include_ses)),'_'];
end
model_timepoints = [model_timepoints, 'MO'];

if time_catergorical == 1
    timepoints = nominal(timepoints);
    tp_name = 'Cater';
else
    tp_name = 'Conti';
end
timepoints_weight = (repmat(actual_time(include_ses_change), length(subjects_list), 1));timepoints_weight = timepoints_weight(:);
subjects = repmat([1: length(subjects_list)]', 1, length(include_ses_change)); subjects =  subjects(:); subjects = nominal(subjects);
age_sex_tl = repmat(age_sex,length(include_ses_change),1);
group_info = Interv_Info; group_info = repmat(group_info, 1, length(include_ses_change)); group_info = nominal(group_info(:));
response_info = 2-resp_ident_long(:, include_ses_change); response_info = response_info(:);
remit_info = 2-remi_ident_long(:, include_ses_change); remit_info = remit_info(:); 
% prepare change
beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));beh_ibeh_change = beh_ibeh_change(:);
beh_ibeh_change_weighted = beh_ibeh_change./timepoints_weight; beh_ibeh_change_weighted(isinf(beh_ibeh_change_weighted)) = nan;
beh_ibeh_baseline = repmat(squeeze(beh_all(:, i_beh, 1)), 1, length(include_ses_change)); beh_ibeh_baseline = beh_ibeh_baseline(:);
beh_ibeh_change_percent = beh_ibeh_change./beh_ibeh_baseline; beh_ibeh_change_percent(isinf(beh_ibeh_change_percent)) = nan;
% prepare baseline BMI
bmi_baseline = repmat(squeeze(beh_all(:, 21, 1)), 1, length(include_ses_change));bmi_baseline = bmi_baseline(:);
 

% model name
model_names = {'mechanism';'predict'};
model_equations = {
    % Identifying cognitive control circuit as a neural mechanisam engaged by the I-CARE intervention
    sprintf('%s ~ %s + %s:%s +%s:%s+ %s:%s:%s+ %s + %s + %s + %s + (%s|%s)', ...
    y_variable,'fMRI', 'fMRI', 'Time','fMRI', 'Group', 'fMRI', 'Time','Group', [beh_name{i_beh},'_baseline'],'BMI_baseline','age','sex','1',char('Subj'));...
    
    % Identifying cognitive control circuit engagement as a predictor of later behavioral improvement
    sprintf('%s ~ %s + %s:%s +%s:%s+ %s:%s:%s+ %s + %s + %s + %s + (%s|%s)', ...
    y_variable,'fMRI_2MO', 'fMRI_2MO', 'Time','fMRI_2MO', 'Group', 'fMRI_2MO', 'Time','Group',[beh_name{i_beh},'_baseline'],'BMI_baseline','age','sex','1',char('Subj'))
    
    }
model_varN = [9;9];



%% define Mask
% grey matter
Mask=['/Users/xuezhang/Downloads/Software/DPABI_V4.3_200401/Templates/GreyMask_02_91x109x91.img'];
[Mask3d,~,~, Header] = y_ReadAll(Mask); Masksize = size(Mask3d); thresh=0; MaskInd = find(Mask3d > thresh);
Header.pinfo = [1;0;0]; Header.dt =[16,0]; Mask_suffix = ['GM_',num2str(thresh)];
% make mask for GRF correction
Mask_GRF= zeros(Masksize); Mask_GRF(MaskInd) = 1; [a, b, c] = fileparts(Mask); Mask_GRF_outname = fullfile(a,[b,'_thresh_',num2str(thresh),c]); y_Write(Mask_GRF, Header, Mask_GRF_outname);

%% prepare fmri data
DATADIR = fullfile('..','engage');
fmri_wb = nan(length(subjects_list), length(MaskInd), length(session_list));
for i_task = 1
    for i_con = 1
        for i_ses = 1: length(session_list)
            for sub = 1:length(subjects_list)
                subject = subjects_list{sub};
%                 tic
                if FD_alltasks(sub, i_ses, i_task) <= frame_N * perc_motion
                    data = [DATADIR,'/',task_list{i_task},'/', session_list{i_ses},'/',subject,'/con_00', task_contrasts_select_con{i_task}{i_con},'_mni.nii'];
                    if exist(data,'file')
                        f3d = y_ReadAll(data);
                        fmri_wb(sub,:,i_ses) = f3d(MaskInd);
                        clear f3d
                    end
                    clear data
                end
%                 toc
            end
        end
        mkdir(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}))
        save(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_wholebrain_',Mask_suffix,'.mat']), 'fmri_wb', '-v7.3');
    end
end
      

%% modeling
for i_task = 1
    for i_con = 1
        OutDir = fullfile('..','engage','SPM','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh},[Mask_suffix,'_',model_timepoints,'_',binary_suffix]);
        mkdir(OutDir)
        % load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_wholebrain_',Mask_suffix,'.mat']));
        
        FD_include_change = FD_alltasks(:, include_ses_change, i_task) - FD_alltasks(:, 1, i_task); FD_include_change = FD_include_change(:);
        FD_include_baseline = repmat(squeeze(FD_alltasks(:, 1, i_task)), 1, length(include_ses_change)); FD_include_baseline = FD_include_baseline(:);
        % output
        model_imodel_F = nan(length(MaskInd), max(model_varN), length(model_names));
        model_df = nan(length(MaskInd),max(model_varN), 2, length(model_names));
        tic
%         ppm = ParforProgMon( 'Progress bar', length(MaskInd));
        parfor i_ROI = 1: length(MaskInd)
            % create change data
            fmri_data_change = prepare_data_change(squeeze(fmri_wb(:,i_ROI, :)), include_ses);
            % fmri baseline
            fmri_baseline = squeeze(fmri_wb(:,i_ROI, 1)); 
            fmri_baseline = repmat(fmri_baseline, 1, length(include_ses_change)); fmri_baseline = fmri_baseline(:);
            % fmri 2MO
            fmri_2mo = squeeze(fmri_wb(:,i_ROI, 2)); 
            fmri_2mo = repmat(fmri_2mo, 1, length(include_ses_change)); fmri_2mo = fmri_2mo(:) - fmri_baseline;
            
            %% mechanism: △SCL-20 ~ △fMRI + △fMRI * Time + △fMRI * Group + △fMRI* Time * Group + SCL-20_baseline + age + gender + BMI
            inter_list = intersect(intersect(find(~isnan(fmri_data_change)), find(~isnan(beh_ibeh_change))), ...
                intersect(find(~isnan(FD_include_change)), intersect(find(~isnan(age_sex_tl(:,1))), find(~isnan(age_sex_tl(:,2))))));
            
            if length(inter_list) > 0.3*length(include_ses_change)*length(subjects_list) && rank([beh_ibeh_change(inter_list), fmri_data_change(inter_list), age_sex_tl(inter_list,:), FD_include_change(inter_list)]) == 5 ...
                    && length(find(fmri_baseline(inter_list)==0))< 0.7*length(include_ses_change)*length(subjects_list) && length(find(fmri_data_change(inter_list)==0))< 0.7*length(include_ses_change)*length(subjects_list)
               
                tbl_data = table(beh_ibeh_change, beh_ibeh_baseline, timepoints,subjects, fmri_baseline, fmri_2mo, fmri_data_change, response_info, remit_info, age_sex_tl(:,1), age_sex_tl(:,2), group_info,bmi_baseline,...
                    'VariableNames',{beh_name{i_beh},[beh_name{i_beh},'_baseline'],'Time','Subj','fMRI_baseline', 'fMRI_2MO', 'fMRI','Responder','Remitter','age','sex', 'Group','BMI_baseline'});
                % loop across all models
                for i_model = 1: length(model_names)
                    mixmodel_imodel = fitlme(tbl_data, model_equations{i_model},...
                        'DummyVarCoding','effects','StartMethod','random');
                    ss_imodel = anova(mixmodel_imodel);
                    model_df(i_ROI, :, :, i_model) = [[ss_imodel.DF1, ss_imodel.DF2];nan(max(model_varN) -length(ss_imodel.DF1), 2)];
                    model_imodel_F(i_ROI,:,i_model) = [ss_imodel.FStat;nan(max(model_varN) -length(ss_imodel.DF1), 1)]; %model_imodel_p(i_ROI, :) = ss_imodel.pValue;
                end
            end
%             ppm.increment();
        end
        toc
        %% save data
        for i_model = 1: length(model_names)
            for i_sub_con = 1: model_varN(i_model)
                data = model_imodel_F(:, i_sub_con, i_model); data_3d = zeros(Masksize); data_3d(MaskInd) = data;
                outname = [OutDir,'/F_model', model_names{i_model},'_var',num2str(i_sub_con), '_Time_',tp_name,'.nii'];
                y_Write(data_3d, Header, outname);
                [~,filename,~] = fileparts(outname);
                OutputName_corr = [OutDir,'/',filename,'_GRF_',num2str(VoxelPThreshold), '.nii'];
                Df1 = max(model_df(~isnan(model_df(:, i_sub_con, 1,i_model)), i_sub_con, 1, i_model));
                Df2 = max(model_df(~isnan(model_df(:, i_sub_con, 2,i_model)), i_sub_con, 2, i_model));
                y_GRF_Threshold(outname,VoxelPThreshold,IsTwoTailed,ClusterPThreshold,OutputName_corr,Mask_GRF_outname,Flag,Df1,Df2);
                
                save([OutDir,'/',model_names{i_model},'_item',num2str(i_sub_con),'_DF.mat'],'Df1','Df2')
                clear Df1; clear Df2; clear data; clear data_3d;
            
            end
        end
        
    end
end

                    
%% plot results

% Volume Rendering
[DPABIPath, fileN, extn] = fileparts(which('DPABI.m'));
UnderlayFileName = [DPABIPath, filesep, 'Templates', filesep, 'ch2.nii'];
for i_task = 1
    for i_con = 1
        DataDir = fullfile('..','engage','SPM','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh},[Mask_suffix,'_',model_timepoints,'_',binary_suffix]);
        Plot_dir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, [Mask_suffix,'_',model_timepoints,'_',binary_suffix]);
        mkdir(Plot_dir);
        for i_model = 1: length(model_names)
            for i_var = 1:model_varN(i_model)
                ImageFile0 = [DataDir,'/Z_ClusterThresholded_F_model', model_names{i_model}, '_var',num2str(i_var), '_Time_',tp_name,'_GRF_',num2str(VoxelPThreshold),'.nii'];
                if exist(ImageFile0,'file')
                    [~,filename,~] = fileparts(ImageFile0);
                    H0=w_Call_DPABI_VIEW(ImageFile0,NMin,PMin,ClusterSize,ConnectivityCriterion,UnderlayFileName,ColorMap,NMax,PMax);
                    [ImageA0,Space0]=w_MontageImage([-46:slice_interval:-22;-22:slice_interval:2;2:slice_interval:26;26:slice_interval:50;50:slice_interval:74],'T',H0);
                    ImageA0=flipdim(ImageA0,1);
                    imwrite(ImageA0,[Plot_dir,'/',filename, '.tif']);
                    close all
                end
            end
            
        end
    end
end

% Surface Rendering
Mask_sur=[];   
Space='fsaverage';%'fsaverage5'; 
surface = 'white';
DPABISurfPath=fileparts(which('DPABISurf.m'));
SurfUnderlay={fullfile(DPABISurfPath,'SurfTemplates',[Space,'_lh_',surface,'.surf.gii']);fullfile(DPABISurfPath,'SurfTemplates',[Space,'_rh_',surface,'.surf.gii'])};
for i_task = 1
    for i_con = 1
        DataDir = fullfile('/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis','engage','SPM','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh},[Mask_suffix,'_',model_timepoints,'_',binary_suffix]);
        Plot_dir = fullfile('/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis','engage','SPM','Figures','Level2','MixedModel_Surface',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, [Mask_suffix,'_',model_timepoints,'_',binary_suffix]);
        mkdir(Plot_dir);
        for i_model = 1: length(model_names)
            for i_var = 1: length(model_names)%[3:4, 6:9]%[9:10]%
                ImageFile0 = [DataDir,'/Z_ClusterThresholded_F_model', model_names{i_model}, '_var',num2str(i_var), '_Time_',tp_name,'_GRF_',num2str(VoxelPThreshold),'.nii'];
                if exist(ImageFile0,'file')
                    [~,filename,~] = fileparts(ImageFile0);
                    tic
                    mkdir([Plot_dir,'/',Space,'_',surface])
                    JPGFile = [Plot_dir,'/',Space,'_',surface,'/',filename,'.tif'];
                    system(['killall Docker && open /Applications/Docker.app'])% restart docker
                    pause(100)
                    y_Call_DPABISurf_VIEW_FromVolume(ImageFile0,JPGFile,NMin,PMin,Mask_sur,'',ClusterSize,ConnectivityCriterion,Space,SurfUnderlay,ColorMap,NMax,PMax);
                    close all
                    toc 
                end
            end
        end
    end
end



function [data_change] = prepare_data_change(data, include_ses)
% data_fil = filloutliers(data,nan,'mean');
data_fil = data;
include_ses_change = setdiff(include_ses, [1]);
data_change = data_fil(:, include_ses_change) - repmat(data_fil(:, 1), 1, length(include_ses_change));% 108 * length(include_ses_change)
data_change = data_change(:);
end











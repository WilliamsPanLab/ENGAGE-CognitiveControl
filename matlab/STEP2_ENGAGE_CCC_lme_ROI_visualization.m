%% if interested in other contrasts than nogovsgo, should modify the code to save other contrasts first
clear, clc, close all
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
task_contrasts_select = {{'NogovsGo'}};
task_contrasts_select_con = {{'03'}};
%% load FD
load(fullfile('..','engage','FD_alltasks.mat'));% FD_alltasks
FD_alltasks = FD_alltasks(:,:,3);
perc_motion = 0.25;
frame_N = 151; 
%% load data
% fmri gonogo_data: task_beh_data,beh_name_itask
load(['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/engage/Behavior/gng/ENGAGE/gng_performance.mat']);
beh_name_itask = cellfun(@(x) ['fmri_', x], beh_name_itask,'UniformOutput',false);
fmri_gonogo_data = task_beh_data;


% load intervention
load(fullfile('..','engage','Interv_Info.mat'));
load(fullfile('..','engage','Responder_list_6mo.mat'))
load(fullfile('..','engage','Remission_list_6mo.mat'))
load(fullfile('..','engage','beh_all_18-Oct-2021.mat'));
age_sex = beh_all(:,19:20, 1);


%% combine with beh_all
beh_all = [beh_all, fmri_gonogo_data];
beh_name = [beh_name; beh_name_itask];

%% ROIs from SCL-20
ROI_dir = ['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Mask/engage/Cor_Beh/MixedModel/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/SCL_20'];
ROI_label = ['MixedModel_SCL20'];
ROI_list = {'lAIn';'lCaudate';'rIPG';'rIFG'; 'rIPG3';'lSPL'};
ROI_plot_list = {'lAIn';'lCaudate';'rIPL';'rDLPFC'; 'rIPL';'lSPL'};

%% ROIs from spsi
% ROI_dir = ['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Mask/engage/Cor_Beh/MixedModel/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/SPSI'];
% ROI_label = ['MixedModel_SPSI'];
% ROI_list = {'lDLPFC1';'rSMG';'Precuneus2';'lFusiform';'rIPL';'lIPL'};
% ROI_plot_list = {'lDLPFC';'rIPL';'Precuneus';'lFusiform';'rIPL';'lIPL'};

%% prepare fmri data
radius = 2;
ROI_list = cellfun(@(x) [x, '_',num2str(radius)], ROI_list, 'UniformOutput', false);
DATADIR = fullfile('..','engage');
for i_task = 1
    for i_con = 1
        task_activation_ROI_icon_itask = nan(length(subjects_list), length(ROI_list), length(session_list));
        for i_ses = 1: length(session_list)
            scan = {};
            for sub = 1:length(subjects_list)
                subject = subjects_list{sub};
                data = [];
                if FD_alltasks(sub, i_ses, i_task) <= frame_N * perc_motion
                    data = [DATADIR,'/',task_list{i_task},'/', session_list{i_ses},'/',subject,'/con_00', task_contrasts_select_con{i_task}{i_con},'_mni.nii'];
                end
                scan{sub,1}=data;
            end
            for i_ROI = 1: length(ROI_list)
                ROI_i_dir = fullfile(ROI_dir, [ROI_list{i_ROI},'.nii']);
                task_activation_ROI_icon_itask(:,i_ROI, i_ses) =ExtValMul3D(scan,ROI_i_dir);
            end
        end
        mkdir(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}))
        save(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']), 'task_activation_ROI_icon_itask', 'ROI_list', 'ROI_dir'); 
    end
end


%% generating statistical results for all LMMs
include_ses = [1,3:5];
actual_time = [0, 2, 6, 12, 24];
covariates = repmat(age_sex,[1 1 5]);
covariates_name = {'age';'gender'};
% outcome
i_beh = 24;% scl-20: 24 spsi 25;

 for i_task = 1
     for i_con = 1
        
        OutDir = fullfile('..','engage','SPM','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
        mkdir(OutDir)
%         load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
        %% mechanism and predict LME - include changes at all timepoints
        % prepare data
        include_ses_change = setdiff(include_ses, [1]);
        timepoints = repmat(include_ses_change, length(subjects_list), 1); timepoints = timepoints(:); timepoints = nominal(timepoints);
        timepoints1 = repmat(include_ses_change, length(subjects_list), 1); timepoints1 = timepoints1(:); timepoints1 = nominal(timepoints1);
        timepoints_weight = (repmat(actual_time(include_ses_change), length(subjects_list), 1));timepoints_weight = timepoints_weight(:);
        subjects = repmat([1: length(subjects_list)]', 1, length(include_ses_change)); subjects =  subjects(:); subjects = nominal(subjects);
        FD_include_change = FD_alltasks(:, include_ses_change, i_task) - FD_alltasks(:, 1, i_task); FD_include_change = FD_include_change(:);
        FD_include_baseline = repmat(squeeze(FD_alltasks(:, 1, i_task)), 1, length(include_ses_change)); FD_include_baseline = FD_include_baseline(:);
        age_sex_tl = repmat(age_sex,length(include_ses_change),1);
        group_info = Interv_Info; group_info = repmat(group_info, 1, length(include_ses_change)); group_info = nominal(group_info(:));
        beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_change_percent = beh_ibeh_change./repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_change = beh_ibeh_change(:);
        beh_ibeh_change_weighted = beh_ibeh_change./timepoints_weight;
        beh_ibeh_change_percent = beh_ibeh_change_percent(:);
        beh_ibeh_baseline = repmat(squeeze(beh_all(:, i_beh, 1)), 1, length(include_ses_change));
        beh_ibeh_baseline = beh_ibeh_baseline(:);
        % prepare baseline BMI
        bmi_baseline = repmat(squeeze(beh_all(:, 21, 1)), 1, length(include_ses_change));bmi_baseline = bmi_baseline(:);

        for i_ROI = 1: length(ROI_list)
            ROI_list{i_ROI}
            % fmri change
            task_activation_ROI_icon_itask_iroi_change = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses_change)) ...
                - repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
            task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change(:);
            % baseline fmri
            task_activation_ROI_icon_itask_iroi_baseline = repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
            task_activation_ROI_icon_itask_iroi_baseline = task_activation_ROI_icon_itask_iroi_baseline(:);
            % 2mo fmri
            task_activation_ROI_icon_itask_iroi_2mo = repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 2) - task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
            task_activation_ROI_icon_itask_iroi_2mo = task_activation_ROI_icon_itask_iroi_2mo(:);
            % data -> table
            tbl2 = table(beh_ibeh_change, beh_ibeh_baseline, timepoints,subjects, FD_include_baseline, task_activation_ROI_icon_itask_iroi_baseline, task_activation_ROI_icon_itask_iroi_2mo, task_activation_ROI_icon_itask_iroi_change, age_sex_tl(:,1), age_sex_tl(:,2), group_info,bmi_baseline,...
                'VariableNames',{beh_name{i_beh},[beh_name{i_beh},'_baseline'],'Time','Subj','fd', 'fMRI_baseline', 'fMRI_2MO', 'fMRI','age','sex', 'Group','BMI_baseline'});
            tbl2.sex = categorical(tbl2.sex); tbl2.Group = nominal(tbl2.Group);

            if i_ROI < 5 % brain-behavior association
                
                tbl2 = rmmissing(tbl2, 'DataVariables',{'fMRI',beh_name{i_beh}, 'Time','Group',[beh_name{i_beh},'_baseline'],'age','sex', 'BMI_baseline'});
                mixmodel2_2 = fitlme(tbl2, sprintf('%s ~ %s + %s:%s +%s:%s+ %s:%s:%s+ %s + %s + %s + %s + (%s|%s)', ...
                    beh_name{i_beh},'fMRI', 'fMRI', 'Time','fMRI', 'Group', 'fMRI', 'Time','Group', [beh_name{i_beh},'_baseline'],'age','sex','BMI_baseline','1',char('Subj')),...
                    'DummyVarCoding','effects','StartMethod','random');
                ss2_2 = anova(mixmodel2_2,'DFMethod','satterthwaite'); ss2_2_table = dataset2table(ss2_2)
               
                % comparing to control model
                mixmodel2_2_baseline = fitlme(tbl2, sprintf('%s ~ %s + %s + %s + (%s|%s)', ...
                    beh_name{i_beh},[beh_name{i_beh},'_baseline'],'age','sex','1',char('Subj')),...
                    'DummyVarCoding','effects','StartMethod','random');
                ss2_2_comparison = dataset2table(compare(mixmodel2_2_baseline, mixmodel2_2));
                writetable(ss2_2_comparison,[OutDir,'/model2_2_likelihood.xls'], 'Sheet', ROI_list{i_ROI});
            else % brain-behavior prediction
                
                tbl2 = rmmissing(tbl2, 'DataVariables',{'fMRI_2MO',beh_name{i_beh}, 'Time','Group',[beh_name{i_beh},'_baseline'],'age','sex'});
                mixmodel6_3 = fitlme(tbl2, sprintf('%s ~ %s + %s:%s +%s:%s+ %s:%s:%s+ %s + %s + %s + %s + (%s|%s)', ...
                    [beh_name{i_beh}],'fMRI_2MO', 'fMRI_2MO', 'Time','fMRI_2MO', 'Group', 'fMRI_2MO', 'Time','Group',[beh_name{i_beh},'_baseline'],'age','sex','BMI_baseline','1',char('Subj')),...
                    'DummyVarCoding','effects','StartMethod','random');
                ss6_3 = anova(mixmodel6_3,'DFMethod','satterthwaite'); ss6_3_table = dataset2table(ss6_3)
                
                
                % comparing to control model
                mixmodel6_3_baseline = fitlme(tbl2, sprintf('%s ~ %s + %s + %s + (%s|%s)', ...
                    [beh_name{i_beh}],[beh_name{i_beh},'_baseline'],'age','sex','1',char('Subj')),...
                    'DummyVarCoding','effects','StartMethod','random');
                ss6_3_comparison = dataset2table(compare(mixmodel6_3_baseline, mixmodel6_3))
                writetable(ss6_3_comparison,[OutDir,'/model6_3_likelihood.xls'], 'Sheet', ROI_list{i_ROI});
            end
            
        end
        
     end
 end



%% creating scatter plots for above significant ROIs
include_ses_change = setdiff(include_ses, [1]);
age_sex_tl = repmat(age_sex,length(include_ses_change),1);
beh_name_plot = beh_name;
beh_name_plot{24} = 'Depression symptoms';
beh_name_plot{25} = 'Problem-solving ability';

%% one group: intervention
group_Info = Interv_Info;%resp_ident;%remi_ident;%
leg = {['U-CARE'],['I-CARE'];};
group_Info_ses_change = repmat(group_Info, 1, length(include_ses_change)); group_Info_ses_change = group_Info_ses_change(:);
group_Info_ses_change_name = cell(length(group_Info_ses_change),1);
group_Info_ses_change_name(group_Info_ses_change == 0) = {leg{1}};
group_Info_ses_change_name(group_Info_ses_change == 1) = {leg{2}};
group_N = 2;

%% one grouop: time
% group_Info = [ones(length(subjects_list),1), 2*ones(length(subjects_list),1), 3*ones(length(subjects_list),1),4*ones(length(subjects_list),1),5*ones(length(subjects_list),1)];%Interv_Info;%resp_ident;%remi_ident;%
% leg = {'Baseline';'2MO';'6MO';'12MO';'24MO'};
% group_Info_ses_change = group_Info(:,include_ses_change); group_Info_ses_change = group_Info_ses_change(:);
% group_uni = unique(group_Info_ses_change);
% leg = leg(include_ses_change);
% group_N = length(include_ses_change);



%% plot the correlation between roi activation and SCL-20 across timepoints -> mechanism markers
for i_task = 1
    for i_con = 1
        % load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
        beh_ibeh_change_raw = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_change = beh_ibeh_change_raw(:);
        
        %% plot fmri*group
        selected_rois = [1:3]
        for i_ROI = selected_rois
            c_f = figure('Position',[538   359   180   120],'Color','White'),
            task_activation_ROI_icon_itask_iroi_change_raw = prepare_data_change(squeeze(task_activation_ROI_icon_itask(:,i_ROI,:)), include_ses);
            task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change_raw(:);
            % overall corr
            [r, p] = partialcorr(beh_ibeh_change, task_activation_ROI_icon_itask_iroi_change, age_sex_tl,'Rows','pairwise','Type','Spearman')
            mdl = fitlm([task_activation_ROI_icon_itask_iroi_change], beh_ibeh_change); 
            % group 1
            [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 0), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 0), age_sex_tl(group_Info_ses_change == 0),'Rows','pairwise')
            % group 2
            [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 1), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 1), age_sex_tl(group_Info_ses_change == 1),'Rows','pairwise')

            x_label = ['\Delta Activation in ', ROI_plot_list{i_ROI}];
            y_label = ['\Delta ', beh_name_plot{i_beh}];
            fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
            [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, 1- group_Info_ses_change,[183, 90, 101; 216-50, 216-50, 222-50]/255.0,1,0.8,x_label,y_label,'',fliplr(leg), fig_p, '', 0);
            if i_beh == 24
                set(gca, 'YDir','reverse')
            end
            legend off
            OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
            mkdir(OutDir);
            savepath = [OutDir,'/model_2_2_var7_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];           
        end
    end
end



%% plot the prediction of 6, 12, and 24 mo behavior outcome using 2mo activation -> preditive markers
include_ses_change = setdiff(include_ses, [1]);
for i_task = 1
    for i_con = 1
        % load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
        beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_change = beh_ibeh_change(:);
        beh_ibeh_baseline = repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_baseline = beh_ibeh_baseline(:);
        timepoints_weight = (repmat(actual_time(include_ses_change), length(subjects_list), 1));timepoints_weight = timepoints_weight(:);
        beh_ibeh_change_weighted = beh_ibeh_change./timepoints_weight;
        % plot fmri
        for i_ROI = [1:length(ROI_list)]
            c_f = figure('Position',[538   359   180   120]),
            %% fmri*group
            task_activation_ROI_icon_itask_iroi_2mo = repmat(filloutliers(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 2) - task_activation_ROI_icon_itask(:,i_ROI, 1)), nan,'mean'), 1, length(include_ses_change));% 108 * length(include_ses_change)
            task_activation_ROI_icon_itask_iroi_2mo = task_activation_ROI_icon_itask_iroi_2mo(:);
            
            for i_group = include_ses_change
                [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_group), task_activation_ROI_icon_itask_iroi_2mo(group_Info_ses_change == i_group), age_sex_tl(group_Info_ses_change == i_group),'Rows','pairwise')
            end
            x_label = ['2MO - BL \DeltaNoGo > Go Activation in ', ROI_plot_list{i_ROI}];
            y_label = ['\Delta', beh_name_plot{i_beh}];
            fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
            [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_2mo, beh_ibeh_change, 1- group_Info_ses_change,[183, 90, 101; 216-50, 216-50, 222-50]/255.0,1,0.8,x_label,y_label,'',fliplr(leg), fig_p,'','','',0)
            legend off; 
            if i_beh == 24
                set(gca, 'YDir','reverse')
            end
            OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
            mkdir(OutDir);
%             savepath = [OutDir,'/model_6_3_var3_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')

        end
    end
end


%% plot linear prediction from R
data = readtable('/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/engage/SPM/Level2/MixedModel/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/SCL_20/MixedModel_SCL20/DataforR/linear_prediction.csv', 'PreserveVariableNames', true);
leg = {'6MO';'12MO';'24MO'};
s = rng('default'); rng(s); ColorMap=rand(400,3);
alpha = 0.5;

x_label = ['Predicted \DeltaSCL-20'];
y_label = ['Actual \DeltaSCL-20'];
fig_p = [2 2 8 7.5];% [2, 2, 6, 4]%
[s_f]= fmri_behavior_plot4(data.SCL_20_change, data.predicted_SCL_fmri, data.Time,'',1,'',x_label,y_label,'',leg, fig_p)
mkdir(OutDir);
savepath = [OutDir,'/ccn_linear_prediction.tiff'];
export_fig(savepath, '-r300', '-transparent')



function [data_change] = prepare_data_change(data, include_ses)
% data_fil = filloutliers(data,nan,'mean');
data_fil = data;
include_ses_change = setdiff(include_ses, [1]);
data_change = data_fil(:, include_ses_change) - repmat(data_fil(:, 1), 1, length(include_ses_change));% 108 * length(include_ses_change)
% data_change = data_change(:);
end


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
task_list = {'115_fMRI_stats_spikesonly_FD_fromFile_FACES-NONCONSCIOUS';'113_fMRI_stats_spikesonly_FD_fromFile_FACES-CONSCIOUS';...
    '111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO';};
task_contrasts_select = {{};{};{'Go','Nogo','NogovsGo','GovsNogo','Main Effect';}}
task_contrasts_select_con = {{};{};{'01','02','03','04','05';}};
%% load FD
load(fullfile('..','engage','FD_alltasks.mat'));% FD_alltasks
perc_motion = 0.25;
frame_N = 151; 
%% load data
load(fullfile('..','engage','wb_cog1.mat'));
load(fullfile('..','engage','wb_emo1.mat'));
load(fullfile('..','engage','wb_att1.mat'));
wb_data = [wn_cog, wn_emo, wn_att];
wb_data_names = [variable_names_cog;variable_names_emo;variables_att];
wb_data_names = cellfun(@(x) ['wb_', x], wb_data_names,'UniformOutput',false);
% fmri gonogo_data: task_beh_data,beh_name_itask
load(['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/engage/Behavior/gng/ENGAGE/gng_performance.mat']);
beh_name_itask = cellfun(@(x) ['fmri_', x], beh_name_itask,'UniformOutput',false);
% concatenate gonogo beh
fmri_gonogo_data = task_beh_data;
webneuro_gonogo_ind = [1:3];
webneuro_gonogo_data = wn_cog(:,webneuro_gonogo_ind,:);
% add FT
webneuro_gonogo_data = [webneuro_gonogo_data, webneuro_gonogo_data(:, 2,:) + webneuro_gonogo_data(:, 3,:)];
webneuro_gonogo_name = {'wb_nRT';'wb_nFA';'wb_nFM';'wb_nFT'};
gonogo_all_data = [fmri_gonogo_data,webneuro_gonogo_data];
gonogo_all_name = [beh_name_itask;webneuro_gonogo_name];
gonogo_all_baseline = gonogo_all_data(:,:,1);

% load intervention
load(fullfile('..','engage','Interv_Info.mat'));
load(fullfile('..','engage','Responder_list_6mo.mat'))
load(fullfile('..','engage','Remission_list_6mo.mat'))
% load beh_all and beh_name
% load(fullfile('..','engage','beh_all.mat')); 
% load(fullfile('..','engage','beh_all_06012021.mat'));
load(fullfile('..','engage','beh_all_18-Oct-2021.mat'));
if size(beh_all,2) == 39
    beh_all = beh_all(:,[1:18,26:27, 29, 28, 30:39,19:25],:);
    beh_name = beh_name([1:18,26:27, 29, 28, 30:39,19:25],1);
end
if size(beh_all,2) == 36
    beh_all = beh_all(:,[1:18,23:24, 26, 25, 27:36,19:22],:);
    beh_name = beh_name([1:18,23:24, 26, 25, 27:36,19:22],1);
end


age_sex = beh_all(:,19:20, 1);
beh_all = [beh_all, gonogo_all_data, wb_data];
beh_name = [beh_name;gonogo_all_name; wb_data_names];
%% define ROIs
% % ROI1
% % ROI_dir=['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Mask/engage/Cor_Beh/6MO-000/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/SCL_20_6MOvsBaseline'];
% % ROI_list = {'lAIn';'lIFG';'rDLPFC';'rTPJ';'ACC';'lIFG2';'rPrecenG';'lparaCin';'rFronPole'};
% % ROI_label = '6MOcorrROIs'
% % ROI 2
% % ROI_dir = ['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Mask/engage/Onesample/000/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/Activation'];
% % ROI_list = {'lAIn';'ldlpfc';'lMTG';'lSMG';'ACC';'rFP';'rdlpfc';'rMFG';'rIFG';'rAIn';'rSMG';'rMTG'};
% % ROI_label = ['BaselineActivationROIs'];

%% ROI 3 SCL-20
% ROI_dir = ['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Mask/engage/Cor_Beh/MixedModel/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/SCL_20'];
% % ROI_list = {'test'}
% % ROI_label = 'test'
% 
% ROI_label = ['MixedModel_SCL20'];
% ROI_list = {'lAIn';'lCaudate';'rIPG';'rIFG'; 'lIFG';'lPHG';'Precuneus'; 'rIPG3';'lSPL'};
% ROI_plot_list = {'lAIn';'lCaudate';'rIPL';'rDLPFC'; 'lDLPFC';'lPHG'; 'Precuneus';'rIPL';'lSPL'};


%% ROI 4 for BMI
% ROI_dir = ['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Mask/engage/Cor_Beh/MixedModel/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/BMI'];
% ROI_list = {'Overlay_model1_var_3_6_7_8'};
% ROI_label = ['MixedModel_BMI'];

%% ROI 5 for spsi
% 
ROI_dir = ['/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/Mask/engage/Cor_Beh/MixedModel/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/SPSI'];
ROI_label = ['MixedModel_SPSI'];
ROI_list = {'lDLPFC2';'rSMG';'Precuneus2';'lFusiform';'rIPL';'lIPL'};
ROI_plot_list = {'lDLPFC';'rIPL';'Precuneus';'lFusiform';'rIPL';'lIPL'};

%% prepare fmri data
radius = 2;
ROI_list = cellfun(@(x) [x, '_',num2str(radius)], ROI_list, 'UniformOutput', false);
DATADIR = fullfile('..','engage');
% for i_task = 3%1: length(task_list)
%     for i_con = 3%1: length(task_contrasts_select{i_task})
%         task_activation_ROI_icon_itask = nan(length(subjects_list), length(ROI_list), length(session_list));
%         for i_ses = 1: length(session_list)
%             scan = {};
%             for sub = 1:length(subjects_list)
%                 subject = subjects_list{sub};
%                 data = [];
%                 if FD_alltasks(sub, i_ses, i_task) <= frame_N * perc_motion
%                     data = [DATADIR,'/',task_list{i_task},'/', session_list{i_ses},'/',subject,'/con_00', task_contrasts_select_con{i_task}{i_con},'_mni.nii'];
%                 end
%                 scan{sub,1}=data;
%             end
%             for i_ROI = 1: length(ROI_list)
%                 ROI_i_dir = fullfile(ROI_dir, [ROI_list{i_ROI},'.nii']);
%                 task_activation_ROI_icon_itask(:,i_ROI, i_ses) =ExtValMul3D(scan,ROI_i_dir);
%             end
%         end
%         mkdir(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}))
%         save(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']), 'task_activation_ROI_icon_itask', 'ROI_list', 'ROI_dir'); 
%     end
% end


%% mixed modeling between roi activation and SCL-20 across timepoints
include_ses = [1,3:5];
actual_time = [0, 2, 6, 12, 24];
covariates = repmat(age_sex,[1 1 5]);
covariates_name = {'age';'gender'};



% for i_beh = [40:43]%[1:18, 24:25,40:79]
% for i_task = 3%1: length(task_list)
%     for i_con = 3%1: length(task_contrasts_select{i_task})
%         
%         OutDir = fullfile('..','engage','SPM','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%         mkdir(OutDir)
%         % load data
%         load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
%         
%         %% baseline corr between fmri and beh
% %         data_test = [task_activation_ROI_icon_itask(:,:,1), beh_all(:,[24:25,40:43,48:71],1)];
% %         data_name = [ROI_list; beh_name([24:25,40:43,48:71])]
% %         xylabel = strrep(data_name,'_','-');
% %         % total subjects
% %         [r, p] = corrcoef(data_test, 'Rows','pairwise');
% %         r(p >= 0.05) = 0;
% %         figure, subplot(2,3,1), imagesc(r), axis xy, colormap('jet')
% %         set(gca, 'XTick', 1: length(xylabel), 'XTickLabel', xylabel), xtickangle(45)
% %         set(gca, 'YTick', 1: length(xylabel), 'YTickLabel', xylabel)
% %         xlim([1,length(ROI_list)]);
% %         title({'All subjects baseline correlation between WebNeuo and treatment outcomes'; ['N = ', num2str(length(Interv_Info))]});
% %         
% %         % I-CARE
% %         [r, p] = corrcoef(data_test(Interv_Info == 1,:), 'Rows','pairwise');
% %         r(p >= 0.05) = 0;
% %         subplot(2,3,2), imagesc(r), axis xy, colormap('jet'),
% %         set(gca, 'XTick', 1: length(xylabel), 'XTickLabel', xylabel), xtickangle(45)
% %         set(gca, 'YTick', 1: length(xylabel), 'YTickLabel', xylabel),
% %         xlim([1,length(ROI_list)]);
% %         title({'I-CARE baseline correlation between WebNeuo and treatment outcomes'; ['N = ', num2str(length(find(Interv_Info == 1)))]});
% %         
% %         % U-CARE
% %         [r, p] = corrcoef(data_test(Interv_Info == 0,:), 'Rows','pairwise');
% %         r(p >= 0.05) = 0;
% %         subplot(2,3,3), imagesc(r), axis xy, colormap('jet'),
% %         set(gca, 'XTick', 1: length(xylabel), 'XTickLabel', xylabel), xtickangle(45)
% %         set(gca, 'YTick', 1: length(xylabel), 'YTickLabel', xylabel),
% %         xlim([1,length(ROI_list)]);
% %         title({'U-CARE baseline correlation between WebNeuo and treatment outcomes'; ['N = ', num2str(length(find(Interv_Info == 0)))]});
% 
%         %% correlation between change at 6 months for wb, fmri performance and treatment outcomes: SPSI and SCL-20
% %         data_test = [task_activation_ROI_icon_itask(:,:,3) - task_activation_ROI_icon_itask(:,:,1), beh_all(:,[24:25,40:43,48:71],3) - beh_all(:,[24:25,40:43,48:71],1)];
% %         data_name = [ROI_list; beh_name([24:25,40:43,48:71])]
% %         xylabel = strrep(data_name,'_','-');
% %         % total subjects
% %         [r, p] = corrcoef(data_test, 'Rows','pairwise');
% %         r(p >= 0.05) = 0;
% %         subplot(2,3,4), imagesc(r), axis xy, colormap('jet'),
% %         set(gca, 'XTick', 1: length(xylabel), 'XTickLabel', xylabel), xtickangle(45)
% %         set(gca, 'YTick', 1: length(xylabel), 'YTickLabel', xylabel)
% %         xlim([1,length(ROI_list)]);
% %         title({'All subjects correlation between 6MO change of WebNeuo and treatment outcomes'; ['N = ', num2str(length(Interv_Info))]});
% %         
% %         % I-CARE
% %         [r, p] = corrcoef(data_test(Interv_Info == 1,:), 'Rows','pairwise');
% %         r(p >= 0.05) = 0;
% %         subplot(2,3,5), imagesc(r), axis xy, colormap('jet'),
% %         set(gca, 'XTick', 1: length(xylabel), 'XTickLabel', xylabel), xtickangle(45)
% %         set(gca, 'YTick', 1: length(xylabel), 'YTickLabel', xylabel),
% %         xlim([1,length(ROI_list)]);
% %         title({'I-CARE correlation between 6MO change of WebNeuo and treatment outcomes'; ['N = ', num2str(length(find(Interv_Info == 1)))]});
% %         
% %         % U-CARE
% %         [r, p] = corrcoef(data_test(Interv_Info == 0,:), 'Rows','pairwise');
% %         r(p >= 0.05) = 0;
% %         subplot(2,3,6), imagesc(r), axis xy, colormap('jet'),
% %         set(gca, 'XTick', 1: length(xylabel), 'XTickLabel', xylabel), xtickangle(45)
% %         set(gca, 'YTick', 1: length(xylabel), 'YTickLabel', xylabel),
% %         xlim([1,length(ROI_list)]);
% %         title({'U-CARE correlation between 6MO change of WebNeuo and treatment outcomes'; ['N = ', num2str(length(find(Interv_Info == 0)))]});
% 
%         %% model2~3 - include changes at all timepoints
%         % prepare for model 2 & 3
%         include_ses_change = setdiff(include_ses, [1]);
%         timepoints = repmat(include_ses_change, length(subjects_list), 1); timepoints = timepoints(:); timepoints = nominal(timepoints);
%         timepoints1 = repmat(include_ses_change, length(subjects_list), 1); timepoints1 = timepoints1(:); timepoints1 = nominal(timepoints1);
%         timepoints_weight = (repmat(actual_time(include_ses_change), length(subjects_list), 1));timepoints_weight = timepoints_weight(:);
%         subjects = repmat([1: length(subjects_list)]', 1, length(include_ses_change)); subjects =  subjects(:); subjects = nominal(subjects);
%         FD_include_change = FD_alltasks(:, include_ses_change, i_task) - FD_alltasks(:, 1, i_task); FD_include_change = FD_include_change(:);
%         FD_include_baseline = repmat(squeeze(FD_alltasks(:, 1, i_task)), 1, length(include_ses_change)); FD_include_baseline = FD_include_baseline(:);
%         age_sex_tl = repmat(age_sex,length(include_ses_change),1);
%         group_info = Interv_Info; group_info = repmat(group_info, 1, length(include_ses_change)); group_info = nominal(group_info(:));
%         beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
%         beh_ibeh_change_percent = beh_ibeh_change./repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
%         beh_ibeh_change = beh_ibeh_change(:);
%         beh_ibeh_change_weighted = beh_ibeh_change./timepoints_weight;
%         beh_ibeh_change_percent = beh_ibeh_change_percent(:);
%         beh_ibeh_baseline = repmat(squeeze(beh_all(:, i_beh, 1)), 1, length(include_ses_change));
%         beh_ibeh_baseline = beh_ibeh_baseline(:);
%         
%         scl_change = squeeze(beh_all(:, 24, include_ses_change)) - repmat(beh_all(:, 24, 1), 1, length(include_ses_change));
%         scl_change = scl_change(:);
%         spsi_change = squeeze(beh_all(:, 25, include_ses_change)) - repmat(beh_all(:, 25, 1), 1, length(include_ses_change));
%         spsi_change = spsi_change(:);
%         
%         for i_ROI = 1: length(ROI_list)
% %             ROI_list{i_ROI}
%             % fmri change
%             task_activation_ROI_icon_itask_iroi_change = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses_change)) ...
%                 - repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
%             %             task_activation_ROI_icon_itask_iroi_change = filloutliers(task_activation_ROI_icon_itask_iroi_change, nan, 'mean');
%             task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change(:);
%             % baseline fmri
%             task_activation_ROI_icon_itask_iroi_baseline = repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
%             task_activation_ROI_icon_itask_iroi_baseline = task_activation_ROI_icon_itask_iroi_baseline(:);
%             % 2mo fmri
%             task_activation_ROI_icon_itask_iroi_2mo = repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 2) - task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
%             task_activation_ROI_icon_itask_iroi_2mo = task_activation_ROI_icon_itask_iroi_2mo(:);
%             % data -> table
%             tbl2 = table(beh_ibeh_change, beh_ibeh_baseline, timepoints,subjects, FD_include_baseline, task_activation_ROI_icon_itask_iroi_baseline, task_activation_ROI_icon_itask_iroi_2mo, task_activation_ROI_icon_itask_iroi_change, age_sex_tl(:,1), age_sex_tl(:,2), group_info,scl_change, spsi_change, ...
%                 'VariableNames',{beh_name{i_beh},[beh_name{i_beh},'_baseline'],'Time','Subj','fd', 'fMRI_baseline', 'fMRI_2MO', 'fMRI','age','sex', 'Group', 'SCL','SPSI'});
%             tbl2.sex = categorical(tbl2.sex); tbl2.Group = nominal(tbl2.Group);
%             %
%             tbl2_2 = table(beh_ibeh_change, beh_ibeh_baseline, timepoints1,subjects, FD_include_baseline, task_activation_ROI_icon_itask_iroi_baseline, task_activation_ROI_icon_itask_iroi_change, age_sex_tl(:,1), age_sex_tl(:,2), group_info,...
%                 'VariableNames',{beh_name{i_beh},[beh_name{i_beh},'_baseline'],'Time','Subj','fd', 'fMRI_baseline', 'fMRI','age','sex', 'Group'});
%             
%             % I-CARE vs. U-CARE behavior
% %             mixmodel_beh = fitlme(tbl2, sprintf('%s ~ %s + %s:%s + %s + (%s|%s)', ...
% %                 beh_name{i_beh},'Group','Group','Time', 'Time', '1',char('Subj')),...
% %                 'DummyVarCoding','effects','StartMethod','random');
% %             ss_beh = anova(mixmodel_beh); ss_beh_table = dataset2table(ss_beh)
%             
%             %% model2: change predicts change
%             mixmodel2 = fitlme(tbl2, sprintf('%s ~ %s + %s:%s + (%s|%s)', ...
%                 beh_name{i_beh},'fMRI', 'fMRI',  'Group', '1',char('Subj')),...
%                 'DummyVarCoding','effects','StartMethod','random');
%             ss2 = anova(mixmodel2); ss2_table = dataset2table(ss2);
%             if ss2.pValue(3)< 0.05 || ss2.pValue(3)< 0.05
%                 ROI_list{i_ROI}
%                 beh_name{i_beh}
%                 ss2
%                 mixmodel_SCL = fitlme(tbl2, sprintf('%s ~ %s + %s:%s + (%s|%s)', ...
%                     'SCL',beh_name{i_beh}, beh_name{i_beh},  'Group', '1',char('Subj')),...
%                     'DummyVarCoding','effects','StartMethod','random');
%                 ss_SCL = anova(mixmodel_SCL)
%                 
%                 mixmodel_SPSI = fitlme(tbl2, sprintf('%s ~ %s + %s:%s + (%s|%s)', ...
%                     'SPSI',beh_name{i_beh}, beh_name{i_beh},  'Group', '1',char('Subj')),...
%                     'DummyVarCoding','effects','StartMethod','random');
%                 ss_SPSI = anova(mixmodel_SPSI)
%             end
%             
% %             %% baseline predicts baseline
% %             
% %             mixmodel3_3 = fitlm(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)),beh_all(:, i_beh, 1),...
% %                 'DummyVarCoding','effects')
% %             ss3_3 = anova(mixmodel3_3);
%         end
%         
%     end
%     
% end
% end

%% Plot Results
include_ses_change = setdiff(include_ses, [1]);
age_sex_tl = repmat(age_sex,length(include_ses_change),1);

i_beh = 41;
beh_name_plot = beh_name;
beh_name_plot{24} = 'Depression symptoms';
beh_name_plot{25} = 'Problem-solving ability';
%% two group: intervention by response
% group_Info1 = Interv_Info;
% group_Info2 = resp_ident;%remi_ident;%
% group_Info = nan(size(group_Info1));
% group_Info(intersect(find(group_Info1 == 0), find(group_Info2 == 1))) = 1;
% group_Info(intersect(find(group_Info1 == 0), find(group_Info2 == 0))) = 2;
% group_Info(intersect(find(group_Info1 == 1), find(group_Info2 == 1))) = 3;
% group_Info(intersect(find(group_Info1 == 1), find(group_Info2 == 0))) = 4;
% leg = {['U-care: Responder'];['U-care: Non'];['I-CARE: Responder'];['I-CARE: Non']};
% group_Info_ses_change = repmat(group_Info, 1, length(include_ses_change)); group_Info_ses_change = group_Info_ses_change(:);
% group_N = 4;

%% two group: intervention by time
% group_Info1 = Interv_Info;
% leg = {['U-CARE'],['I-CARE'];};
% group_Info1_ses_change = repmat(group_Info1, 1, length(include_ses_change)); group_Info1_ses_change = group_Info1_ses_change(:);
% group_Info2 = [ones(length(subjects_list),1), 2*ones(length(subjects_list),1), 3*ones(length(subjects_list),1),4*ones(length(subjects_list),1),5*ones(length(subjects_list),1)];
% group_Info2_ses_change = group_Info2(:,include_ses_change); group_Info2_ses_change = group_Info2_ses_change(:);
% group_Info_ses_change = nan(size(group_Info1_ses_change));
% i_group= 0;
% for i_group1 = 0:1
%     for i_group2 = include_ses_change
%         i_group = i_group +1;
%         group_Info_ses_change(intersect(find(group_Info1_ses_change == i_group1), find(group_Info2_ses_change == i_group2))) = i_group;
%     end
% end
% group_N = length(include_ses_change)*2;

%% one group: time invariant
group_Info = Interv_Info;%resp_ident;%remi_ident;%
leg = {['U-CARE'],['I-CARE'];};
group_Info_ses_change = repmat(group_Info, 1, length(include_ses_change)); group_Info_ses_change = group_Info_ses_change(:);
group_Info_ses_change_name = cell(length(group_Info_ses_change),1);
group_Info_ses_change_name(group_Info_ses_change == 0) = {leg{1}};
group_Info_ses_change_name(group_Info_ses_change == 1) = {leg{2}};
group_N = 2;
%% one grouop: time variant
% group_Info = [ones(length(subjects_list),1), 2*ones(length(subjects_list),1), 3*ones(length(subjects_list),1),4*ones(length(subjects_list),1),5*ones(length(subjects_list),1)];%Interv_Info;%resp_ident;%remi_ident;%
% leg = {'Baseline';'2MO';'6MO';'12MO';'24MO'};
% group_Info_ses_change = group_Info(:,include_ses_change); group_Info_ses_change = group_Info_ses_change(:);
% group_uni = unique(group_Info_ses_change);
% leg = leg(include_ses_change);
% group_N = length(include_ses_change);

%% plot the correlation between roi activation and SCL-20 across timepoints -> 2mo - baseline
% include_ses_change = setdiff(include_ses, [1]);
% for i_task = 3%1: length(task_list)
%     for i_con = 3%1: length(task_contrasts_select{i_task})
%         % load data
%         load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
%         beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
%         beh_ibeh_change = beh_ibeh_change(:);
%         beh_ibeh_baseline = repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
%         beh_ibeh_baseline = beh_ibeh_baseline(:);
%         timepoints_weight = (repmat(actual_time(include_ses_change), length(subjects_list), 1));timepoints_weight = timepoints_weight(:);
%         beh_ibeh_change_weighted = beh_ibeh_change./timepoints_weight;
%         % plot fmri
%         for i_ROI = [1:length(ROI_list)]
%             c_f = figure('Position',[538   359   180   120]),
%             x
%             %% fmri main effect, only one color
%             task_activation_ROI_icon_itask_iroi_2mo = repmat(filloutliers(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 2) - task_activation_ROI_icon_itask(:,i_ROI, 1)), nan,'mean'), 1, length(include_ses_change));% 108 * length(include_ses_change)
%             task_activation_ROI_icon_itask_iroi_2mo = task_activation_ROI_icon_itask_iroi_2mo(:);
%             
%             for i_group = include_ses_change
%                 [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_group), task_activation_ROI_icon_itask_iroi_2mo(group_Info_ses_change == i_group), age_sex_tl(group_Info_ses_change == i_group),'Rows','pairwise')
%             end
%             x_label = ['2MO - BL \DeltaNoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%             y_label = ['\Delta', beh_name_plot{i_beh}];
%             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%             [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_2mo, beh_ibeh_change, ones(length(beh_ibeh_change), 1),[0 0 0],1,0.5,x_label,y_label,'','', fig_p,'','','',0)
%             legend off; 
%             if i_beh == 24
%                 set(gca, 'YDir','reverse')
%             end
%             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%             mkdir(OutDir);
% %             savepath = [OutDir,'/model_6_3_var3_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
% %             export_fig(savepath, '-r300', '-transparent')
%             
%             %% fmri by time, all timepoints included and color coded
% %             task_activation_ROI_icon_itask_iroi_2mo = repmat(filloutliers(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 2) - task_activation_ROI_icon_itask(:,i_ROI, 1)), nan,'mean'), 1, length(include_ses_change));% 108 * length(include_ses_change)
% %             task_activation_ROI_icon_itask_iroi_2mo = task_activation_ROI_icon_itask_iroi_2mo(:);
% %             
% %             for i_group = include_ses_change
% %                 [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_group), task_activation_ROI_icon_itask_iroi_2mo(group_Info_ses_change == i_group), age_sex_tl(group_Info_ses_change == i_group),'Rows','pairwise')
% %             end
% %             x_label = ['2MO - BL \DeltaNoGo > Go Activation in ', ROI_plot_list{i_ROI}];
% %             y_label = ['\Delta', beh_name_plot{i_beh}];
% %             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
% %             [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_2mo, beh_ibeh_change, group_Info_ses_change,'',1,'',x_label,y_label,'',leg, fig_p)
% %             legend off
% %             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
% %             mkdir(OutDir);
% % %             savepath = [OutDir,'/model_6_3_var3_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
% % %             export_fig(savepath, '-r300', '-transparent')
%             
% 
%             %% overall prediction
% %             mixmodel6_3_sim = fitlm([task_activation_ROI_icon_itask_iroi_2mo, beh_ibeh_baseline, age_sex_tl], beh_ibeh_change);
% %             [r,p] = corrcoef(beh_ibeh_change, mixmodel6_3_sim.Fitted,'Rows','pairwise')
% %             x_label = ['Actual \Delta', beh_name_plot{i_beh}];
% %             y_label = ['Predicted \Delta', beh_name_plot{i_beh}];
% %             figure,
% %             fmri_behavior_plot4(beh_ibeh_change, mixmodel6_3_sim.Fitted, group_Info_ses_change,'',1,'',x_label,y_label,'',leg, fig_p);
% % %             ylim([-1.5,1])
%         end
%     end
% end

%% plot the correlation between roi activation and SCL-20 across timepoints -> △

for i_task = 3%1: length(task_list)
    for i_con = 3%1: length(task_contrasts_select{i_task})
        % load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
        beh_ibeh_change_raw = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_change = beh_ibeh_change_raw(:);
        % plot for model 1
        s = rng('default'); rng(s); ColorMap=rand(400,3);
        alpha = 0.5;
        fp_dynamic4_1group = [538, 359, 270, 180];
        %% plot fmri main effect, one color for all subjects
%         for i_ROI = 1:length(ROI_list)
%             ROI_plot_list{i_ROI}
%             c_f = figure('Position',fp_dynamic4_1group);
%             task_activation_ROI_icon_itask_iroi_change_raw = prepare_data_change(squeeze(task_activation_ROI_icon_itask(:,i_ROI,:)), include_ses);
%             task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change_raw(:);
%             [r, p] = partialcorr(beh_ibeh_change, task_activation_ROI_icon_itask_iroi_change, age_sex_tl,'Rows','pairwise')
%             for i_ses = 3:5
%                 [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_ses), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == i_ses), age_sex_tl(group_Info_ses_change == i_ses),'Rows','pairwise')
%             end
% % %             % group 1
% % %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 0), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 0), age_sex_tl(group_Info_ses_change == 0),'Rows','pairwise')
% % %             % group 2
% % %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 1), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 1), age_sex_tl(group_Info_ses_change == 1),'Rows','pairwise')
% 
%             x_label = ['\DeltaNoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%             y_label = ['\Delta', beh_name_plot{i_beh}];
%             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%             [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, ones(length(beh_ibeh_change), 1),[0 0 0],1,0.5,x_label,y_label,'','', fig_p,'',1,'',0)
% %             ylim([-2.5 2]);
%             if i_beh == 24
%                 set(gca, 'YDir','reverse')
%             end
%             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%             mkdir(OutDir);
%             savepath = [OutDir,'/model_2_2_var5_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')
%         end
        
        
        %% plot fmri*time, all timepoints inlcuded and color coded
%         for i_ROI = [2]%[1:2,19:20]%[3:6]%length(ROI_list)
%             ROI_plot_list{i_ROI}
%             c_f = figure('Position',fp_dynamic4_1group);
%             task_activation_ROI_icon_itask_iroi_change_raw = prepare_data_change(squeeze(task_activation_ROI_icon_itask(:,i_ROI,:)), include_ses);
%             task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change_raw(:);
%             [r, p] = partialcorr(beh_ibeh_change, task_activation_ROI_icon_itask_iroi_change, age_sex_tl,'Rows','pairwise')
%             for i_ses = 3:5
%                 [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_ses), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == i_ses), age_sex_tl(group_Info_ses_change == i_ses),'Rows','pairwise')
%             end
% % %             % group 1
% % %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 0), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 0), age_sex_tl(group_Info_ses_change == 0),'Rows','pairwise')
% % %             % group 2
% % %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 1), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 1), age_sex_tl(group_Info_ses_change == 1),'Rows','pairwise')
% 
%             x_label = ['\DeltaNoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%             y_label = ['\Delta', beh_name_plot{i_beh}];
%             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%             [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, group_Info_ses_change,'',1,'',x_label,y_label,'',leg, fig_p)
%             ylim([-2.5 2])
%             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%             mkdir(OutDir);
%             savepath = [OutDir,'/model_2_2_var6_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')
%         end
        %% plot fmri*time, individual plot for each timepoint
%         selected_rois = [6:9]%[3:6]; %
%         for i_ROI = selected_rois
%                 c_f = figure('Position',[538   359   270*3   180],'Color','White'),
% %                 suptitle(ROI_list{i_ROI})
%                 for i_sess_change = 1: length(include_ses_change)
%                     include_ses_change(i_sess_change)
%                     beh_each = squeeze(beh_all(:, i_beh, include_ses_change(i_sess_change))) - beh_all(:, i_beh, 1);
%                     task_activation_ROI_each = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses_change(i_sess_change))) ...
%                     - squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1));
%                     task_activation_ROI_each = filloutliers(task_activation_ROI_each,nan,'mean');
%                     [r, p] = partialcorr(beh_each, task_activation_ROI_each, age_sex,'Rows','pairwise')
%                     mdl = fitlm([task_activation_ROI_each], beh_each);
%                     figure(c_f), subplot(1, length(include_ses_change), i_sess_change)
% %                     plot(task_activation_ROI_each, beh_each,'o');
%                     hold on; h = plot(mdl,'Color',cl(i_ROI-selected_rois(1)+1,:), 'MarkerSize',2,'Marker','o','MarkerEdgeColor', cl(i_ROI-selected_rois(1)+1,:),'MarkerFaceColor', cl(i_ROI-selected_rois(1)+1,:)); hold off
%                     for i_h = 2:4
%                         set(h(i_h),'Color', cl(i_ROI-selected_rois(1)+1,:));
%                     end
%                     title(session_list{include_ses_change(i_sess_change)})% newline 'r=', num2str(r),'; p=',num2str(p)]);
%                     y_label = strrep(beh_name{i_beh}, '_','\_');
%                     xlabel(['\Delta','Nogo > Go'],'Interpreter','tex'); ylabel(['\Delta',y_label],'Interpreter','tex'); 
% %                     xlim([-1.8 1.8]); ylim([-2 1.1])
%                     box on
%                     legend off
%                 end
%                 
%         end
        %% plot fmri*group
        
        
        selected_rois = 1:length(ROI_list)
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
            fig_p = [2 2 8 6.7];% 
            [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, 1- group_Info_ses_change,[183, 90, 101; 216-50, 216-50, 222-50]/255.0,1,0.8,x_label,y_label,'',fliplr(leg), fig_p, '', 0);
            legend off
            if i_beh == 41
                set(gca, 'YDir','reverse')
            end
            OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
            mkdir(OutDir);
            savepath = [OutDir,'/model_2_2_var7_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')
            
            
%             h = gscatter(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, group_Info_ses_change, cl(1: group_N,:),'.',[10])
%             lsline
%             if i_ROI == selected_rois(1)
%                 legend([h(1),h(2)],leg);
%             else
%                 legend off
%             end
% %             hold on; plot(mdl,'Color',[0 0 0], 'Marker','o','MarkerSize',2,'MarkerEdgeColor', [0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5]); hold off
%             title([ROI_list{i_ROI}],'Interpreter','none');% newline 'r=', num2str(r),'; p=',num2str(p)
%             y_label = strrep(beh_name{i_beh}, '_','\_');
%             xlabel(['\Delta','Nogo > Go'],'Interpreter','tex'); ylabel(['\Delta', y_label],'Interpreter','tex');
%             %                     xlim([-1.8 1.8]); ylim([-2 1.1])
            
            
        end

                  %% plot fmri*group*time
%         for i_ROI = [3,5]%1:6
%             ROI_plot_list{i_ROI}
%             i_plot = 0;
%             
%             for i_sess_change = 1: length(include_ses_change)
%                 i_sess_change
%                 figure('Position',[538   359   180*3*1.5   120*2],'Color','White');
%                 beh_each = squeeze(beh_all(:, i_beh, include_ses_change(i_sess_change))) - beh_all(:, i_beh, 1);
%                 task_activation_ROI_icon_itask_iroi_change_raw = squeeze(task_activation_ROI_icon_itask(:,i_ROI,include_ses_change(i_sess_change))) - squeeze(task_activation_ROI_icon_itask(:,i_ROI,1));
%                 [r, p] = corrcoef(task_activation_ROI_icon_itask_iroi_change_raw, beh_each,'Rows','pairwise');
%                 'U-CARE'
%                 [r, p] = corrcoef(task_activation_ROI_icon_itask_iroi_change_raw(group_Info1 == 0), beh_each(group_Info1 == 0),'Rows','pairwise')
%                 'I_care'
%                 [r, p] = corrcoef(task_activation_ROI_icon_itask_iroi_change_raw(group_Info1 == 1), beh_each(group_Info1 == 1),'Rows','pairwise')
%                 x_label = ['\Delta NoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%                 y_label = ['\Delta', beh_name_plot{i_beh}];
%                 fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%                 subplot(1,3,i_sess_change)
%                 hold on
%                 [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change_raw, beh_each, 1- group_Info1,[cl(i_sess_change,:)/2;cl(i_sess_change,:)],1,'',x_label,y_label,'',fliplr(leg), fig_p);
%                 OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%                 mkdir(OutDir);
%                 
%             end
%             axes_p = [0.2    0.2    0.75    0.75];
%             fontsize = 10;
%             fontname = 'Helvetica';
%             set(get(gca,'title'),'FontSize',fontsize,'FontName',fontname);%设置标题字体大小，字型
%             set(get(gca,'XLabel'),'FontSize',fontsize,'FontName',fontname);%设置X坐标标题字体大小，字型
%             set(get(gca,'YLabel'),'FontSize',fontsize,'FontName',fontname);%设置Y坐标标题字体大小，字型
%             set(gca,'FontName',fontname,'FontSize',fontsize)%设置坐标轴字体大小，字型
%             set(get(gca,'Legend'), 'FontSize',fontsize,'FontName',fontname);
% %             savepath = [OutDir,'/model_3_2_var8_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
% %             export_fig(savepath, '-r300', '-transparent')
%         end  
        
    end
end



%% plot the correlation between behavior and SCL-20/SPSI across timepoints -> △
% outcome
i_beh = 41% bmi: 22; scl-20: 24 spsi 25; fmri RT: 26
for i_task = 3%1: length(task_list)
    for i_con = 3%1: length(task_contrasts_select{i_task})
        % load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
        beh_ibeh_change_raw = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_change = beh_ibeh_change_raw(:);
        % plot for model 1
        s = rng('default'); rng(s); ColorMap=rand(400,3);
        alpha = 0.5;
        %% plot fmri main effect, one color for all subjects
        for i_ROI = 1:length(ROI_list)
            ROI_plot_list{i_ROI}
            c_f = figure;
            task_activation_ROI_icon_itask_iroi_change_raw = prepare_data_change(squeeze(task_activation_ROI_icon_itask(:,i_ROI,:)), include_ses);
            task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change_raw(:);
            [r, p] = partialcorr(beh_ibeh_change, task_activation_ROI_icon_itask_iroi_change, age_sex_tl,'Rows','pairwise')
            for i_ses = 3:5
                [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_ses), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == i_ses), age_sex_tl(group_Info_ses_change == i_ses),'Rows','pairwise')
            end
% %             % group 1
% %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 0), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 0), age_sex_tl(group_Info_ses_change == 0),'Rows','pairwise')
% %             % group 2
% %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 1), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 1), age_sex_tl(group_Info_ses_change == 1),'Rows','pairwise')

            x_label = ['\DeltaNoGo > Go Activation in ', ROI_plot_list{i_ROI}];
            y_label = ['\Delta', beh_name_plot{i_beh}];
            fig_p = [2 2 7.5 6.35];% [2, 2, 6, 4]%
            [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, ones(length(beh_ibeh_change), 1),[0 0 0],1,0.5,x_label,y_label,'','', fig_p,'',0,'',0)
%             ylim([-2.5 2]);
            if i_beh == 24 || i_beh == 41
                set(gca, 'YDir','reverse')
            end
            OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
            mkdir(OutDir);
            savepath = [OutDir,'/model_2_2_var5_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
            export_fig(savepath, '-r300', '-transparent')
        end
        
        
        %% plot fmri*time, all timepoints inlcuded and color coded
%         for i_ROI = [2]%[1:2,19:20]%[3:6]%length(ROI_list)
%             ROI_plot_list{i_ROI}
%             c_f = figure('Position',fp_dynamic4_1group);
%             task_activation_ROI_icon_itask_iroi_change_raw = prepare_data_change(squeeze(task_activation_ROI_icon_itask(:,i_ROI,:)), include_ses);
%             task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change_raw(:);
%             [r, p] = partialcorr(beh_ibeh_change, task_activation_ROI_icon_itask_iroi_change, age_sex_tl,'Rows','pairwise')
%             for i_ses = 3:5
%                 [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_ses), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == i_ses), age_sex_tl(group_Info_ses_change == i_ses),'Rows','pairwise')
%             end
% % %             % group 1
% % %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 0), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 0), age_sex_tl(group_Info_ses_change == 0),'Rows','pairwise')
% % %             % group 2
% % %             [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == 1), task_activation_ROI_icon_itask_iroi_change(group_Info_ses_change == 1), age_sex_tl(group_Info_ses_change == 1),'Rows','pairwise')
% 
%             x_label = ['\DeltaNoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%             y_label = ['\Delta', beh_name_plot{i_beh}];
%             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%             [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, group_Info_ses_change,'',1,'',x_label,y_label,'',leg, fig_p)
%             ylim([-2.5 2])
%             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%             mkdir(OutDir);
%             savepath = [OutDir,'/model_2_2_var6_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')
%         end
        %% plot fmri*time, individual plot for each timepoint
%         selected_rois = [6:9]%[3:6]; %
%         for i_ROI = selected_rois
%                 c_f = figure('Position',[538   359   270*3   180],'Color','White'),
% %                 suptitle(ROI_list{i_ROI})
%                 for i_sess_change = 1: length(include_ses_change)
%                     include_ses_change(i_sess_change)
%                     beh_each = squeeze(beh_all(:, i_beh, include_ses_change(i_sess_change))) - beh_all(:, i_beh, 1);
%                     task_activation_ROI_each = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses_change(i_sess_change))) ...
%                     - squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1));
%                     task_activation_ROI_each = filloutliers(task_activation_ROI_each,nan,'mean');
%                     [r, p] = partialcorr(beh_each, task_activation_ROI_each, age_sex,'Rows','pairwise')
%                     mdl = fitlm([task_activation_ROI_each], beh_each);
%                     figure(c_f), subplot(1, length(include_ses_change), i_sess_change)
% %                     plot(task_activation_ROI_each, beh_each,'o');
%                     hold on; h = plot(mdl,'Color',cl(i_ROI-selected_rois(1)+1,:), 'MarkerSize',2,'Marker','o','MarkerEdgeColor', cl(i_ROI-selected_rois(1)+1,:),'MarkerFaceColor', cl(i_ROI-selected_rois(1)+1,:)); hold off
%                     for i_h = 2:4
%                         set(h(i_h),'Color', cl(i_ROI-selected_rois(1)+1,:));
%                     end
%                     title(session_list{include_ses_change(i_sess_change)})% newline 'r=', num2str(r),'; p=',num2str(p)]);
%                     y_label = strrep(beh_name{i_beh}, '_','\_');
%                     xlabel(['\Delta','Nogo > Go'],'Interpreter','tex'); ylabel(['\Delta',y_label],'Interpreter','tex'); 
% %                     xlim([-1.8 1.8]); ylim([-2 1.1])
%                     box on
%                     legend off
%                 end
%                 
%         end
        %% plot fmri*group
        
        
        selected_rois = 1:length(ROI_list)
        for i_ROI = selected_rois
            c_f = figure('Color','White'),
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
            fig_p = [2 2 7.5 6.35];% [2, 2, 6, 4]%
            [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, 1- group_Info_ses_change,[0, 0, 0; 0.7, 0.7, 0.7],1,0.5,x_label,y_label,'',fliplr(leg), fig_p, '', 0);
            legend off
            OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
            mkdir(OutDir);
            savepath = [OutDir,'/model_2_2_var7_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')
            
            
%             h = gscatter(task_activation_ROI_icon_itask_iroi_change, beh_ibeh_change, group_Info_ses_change, cl(1: group_N,:),'.',[10])
%             lsline
%             if i_ROI == selected_rois(1)
%                 legend([h(1),h(2)],leg);
%             else
%                 legend off
%             end
% %             hold on; plot(mdl,'Color',[0 0 0], 'Marker','o','MarkerSize',2,'MarkerEdgeColor', [0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5]); hold off
%             title([ROI_list{i_ROI}],'Interpreter','none');% newline 'r=', num2str(r),'; p=',num2str(p)
%             y_label = strrep(beh_name{i_beh}, '_','\_');
%             xlabel(['\Delta','Nogo > Go'],'Interpreter','tex'); ylabel(['\Delta', y_label],'Interpreter','tex');
%             %                     xlim([-1.8 1.8]); ylim([-2 1.1])
            
            
        end

                  %% plot fmri*group*time
%         for i_ROI = [3,5]%1:6
%             ROI_plot_list{i_ROI}
%             i_plot = 0;
%             
%             for i_sess_change = 1: length(include_ses_change)
%                 i_sess_change
%                 figure('Position',[538   359   180*3*1.5   120*2],'Color','White');
%                 beh_each = squeeze(beh_all(:, i_beh, include_ses_change(i_sess_change))) - beh_all(:, i_beh, 1);
%                 task_activation_ROI_icon_itask_iroi_change_raw = squeeze(task_activation_ROI_icon_itask(:,i_ROI,include_ses_change(i_sess_change))) - squeeze(task_activation_ROI_icon_itask(:,i_ROI,1));
%                 [r, p] = corrcoef(task_activation_ROI_icon_itask_iroi_change_raw, beh_each,'Rows','pairwise');
%                 'U-CARE'
%                 [r, p] = corrcoef(task_activation_ROI_icon_itask_iroi_change_raw(group_Info1 == 0), beh_each(group_Info1 == 0),'Rows','pairwise')
%                 'I_care'
%                 [r, p] = corrcoef(task_activation_ROI_icon_itask_iroi_change_raw(group_Info1 == 1), beh_each(group_Info1 == 1),'Rows','pairwise')
%                 x_label = ['\Delta NoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%                 y_label = ['\Delta', beh_name_plot{i_beh}];
%                 fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%                 subplot(1,3,i_sess_change)
%                 hold on
%                 [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_change_raw, beh_each, 1- group_Info1,[cl(i_sess_change,:)/2;cl(i_sess_change,:)],1,'',x_label,y_label,'',fliplr(leg), fig_p);
%                 OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%                 mkdir(OutDir);
%                 
%             end
%             axes_p = [0.2    0.2    0.75    0.75];
%             fontsize = 10;
%             fontname = 'Helvetica';
%             set(get(gca,'title'),'FontSize',fontsize,'FontName',fontname);%设置标题字体大小，字型
%             set(get(gca,'XLabel'),'FontSize',fontsize,'FontName',fontname);%设置X坐标标题字体大小，字型
%             set(get(gca,'YLabel'),'FontSize',fontsize,'FontName',fontname);%设置Y坐标标题字体大小，字型
%             set(gca,'FontName',fontname,'FontSize',fontsize)%设置坐标轴字体大小，字型
%             set(get(gca,'Legend'), 'FontSize',fontsize,'FontName',fontname);
% %             savepath = [OutDir,'/model_3_2_var8_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
% %             export_fig(savepath, '-r300', '-transparent')
%         end  
        
    end
end








%% plot the correlation between roi activation and SCL-20 across timepoints -> baseline
include_ses_change = setdiff(include_ses, [1]);
% for i_task = 3%1: length(task_list)
%     for i_con = 3%1: length(task_contrasts_select{i_task})
%         % load data
%         load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
%         beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
%         beh_ibeh_change = beh_ibeh_change(:);
%         timepoints_weight = (repmat(actual_time(include_ses_change), length(subjects_list), 1));timepoints_weight = timepoints_weight(:);
%         beh_ibeh_change_weighted = beh_ibeh_change./timepoints_weight;
%         % plot fmri
%         for i_ROI = 1: length(ROI_list)
%             ROI_plot_list{i_ROI}
%             c_f = figure('Position',[538   359   180   120]),
%             task_activation_ROI_icon_itask_iroi_baseline = filloutliers(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)),nan,'mean');
%             task_activation_ROI_icon_itask_iroi_baseline = repmat(task_activation_ROI_icon_itask_iroi_baseline, 1, length(include_ses_change));% 108 * length(include_ses_change)
%             task_activation_ROI_icon_itask_iroi_baseline = task_activation_ROI_icon_itask_iroi_baseline(:);
%             
%             for i_group = include_ses_change
%                 [r, p] = partialcorr(beh_ibeh_change(group_Info_ses_change == i_group), task_activation_ROI_icon_itask_iroi_baseline(group_Info_ses_change == i_group), age_sex_tl(group_Info_ses_change == i_group),'Rows','pairwise')
%             end
%             % overall correlation 
%             'Overall'
%             [r, p] = partialcorr(beh_ibeh_change, task_activation_ROI_icon_itask_iroi_baseline, age_sex_tl,'Rows','pairwise')
%             x_label = ['Baseline NoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%             y_label = ['\Delta', beh_name_plot{i_beh}];
%             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%             [s_f]= fmri_behavior_plot4(task_activation_ROI_icon_itask_iroi_baseline, beh_ibeh_change, group_Info_ses_change,'',1,'',x_label,y_label,'',leg, fig_p)
%             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%             mkdir(OutDir);
% %             savepath = [OutDir,'/model_3_2_var6_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
% %             export_fig(savepath, '-r300', '-transparent')
%             
%         end
%         %% plot fmri*time
% %         for i_ROI =11: 12
% %                 c_f = figure('Position',[538   359   180*3   120],'Color','White'),
% % %                 suptitle(ROI_list{i_ROI})
% %                 for i_sess_change = 1: length(include_ses_change)
% %                     include_ses_change(i_sess_change)
% %                     beh_each = squeeze(beh_all(:, i_beh, include_ses_change(i_sess_change))) - beh_all(:, i_beh, 1);
% %                     task_activation_ROI_each = squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1));
% %                     task_activation_ROI_each = filloutliers(task_activation_ROI_each,nan,'mean');
% % %                     [r, p] = partialcorr(beh_each, task_activation_ROI_each, age_sex,'Rows','pairwise')
% %                     mdl = fitlm([task_activation_ROI_each], beh_each);
% %                     figure(c_f), subplot(1, length(include_ses_change), i_sess_change)
% % %                     plot(task_activation_ROI_each, beh_each,'o');
% %                     hold on; h = plot(mdl,'Color',cl(i_ROI-10,:), 'MarkerSize',2,'Marker','o','MarkerEdgeColor', cl(i_ROI-10,:),'MarkerFaceColor', cl(i_ROI-10,:)); hold off
% %                     for i_h = 2:4
% %                         set(h(i_h),'Color', cl(i_ROI-10,:));
% %                     end
% %                     title(session_list{include_ses_change(i_sess_change)})% newline 'r=', num2str(r),'; p=',num2str(p)]);
% %                     y_label = strrep(beh_name{i_beh}, '_','\_');
% %                     xlabel(['Baseline','Nogo > Go Activation in '],'Interpreter','tex'); ylabel(['\Delta',y_label],'Interpreter','tex'); 
% % %                     xlim([-1.8 1.8]); ylim([-2 1.1])
% %                     box on
% %                     legend off
% %                 end
% %                 
% %         end
%         %% plot fmri*group*time
%         % newer
%         for i_ROI = 14:18
%             c_f = figure('Position',[538   359   180*3*1.5   120*2],'Color','White'),
%             i_plot = 0;
%             
%             for i_sess_change = 1: length(include_ses_change)
%                 task_activation_ROI_each = squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1));
%                 task_activation_ROI_each = filloutliers(task_activation_ROI_each,nan,'mean');
%                 beh_each = squeeze(beh_all(:, i_beh, include_ses_change(i_sess_change))) - beh_all(:, i_beh, 1);
%                 
%                 
%                 x_label = ['Baseline NoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%             y_label = ['\Delta', beh_name_plot{i_beh}];
%             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%             subplot(1,3,i_sess_change)
%             hold on
%             [s_f]= fmri_behavior_plot4(task_activation_ROI_each, beh_each, 1- group_Info1,[cl(i_sess_change,:)/2;cl(i_sess_change,:)],1,'',x_label,y_label,'',fliplr(leg), fig_p)
%             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%             mkdir(OutDir);
%             
%             end
%             axes_p = [0.2    0.2    0.75    0.75];
%             fontsize = 10;
%             fontname = 'Helvetica';
%             set(get(gca,'title'),'FontSize',fontsize,'FontName',fontname);%设置标题字体大小，字型
%             set(get(gca,'XLabel'),'FontSize',fontsize,'FontName',fontname);%设置X坐标标题字体大小，字型
%             set(get(gca,'YLabel'),'FontSize',fontsize,'FontName',fontname);%设置Y坐标标题字体大小，字型
%             set(gca,'FontName',fontname,'FontSize',fontsize)%设置坐标轴字体大小，字型
%             set(get(gca,'Legend'), 'FontSize',fontsize,'FontName',fontname);
%             savepath = [OutDir,'/model_3_2_var8_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')
%         end
%         % new
%         for i_ROI =14:18
%             c_f = figure('Position',[538   359   180*3*1.5   120*2],'Color','White'),
%             i_plot = 0;
%             for i_sess_change = 1: length(include_ses_change)
%                 task_activation_ROI_each = squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1));
%                 task_activation_ROI_each = filloutliers(task_activation_ROI_each,nan,'mean');
%                 beh_each = squeeze(beh_all(:, i_beh, include_ses_change(i_sess_change))) - beh_all(:, i_beh, 1);
%                 
%                 
%                 x_label = ['Baseline NoGo > Go Activation in ', ROI_plot_list{i_ROI}];
%             y_label = ['\Delta', beh_name{i_beh}];
%             fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
%             subplot(1,3,i_sess_change)
%             hold on
%             [s_f]= fmri_behavior_plot4(task_activation_ROI_each, beh_each, 1- group_Info1,[cl(i_sess_change,:)/2;cl(i_sess_change,:)],1,'',x_label,y_label,'',fliplr(leg), fig_p)
%             OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
%             mkdir(OutDir);
%             
%             end
%             axes_p = [0.2    0.2    0.75    0.75];
%             fontsize = 10;
%             fontname = 'Helvetica';
%             set(get(gca,'title'),'FontSize',fontsize,'FontName',fontname);%设置标题字体大小，字型
%             set(get(gca,'XLabel'),'FontSize',fontsize,'FontName',fontname);%设置X坐标标题字体大小，字型
%             set(get(gca,'YLabel'),'FontSize',fontsize,'FontName',fontname);%设置Y坐标标题字体大小，字型
%             set(gca,'FontName',fontname,'FontSize',fontsize)%设置坐标轴字体大小，字型
%             set(get(gca,'Legend'), 'FontSize',fontsize,'FontName',fontname);
%             savepath = [OutDir,'/model_3_2_var8_',ROI_list{i_ROI},'_',beh_name{i_beh},'.tiff'];
%             export_fig(savepath, '-r300', '-transparent')
%         end
%         % old - individual group * time
% %         for i_ROI =14:18
% %             c_f = figure('Position',[538   359   180*3   120*2],'Color','White'),
% %             i_plot = 0;
% %             for i_group = 1: 2
% %                 for i_sess_change = 1: length(include_ses_change)
% %                     
% %                     i_plot = i_plot +1;
% %                     task_activation_ROI_each = squeeze(task_activation_ROI_icon_itask((group_Info1 == i_group - 1),i_ROI, 1));
% %                     task_activation_ROI_each = filloutliers(task_activation_ROI_each,nan,'mean');
% %                     beh_each = squeeze(beh_all((group_Info1 == i_group - 1), i_beh, include_ses_change(i_sess_change))) - beh_all((group_Info1 == i_group - 1), i_beh, 1);
% %                     [r, p] = corr(beh_each, task_activation_ROI_each,'Rows','pairwise')
% %                     mdl = fitlm([task_activation_ROI_each], beh_each);
% %                     figure(c_f), subplot(2, length(include_ses_change), i_plot)
% %                     %                     plot(task_activation_ROI_each, beh_each,'o');
% %                     hold on; h = plot(mdl,'Color',cl(i_ROI-12,:), 'MarkerSize',2,'Marker','o','MarkerEdgeColor', cl(i_ROI-12,:),'MarkerFaceColor', cl(i_ROI-12,:)); hold off
% %                     for i_h = 2:4
% %                         set(h(i_h),'Color', cl(i_ROI-12,:));
% %                     end
% %                     title([session_list{include_ses_change(i_sess_change)},' ', leg{i_group} newline 'r=', num2str(r),'; p=',num2str(p)])% 
% %                     y_label = strrep(beh_name{i_beh}, '_','\_');
% %                     xlabel(['Baseline','Nogo > Go'],'Interpreter','tex'); ylabel(['\Delta',y_label],'Interpreter','tex');
% %                     %                     xlim([-1.8 1.8]); ylim([-2 1.1])
% %                     box on
% %                     legend off
% %                 end
% %             end
% %             
% %         end
%     end
% end


%% plot two behavior data
% i_beh1 = 24;% bmi: 21; scl-20: 24
% for i_beh2 = [4, 11, 17, 18, 23, 25, 26:58]
% beh_ibeh_change = squeeze(beh_all(:, i_beh1, include_ses_change)) - repmat(beh_all(:, i_beh1, 1), 1, length(include_ses_change));
% beh_ibeh_change_2 = squeeze(beh_all(:, i_beh2, include_ses_change)) - repmat(beh_all(:, i_beh2, 1), 1, length(include_ses_change));
% beh_ibeh_change = beh_ibeh_change(:);
% beh_ibeh_change_2 = beh_ibeh_change_2(:);
% [cl] = cbrewer('qual', 'Set1', 9);
% x_label = ['\Delta', beh_name{i_beh2}];
% y_label = ['\Delta', beh_name{i_beh1}];
% leg = {'Baseline';'2MO';'6MO';'12MO';'24MO'};
% fig_p = [2 2 9.5 6.35];% [2, 2, 6, 4]%
% mdl = fitlm(beh_ibeh_change_2, beh_ibeh_change);
% if mdl.Coefficients.pValue(2) < 0.05
%     i_beh2
%     x_label
%     mdl.Coefficients.pValue(2)
% end
% % for i_ses = include_ses_change
% %     fitlm(beh_ibeh_change_2(group_Info_ses_change == i_ses), beh_ibeh_change(group_Info_ses_change == i_ses))
% %     figure, fmri_behavior_plot4(beh_ibeh_change_2(group_Info_ses_change == i_ses), beh_ibeh_change(group_Info_ses_change == i_ses), ones(108,1),cl(i_ses-2,:),1,'',x_label,y_label,'',leg{i_ses}, fig_p)
% % end
% end
%% plot the dynamic pattern of fmri activation
i_beh = 59
fp_dynamic4_1group = [538, 359, 270, 180];
x_lim_scale = 0.05;
include_ses = [1,3];
measureNames = {'BL','2MO','6MO','12MO','24MO'};
include_ses_change = setdiff(include_ses, [1]);
% beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
% c_f = figure('Position',fp_dynamic4_1group,'Color','White'), 
% set(gca, 'InnerPosition', [0.2, 0.2, 0.75, 0.75]);%  [0.13, 0.11, 0.775, 0.8150]
beh_ibeh = squeeze(beh_all(:, i_beh, include_ses));
% % beh_ibeh = (beh_ibeh - repmat(beh_ibeh(:, 1), 1, length(include_ses)))./repmat(beh_ibeh(:, 1), 1, length(include_ses));
% beh_ibeh = filloutliers(beh_ibeh, nan, 'mean');
% h = run_rm_raincloud(beh_ibeh, ones(size(beh_ibeh,1),1),cl(9,:),1,0); % cl(1,:)

% I-CARE vs. U-CARE
c_f = figure('Position',[538, 359, 1400, 360],'Color','White'), 
select_beh = [1,16:17];
for beh_i = 1: length(select_beh)
    i_beh = select_beh(beh_i);
    subplot(1,3,beh_i)
    beh_ibeh = squeeze(beh_all(:, i_beh, include_ses));
    h = run_rm_raincloud(beh_ibeh, Interv_Info,cl(1:2,:),0,0); % cl(1,:) %ones(size(beh_ibeh,1),1)
    xlabel(strrep(beh_name{i_beh},'_','-'),'Interpreter','none')
    leg = {['U-CARE'],['I-CARE']};%{['Non'],['Resp']};%{['Non'],['Remission']};%
    interventionNames = leg;
    legend([h.s{1,1} h.s{1,2}], interventionNames,'Location','northeast');
    set(gca, 'yTickLabel', fliplr(measureNames(include_ses)),'FontSize',fontsize);
    x_lim = get(gca, 'xlim');x_lim = x_lim(2) - x_lim(1);
    xlim([min(beh_ibeh(:)) - x_lim_scale*x_lim, max(beh_ibeh(:)) + x_lim_scale*x_lim]);
   
end

%% I-CARE vs. U-CARE behavior statistics
for i_beh = [1:18]
    
    timepoints = repmat(include_ses, length(subjects_list), 1); timepoints = timepoints(:); timepoints = nominal(timepoints);
    subjects = repmat([1: length(subjects_list)]', 1, length(include_ses)); subjects =  subjects(:); subjects = nominal(subjects);
    age_sex_tl = repmat(age_sex,length(include_ses),1);
    group_info = Interv_Info; group_info = repmat(group_info, 1, length(include_ses)); group_info = nominal(group_info(:));
    beh_ibeh = squeeze(beh_all(:, i_beh, include_ses)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses));
    beh_ibeh = beh_ibeh(:);
    beh_ibeh_baseline = repmat(squeeze(beh_all(:, i_beh, 1)), 1, length(include_ses));
    beh_ibeh_baseline = beh_ibeh_baseline(:);
    % data -> table
    tbl2 = table(beh_ibeh, beh_ibeh_baseline, timepoints,subjects, age_sex_tl(:,1), age_sex_tl(:,2), group_info,...
        'VariableNames',{beh_name{i_beh},[beh_name{i_beh},'_baseline'],'Time','Subj','age','sex', 'Group'});
    tbl2.sex = categorical(tbl2.sex); tbl2.Group = nominal(tbl2.Group);
    
    % I-CARE vs. U-CARE behavior
    mixmodel_beh = fitlme(tbl2, sprintf('%s ~ %s*%s + (%s|%s)', ...
        beh_name{i_beh},'Group','Time', '1',char('Subj')),...
        'DummyVarCoding','effects','StartMethod','random');
    beh_name{i_beh}
    ss_beh = anova(mixmodel_beh); ss_beh_table = dataset2table(ss_beh)
    
end

% ICC SCL-20
%% calculating ICC (3,1)
% including BL
ICC_type = 'A-1'%'C-1';
[ICC_wBL,ICC_woBL] = ICC_data(beh_ibeh, ICC_type)
figure, corrplot(beh_ibeh,'testR','on','varNames', measureNames(include_ses));
beh_ibeh_table = array2table(beh_ibeh,'VariableNames',measureNames(include_ses));
i_task = 3; i_con = 3;
OutDir = fullfile('..','engage','SPM','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label,'DataforR');
writetable(beh_ibeh_table, [OutDir,'/',beh_name{i_beh},'_short.csv'])

% fmri_short_table = table();
% for i_task = 3%1: length(task_list)
%     for i_con = 3%1: length(task_contrasts_select{i_task})
%         % load data
%         load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
%         for i_ROI = [1:2, 19:20]%1:length(ROI_list)%1:2%[19:20]%
% %             ROI_list{i_ROI}
%             % fmri change
%             task_activation_ROI_icon_itask_iroi = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses));
% %             task_activation_ROI_icon_itask_iroi = (squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses)) ...
% %                 - repmat(squeeze(task_activation_ROI_icon_itask_iroi(:, 1)), 1, length(include_ses)))./repmat(squeeze(task_activation_ROI_icon_itask_iroi(:, 1)), 1, length(include_ses));% 108 * length(include_ses_change)
%             task_activation_ROI_icon_itask_iroi = filloutliers(task_activation_ROI_icon_itask_iroi,nan,'mean');
%             % ICC(3,1)
% %             [ICC_wBL,ICC_woBL] = ICC_data(task_activation_ROI_icon_itask_iroi, ICC_type)
% %             figure, corrplot(task_activation_ROI_icon_itask_iroi,'testR','on','varNames', measureNames(include_ses));
%             % save data
%             for i_ses = 1: length(include_ses)
%                 varname = [ROI_plot_list{i_ROI},'_', measureNames{include_ses(i_ses)}];
%                 fmri_short_table.(char(varname)) = task_activation_ROI_icon_itask_iroi(:, i_ses);
%             end 
% %             % correlation with beh
% % %             [r,p] = corrcoef([task_activation_ROI_icon_itask_iroi', beh_ibeh'],'Rows','pairwise');
% % %             r_subj = diag(r((size(task_activation_ROI_icon_itask_iroi, 1)+1:end), 1: size(beh_ibeh)));
% % %             p_subj = diag(p((size(task_activation_ROI_icon_itask_iroi, 1)+1:end), 1: size(beh_ibeh)));
%             figure('Position',fp_dynamic4_1group,'Color','White'),
%             set(gca, 'InnerPosition', [0.2, 0.2, 0.75, 0.75]);%
% %             plot(task_activation_ROI_icon_itask_iroi','o-')
%             h = run_rm_raincloud(task_activation_ROI_icon_itask_iroi, ones(size(task_activation_ROI_icon_itask_iroi,1),1), cl(10,:), 1, 1);% 
%             set(gca, 'yTickLabel', fliplr(measureNames(include_ses)),'FontSize',fontsize);
%             xlabel(['NoGo > Go Activation in ', ROI_plot_list{i_ROI}],'Interpreter','none','FontSize',fontsize); %task_list{i_task} newline 
%             x_lim = get(gca, 'xlim'); x_lim = x_lim(2) - x_lim(1);
%             xlim([min(task_activation_ROI_icon_itask_iroi(:)) - x_lim_scale*x_lim, max(task_activation_ROI_icon_itask_iroi(:)) + x_lim_scale*x_lim])
%         end
%     end
% end
% % 
% % OutDir = fullfile('..','engage','SPM','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label,'DataforR');
% % writetable(fmri_short_table, [OutDir,'/fmri_short.csv'])

            
%% plot results and do statistical analysis - only 1 group info - new
include_ses = [1,3:5];
include_ses_change = setdiff(include_ses, [1]);
covariates = repmat(age_sex,[1 1 5]);
covariates_name = {'age';'gender'};
leg = {['U-CARE'],['I-CARE']};%{['Non'],['Resp']};%{['Non'],['Remission']};% 
interventionNames = leg;
measureNames = {'BL','2MO','6MO','12MO','24MO'};
timeNames = {'0','2','6','12','24'};
group_info =  Interv_Info;%resp_ident;%Interv_Info;%resp_ident;%%remi_ident; 
group_info = repmat(group_info, 1, length(include_ses_change));  group_info = nominal(group_info(:));%group_Info = group_info; %
load(fullfile('..','engage','Responder_list_long.mat'));
load(fullfile('..','engage','Remission_list_long.mat'));
% group_info = resp_ident_long(:, include_ses_change); group_info = nominal(group_info(:));
for i_task = 3%1: length(task_list)
    for i_con = 3%1: length(task_contrasts_select{i_task})
        
        % load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
        %% prepare for model 4
        
        timepoints = repmat(include_ses_change, length(subjects_list), 1); timepoints = timepoints(:); timepoints = nominal(timepoints);
        subjects = repmat([1: length(subjects_list)]', 1, length(include_ses_change)); subjects =  subjects(:); subjects = nominal(subjects);
        FD_include_change = FD_alltasks(:, include_ses_change, i_task) - FD_alltasks(:, 1, i_task); FD_include_change = FD_include_change(:);
        age_sex_tl = repmat(age_sex,length(include_ses_change),1);
        
        beh_ibeh_change = squeeze(beh_all(:, i_beh, include_ses_change)) - repmat(beh_all(:, i_beh, 1), 1, length(include_ses_change));
        beh_ibeh_change = beh_ibeh_change(:);
        beh_ibeh_baseline = repmat(squeeze(beh_all(:, i_beh, 1)), 1, length(include_ses_change));
        beh_ibeh_baseline = beh_ibeh_baseline(:);
        for i_ROI = 1:9%3:6%1: 2%length(ROI_list)
            ROI_list{i_ROI}
            % fmri change
            task_activation_ROI_icon_itask_iroi_change = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses_change)) ...
                - repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
            task_activation_ROI_icon_itask_iroi_change = task_activation_ROI_icon_itask_iroi_change(:);
            % baseline fmri
            task_activation_ROI_icon_itask_iroi_baseline = repmat(squeeze(task_activation_ROI_icon_itask(:,i_ROI, 1)), 1, length(include_ses_change));% 108 * length(include_ses_change)
            task_activation_ROI_icon_itask_iroi_baseline = task_activation_ROI_icon_itask_iroi_baseline(:);
            % data -> table
            tbl2 = table(beh_ibeh_change, beh_ibeh_baseline, timepoints,subjects, FD_include_change, task_activation_ROI_icon_itask_iroi_baseline, task_activation_ROI_icon_itask_iroi_change, age_sex_tl(:,1), age_sex_tl(:,2), group_info,...
                'VariableNames',{beh_name{i_beh},[beh_name{i_beh},'_baseline'],'Time','Subj','fd', 'fMRI_baseline', 'fMRI','age','sex', 'Group'});
            tbl2.Group = nominal(tbl2.Group); tbl2.sex = nominal(tbl2.sex);
            % model2: change predicts change
             task_activation_ROI_icon_itask(:,i_ROI, :) = filloutliers(squeeze(task_activation_ROI_icon_itask(:,i_ROI, :)),nan,'mean');
            if length(include_ses_change)> 1
                % behavior 
                % excluding baseline
                mixmodel7 = fitlme(tbl2, sprintf('%s ~ %s*%s + %s + %s + %s + (%s|%s)', ...
                    beh_name{i_beh}, 'Time', 'Group',[beh_name{i_beh},'_baseline'],'age','sex','1',char('Subj')),...
                    'DummyVarCoding','effects','StartMethod','random');
                ss7 = anova(mixmodel7);
                % including baseline
                
                % fmri
                % excluding baseline
                mixmodel4 = fitlme(tbl2, sprintf('%s ~ %s*%s + %s + %s + %s + (%s|%s)', ...
                    'fMRI', 'Time', 'Group','fd','age','sex','1',char('Subj')),...
                    'DummyVarCoding','effects','StartMethod','random');
                ss4 = anova(mixmodel4);
                
                % two sample t-test
                '6-month'
                figure('Position', [538   359   180*1.5   180*1.5]);
                h = run_rm_raincloud([task_activation_ROI_icon_itask(:,i_ROI,3) - task_activation_ROI_icon_itask(:,i_ROI,1)], Interv_Info, cl(1:2,:),0,0);%resp_ident_long(:, include_ses)
                [h,p] = ttest2([task_activation_ROI_icon_itask(Interv_Info==0,i_ROI,3) - task_activation_ROI_icon_itask(Interv_Info==0,i_ROI,1)],...
                    [task_activation_ROI_icon_itask(Interv_Info==1,i_ROI,3) - task_activation_ROI_icon_itask(Interv_Info==1,i_ROI,1)])
%                 '12-month'
%                 [h,p] = ttest2([task_activation_ROI_icon_itask(Interv_Info==0,i_ROI,4) - task_activation_ROI_icon_itask(Interv_Info==0,i_ROI,1)],...
%                     [task_activation_ROI_icon_itask(Interv_Info==1,i_ROI,4) - task_activation_ROI_icon_itask(Interv_Info==1,i_ROI,1)])
            else
                mdl4 = fitlm(tbl2, sprintf('%s ~ %s + %s + %s + %s', ...
                    'fMRI', 'Group','fd','age','sex'))
            end
% %             if ss4.pValue(7) < 0.05 || ss4.pValue(6) < 0.05
% %                 ROI_list{i_ROI}
% %                 ss4
%                 task_activation_ROI_icon_itask(:,i_ROI, :) = filloutliers(squeeze(task_activation_ROI_icon_itask(:,i_ROI, :)),nan,'mean');
%             task_activation_ROI_icon_itask_iroi = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses));
%                 %% baseline difference
% %                 mdl = fitlm([resp_ident_long(:,1), age_sex], task_activation_ROI_icon_itask(:,i_ROI, 1))
%                 % only those with all timepoints rcplot
%                 task_activation_ROI_icon_itask_iroi_alltime = task_activation_ROI_icon_itask_iroi;
%                 task_activation_ROI_icon_itask_iroi_alltime(find(isnan(sum(task_activation_ROI_icon_itask_iroi_alltime, 2))), :) = nan;
%                 figure('Position', [538   359   180*1.5   180*1.5]);
%                 h = run_rm_raincloud(task_activation_ROI_icon_itask_iroi_alltime, Interv_Info, cl(1:2,:),0,0);%resp_ident_long(:, include_ses)
%                 legend([h.s{1,1} h.s{1,2}], interventionNames,'Location','southeast');
%                 set(gca, 'yTickLabel', fliplr(measureNames(include_ses)),'FontSize',fontsize);
%                 title([task_list{i_task} newline ROI_list{i_ROI}],'Interpreter','none','FontSize',fontsize);
%                 xlim([min(task_activation_ROI_icon_itask_iroi_alltime(:)) - x_lim_scale*(max(task_activation_ROI_icon_itask_iroi_alltime(:)) - min(task_activation_ROI_icon_itask_iroi_alltime(:))), max(task_activation_ROI_icon_itask_iroi_alltime(:)) + x_lim_scale*(max(task_activation_ROI_icon_itask_iroi_alltime(:)) - min(task_activation_ROI_icon_itask_iroi_alltime(:)))])
% %             end
        end
    end
end


%% plot results and do statistical analysis - only 1 group info -old
include_ses = [1,2];
group_Info = resp_ident;%remi_ident;%Interv_Info;%resp_ident;%
leg = {['Non'],['Resp']};% {['Non'],['Remission']};%leg = {['uCARE'],['iCARE']};%
interventionNames = leg;
measureNames = {'BL','2MO','6MO','12MO','24MO'};
timeNames = {'0','2','6','12','24'};
withinsubNames = {'time'};
betweensubNames = {'groups'};
subjNames = {'subj'};
covariates = repmat(age_sex,[1 1 5]);
covariates_name = {'age';'gender'};
for i_task = 3%1: length(task_list)
    for i_con = 3%1: length(task_contrasts_select{i_task})
        % load data
        load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
        for i_ROI = 1:2% length(ROI_list)
            task_activation_ROI_icon_itask(:,i_ROI, :) = filloutliers(squeeze(task_activation_ROI_icon_itask(:,i_ROI, :)),nan,'mean');
            task_activation_ROI_icon_itask_iroi = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses));
            [ranovatbl, c, mixmodel2, ss2] = statistics_longtidutinal_data(squeeze(task_activation_ROI_icon_itask(:,i_ROI, :)), include_ses, group_Info, FD_alltasks(:, :, i_task), measureNames, betweensubNames, withinsubNames,timeNames);
%             if ss2.pValue(5) < 0.05  || ranovatbl.pValue(5) < 0.05 %|| ss2.pValue(3) < 0.05
                ROI_list{i_ROI}
                ss2
                %% baseline difference
                mdl = fitlm([group_Info, age_sex], task_activation_ROI_icon_itask(:,i_ROI, 1))
                % only those with all timepoints rcplot
                task_activation_ROI_icon_itask_iroi_alltime = task_activation_ROI_icon_itask_iroi;
                task_activation_ROI_icon_itask_iroi_alltime(find(isnan(sum(task_activation_ROI_icon_itask_iroi_alltime, 2))), :) = nan;
                figure('Position', fig_position1);
                h = run_rm_raincloud(task_activation_ROI_icon_itask_iroi_alltime, group_Info, cl(1:2,:),1,0);
                legend([h.p{1,1} h.p{1,2}], interventionNames,'Location','southeast');
                set(gca, 'yTickLabel', fliplr(measureNames(include_ses)),'FontSize',fontsize);
                title([task_list{i_task} newline ROI_list{i_ROI}],'Interpreter','none','FontSize',fontsize);
                xlim([min(task_activation_ROI_icon_itask_iroi_alltime(:)) - x_lim_scale*(max(task_activation_ROI_icon_itask_iroi_alltime(:)) - min(task_activation_ROI_icon_itask_iroi_alltime(:))), max(task_activation_ROI_icon_itask_iroi_alltime(:)) + x_lim_scale*(max(task_activation_ROI_icon_itask_iroi_alltime(:)) - min(task_activation_ROI_icon_itask_iroi_alltime(:)))])
%             end
        end
    end
end

%% plot results and do statistical analysis - 2 group infos

% include_ses = [1: 5];
% group_Info1 = Interv_Info;
% group_Info2 = resp_ident;%remi_ident;%
% group_Info = nan(size(group_Info1));
% group_Info(intersect(find(group_Info1 == 0), find(group_Info2 == 1))) = 1;
% group_Info(intersect(find(group_Info1 == 0), find(group_Info2 == 0))) = 2;
% group_Info(intersect(find(group_Info1 == 1), find(group_Info2 == 1))) = 3;
% group_Info(intersect(find(group_Info1 == 1), find(group_Info2 == 0))) = 4;
% leg = {['U-CARE: Responder'];['U-CARE: Non'];['I-CARE: Responder'];['I-CARE: Non']};
% interventionNames = leg;
% measureNames = {'BL','2MO','6MO','12MO','24MO'};
% timeNames = {'0','2','6','12','24'};
% withinsubNames = {'time'};
% betweensubNames = {'group1','group2'};
% subjNames = {'subj'};
% covariates = repmat(age_sex,[1 1 5]);
% covariates_name = {'age';'gender'};
% for i_task = 3%1: length(task_list)
%     for i_con = 3%1: length(task_contrasts_select{i_task})
%         % load data
%         load(fullfile(DATADIR, task_list{i_task},task_contrasts_select{i_task}{i_con}, ['task_activation_', ROI_label, '.mat']));
%         for i_ROI = 1 : length(ROI_list)
%             task_activation_ROI_icon_itask_iroi = squeeze(task_activation_ROI_icon_itask(:,i_ROI, include_ses));
%             % rm outliers
%             task_activation_ROI_icon_itask_iroi = filloutliers(task_activation_ROI_icon_itask_iroi,nan,'mean');
%             task_activation_ROI_icon_itask(:,i_ROI, :) = filloutliers(task_activation_ROI_icon_itask(:,i_ROI, :),nan,'mean');
%             [ranovatbl, mixmodel2, ss2] = statistics_longtidutinal_data_2group(squeeze(task_activation_ROI_icon_itask(:,i_ROI, :)), include_ses, group_Info1, group_Info2, FD_alltasks(:, :, i_task), measureNames, betweensubNames, withinsubNames,timeNames);
% %             if ss2.pValue(end) < 0.05  || ranovatbl.pValue(9) < 0.05 %|| ss2.pValue(3) < 0.05
%                 ROI_list{i_ROI}
%                 ss2
%                 %% baseline difference
% %                 mdl = fitlm([group_Info1, age_sex], task_activation_ROI_icon_itask(:,i_ROI, 1))
% %                 mdl1 = fitlm([group_Info1, age_sex, task_activation_ROI_icon_itask(:,i_ROI, 1)], task_activation_ROI_icon_itask(:,i_ROI, 3) - task_activation_ROI_icon_itask(:,i_ROI, 1))
%                 % only those with all timepoints rcplot
%                 task_activation_ROI_icon_itask_iroi_alltime = task_activation_ROI_icon_itask_iroi;
%                 task_activation_ROI_icon_itask_iroi_alltime(find(isnan(sum(task_activation_ROI_icon_itask_iroi_alltime, 2))), :) = nan;
%                 task_activation_ROI_icon_itask_iroi_alltime = task_activation_ROI_icon_itask_iroi;
%                 figure('Position', fig_position1);
%                 h7 = run_rm_raincloud(task_activation_ROI_icon_itask_iroi_alltime, group_Info, cb([1,3:4,5],:));
%                 legend([h7.p{1,1} h7.p{1,2} h7.p{1,3} h7.p{1,4}], leg, 'Location','southeast');
%                 set(gca, 'yTickLabel', fliplr(measureNames(include_ses)),'FontSize',fontsize);
%                 title([task_list{i_task} newline ROI_list{i_ROI}],'Interpreter','none','FontSize',fontsize);
%                 xlim([min(task_activation_ROI_icon_itask_iroi_alltime(:)) - x_lim_scale*(max(task_activation_ROI_icon_itask_iroi_alltime(:)) - min(task_activation_ROI_icon_itask_iroi_alltime(:))), max(task_activation_ROI_icon_itask_iroi_alltime(:)) + x_lim_scale*(max(task_activation_ROI_icon_itask_iroi_alltime(:)) - min(task_activation_ROI_icon_itask_iroi_alltime(:)))])
% %             end
%         end
%     end
% end



%% plot linear prediction from R
data = readtable('/Users/xuezhang/Documents/Stanford/Projects/ENGAGE/Analysis/engage/SPM/Level2/MixedModel/111_fMRI_stats_spikesonly_FD_fromFile_GO_NO_GO/NogovsGo/SCL_20/MixedModel_SCL20/DataforR/linear_prediction.csv', 'PreserveVariableNames', true);
leg = {'6MO';'12MO';'24MO'};
% plot for model 1
s = rng('default'); rng(s); ColorMap=rand(400,3);
alpha = 0.5;

x_label = ['Predicted \DeltaSCL-20'];
y_label = ['Actual \DeltaSCL-20'];
fig_p = [2 2 8 7.5];% [2, 2, 6, 4]%
[s_f]= fmri_behavior_plot4(data.SCL_20_change, data.predicted_SCL_fmri, data.Time,'',1,'',x_label,y_label,'',leg, fig_p)
%             ylim([-2.5 2])
% OutDir = fullfile('..','engage','SPM','Figures','Level2','MixedModel',task_list{i_task},task_contrasts_select{i_task}{i_con},beh_name{i_beh}, ROI_label);
mkdir(OutDir);
savepath = [OutDir,'/ccn_linear_prediction.tiff'];
export_fig(savepath, '-r300', '-transparent')



function [data_change] = prepare_data_change(data, include_ses)
data_fil = filloutliers(data,nan,'mean');
% data_fil = data;
include_ses_change = setdiff(include_ses, [1]);
data_change = data_fil(:, include_ses_change) - repmat(data_fil(:, 1), 1, length(include_ses_change));% 108 * length(include_ses_change)
% data_change = data_change(:);
end
function [data_fil_in] = prepare_data_all(data, include_ses)
% data_fil = filloutliers(data,nan,'mean');
data_fil = data;
data_fil_in = squeeze(data_fil(:, include_ses));
data_fil_in = data_fil_in(:);
end

function [ICC_wBL,ICC_woBL] = ICC_data(data, type)
% including BL
data_index = find(~isnan(sum(data, 2)));
data_new = data(data_index,:); % remove nans
[ICC_wBL] = ICC(data_new, type);
clear data_new, clear data_index
% excluding BL
data_index = find(~isnan(sum(data(:,2:end), 2)));
data_new = data(data_index,2:end); % remove nans
[ICC_woBL] = ICC(data_new, type);
end
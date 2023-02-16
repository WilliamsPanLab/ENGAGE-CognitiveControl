clear, clc, close all

%% Main effect
fig_p = [10,10,2.8,3];%[10,10,7,7];
xx = rand(108,1);
y = 1*xx + rand(108,1)*0.5;
group = [ones(length(y),1);];
figure,fmri_behavior_plot4(xx, y, group,[0 0 0],0,1,[],[],[],[],fig_p);
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);

leg = {['I-CARE 6MO'], ['U-CARE 6MO'],['I-CARE 12MO'],['U-CARE 12MO'],['I-CARE 24MO'],['U-CARE 24MO'];};
cl = [[227, 136, 47];  [0, 157, 118]; [73, 84, 138];]./255.0; c_m = [0.7*cl(1,:); 1.1*cl(1,:);0.7*cl(2,:); 1.2*cl(2,:); 0.7*cl(3,:); 1.1*cl(3,:)]


%% Circuit x Intervention
leg = {['I-CARE'],['U-CARE'];};
group_N = 2;
fig_p = [10,10,2.8,3];
x = rand(100,1);
y = 1*x + rand(100,1)*0.5;
y1 = -1*x + rand(100,1)*1;
y2 = 1.5*x + rand(100,1)*1 - 1.2;
group = [ones(length(y1),1);2*ones(length(y1),1);];
x_label = [];%['2MO-BL \DeltaNoGo > Go Activation'];
y_label = [];%['\Delta Symptom-behavior'];
figure, fmri_behavior_plot4([x;x], [y1;y2], group,[183, 90, 101; 216, 216, 222]/255.0,0,1,x_label,y_label,'',leg, fig_p)
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
hLegend = findobj(gcf, 'Type', 'Legend');
hLegend.FontSize = 8;
legend off 

%% Circuit x Time
include_ses = [1,3:5];
include_ses_change = setdiff(include_ses, [1]);
leg = {'Baseline';'2MO';'6MO';'12MO';'24MO'}; leg = leg(include_ses_change);
fig_p = [10,10,2.8,3];
x = rand(100,1);
y = 1*x + rand(100,1)*0.5;
y1 = 1.5*x + rand(100,1)*1;
y2 = 1*x + rand(100,1)*1;
y3 = 0.3*x + rand(100,1)*1;
group = [ones(length(y1),1);2*ones(length(y1),1);3*ones(length(y1),1)];
x_label = [];%['2MO-BL \DeltaNoGo > Go Activation'];
y_label = [];%['\Delta Symptom-behavior'];
figure, fmri_behavior_plot4([x;x;x], [y1;y2;y3], group,cl,0,1,x_label,y_label,'','', fig_p)
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
hLegend = findobj(gcf, 'Type', 'Legend');
if ~isempty(hLegend)
    hLegend.FontSize = 8;
    legend boxoff
end

%% Circuit x Time x Intervention 
group_N = 2;
fig_p = [10,10,2.8,3];%[10,10,7,7];
x = rand(54,1);
xx = rand(108*3,1);
y = 1*xx + rand(108*3,1)*0.5;
y1 = 2*x + rand(54,1)*1;
y2 = 1*x + rand(54,1)*1;
y3 = 0.3*x + rand(54,1)*1;
y4 = 1.5*x + rand(54,1)*1;
y5 = 0.6*x + rand(54,1)*1;
y6 = 0.1*x + rand(54,1)*1;
group = [ones(length(y1),1);2*ones(length(y1),1);3*ones(length(y1),1);4*ones(length(y1),1);5*ones(length(y1),1);6*ones(length(y1),1);];
x_label = []; %['2MO-BL \DeltaNoGo > Go Activation'];
y_label = []; %['\Delta Symptom-behavior'];
figure, fmri_behavior_plot4([x;x;x;x;x;x], [y1;y2;y3;y4;y5;y6], group,c_m,0,1,x_label,y_label,'',leg, fig_p);
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
legend off

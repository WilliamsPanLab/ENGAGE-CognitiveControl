%% plot
% fig_p = [2   2    12    10];
fig_p = [2   2  9.5 7];
axes_p = [0.2    0.2    0.70    0.70];
fontsize= 10;
fontname='Times New Roman';
linewidth = 2;

% parameters for raincloud
% [cb] = cbrewer('qual', 'Set3', 10, 'pchip');
% cl(1, :) = cb(4, :);
% cl(2, :) = cb(5, :);
% cl(3, :) = cb(1, :);
% cl(4, :) = cb(3, :);
% cl(5, :) = cb(2, :);
% cl(6, :) = cb(6, :);
addpath(genpath('/Users/xuezhang/Downloads/Software/RainCloudPlots'))
[cb] = cbrewer('qual', 'Set1', 10);  cl = cb;
fig_position01 = [200 200 250 200];
fig_position0 = [200   247   138   128];%[200 200 200 175];
fig_position1 = [200 200 400 350];
fig_position12 = [200 200 400 450];
fig_position2 = [200 200 500 350]; % coordinates for figures
fig_position3 = [200 200 1400 1000]; % coordinates for figures
fig_position34 = [200 200 1000 800]; % coordinates for figures
fig_position4 = [60 244 1400 300]; % wide
xlabel_angle = -45;
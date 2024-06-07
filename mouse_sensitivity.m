% Sensitivity analysis for a simple mouse model
clear; clc; close all;
%%
% load m2_gauss_final.mat P_LV_avg P_SA_avg V_LV_avg
% load m2_timeseries_final.mat PLV_stack VLV_stack PSA_stack
load AJP_test_data_mouse.mat
T = 0.11;
param = mouse_parameters_newm1(V_LV_avg,P_LV_avg,P_SA_avg,T);
Names = {'$P_{LA}$','$P_{Sys}$','$R_{mv}$','$R_{av}$',...
        '$R_{art}$','$E_{es}$','$E_{ed}$','$T_{s}$',...
        '$T_{e}$','$T$','$V_{lv,d}$','$C_{Ao}$'};

%%
% Definie initial conditions, time variables, and ODE options
Vlv_init = 80; Vao_init = 20;
tstart = 0; tend   = 30.*T; % 30 cycles
tspace = tstart:1e-3:tend;
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Necessary for the ODE solver
IC = [Vlv_init; Vao_init]; % Our initial conditions

%%
% Calculate the model sensitivity
sens_func = @(q) call_model(q,IC,tspace);
ynom= sens_func(param);
% ids = 1:12;
ids = [3:9 11:12];
% ids = [3:9 11:12];

% ids = [3 5:6 7:9 11:12];
% ids = [3 6:9];
% ids = [7 8 9];
[sens,F]= local_sensitivity(param,sens_func,ynom,1e-1,ids);

%% Formulate the residual vector
ids2 = 1:9;
% ids2 = [1 3 4 5 6 7 8 9]; % o Av valve
% ids2 = [1 3 4 5 6 7 9];   % no dead volume
% ids2 = [3 4 5 6 7 9];   % no mitral valve


s2 = [1 1 1];%[3 2 0.2];
ynom = ynom./s2';
y_nom_stack = [ynom(1,:) ynom(2,:) ynom(3,:)];
sens_R = [squeeze(sens(ids2,1,:)) squeeze(sens(ids2,2,:)) squeeze(sens(ids2,3,:))].*y_nom_stack;
% sens_R(:,1:51)    = sens_R(:,1:51)./mean(ynom(1,:));
% sens_R(:,52:102)  = sens_R(:,52:102)./mean(ynom(2,:));
% sens_R(:,103:end) = sens_R(:,103:end)./mean(ynom(3,:));

sens_R(:,1:51)    = sens_R(:,1:51)./mean(V_LV_avg);
sens_R(:,52:102)  = sens_R(:,52:102)./mean(P_LV_avg);
sens_R(:,103:end) = sens_R(:,103:end)./mean(P_SA_avg);



%% We investigate the sensitivity of three ouputs: LV pressure, LV volume,
% and systemic arterial pressure
sens_LVV = squeeze(sens(ids2,1,:));
sens_LVP = squeeze(sens(ids2,2,:));
sens_SAP = squeeze(sens(ids2,3,:));

rank_sensLVV = sqrt(sum(sens_LVV.^2,2));
rank_sensLVP = sqrt(sum(sens_LVP.^2,2));
rank_sensSAP = sqrt(sum(sens_SAP.^2,2));
rank_sensR   = sqrt(sum(sens_R.^2,2));
%% Dots
figure; clf;
subplot(4,1,1); 
plot(rank_sensLVV,'*','MarkerSize',8,'LineWidth',2);
ylabel('$V_{LV}$ Sens'); set(gca,'FontSize',20);xticks(1:length(ids)); 
xticklabels(''); grid on;
% set(gca,'YScale','log')

subplot(4,1,2);
plot(rank_sensLVP,'*','MarkerSize',8,'LineWidth',2);
ylabel('$p_{LV}$ Sens'); set(gca,'FontSize',20);xticks(1:length(ids)); 
xticklabels(''); grid on;
% set(gca,'YScale','log')


subplot(4,1,3); 
plot(rank_sensSAP,'*','MarkerSize',8,'LineWidth',2);
ylabel('$p_{SA}$ Sens'); set(gca,'FontSize',20);xticks(1:length(ids));
xticklabels(''); grid on;
% set(gca,'YScale','log')

subplot(4,1,4); 
plot(rank_sensR,'*','MarkerSize',8,'LineWidth',2);
ylabel('Combined'); set(gca,'FontSize',20);xticks(1:length(ids));
xticklabels(Names(ids)); grid on;
% set(gca,'YScale','log')
%% Bars
figure; clf;
subplot(4,1,1); 
bar(rank_sensLVV);
ylabel('$V_{LV}$ Sens'); set(gca,'FontSize',20);xticks(1:length(ids2)); 
xticklabels(''); grid on;
% set(gca,'YScale','log')

subplot(4,1,2);
bar(rank_sensLVP);
ylabel('$p_{LV}$ Sens'); set(gca,'FontSize',20);xticks(1:length(ids2)); 
xticklabels(''); grid on;
% set(gca,'YScale','log')


subplot(4,1,3); 
bar(rank_sensSAP);
ylabel('$p_{SA}$ Sens'); set(gca,'FontSize',20);xticks(1:length(ids2));
xticklabels(''); grid on;
% set(gca,'YScale','log')

subplot(4,1,4); 
bar(rank_sensR);
ylabel('Combined'); set(gca,'FontSize',20);xticks(1:length(ids2));
xticklabels(Names(ids(ids2))); grid on;
% set(gca,'YScale','log')
%% Identifiability analysis
F = sens_R*sens_R';
cov_mat = inv(F);
corr_mat = corrcov(cov_mat);

figure; h = heatmap(corr_mat); colormap(jet)
h.XDisplayLabels = Names(ids(ids2));
h.YDisplayLabels = Names(ids(ids2));
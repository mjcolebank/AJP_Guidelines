% Mouse tornado plot
% Sensitivity analysis for a simple mouse model
clear; clc; close all;
%%

load AJP_test_data_m1.mat
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

ids = [3:9 11:12];
delta_VLV = zeros(length(ids),2,2);
delta_PLV = zeros(length(ids),2,2);
delta_PSA = zeros(length(ids),2,2);

param0 = param;
ynom = sens_func(param0);
ynom_max = max(ynom,[],2);
ynom_min = min(ynom,[],2);
for i=1:length(ids)
    param_plus = param0;
    param_plus(ids(i)) = param0(ids(i)).*1.1;
    yplus = sens_func(param_plus);

    param_minus = param0;
    param_minus(ids(i)) = param0(ids(i)).*0.9;
    yminus = sens_func(param_minus);
    
    delta_VLV(i,1,1) = max(yplus(1,:),[],2)-ynom_max(1);
    delta_PLV(i,1,1) = max(yplus(2,:),[],2)-ynom_max(2);
    delta_PSA(i,1,1) = max(yplus(3,:),[],2)-ynom_max(3);

    delta_VLV(i,1,2) = max(yminus(1,:),[],2)-ynom_max(1);
    delta_PLV(i,1,2) = max(yminus(2,:),[],2)-ynom_max(2);
    delta_PSA(i,1,2) = max(yminus(3,:),[],2)-ynom_max(3);

    delta_VLV(i,2,1) = min(yplus(1,:),[],2)-ynom_min(1);
    delta_PLV(i,2,1) = min(yplus(2,:),[],2)-ynom_min(2);
    delta_PSA(i,2,1) = min(yplus(3,:),[],2)-ynom_min(3);

    delta_VLV(i,2,2) = min(yminus(1,:),[],2)-ynom_min(1);
    delta_PLV(i,2,2) = min(yminus(2,:),[],2)-ynom_min(2);
    delta_PSA(i,2,2) = min(yminus(3,:),[],2)-ynom_min(3);


    
end

% Relative
    delta_VLV(:,1,:) = 100*delta_VLV(:,1,:)./ynom_max(1);
    delta_VLV(:,2,:) = 100*delta_VLV(:,2,:)./ynom_min(1);

    delta_PLV(:,1,:) = 100*delta_PLV(:,1,:)./ynom_max(2);
    delta_PLV(:,2,:) = 100*delta_PLV(:,2,:)./ynom_min(2);

    delta_PSA(:,1,:) = 100*delta_PSA(:,1,:)./ynom_max(3);
    delta_PSA(:,2,:) = 100*delta_PSA(:,2,:)./ynom_min(3);

%%
% VLV_P = mean(squeeze(delta_VLV(:,:,1))')

figure(100); clf; %
tiledlayout(2, 3, "TileSpacing", "compact");

nexttile; hold on;
barh(squeeze(delta_VLV(:,1,1))')
barh(squeeze(delta_VLV(:,1,2))')
yticks(1:11); yticklabels(Names(ids));
grid on; set(gca,'FontSize',20); 
xlabel('\% Change in max $V_{LV}$')
xlim([-6.5 6.5])

nexttile; hold on;
barh(squeeze(delta_PLV(:,1,1))')
barh(squeeze(delta_PLV(:,1,2))')
grid on; set(gca,'FontSize',20);
yticks(1:10); yticklabels('');
xlabel('\% Change in max $p_{LV}$')
xlim([-6.5 6.5])

nexttile; hold on;
barh(squeeze(delta_PSA(:,1,1))')
barh(squeeze(delta_PSA(:,1,2))')
grid on; set(gca,'FontSize',20);
yticks(1:10); yticklabels('');
xlabel('\% Change in max $p_{Ao}$')
xlim([-6.5 6.5])

nexttile; hold on;
barh(squeeze(delta_VLV(:,2,1))')
barh(squeeze(delta_VLV(:,2,2))')
yticks(1:11); yticklabels(Names(ids));
grid on; set(gca,'FontSize',20);
xlabel('\% Change in min $V_{LV}$')
xlim([-6.5 6.5])



nexttile; hold on;
barh(squeeze(delta_PLV(:,2,1))')
barh(squeeze(delta_PLV(:,2,2))')
grid on; set(gca,'FontSize',20); 
yticks(1:10); yticklabels('');
xlabel('\% Change in min $p_{LV}$')
xlim([-6.5 6.5])



nexttile; hold on;
barh(squeeze(delta_PSA(:,2,1))')
barh(squeeze(delta_PSA(:,2,2))')
grid on; set(gca,'FontSize',20);
yticks(1:10); yticklabels('');
xlabel('\% Change in min $p_{Ao}$')
xlim([-6.5 6.5])


legend('$+10\%$','$-10\%$')

%%
for i=1:6
    figure(i);xlim([-0.1 0.1]);
end
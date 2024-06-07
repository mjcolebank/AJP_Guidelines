% Driver for AJP Guidelines Model
%
% Author: Mitchel J. Colebank
% Date: 4/17/2024
%%
clear; clc; close all;
%%
load AJP_test_data_mouse.mat PLV_stack PSA_stack P_LV_avg P_SA_avg VLV_stack V_LV_avg
data_shift = 0;
model_shift = -2;
T = 0.11;

% Get parameters
param = mouse_parameters_newm1(V_LV_avg,P_LV_avg,P_SA_avg,T);

% First, define some initial volumes for the differential equation
% Volumes are in microliters
Vlv_init = 80; % LV Volume
Vao_init = 20;  % Aortic Volume
% Starting and ending time
tstart = 0;
tend   = 30.*T; % 30 cycles
dt = 1e-3;
tspace = tstart:dt:tend;
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Necessary for the ODE solver


Initial_Conditions = [Vlv_init; Vao_init]; % Our initial conditions
% Now solve our system of ODEs
y = ode45(@LV_Ao,[tstart, tend],Initial_Conditions,options,param);


% We solved the model for 30 heartbeats; we want to extract the last two
% for plotting
tplot = linspace(28*T,30*T,99);
yout = deval(y,tplot);
tplot = tplot-tplot(1);


% Extract the solutions to the two differential equations
Vlv = yout(1,:);
Vao = yout(2,:);
% param = [P_LA,P_SysCap,Rmv,Rav,Rart,Emax,Emin,T_peak,T_relax,T,Vlv_d,Cao];
P_LA     = param(1);
P_SysCap = param(2);
Rmv      = param(3);
Rav      = param(4);
Rart     = param(5);
Cao      = param(12);
% Now, recompute the pressures and flows
plv = LinearElastance(Vlv,tplot,param([6:7 11 8:10]));%[Emax,Emin,Vlv_d,T_peak,T_relax,T]);
pao = Vao./Cao;

% Use the 'max' operator to keep positive flows for valves
QMV = max((P_LA-plv)./Rmv,0);
QAV = max((plv-pao)./Rav,0);
Qart_sys = (pao-P_SysCap)./Rart;

tdata = linspace(0,T,50);

%% Shift output that is for comparison
plv = circshift(plv,model_shift);
Vlv = circshift(Vlv,model_shift);
pao = circshift(pao,model_shift);

plv = plv(1:50);
Vlv = Vlv(1:50);
pao = pao(1:50);
tplot = tplot(1:50);

%% Plotting
% Pressure
figure(1); clf; hold on;
plot(tdata,P_LV_avg,'k','LineWidth',3);
plot(tplot,plv,':c','LineWidth',2);
plot(tdata,P_SA_avg,'r','LineWidth',3)
plot(tplot,pao,':m','LineWidth',2);
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
grid on;
set(gca,'FontSize',20)

%% PV loop
figure(2); hold on;
plot(Vlv,plv,'r','LineWidth',3);
plot(V_LV_avg,P_LV_avg,'k','LineWidth',3)
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
grid on;
set(gca,'FontSize',20);

%% Individual pressures
tdata2 = linspace(0,T,50);
figure(3);clf; 
subplot(2,2,1); hold on;
plot(tdata2,VLV_stack,'Color',[0.8 0.8 0.8]);
plot(tplot,Vlv,'r','LineWidth',2);
subplot(2,2,2); hold on;
plot(tdata2,PLV_stack,'Color',[0.8 0.8 0.8]);
plot(tplot,plv,'r','LineWidth',2);
subplot(2,2,3); hold on;
plot(VLV_stack,PLV_stack,'Color',[0.8 0.8 0.8]);
plot(Vlv,plv,'r','LineWidth',2);
subplot(2,2,4); hold on;
plot(tdata2,PSA_stack,'Color',[0.8 0.8 0.8]);
plot(tplot,pao,'r','LineWidth',2);

tdata2 = linspace(0,T,50);
figure;clf; 
subplot(1,3,1); hold on;
plot(tdata2,VLV_stack,'Color',[0.3 0.3 0.3]);
plot(tplot,Vlv,'r','LineWidth',3);
xlabel('Time (s)'); ylabel('$V_{LV}$ ($\mu$L)'); 
grid on; set(gca,'FontSize',20);

subplot(1,3,2); hold on;
plot(tdata2,PLV_stack,'Color',[0.3 0.3 0.3]);
plot(tplot,plv,'r','LineWidth',3);
xlabel('Time (s)'); ylabel('$p_{LV}$ (mmHg)'); 
grid on; set(gca,'FontSize',20);

subplot(1,3,3); hold on;
h1 = plot(tdata2,PSA_stack(:,2:end),'Color',[0.3 0.3 0.3]);
h2 = plot(tplot,pao,'r','LineWidth',3);
xlabel('Time (s)'); ylabel('$p_{Ao}$ (mmHg)'); 
grid on; set(gca,'FontSize',20);
legend([h1(1) h2],{'Data','Model'})
set(gcf,'Position',[153.00        233.80       1153.60        539.20])

res_VLV = V_LV_avg - Vlv';
res_PLV = P_LV_avg - plv';
res_PSA = P_SA_avg - pao';

figure;
subplot(1,3,1); qqplot(res_VLV);
subplot(1,3,2); qqplot(res_PLV);
subplot(1,3,3); qqplot(res_PSA);


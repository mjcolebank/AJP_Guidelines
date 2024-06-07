% Run optimization with the mice data
clear; clc; close all;
%%
load AJP_test_data_m1.mat
data_shift = 0;
model_shift = -2;
VLV_stack_shift = circshift(VLV_stack,data_shift);
PLV_stack_shift = circshift(PLV_stack,data_shift);
PSA_stack = circshift(PSA_stack,data_shift);
V_LV_avg = circshift(V_LV_avg,data_shift);
P_LV_avg = circshift(P_LV_avg,data_shift);
P_SA_avg = circshift(P_SA_avg,data_shift);
T = 0.11;


param = mouse_parameters_newm1(V_LV_avg,P_LV_avg,P_SA_avg,T);
Names = {'$P_{LA}$','$P_{Sys}$','$R_{mv}$','$R_{av}$',...
        '$R_{art}$','$E_{es}$','$E_{ed}$','$T_{max}$',...
        '$T_{min}$','$T$','$V_{lv,d}$','$C_{Ao}$'};

% Definie initial conditions, time variables, and ODE options
Vlv_init = 80; Vao_init = 20;
tstart = 0; tend   = 30.*T; % 30 cycles
tspace = tstart:1e-3:tend;
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Necessary for the ODE solver
IC = [Vlv_init; Vao_init]; % Our initial conditions

UB = [10; 45; 1e-2; 1e-2; 1.0; 6; 1; 0.6.*T; 0.8.*T; 0.11; 10; 1.5];
LB = [2; 10; 1e-4; 1e-4;  0.1; 0.5; 1e-2; 0.3.*T; 0.6.*T; 0.11; 1; 0.4];


%%
% Run an optimization
% ids = [1:9 11:12];%[3 5:9];
% ids = [3:9 12]
% ids = [3 5 7:9 11:12];
% ids = [3 6:9];
% ids = [7 8 9];
% ids = [8 9];

% ids = [3 5 6 7 8 9 11 12];


% No Vd
% ids = [3 5 6 7 8 9 12];

% No RMV
% ids = [5 6 7 8 9 12];

% No Cao, but with Rmv
% ids = [3 5 6 7 8 9];

% No Ees
ids = [3 5 7 8 9 12]

% ids = [12];
par0_OG = param(ids);
n_samp = 100;
Xopt_save = zeros(n_samp,length(ids));
par0_save = zeros(n_samp,length(ids));
J_save    = zeros(n_samp,1);
jacobian_save = cell(n_samp,1);
    data = [V_LV_avg(:); P_LV_avg(:); P_SA_avg(:)];

    s2vec = [1 1 1]; % [mean(V_LV_avg); mean(P_LV_avg); mean(P_SA_avg)]'; % [6 2 1];
    J = @(q) mouse_residual(q,data,ids,param,IC,tspace,s2vec);
    LB = LB(ids);
    UB = UB(ids);

%     LB = log(LB);
%     UB = log(UB);
% tpred = linspace(29*T,30*T,50);
% model_func = @(q,dummy) mouse_fwd_UQ(exp(q),ids,param,IC,tspace);

parfor i=1:n_samp
    disp(i)
%     par0 = unifrnd(exp(LB),exp(UB));
    par0 = par0_OG.*unifrnd(0.1,3,1,length(ids));
    options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','FiniteDifferenceType','central');
        [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(J,par0,LB,UB,options);

        % Nlinfit
%         options = statset('Display','iter');
%         [X,RESIDUAL,JACOBIAN,CovB,RESNORM,ErrorModelInfo] = nlinfit(tpred,data,model_func,log(par0),options)

    par0_save(i,:) = par0;
    Xopt_save(i,:) = X;
    J_save(i) = RESNORM;
    jacobian_save{i} = JACOBIAN;
    %%
   
end
save('mouseNEW_opt_100samp_6par_noEes_shift2','par0_OG','Xopt_save','par0_save','J_save','jacobian_save','ids');



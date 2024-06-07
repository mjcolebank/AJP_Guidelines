%% Construct parameter and output confidence intervals
load AJP_test_data_m1.mat
data_shift = 0;
model_shift = -2;
VLV_stack_shift = circshift(VLV_stack,data_shift); PLV_stack_shift = circshift(PLV_stack,data_shift); PSA_stack = circshift(PSA_stack,data_shift);
V_LV_avg = circshift(V_LV_avg,data_shift); P_LV_avg = circshift(P_LV_avg,data_shift); P_SA_avg = circshift(P_SA_avg,data_shift);
T = 0.11;
data      = [V_LV_avg(:); P_LV_avg(:); P_SA_avg(:)];
data_mean = [ones(50,1).*mean(V_LV_avg); ones(50,1).*mean(P_LV_avg); ones(50,1).*mean(P_SA_avg)];

param = mouse_parameters_newm1(V_LV_avg,P_LV_avg,P_SA_avg,T);

% load mouseNEW_opt_100samp_8par_withVd_shift2.mat J_save Xopt_save ids jacobian_save par0_save 
load mouseNEW_opt_100samp_7par_shift2 J_save Xopt_save ids jacobian_save par0_save 
% load mouseNEW_opt_100samp_6par_noRMV_shift2.mat J_save Xopt_save ids jacobian_save par0_save 


[Jbest,id_best] = sort(J_save);
par_opt = Xopt_save(id_best(1),:);
param(ids) = par_opt;
Jacob_opt  = jacobian_save{id_best(1)};

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

out_func = @(q,dummyT) mouse_fwd_UQ(q,ids,param,Initial_Conditions,tspace);
yopt = out_func(par_opt,[]);
res_opt = yopt(:) - data(:);
RSS_opt = sum (res_opt).^2 ;
s2_estimate = RSS_opt./(length(data(:))-length(ids));
%%
tpred = linspace(0,T,50);
[Ypred,delta_CI] = nlpredci(out_func,tpred,par_opt,res_opt,'Jacobian',Jacob_opt);
[Ypred,delta_PI] = nlpredci(out_func,tpred,par_opt,res_opt,'Jacobian',Jacob_opt,'PredOpt','observation');

par_CI = nlparci(par_opt,res_opt,'jacobian',Jacob_opt);
[par_CI(:,1) par_opt' par_CI(:,2)]
[par_CI(:,1)./par_opt' par_CI(:,2)./par_opt']-1

%%
figure(99); hold on;
plot(delta_CI,'--');
plot(delta_PI,'-');

figure;
subplot(1,3,1);hold on;
fillyy(tpred,Ypred(1:50)+delta_PI(1:50),Ypred(1:50)-delta_PI(1:50),[0.9 0.9 0.9])
fillyy(tpred,Ypred(1:50)+delta_CI(1:50),Ypred(1:50)-delta_CI(1:50),[0.6 0.6 0.6])
plot(tpred,Ypred(1:50),'r','LineWidth',3);
plot(tpred,data(1:50),'k','LineWidth',3);
set(gca,'FontSize',20); grid on; ylabel('$V_{LV}$ ($\mu$L)'); 
xlabel('Time (s)');

subplot(1,3,2);hold on;
fillyy(tpred,Ypred(51:100)+delta_PI(51:100),Ypred(51:100)-delta_PI(51:100),[0.9 0.9 0.9])
fillyy(tpred,Ypred(51:100)+delta_CI(51:100),Ypred(51:100)-delta_CI(51:100),[0.6 0.6 0.6])
plot(tpred,Ypred(51:100),'r','LineWidth',3);
plot(tpred,data(51:100),'k','LineWidth',3);
set(gca,'FontSize',20); grid on; ylabel('$p_{LV}$ (mmHg)'); 
xlabel('Time (s)');

subplot(1,3,3);hold on;
fillyy(tpred,Ypred(101:150)+delta_PI(101:150),Ypred(101:150)-delta_PI(101:150),[0.9 0.9 0.9])
fillyy(tpred,Ypred(101:150)+delta_CI(101:150),Ypred(101:150)-delta_CI(101:150),[0.6 0.6 0.6])
plot(tpred,Ypred(101:150),'r','LineWidth',3);
plot(tpred,data(101:150),'k','LineWidth',3);
set(gca,'FontSize',20); grid on; ylabel('$p_{Ao}$ (mmHg)'); 
xlabel('Time (s)');

h1 = plot(nan,nan,'Color',[0.9 0.9 0.9],'LineWidth',3);
h2 = plot(nan,nan,'Color',[0.6 0.6 0.6],'LineWidth',3);
h3 = plot(nan,nan,'Color','r','LineWidth',3);
h4 = plot(nan,nan,'Color','k','LineWidth',3);

legend([h4 h3 h2 h1],{'Data','Optimum','CI','PI'});
set(gcf,'Position',[153.00        233.80       1153.60        539.20])
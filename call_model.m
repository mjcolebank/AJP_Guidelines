% Call mouse model and return values

function yout = call_model(param,IC,tspace)
model_shift = -2; % This is an offset between the model and the data. It 
% corresponds to a slight delay in systolic contraction inherent in the model versus data

% Definie initial conditions, time variables, and ODE options
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Necessary for the ODE solver
T = param(10); % Cardiac Cycle Length

% Now solve our system of ODEs
y = ode45(@LV_Ao,[tspace(1), tspace(end)],IC,options,param);
tplot = linspace(28*T,30*T,99);
yout = deval(y,tplot);


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


plv = circshift(plv,model_shift);
Vlv = circshift(Vlv,model_shift);
pao = circshift(pao,model_shift);

plv = plv(1:50);
Vlv = Vlv(1:50);
pao = pao(1:50);

yout = [Vlv;plv;pao];


end
% Convert data into a parameter vector for a lumped parameter model

function param = mouse_parameters_newm1(VLV,PLV,PSA,T)
SV = max(VLV) - min(VLV);
HR = 60./T;
CO = SV.*HR./60;
PP_SA = max(PSA) - min(PSA);
%%
% The mathematical model consists of a fixed, left atrial pressure source,
% a dynamics left ventricle (modeled by an elastance function), and a
% compliance system artery, with a constant, systemic capillary pressure
% source.



Vlv_d    = 5;   % LV dead volume (ml)

% Now define the parameters for the system
P_LA     = 5;%min(PLV).*1.2;    % left atrial pressure (mmHg)
P_SysCap = 20;%mean(PSA).*0.8;   % systemic capillary pressure (mmHg)
Rmv      = 5e-3; % resistance in the mitral valve (mmHg s/micro l)
Rav      = 1e-3; % resistance in the aortic valve (mmHg s/micro l)
Rart     = (mean(PSA)-P_SysCap)./CO;  % resistance of the systemic arteries/arterioles (mmHg s/ml)
Emax     = max(PLV)./(min(VLV)-Vlv_d);  % End systolic elastance (mmHg/ml)
Emin     = min(PLV)./max(VLV);% Nonlinear elastance term (1/ml) 
T_peak   = 0.4.*T; % peak elastance (s)
T_relax  = 0.7.*T;  % end of systole (s)
Cao      = SV./PP_SA;  % aortic compliance (ml s / mmHg)

% Stack all the parameters into a vector
param = [P_LA,P_SysCap,Rmv,Rav,Rart,Emax,Emin,T_peak,T_relax,T,Vlv_d,Cao];

end
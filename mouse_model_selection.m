% Compare models for model selection and parameter variability
clear; clc; close all;
AIC_vals = zeros(5,1);
for i=1:5
    if i==1
        load mouseNEW_opt_100samp_8par_withVd_shift2.mat J_save Xopt_save ids jacobian_save par0_save
    elseif i==2
        load mouseNEW_opt_100samp_7par_shift2 J_save Xopt_save ids jacobian_save par0_save
    elseif i==3
        load mouseNEW_opt_100samp_6par_noRMV_shift2.mat J_save Xopt_save ids jacobian_save par0_save
    elseif i==4
        load mouseNEW_opt_100samp_6par_noCao_shift2.mat J_save Xopt_save ids jacobian_save par0_save
    else
        load mouseNEW_opt_100samp_6par_noEes_shift2.mat J_save Xopt_save ids jacobian_save par0_save
    end
    num_par = length(ids);
    [Jbest,id_best] = sort(J_save);
    AIC_vals(i) = 2*length(ids) + Jbest(1);
    mu_par = mean(Xopt_save(id_best(1:20),:));
    std_par = std(Xopt_save(id_best(1:20),:));
    CoV = std_par./mu_par;
    CoV0 = std(par0_save)./mean(par0_save);
    disp([mu_par; std_par; CoV; CoV0]')
end










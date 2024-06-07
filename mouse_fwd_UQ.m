% Residual vector

function yout = mouse_fwd_UQ(q,ids,param,IC,tspace)
    param(ids) = q;
    yout = call_model(param,IC,tspace)';
    yout = yout(:);
end
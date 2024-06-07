% Residual vector

function res = mouse_residual(q,data,ids,param,IC,tspace,s2vec)
% q = exp(q);
% q
param(ids) = q;
yout = call_model(param,IC,tspace);
yout = yout';
 
% res = (data(:) - yout(:))./sqrt(length(data(:)));

% data = [data(1:51) data(52:102) data(103:end)]./max(yout); yout = yout./max(yout);
% res = (data(:) - yout(:));%./sqrt(length(data(:)));

% if length(s2vec)==3
%     data = [data(1:50) data(51:100) data(101:end)]./s2vec; yout = yout./s2vec;
% else
%     data = data(:)./s2vec(:); yout = yout(:)./s2vec(:);
% ends
res = (data(:) - yout(:));%./sqrt(length(data(:)));

end
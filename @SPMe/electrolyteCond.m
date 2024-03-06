
function [kappa,varargout] = electrolyteCond(mean_electrolyte_concentration)
%% Electrolyte Conductivity Function: kappa(c_e) [1/Ohms*m]
%   Created July 12, 2011 by Scott Moura

% Identified function from Joel Forman's JPS ParamID paper
% xx = 0:1000:4000;
% yy = [1.050e-1, 1.760e-1, 2.190e-1, 8.166e-2, 3.014e-2];

% kappa = spline(xx*1e-6,yy,c_e) * 1e-2;
% kappa = ones(size(c_e)) * 0.1330 * 1e-2;
% kappa = 0.0370 * ones(size(c_e));

% From DUALFOIL LiPF6 in EC:DMC, Capiaglia et al. 1999
kappa = 0.0911+1.9101*mean_electrolyte_concentration/1e3 - 1.052*(mean_electrolyte_concentration/1e3).^2 + 0.1554*(mean_electrolyte_concentration/1e3).^3;

if(nargout == 2)
    dkappa = 1.9101/1e3 - 2*1.052*mean_electrolyte_concentration/1e3/1e3 + 0.1554*3*(mean_electrolyte_concentration/1e3)^2/1e3;
    varargout{1} = dkappa;
end
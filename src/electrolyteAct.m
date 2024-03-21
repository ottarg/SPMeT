function [dActivity,varargout] = electrolyteAct(model,electrolyte_concentration,temperature)
%% Electrolyte Activity Coefficient Function: dlnf/dln(electrolyte_concentration)
%   Created Oct 25, 2016 by Saehong Park

% From LiPF6, Valoen et al. 2005
% Fig.6 in the paper

% DataFitting Coefficients

v00 = 0.601;
v01 = 0;
v10 = -0.24;
v11 = 0;
v20 = 0;
v21 = 0;
v30 = 0.982;
v31 = -0.0052;

electrolyte_concentration = electrolyte_concentration/1000; % UnitConversion: 1 mol/L -> 1000 mol/m^3

dActivity = ((v00 + v10.*(electrolyte_concentration).^(0.5) + v30*(1+v31*(temperature - model.nominal_temperature)) .* (electrolyte_concentration).^(1.5))./(1-model.transference_number))-1;

if(nargout == 2)
    d_dactivity = (0.5 * v10 * (electrolyte_concentration).^(-0.5) + (1.5)*v30*(1+v31*(temperature - model.nominal_temperature))*(electrolyte_concentration).^(0.5))/(1-model.transference_number);
    varargout{1} = d_dactivity;
end
end
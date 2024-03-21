function [diffusion_coefficient,varargout] = electrolyte_diffusion_coefficient(electrolyte_concentration)
%% Electrolyte Diffusion Coefficient Function: D_e(c_e) [m^2/s]
%   Created July 12, 2011 by Scott Moura

% From DUALFOIL LiPF6 in EC:DMC, Capiglia et al. 1999
diffusion_coefficient = 5.34e-10*exp(-0.65*electrolyte_concentration/1e3);

if(nargout == 2)
    diffusion_coefficient_dot = -0.65*diffusion_coefficient/1e3;
    varargout{1} = diffusion_coefficient_dot;
end
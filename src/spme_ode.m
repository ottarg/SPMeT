function [x_dot,varargout] = spme_ode(obj,t,x,data)

%% Shorthand
% activation_energy = obj.activation_energy;
anode = obj.anode;
cathode = obj.cathode;
%% Parse Input Data
% Parse and interpolate current
current     = interp1(data.time,data.current,t,[]);
temperature = interp1(data.time,data.temperature,t,[])+273.15;
% Parse states
anode_solid_concentration   = x(1:(obj.discretization.radial_divisions-1));
cathode_solid_concentration = x(obj.discretization.radial_divisions : 2*(obj.discretization.radial_divisions-1));
electrolyte_concentration   = x(2*obj.discretization.radial_divisions-1 : 2*obj.discretization.radial_divisions-1+obj.discretization.Nx-4);
delta_sei = x(end);

%% Pre-calculations with current states

%%% MOLAR FLUXES
% Compute total molar flux
jn_tot =  current/(F*  anode.specific_interfacial_area*obj.electrode_area*  anode.electrode_thickness);
jp_tot = -current/(F*cathode.specific_interfacial_area*obj.electrode_area*cathode.electrode_thickness);

%%% SOLID PHASE DYNAMICS
% Solid phase diffusivity temperature dependence
anode_diffusion_coefficient   = anode.diffusion_coefficient * exp(anode.diffusion_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));
cathode_diffusion_coefficient = cathode.diffusion_coefficient * exp(cathode.diffusion_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));

% Construct (A,B) matrices for solid-phase Li diffusion
obj = initialize_solid_phase_matrices(obj,anode_diffusion_coefficient,cathode_diffusion_coefficient);
% Compute surface concentrations
anode_solid_surface_concentration = obj.solid_phase_matrices.C_n*anode_solid_concentration + obj.solid_phase_matrices.D_n*jn_tot;
cathode_solid_surface_concentration = obj.solid_phase_matrices.C_p*cathode_solid_concentration + obj.solid_phase_matrices.D_p*jp_tot;
% Remark: I am cheating slightly here. jn_tot should be jn, but doing so
% imposes an algebraic equation. This forms a DAE. I am going to
% approximate jn by jn_tot, which should be ok, since jn and jn_tot have a
% error on the order of 0.001%


%%% ELECTROLYTE PHASE DYNAMICS
% Compute electrolyte Boundary Conditions
c_e_bcs = obj.electrolyte_matrices.C * electrolyte_concentration;

ce0n = c_e_bcs(1);
cens = c_e_bcs(2);
cesp = c_e_bcs(3);
ce0p = c_e_bcs(4);

% Separate and aggregate electrolyte concentration
anode_electrolyte_concentration = electrolyte_concentration(1:(obj.discretization.Nxn-1));
separator_electrolyte_concentration = electrolyte_concentration((obj.discretization.Nxn-1)+1:(obj.discretization.Nxn-1)+(obj.discretization.Nxs-1));
cathode_electrolyte_concentration = electrolyte_concentration((obj.discretization.Nxn-1)+obj.discretization.Nxs : end);
electrolyte_concentrations = [ce0n; anode_electrolyte_concentration; cens; separator_electrolyte_concentration; cesp; cathode_electrolyte_concentration; ce0p];


%% Voltage output

% Average electrolyte concentrations
mean_electrolyte_concentration_anode = mean(electrolyte_concentrations(1:obj.discretization.Nxn+1,:));
mean_electrolyte_concentration_separator = mean(electrolyte_concentrations((obj.discretization.Nxn+1):(obj.discretization.Nxn+obj.discretization.Nxs+1),:));
mean_electrolyte_concentration_cathode = mean(electrolyte_concentrations((obj.discretization.Nxn+obj.discretization.Nxs+1):(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp+1),:));

% Overpotentials due to electrolyte subsystem
nominal_electrolyte_conductivity_anode = electrolyteCond(mean_electrolyte_concentration_anode);
nominal_electrolyte_conductivity_separator = electrolyteCond(mean_electrolyte_concentration_separator);
nominal_electrolyte_conductivity_cathode = electrolyteCond(mean_electrolyte_concentration_cathode);

% Adjustment for Arrhenius temperature dependence
electrolyte_conductivity_anode = nominal_electrolyte_conductivity_anode * exp(obj.electrolyte.conductivity_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));
electrolyte_conductivity_separator = nominal_electrolyte_conductivity_separator * exp(obj.electrolyte.conductivity_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));
electrolyte_conductivity_cathode = nominal_electrolyte_conductivity_cathode * exp(obj.electrolyte.conductivity_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));

% Bruggeman relationships
kap_n_eff = electrolyte_conductivity_anode * anode.volume_fraction_electrolyte.^(obj.bruggemann_porosity);
kap_s_eff = electrolyte_conductivity_separator * obj.separator.volume_fraction_electrolyte.^(obj.bruggemann_porosity);
kap_p_eff = electrolyte_conductivity_cathode * cathode.volume_fraction_electrolyte.^(obj.bruggemann_porosity);

% Activity coefficient
dfca_n = electrolyteAct(obj,mean_electrolyte_concentration_anode,temperature);
dfca_s = electrolyteAct(obj,mean_electrolyte_concentration_separator,temperature);
dfca_p = electrolyteAct(obj,mean_electrolyte_concentration_cathode,temperature);

% Kinetic reaction rate, adjusted for Arrhenius temperature dependence
anode_reaction_rate = anode.reaction_rate * exp(anode.intercalation_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));
cathode_reaction_rate = cathode.reaction_rate * exp(cathode.intercalation_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));

% Stochiometric Concentration Ratio
theta_n = anode_solid_surface_concentration / anode.maximum_concentration;
theta_p = cathode_solid_surface_concentration / cathode.maximum_concentration;

% Equilibrium Potential
anode_reference_potential = refPotentialAnode(theta_n);
cathode_reference_potential = refPotentialCathode(theta_p);

% Exchange current density
mean_electrolyte_concentrations = [mean_electrolyte_concentration_anode; mean_electrolyte_concentration_separator; mean_electrolyte_concentration_cathode];
[i_0n,i_0p] = exch_cur_dens(obj,cathode_reaction_rate,anode_reaction_rate,anode_solid_surface_concentration,cathode_solid_surface_concentration,mean_electrolyte_concentrations);
% Overpotentials
RTaF=(R*temperature)/(obj.charge_transfer_coefficient*F);
anode_overpotential = RTaF * asinh(current / (2*anode.specific_interfacial_area*obj.electrode_area*anode.electrode_thickness*i_0n(1)));
cathode_overpotential = RTaF * asinh(-current / (2*cathode.specific_interfacial_area*obj.electrode_area*cathode.electrode_thickness*i_0p(end)));

% Total resistance (film + growing SEI layer)
anode_resistivity = anode.sei_resistivity + delta_sei/obj.side_reaction_product.conductivity;
anode_resistance = anode_resistivity/(anode.specific_interfacial_area*anode.electrode_thickness*obj.electrode_area);
cathode_resistivity = cathode.sei_resistivity + 0;
cathode_resistance = cathode_resistivity/(cathode.specific_interfacial_area*cathode.electrode_thickness*obj.electrode_area);
electrode_resistance = anode_resistance + cathode_resistance;
% SPM Voltage (i.e. w/o electrolyte concentration terms)
anode_potential = anode_overpotential + anode_reference_potential;
cathode_potential = cathode_overpotential + cathode_reference_potential;

%Voltage excluding electrolyte
V_noVCE = cathode_potential - anode_potential - electrode_resistance*current;

% Overpotential due to electrolyte conductivity
V_electrolyteCond = (anode.electrode_thickness/(2*kap_n_eff)...
    + 2*obj.separator.thickness/(2*kap_s_eff)...
    + cathode.electrode_thickness/(2*kap_p_eff))*current;

% Overpotential due to electrolyte polarization
V_electrolytePolar = (2*R*temperature)/(F) * (1-obj.transference_number)* ...
    ( (1+dfca_n) * (log(cens) - log(ce0n)) ...
    +(1+dfca_s) * (log(cesp) - log(cens)) ...
    +(1+dfca_p) * (log(ce0p) - log(cesp)));

% Add 'em up!
V = V_noVCE + V_electrolyteCond + V_electrolytePolar;


%% Aging Dynamics

%   SEI Layer Growth model
%   Eqns Adopted from Ramadass et al (2004) [Univ of South Carolina]
%   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
%   DOI: 10.1149/1.1634273
%   NOTE1: This model has NOT been validated experimentally by eCAL
%   NOTE2: We assume this submodel only applies to anode

% Difference btw solid and electrolyte overpotential [V]
phi_se = anode_overpotential + anode_reference_potential + F*anode_resistivity*jn_tot;

% Side rxn overpotential [V]
eta_s = phi_se - obj.side_reaction_product.reference_potential - F*anode_resistivity * jn_tot;

% Molar flux of side rxn [mol/s-m^2]
j_s = -obj.side_reaction_product.exchange_current_density/...
    F * exp((-obj.charge_transfer_coefficient*F)/(R*temperature)*eta_s);

% SEI layer growth model [m/s]
delta_sei_dot = -obj.side_reaction_product.molecular_weight/...
    (obj.side_reaction_product.mass_density) * j_s;

% Molar flux of intercalation
jn = (abs(jn_tot) - abs(j_s)) * sign(jn_tot);
jp = jp_tot;

%% Solid Phase Dynamics
c_s_n_dot = obj.solid_phase_matrices.A_n*anode_solid_concentration + obj.solid_phase_matrices.B_n*jn;
c_s_p_dot = obj.solid_phase_matrices.A_p*cathode_solid_concentration + obj.solid_phase_matrices.B_p*jp;

%% Electrolyte Dynamics

% Compute Electrolyte Diffusion Coefficient and Derivative
[D_en0,dD_en0] = electrolyteDe(anode_electrolyte_concentration);
[D_es0,dD_es0] = electrolyteDe(separator_electrolyte_concentration);
[D_ep0,dD_ep0] = electrolyteDe(cathode_electrolyte_concentration);

% Adjustment for Arrhenius temperature dependence
Arrh_De = exp(obj.electrolyte.diffusion_activation_energy/R*(1/obj.nominal_temperature - 1/temperature));
D_en = D_en0 * Arrh_De;
D_es = D_es0 * Arrh_De;
D_ep = D_ep0 * Arrh_De;
dD_en = dD_en0 * Arrh_De;
dD_es = dD_es0 * Arrh_De;
dD_ep = dD_ep0 * Arrh_De;

% Apply BRUGGEMAN RELATION
D_en_eff = D_en .* anode.volume_fraction_electrolyte.^(obj.bruggemann_porosity-1);
dD_en_eff = dD_en .* anode.volume_fraction_electrolyte.^(obj.bruggemann_porosity-1);

D_es_eff = D_es .* obj.separator.volume_fraction_electrolyte.^(obj.bruggemann_porosity-1);
dD_es_eff = dD_es .* obj.separator.volume_fraction_electrolyte.^(obj.bruggemann_porosity-1);

D_ep_eff = D_ep .* cathode.volume_fraction_electrolyte.^(obj.bruggemann_porosity-1);
dD_ep_eff = dD_ep .* cathode.volume_fraction_electrolyte.^(obj.bruggemann_porosity-1);

% Compute derivative
anode_electrolyte_concentration_dot = dD_en_eff.*(obj.electrolyte_matrices.M1n*anode_electrolyte_concentration + obj.electrolyte_matrices.M2n*c_e_bcs(1:2)).^2 ...
    + D_en_eff.*(obj.electrolyte_matrices.M3n*anode_electrolyte_concentration + obj.electrolyte_matrices.M4n*c_e_bcs(1:2)) + diag(obj.electrolyte_matrices.M5n)*jn;

separator_electrolyte_concentration_dot = dD_es_eff.*(obj.electrolyte_matrices.M1s*separator_electrolyte_concentration + obj.electrolyte_matrices.M2s*c_e_bcs(2:3)).^2 ...
    + D_es_eff.*(obj.electrolyte_matrices.M3s*separator_electrolyte_concentration + obj.electrolyte_matrices.M4s*c_e_bcs(2:3));

cathode_electrolyte_concentration_dot = dD_ep_eff.*(obj.electrolyte_matrices.M1p*cathode_electrolyte_concentration + obj.electrolyte_matrices.M2p*c_e_bcs(3:4)).^2 ...
    + D_ep_eff.*(obj.electrolyte_matrices.M3p*cathode_electrolyte_concentration + obj.electrolyte_matrices.M4p*c_e_bcs(3:4)) + diag(obj.electrolyte_matrices.M5p)*jp;

% Assemble c_e_dot
c_e_dot = [anode_electrolyte_concentration_dot; separator_electrolyte_concentration_dot; cathode_electrolyte_concentration_dot];

%% Thermal Dynamics

% State-of-Charge (Bulk)
r_vec = (0:obj.discretization.delta_r_n:1)';
c_n = [anode_solid_concentration(1); anode_solid_concentration; anode_solid_surface_concentration];
c_p = [cathode_solid_concentration(1); cathode_solid_concentration; cathode_solid_surface_concentration];
anode_SOC = 3/anode.maximum_concentration * trapz(r_vec,r_vec.^2.*c_n);
cathode_SOC = 3/cathode.maximum_concentration * trapz(r_vec,r_vec.^2.*c_p);

% Equilibrium potentials
[anode_OCV] = refPotentialAnode(anode_SOC);
[cathode_OCV] = refPotentialCathode(cathode_SOC);
OCV = cathode_OCV - anode_OCV;

%% Concatenate time derivatives
x_dot = [c_s_n_dot; c_s_p_dot; c_e_dot; delta_sei_dot];

%% Concatenate outputs
varargout{1} = V;
varargout{2} = V_noVCE;
varargout{3} = anode_SOC;
varargout{4} = cathode_SOC;
varargout{5} = anode_solid_surface_concentration;
varargout{6} = cathode_solid_surface_concentration;
varargout{7} = electrolyte_concentrations';
varargout{8} = OCV;
varargout{9} = phi_se;
varargout{10} = cathode_potential;

end
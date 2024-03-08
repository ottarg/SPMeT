function [value, isterminal, direction] = detectImagSolution(obj, t, x, data)

%% This is an "Event" method for the ODE23s solver, which detects if any of these calculated variables has become imaginary and if so, exits the simulation. 
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
jn_tot =  current/(SPMe().F*  anode.specific_interfacial_area*obj.electrode_area*  anode.electrode_thickness);
jp_tot = -current/(SPMe().F*cathode.specific_interfacial_area*obj.electrode_area*cathode.electrode_thickness);

%%% SOLID PHASE DYNAMICS
% Solid phase diffusivity temperature dependence
anode_diffusion_coefficient   = anode.diffusion_coefficient * exp(anode.diffusion_activation_energy/SPMe().R*(1/obj.nominal_temperature - 1/temperature));
cathode_diffusion_coefficient = cathode.diffusion_coefficient * exp(cathode.diffusion_activation_energy/SPMe().R*(1/obj.nominal_temperature - 1/temperature));

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
nominal_electrolyte_conductivity_anode = SPMe().electrolyteCond(mean_electrolyte_concentration_anode);
nominal_electrolyte_conductivity_separator = SPMe().electrolyteCond(mean_electrolyte_concentration_separator);
nominal_electrolyte_conductivity_cathode = SPMe().electrolyteCond(mean_electrolyte_concentration_cathode);

% Adjustment for Arrhenius temperature dependence
electrolyte_conductivity_anode = nominal_electrolyte_conductivity_anode * exp(obj.electrolyte.conductivity_activation_energy/SPMe().R*(1/obj.nominal_temperature - 1/temperature));
electrolyte_conductivity_separator = nominal_electrolyte_conductivity_separator * exp(obj.electrolyte.conductivity_activation_energy/SPMe().R*(1/obj.nominal_temperature - 1/temperature));
electrolyte_conductivity_cathode = nominal_electrolyte_conductivity_cathode * exp(obj.electrolyte.conductivity_activation_energy/SPMe().R*(1/obj.nominal_temperature - 1/temperature));

% Bruggeman relationships
kap_n_eff = electrolyte_conductivity_anode * anode.volume_fraction_electrolyte.^(obj.bruggemann_porosity);
kap_s_eff = electrolyte_conductivity_separator * obj.separator.volume_fraction_electrolyte.^(obj.bruggemann_porosity);
kap_p_eff = electrolyte_conductivity_cathode * cathode.volume_fraction_electrolyte.^(obj.bruggemann_porosity);

% Activity coefficient
dfca_n = SPMe().electrolyteAct(obj,mean_electrolyte_concentration_anode,temperature);
dfca_s = SPMe().electrolyteAct(obj,mean_electrolyte_concentration_separator,temperature);
dfca_p = SPMe().electrolyteAct(obj,mean_electrolyte_concentration_cathode,temperature);

% Kinetic reaction rate, adjusted for Arrhenius temperature dependence
anode_reaction_rate = anode.reaction_rate * exp(anode.intercalation_activation_energy/SPMe().R*(1/obj.nominal_temperature - 1/temperature));
cathode_reaction_rate = cathode.reaction_rate * exp(cathode.intercalation_activation_energy/SPMe().R*(1/obj.nominal_temperature - 1/temperature));

% Stochiometric Concentration Ratio
theta_n = anode_solid_surface_concentration / anode.maximum_concentration;
theta_p = cathode_solid_surface_concentration / cathode.maximum_concentration;

value = ~isreal([cathode_solid_surface_concentration,anode_solid_surface_concentration,jn_tot,jp_tot,anode_diffusion_coefficient,cathode_diffusion_coefficient,electrolyte_concentrations',theta_n,theta_p,cathode_reaction_rate,anode_reaction_rate,dfca_p,dfca_s,dfca_n,kap_p_eff,kap_s_eff,kap_n_eff,electrolyte_conductivity_cathode,electrolyte_conductivity_separator,electrolyte_conductivity_anode,nominal_electrolyte_conductivity_cathode,nominal_electrolyte_conductivity_separator,nominal_electrolyte_conductivity_anode,mean_electrolyte_concentration_cathode,mean_electrolyte_concentration_separator,mean_electrolyte_concentration_anode]);
isterminal = 1;   % Stop the integration
direction  = 0;
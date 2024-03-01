function [value, isterminal, direction] = myEvent(t, x,obj, data)

% Parse and interpolate current
cur = interp1(data.time,data.cur(1:length(data.time)),t,[]);

% Parse states
anode_solid_concentration = x(1:(obj.discretization.Nr-1));
cathode_solid_concentration = x(obj.discretization.Nr : 2*(obj.discretization.Nr-1));
electrolyte_concentration = x(2*obj.discretization.Nr-1 : 2*obj.discretization.Nr-1+obj.discretization.Nx-4);
T1 = x(end-2);
T2 = x(end-1);
delta_sei = x(end);

%% Pre-calculations with current states

%%% MOLAR FLUXES
% Compute total molar flux
jn_tot = cur/(faraday*obj.cell_properties.anode.specific_interfacial_area*obj.cell_properties.electrode_area*obj.cell_properties.anode.electrode_thickness);
jp_tot = -cur/(faraday*obj.cell_properties.cathode.specific_interfacial_area*obj.cell_properties.electrode_area*obj.cell_properties.cathode.electrode_thickness);

%%% SOLID PHASE DYNAMICS
% Solid phase diffusivity temperature dependence
anode_diffusion_coefficient = obj.cell_properties.anode.diffusion_coefficient * exp(obj.cell_properties.E.Dsn/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));
cathode_diffusion_coefficient = obj.cell_properties.cathode.diffusion_coefficient * exp(obj.cell_properties.E.Dsp/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));

% Construct (A,B) matrices for solid-phase Li diffusion
initialize_solid_phase_matrices(obj,anode_diffusion_coefficient,cathode_diffusion_coefficient)
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
c_en = electrolyte_concentration(1:(obj.discretization.Nxn-1));
c_es = electrolyte_concentration((obj.discretization.Nxn-1)+1:(obj.discretization.Nxn-1)+(obj.discretization.Nxs-1));
c_ep = electrolyte_concentration((obj.discretization.Nxn-1)+obj.discretization.Nxs : end);
c_ex = [ce0n; c_en; cens; c_es; cesp; c_ep; ce0p];


%% Voltage output

% Average electrolyte concentrations
mean_electrolyte_concentration_anode = mean(c_ex(1:obj.discretization.Nxn+1,:));
mean_electrolyte_concentration_separator = mean(c_ex((obj.discretization.Nxn+1):(obj.discretization.Nxn+obj.discretization.Nxs+1),:));
mean_electrolyte_concentration_cathode = mean(c_ex((obj.discretization.Nxn+obj.discretization.Nxs+1):(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp+1),:));

% Overpotentials due to electrolyte subsystem
nominal_electrolyte_conductivity_anode = SPMe().electrolyteCond(mean_electrolyte_concentration_anode);
nominal_electrolyte_conductivity_separator = SPMe().electrolyteCond(mean_electrolyte_concentration_separator);
nominal_electrolyte_conductivity_cathode = SPMe().electrolyteCond(mean_electrolyte_concentration_cathode);

% Adjustment for Arrhenius temperature dependence
electrolyte_conductivity_anode = nominal_electrolyte_conductivity_anode * exp(obj.cell_properties.E.kappa_e/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));
electrolyte_conductivity_separator = nominal_electrolyte_conductivity_separator * exp(obj.cell_properties.E.kappa_e/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));
electrolyte_conductivity_cathode = nominal_electrolyte_conductivity_cathode * exp(obj.cell_properties.E.kappa_e/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));

% Bruggeman relationships
kap_n_eff = electrolyte_conductivity_anode * obj.cell_properties.anode.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity);
kap_s_eff = electrolyte_conductivity_separator * obj.cell_properties.separator.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity);
kap_p_eff = electrolyte_conductivity_cathode * obj.cell_properties.cathode.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity);

% Activity coefficient
dfca_n = electrolyteAct(obj,mean_electrolyte_concentration_anode,T1);
dfca_s = electrolyteAct(obj,mean_electrolyte_concentration_separator,T1);
dfca_p = electrolyteAct(obj,mean_electrolyte_concentration_cathode,T1);

% Kinetic reaction rate, adjusted for Arrhenius temperature dependence
anode_reaction_rate = obj.cell_properties.anode.reaction_rate * exp(obj.cell_properties.E.kn/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));
cathode_reaction_rate = obj.cell_properties.cathode.reaction_rate * exp(obj.cell_properties.E.kp/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));

% Stochiometric Concentration Ratio
theta_n = anode_solid_surface_concentration / obj.cell_properties.anode.maximum_concentration;
theta_p = cathode_solid_surface_concentration / obj.cell_properties.cathode.maximum_concentration;
value = ~isreal([cathode_solid_surface_concentration,anode_solid_surface_concentration,jn_tot,jp_tot,anode_diffusion_coefficient,cathode_diffusion_coefficient,c_ex',theta_n,theta_p,cathode_reaction_rate,anode_reaction_rate,dfca_p,dfca_s,dfca_n,kap_p_eff,kap_s_eff,kap_n_eff,electrolyte_conductivity_cathode,electrolyte_conductivity_separator,electrolyte_conductivity_anode,nominal_electrolyte_conductivity_cathode,nominal_electrolyte_conductivity_separator,nominal_electrolyte_conductivity_anode,mean_electrolyte_concentration_cathode,mean_electrolyte_concentration_separator,mean_electrolyte_concentration_anode]);
isterminal = 1;   % Stop the integration
direction  = 0;
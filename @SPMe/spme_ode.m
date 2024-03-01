function [x_dot,varargout] = spme_ode(obj,t,x,data)

%% Parse Input Data

% Parse and interpolate current
current = interp1(data.time,data.current(1:length(data.time)),t,[]);

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
jn_tot = current/(faraday*obj.cell_properties.anode.specific_interfacial_area*obj.cell_properties.electrode_area*obj.cell_properties.anode.electrode_thickness);
jp_tot = -current/(faraday*obj.cell_properties.cathode.specific_interfacial_area*obj.cell_properties.electrode_area*obj.cell_properties.cathode.electrode_thickness);

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

% Equilibrium Potential
anode_reference_potential = SPMe().refPotentialAnode(theta_n);
cathode_reference_potential = SPMe().refPotentialCathode(theta_p);
OCV = cathode_reference_potential - anode_reference_potential;
% Exchange current density
c_e_bar = [mean_electrolyte_concentration_anode; mean_electrolyte_concentration_separator; mean_electrolyte_concentration_cathode];
[i_0n,i_0p] = exch_cur_dens(obj,cathode_reaction_rate,anode_reaction_rate,anode_solid_surface_concentration,cathode_solid_surface_concentration,c_e_bar);
% Overpotentials
RTaF=(gas_constant*T1)/(obj.cell_properties.charge_transfer_coefficient*faraday);
anode_overpotential = RTaF * asinh(current / (2*obj.cell_properties.anode.specific_interfacial_area*obj.cell_properties.electrode_area*obj.cell_properties.anode.electrode_thickness*i_0n(1)));
cathode_overpotential = RTaF * asinh(-current / (2*obj.cell_properties.cathode.specific_interfacial_area*obj.cell_properties.electrode_area*obj.cell_properties.cathode.electrode_thickness*i_0p(end)));

% Total resistance (film + growing SEI layer)
anode_resistivity = obj.cell_properties.anode.sei_resistivity + delta_sei/obj.cell_properties.side_reaction_product.conductivity;
anode_resistance = anode_resistivity/(obj.cell_properties.anode.specific_interfacial_area*obj.cell_properties.anode.electrode_thickness*obj.cell_properties.electrode_area);
cathode_resistivity = obj.cell_properties.cathode.sei_resistivity + 0;
cathode_resistance = cathode_resistivity/(obj.cell_properties.cathode.specific_interfacial_area*obj.cell_properties.cathode.electrode_thickness*obj.cell_properties.electrode_area);
electrode_resistance = anode_resistance + cathode_resistance;
% SPM Voltage (i.e. w/o electrolyte concentration terms)
anode_potential = anode_overpotential + anode_reference_potential;
cathode_potential = cathode_overpotential + cathode_reference_potential;
V_noVCE = cathode_potential - anode_potential - electrode_resistance*current;
% cathode_overpotential - anode_overpotential + cathode_reference_potential - anode_reference_potential - (anode_resistance + cathode_resistance)*cur;

% Overpotential due to electrolyte conductivity
V_electrolyteCond = (obj.cell_properties.anode.electrode_thickness/(2*kap_n_eff) + 2*obj.cell_properties.separator.thickness/(2*kap_s_eff) + obj.cell_properties.cathode.electrode_thickness/(2*kap_p_eff))*current;

% Overpotential due to electrolyte polarization
V_electrolytePolar = (2*gas_constant*T1)/(faraday) * (1-obj.cell_properties.t_plus)* ...
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
phi_se = anode_overpotential + anode_reference_potential + faraday*anode_resistivity*jn_tot;

% Side rxn overpotential [V]
eta_s = phi_se - obj.cell_properties.side_reaction_product.reference_potential - faraday*anode_resistivity * jn_tot;

% Molar flux of side rxn [mol/s-m^2]
j_s = -obj.cell_properties.side_reaction_product.exchange_current_density/faraday * exp((-obj.cell_properties.charge_transfer_coefficient*faraday)/(gas_constant*T1)*eta_s);

% SEI layer growth model [m/s]
delta_sei_dot = -obj.cell_properties.side_reaction_product.molecular_weight/(obj.cell_properties.side_reaction_product.mass_density) * j_s;

% Molar flux of intercalation
jn = (abs(jn_tot) - abs(j_s)) * sign(jn_tot);
jp = jp_tot;

%% Solid Phase Dynamics

% ODE for c_s
c_s_n_dot = obj.solid_phase_matrices.A_n*anode_solid_concentration + obj.solid_phase_matrices.B_n*jn;
c_s_p_dot = obj.solid_phase_matrices.A_p*cathode_solid_concentration + obj.solid_phase_matrices.B_p*jp;

%% Electrolyte Dynamics

% Compute Electrolyte Diffusion Coefficient and Derivative
[D_en0,dD_en0] = SPMe().electrolyteDe(c_en);
[D_es0,dD_es0] = SPMe().electrolyteDe(c_es);
[D_ep0,dD_ep0] = SPMe().electrolyteDe(c_ep);

% Adjustment for Arrhenius temperature dependence
Arrh_De = exp(obj.cell_properties.E.De/gas_constant*(1/obj.cell_properties.nominal_temperature - 1/T1));
D_en = D_en0 * Arrh_De;
D_es = D_es0 * Arrh_De;
D_ep = D_ep0 * Arrh_De;
dD_en = dD_en0 * Arrh_De;
dD_es = dD_es0 * Arrh_De;
dD_ep = dD_ep0 * Arrh_De;

% Apply BRUGGEMAN RELATION
D_en_eff = D_en .* obj.cell_properties.anode.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity-1);
dD_en_eff = dD_en .* obj.cell_properties.anode.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity-1);

D_es_eff = D_es .* obj.cell_properties.separator.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity-1);
dD_es_eff = dD_es .* obj.cell_properties.separator.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity-1);

D_ep_eff = D_ep .* obj.cell_properties.cathode.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity-1);
dD_ep_eff = dD_ep .* obj.cell_properties.cathode.volume_fraction_electrolyte.^(obj.cell_properties.bruggemann_porosity-1);

% Compute derivative
c_en_dot = dD_en_eff.*(obj.electrolyte_matrices.M1n*c_en + obj.electrolyte_matrices.M2n*c_e_bcs(1:2)).^2 ...
    + D_en_eff.*(obj.electrolyte_matrices.M3n*c_en + obj.electrolyte_matrices.M4n*c_e_bcs(1:2)) + diag(obj.electrolyte_matrices.M5n)*jn;

c_es_dot = dD_es_eff.*(obj.electrolyte_matrices.M1s*c_es + obj.electrolyte_matrices.M2s*c_e_bcs(2:3)).^2 ...
    + D_es_eff.*(obj.electrolyte_matrices.M3s*c_es + obj.electrolyte_matrices.M4s*c_e_bcs(2:3));

c_ep_dot = dD_ep_eff.*(obj.electrolyte_matrices.M1p*c_ep + obj.electrolyte_matrices.M2p*c_e_bcs(3:4)).^2 ...
    + D_ep_eff.*(obj.electrolyte_matrices.M3p*c_ep + obj.electrolyte_matrices.M4p*c_e_bcs(3:4)) + diag(obj.electrolyte_matrices.M5p)*jp;

% Assemble c_e_dot
c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];

%% Thermal Dynamics

% State-of-Charge (Bulk)
r_vec = (0:obj.discretization.delta_r_n:1)';
c_n = [anode_solid_concentration(1); anode_solid_concentration; anode_solid_surface_concentration];
c_p = [cathode_solid_concentration(1); cathode_solid_concentration; cathode_solid_surface_concentration];
SOC_n = 3/obj.cell_properties.anode.maximum_concentration * trapz(r_vec,r_vec.^2.*c_n);
SOC_p = 3/obj.cell_properties.cathode.maximum_concentration * trapz(r_vec,r_vec.^2.*c_p);

% Equilibrium potentials
[Unb] = SPMe().refPotentialAnode(SOC_n);
[Upb] = SPMe().refPotentialCathode(SOC_p);
OCV = Upb - Unb;
% Heat generation
Qdot = -current*(V - OCV);

% Differential equations
T1_dot = (obj.cell_properties.thermal.h12 * (T2-T1) + Qdot) / obj.cell_properties.thermal.C1;
T2_dot = (obj.cell_properties.thermal.h12 * (T1-T2) + obj.cell_properties.thermal.h2a*(obj.cell_properties.ambient_temperature - T2)) / obj.cell_properties.thermal.C2;

%% Concatenate time derivatives
x_dot = [c_s_n_dot; c_s_p_dot; c_e_dot; T1_dot; T2_dot; delta_sei_dot];

%% Concatenate outputs
varargout{1} = V;
varargout{2} = V_noVCE;
varargout{3} = SOC_n;
varargout{4} = SOC_p;
varargout{5} = anode_solid_surface_concentration;
varargout{6} = cathode_solid_surface_concentration;
varargout{7} = c_ex';
varargout{8} = OCV;
varargout{9} = anode_potential;
varargout{10} = cathode_potential;

end
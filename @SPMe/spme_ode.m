function [x_dot,varargout] = spme_ode(obj,t,x,data)

%% Parse Input Data

% Parse and interpolate current
cur = interp1(data.time,data.cur,t,[]);

% Parse states
c_s_n = x(1:(obj.discretization.Nr-1));
c_s_p = x(obj.discretization.Nr : 2*(obj.discretization.Nr-1));
c_e = x(2*obj.discretization.Nr-1 : 2*obj.discretization.Nr-1+obj.discretization.Nx-4);
T1 = x(end-2);
T2 = x(end-1);
delta_sei = x(end);


%% Pre-calculations with current states

%%% MOLAR FLUXES
% Compute total molar flux
jn_tot = cur/(faraday*obj.cell_properties.a_s_n*obj.cell_properties.Area*obj.cell_properties.L_n);
jp_tot = -cur/(faraday*obj.cell_properties.a_s_p*obj.cell_properties.Area*obj.cell_properties.L_p);


%%% SOLID PHASE DYNAMICS
% Solid phase diffusivity temperature dependence
D_s_n = obj.cell_properties.D_s_n0 * exp(obj.cell_properties.E.Dsn/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
D_s_p = obj.cell_properties.D_s_n0 * exp(obj.cell_properties.E.Dsp/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));

% Construct (A,B) matrices for solid-phase Li diffusion
initialize_solid_phase_matrices(obj,D_s_n,D_s_p)
% Compute surface concentrations
c_ss_n = obj.solid_phase_matrices.C_n*c_s_n + obj.solid_phase_matrices.D_n*jn_tot;
c_ss_p = obj.solid_phase_matrices.C_p*c_s_p + obj.solid_phase_matrices.D_p*jp_tot;
% Remark: I am cheating slightly here. jn_tot should be jn, but doing so
% imposes an algebraic equation. This forms a DAE. I am going to
% approximate jn by jn_tot, which should be ok, since jn and jn_tot have a
% error on the order of 0.001%


%%% ELECTROLYTE PHASE DYNAMICS
% Compute electrolyte Boundary Conditions
c_e_bcs = obj.electrolyte_matrices.C * c_e;

ce0n = c_e_bcs(1);
cens = c_e_bcs(2);
cesp = c_e_bcs(3);
ce0p = c_e_bcs(4);

% Separate and aggregate electrolyte concentration
c_en = c_e(1:(obj.discretization.Nxn-1));
c_es = c_e((obj.discretization.Nxn-1)+1:(obj.discretization.Nxn-1)+(obj.discretization.Nxs-1));
c_ep = c_e((obj.discretization.Nxn-1)+obj.discretization.Nxs : end);
c_ex = [ce0n; c_en; cens; c_es; cesp; c_ep; ce0p];


%% Voltage output

% Average electrolyte concentrations
cen_bar = mean(c_ex(1:obj.discretization.Nxn+1,:));
ces_bar = mean(c_ex((obj.discretization.Nxn+1):(obj.discretization.Nxn+obj.discretization.Nxs+1),:));
cep_bar = mean(c_ex((obj.discretization.Nxn+obj.discretization.Nxs+1):(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp+1),:));

% Overpotentials due to electrolyte subsystem
kap_n_ref = electrolyteCond(cen_bar);
kap_s_ref = electrolyteCond(ces_bar);
kap_p_ref = electrolyteCond(cep_bar);

% Adjustment for Arrhenius temperature dependence
kap_n = kap_n_ref * exp(obj.cell_properties.E.kappa_e/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
kap_s = kap_s_ref * exp(obj.cell_properties.E.kappa_e/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
kap_p = kap_p_ref * exp(obj.cell_properties.E.kappa_e/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));

% Bruggeman relationships
kap_n_eff = kap_n * obj.cell_properties.epsilon_e_n.^(obj.cell_properties.brug);
kap_s_eff = kap_s * obj.cell_properties.epsilon_e_s.^(obj.cell_properties.brug);
kap_p_eff = kap_p * obj.cell_properties.epsilon_e_p.^(obj.cell_properties.brug);

% Activity coefficient
dfca_n = electrolyteAct(obj,cen_bar,T1);
dfca_s = electrolyteAct(obj,ces_bar,T1);
dfca_p = electrolyteAct(obj,cep_bar,T1);

% Kinetic reaction rate, adjusted for Arrhenius temperature dependence
k_n = obj.cell_properties.k_n0 * exp(obj.cell_properties.E.kn/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
k_p = obj.cell_properties.k_p0 * exp(obj.cell_properties.E.kp/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));

% Stochiometric Concentration Ratio
theta_n = c_ss_n / obj.cell_properties.c_s_n_max;
theta_p = c_ss_p / obj.cell_properties.c_s_p_max;

% Equilibrium Potential
Unref = refPotentialAnode(obj,theta_n);
Upref = refPotentialCathode(obj,theta_p);

% Exchange current density
c_e_bar = [cen_bar; ces_bar; cep_bar];
[i_0n,i_0p] = exch_cur_dens(obj,k_p,k_n,c_ss_n,c_ss_p,c_e_bar);
% Overpotentials
RTaF=(obj.cell_properties.R*T1)/(obj.cell_properties.alph*faraday);
eta_n = RTaF * asinh(cur / (2*obj.cell_properties.a_s_n*obj.cell_properties.Area*obj.cell_properties.L_n*i_0n(1)));
eta_p = RTaF * asinh(-cur / (2*obj.cell_properties.a_s_p*obj.cell_properties.Area*obj.cell_properties.L_p*i_0p(end)));

% Total resistance (film + growing SEI layer)
R_tot_n = obj.cell_properties.R_f_n + delta_sei/obj.cell_properties.kappa_P;
R_tot_p = obj.cell_properties.R_f_p + 0;

% SPM Voltage (i.e. w/o electrolyte concentration terms)
V_noVCE = eta_p - eta_n + Upref - Unref ...
    - (R_tot_n/(obj.cell_properties.a_s_n*obj.cell_properties.L_n*obj.cell_properties.Area) + R_tot_p/(obj.cell_properties.a_s_p*obj.cell_properties.L_p*obj.cell_properties.Area))*cur;

% Overpotential due to electrolyte conductivity
V_electrolyteCond = (obj.cell_properties.L_n/(2*kap_n_eff) + 2*obj.cell_properties.L_s/(2*kap_s_eff) + obj.cell_properties.L_p/(2*kap_p_eff))*cur;

% Overpotential due to electrolyte polarization
V_electrolytePolar = (2*obj.cell_properties.R*T1)/(faraday) * (1-obj.cell_properties.t_plus)* ...
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
phi_se = eta_n + Unref + faraday*R_tot_n*jn_tot;

% Side exn overpotential [V]
eta_s = phi_se - obj.cell_properties.Us - faraday*R_tot_n * jn_tot;

% Molar flux of side rxn [mol/s-m^2]
j_s = -obj.cell_properties.i0s/faraday * exp((-obj.cell_properties.alph*faraday)/(obj.cell_properties.R*T1)*eta_s);

% SEI layer growth model [m/s]
delta_sei_dot = -obj.cell_properties.M_P/(obj.cell_properties.rho_P) * j_s;

% Molar flux of intercalation
jn = (abs(jn_tot) - abs(j_s)) * sign(jn_tot);
jp = jp_tot;

%% Solid Phase Dynamics

% ODE for c_s
c_s_n_dot = obj.solid_phase_matrices.A_n*c_s_n + obj.solid_phase_matrices.B_n*jn;
c_s_p_dot = obj.solid_phase_matrices.A_p*c_s_p + obj.solid_phase_matrices.B_p*jp;


%% Electrolyte Dynamics

% Compute Electrolyte Diffusion Coefficient and Derivative
[D_en0,dD_en0] = electrolyteDe(c_en);
[D_es0,dD_es0] = electrolyteDe(c_es);
[D_ep0,dD_ep0] = electrolyteDe(c_ep);

% Adjustment for Arrhenius temperature dependence
Arrh_De = exp(obj.cell_properties.E.De/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
D_en = D_en0 * Arrh_De;
D_es = D_es0 * Arrh_De;
D_ep = D_ep0 * Arrh_De;
dD_en = dD_en0 * Arrh_De;
dD_es = dD_es0 * Arrh_De;
dD_ep = dD_ep0 * Arrh_De;

% Apply BRUGGEMAN RELATION
D_en_eff = D_en .* obj.cell_properties.epsilon_e_n.^(obj.cell_properties.brug-1);
dD_en_eff = dD_en .* obj.cell_properties.epsilon_e_n.^(obj.cell_properties.brug-1);

D_es_eff = D_es .* obj.cell_properties.epsilon_e_s.^(obj.cell_properties.brug-1);
dD_es_eff = dD_es .* obj.cell_properties.epsilon_e_s.^(obj.cell_properties.brug-1);

D_ep_eff = D_ep .* obj.cell_properties.epsilon_e_p.^(obj.cell_properties.brug-1);
dD_ep_eff = dD_ep .* obj.cell_properties.epsilon_e_p.^(obj.cell_properties.brug-1);

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
c_n = [c_s_n(1); c_s_n; c_ss_n];
c_p = [c_s_p(1); c_s_p; c_ss_p];
SOC_n = 3/obj.cell_properties.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n);
SOC_p = 3/obj.cell_properties.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p);

% Equilibrium potentials
[Unb] = refPotentialAnode(obj, SOC_n);
[Upb] = refPotentialCathode(obj, SOC_p);

% Heat generation
% disp(cur)
% disp(V)
% disp((Upb - Unb))
% pause;
Qdot = -cur*(V - (Upb - Unb));

% Differential equations
T1_dot = (obj.cell_properties.thermal.h12 * (T2-T1) + Qdot) / obj.cell_properties.thermal.C1;
T2_dot = (obj.cell_properties.thermal.h12 * (T1-T2) + obj.cell_properties.thermal.h2a*(obj.cell_properties.T_amb - T2)) / obj.cell_properties.thermal.C2;


%% Concatenate time derivatives
x_dot = [c_s_n_dot; c_s_p_dot; c_e_dot; T1_dot; T2_dot; delta_sei_dot];

%% Concatenate outputs
varargout{1} = V;
varargout{2} = V_noVCE;
varargout{3} = SOC_n;
varargout{4} = SOC_p;
varargout{5} = c_ss_n;
varargout{6} = c_ss_p;
varargout{7} = c_ex';


end
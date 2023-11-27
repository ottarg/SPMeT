%% Params for Electrochemical Model
%   Created May 24, 2012 by Scott Moura
%
%   Most parameters are from DUALFOIL model
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3
%   Equilibrium potentials are from Bosch Klein TCST 2011
%   Electrode area chosen to correspond to 2.3 Ah cell
%   Thermal parameters are adopted from L. Zhang (HIT) 2014 10.1016/j.jpowsour.2014.07.110
%   NOTE: Paramters represent an illustrative cell

function p = params_LCO()
    %% Geometric Params
    % Thickness of each layer
    p.anode.electrode_thickness = 100e-6;
    p.separator.thickness = 25e-6;
    p.cathode.electrode_thickness = 100e-6;

    % Particle Radii
    p.anode.particle_radius = 10e-6;   % Radius of solid particles in negative electrode [m]
    p.cathode.particle_radius = 10e-6;   % Radius of solid particles in positive electrode [m]
    
    % Volume fractions
    p.anode.volume_fraction_solid = 0.6;      % Volume fraction in solid for neg. electrode
    p.cathode.volume_fraction_solid = 0.5;      % Volume fraction in solid for pos. electrode
    
    p.anode.volume_fraction_electrolyte = 0.3;   % Volume fraction in electrolyte for neg. electrode
    p.separator.volume_fraction_electrolyte = 1.0;   % Volume fraction in electrolyte for separator
    p.cathode.volume_fraction_electrolyte = 0.3;   % Volume fraction in electrolyte for pos. electrode
    
    % make element to calculate phi_{s} by Saehong Park 
    p.anode.volume_fraction_filler = 1-p.anode.volume_fraction_solid-p.anode.volume_fraction_electrolyte;  % Volume fraction of filler in neg. electrode
    p.cathode.volume_fraction_filler = 1-p.cathode.volume_fraction_solid-p.cathode.volume_fraction_electrolyte;  % Volume fraction of filler in pos. electrode
        
    % Specific interfacial surface area
    p.anode.specific_interfacial_area = 3*p.anode.volume_fraction_solid / p.anode.particle_radius;  % Negative electrode [m^2/m^3]
    p.cathode.specific_interfacial_area = 3*p.cathode.volume_fraction_solid / p.cathode.particle_radius;  % Positive electrode [m^2/m^3]
    
    %% Transport Params
    % Diffusion coefficient in solid
    p.anode.diffusion_coefficient = 3.9e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
    p.cathode.diffusion_coefficient = 1e-13;  % Diffusion coeff for solid in pos. electrode, [m^2/s]
    
    p.bruggemann_porosity = 1.5;       % Bruggeman porosity
    
%     % Conductivity of solid
%     p.anode.solid_conductivity = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
%     p.cathode.solid_conductivity = 10;    % Conductivity of solid in pos. electrode, [1/Ohms*m]
%     
    % Miscellaneous
    p.t_plus = 0.4;       % Transference number
    p.electrode_area = 1;           % Electrode current collector area [m^2]
    
    %% Kinetic Params
%     p.gas_constant = 8.314472;       % Gas constant, [J/mol-K]
    
    p.charge_transfer_coefficient = 0.5; % Charge transfer coefficients
    
    p.anode.sei_resistivity = 1e-3;    % Resistivity of SEI layer, [Ohms*m^2]
    p.cathode.sei_resistivity = 0;     % Resistivity of SEI layer, [Ohms*m^2]
    p.non_jellyroll_resistance = 0;         % Contact Resistance/Current Collector Resistance, [Ohms-m^2]
    
    % Nominal Reaction rates
    p.anode.reaction_rate = 1e-5;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
    p.cathode.reaction_rate = 3e-7; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
    
    %% Thermodynamic Params
    
    % Ambient Temperature
    p.ambient_temperature = 298.15; % [K]
    
    % Activation Energies
    % Taken from Zhang et al (2014) [Harbin]
    % http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
    % All units are [J/mol]
    p.E.kn = 37.48e3;
    p.E.kp = 39.57e3;
    p.E.Dsn = 42.77e3;
    p.E.Dsp = 18.55e3;
    p.E.De = 37.04e3;
    p.E.kappa_e = 34.70e3;
    
    % Reference temperature
    p.nominal_temperature = 298.15; %[K]
    
    % Heat transfer parameters
    % Taken from Zhang et al (2014) [Harbin]
    % http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
    p.thermal.C1 = 62.7;    % [J/K]
    p.thermal.C2 = 4.5;     % [J/K]
    p.thermal.h12 = 10; %1.9386; % [W/K]
    p.thermal.h2a = 21.45;  % [W/K]
    
    %% Aging submodel parameters
    
    %   SEI Layer Growth model
    %   Adopted from Ramadass et al (2004) [Univ of South Carolina]
    %   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
    %   DOI: 10.1149/1.1634273
    %   NOTE: These parameters have NOT been experimentally validated by eCAL
    
    p.kappa_P = 1;      % [S/m] conductivity of side rxn product
    p.M_P = 7.3e1;      % [kg/mol] molecular weight of side rxn product
    p.rho_P = 2.1e3;    % [kg/m^3] mass density of side rxn product
    p.i0s = 0; %1.5e-6;     % [A/m^2] exchange current density of side rxn
    p.Us = 0.4;         % [V] reference potential of side rxn
    
    %% Concentrations
    % Maxima based on DUALFOIL 
    % line 588 in DUALFOIL Fortran code
    
    p.c_s_n_max = 3.6e3 * 372 * 1800 / faraday;   % Max concentration in anode, [mol/m^3]
    %p.c_s_n_max = 3.6e3 * 372 * 2260 / p.Faraday;   % Max concentration in anode, [mol/m^3]
    
    %p.c_s_p_max = 3.6e3 * 247 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]
    p.c_s_p_max = 3.6e3 * 274 * 5010 / faraday;    % Max concentration in cathode, [mol/m^3]
    
    p.n_Li_s = 2.5; %2.781;        % Total moles of lithium in solid phase [mol]
    p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]
    
    %% Cutoff voltages
    p.volt_max = 4.7; %4.1113; %4.7;
    p.volt_min = 3.105; %2.6;

end


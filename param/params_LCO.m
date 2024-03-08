%% Params for Electrochemical Model
%   Created May 24, 2012 by Scott Moura
%
%   Most parameters are from DUALFOIL model
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3
%   Equilibrium potentials are from Bosch Klein TCST 2011
%   Electrode area chosen to correspond to 2.3 Ah cell
%   NOTE: Parameters represent an illustrative cell

function p = params_LCO()
    %% Geometric Params
    % Electrode current collector area [m^2]
    p.electrode_area = 1; 
    
    % Thickness of each layer
    p.anode.electrode_thickness   = 100e-6; % [m]
    p.separator.thickness         = 25e-6;  % [m]
    p.cathode.electrode_thickness = 100e-6; % [m]

    % Particle Radii
    p.anode.particle_radius   = 10e-6; % Radius of solid particles in negative electrode [m]
    p.cathode.particle_radius = 10e-6; % Radius of solid particles in positive electrode [m]
    
    % Volume fractions
    p.anode.volume_fraction_solid   = 0.6; % Volume fraction in solid for neg. electrode
    p.cathode.volume_fraction_solid = 0.5; % Volume fraction in solid for pos. electrode
    
    p.anode.volume_fraction_electrolyte     = 0.3; % Volume fraction in electrolyte for neg. electrode
    p.separator.volume_fraction_electrolyte = 1.0; % Volume fraction in electrolyte for separator
    p.cathode.volume_fraction_electrolyte   = 0.3; % Volume fraction in electrolyte for pos. electrode
    
    % make element to calculate phi_{s} by Saehong Park 
    p.anode.volume_fraction_filler   = 1-p.anode.volume_fraction_solid  -p.anode.volume_fraction_electrolyte;   % Volume fraction of filler in neg. electrode
    p.cathode.volume_fraction_filler = 1-p.cathode.volume_fraction_solid-p.cathode.volume_fraction_electrolyte; % Volume fraction of filler in pos. electrode
        
    % Specific interfacial surface area
    p.anode.specific_interfacial_area   = 3*p.anode.volume_fraction_solid   / p.anode.particle_radius;   % Negative electrode [m^2/m^3]
    p.cathode.specific_interfacial_area = 3*p.cathode.volume_fraction_solid / p.cathode.particle_radius; % Positive electrode [m^2/m^3]
    
    %% Transport Params
    % Diffusion coefficient in solid
    p.anode.diffusion_coefficient   = 3.9e-14; % Diffusion coeff for solid in neg. electrode, [m^2/s]
    p.cathode.diffusion_coefficient = 1e-13;   % Diffusion coeff for solid in pos. electrode, [m^2/s]

    % Bruggeman porosity
    p.bruggemann_porosity = 1.5;      
    
    % Miscellaneous
    p.transference_number = 0.4;

    %% Kinetic Params
    p.charge_transfer_coefficient = 0.5; % Charge transfer coefficients
    
    p.anode.sei_resistivity   = 1e-3; % Resistivity of SEI layer, [Ohms*m^2]
    p.cathode.sei_resistivity = 0;    % Resistivity of SEI layer, [Ohms*m^2]
    % Current Collector Resistance, [Ohms] -- NOT IMPLEMENTED
    % p.non_jellyroll_resistance = 2e-3;         
    
    % Nominal Reaction rates
    p.anode.reaction_rate   = 1e-5; % [(A/m^2)*(mol^3/mol)^(1+alpha)]
    p.cathode.reaction_rate = 3e-7; % [(A/m^2)*(mol^3/mol)^(1+alpha)]
    
    % Activation Energies
    % Taken from Zhang et al (2014) [Harbin]
    % http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
    p.activation_energy.anode_intercalation      = 37.48e3; % [J/mol] (kn)
    p.activation_energy.cathode_intercalation    = 39.57e3; % [J/mol] (kp)
    p.activation_energy.anode_diffusion          = 42.77e3; % [J/mol] (Dsn)
    p.activation_energy.cathode_diffusion        = 18.55e3; % [J/mol] (Dsp)
    p.activation_energy.electrolyte_diffusion    = 37.04e3; % [J/mol] (De)
    p.activation_energy.electrolyte_conductivity = 34.70e3; % [J/mol] (kappa_e)
 
    %% Aging submodel parameters
    % SEI Layer Growth model
    % Adopted from Ramadass et al (2004) [Univ of South Carolina]
    % "Development of First Principles Capacity Fade Model for Li-Ion Cells"
    % DOI: 10.1149/1.1634273
    % NOTE: These parameters have NOT been experimentally validated by eCAL
    p.side_reaction_product.conductivity             = 1;      % [S/m] conductivity of side rxn product
    p.side_reaction_product.molecular_weight         = 7.3e1;  % [kg/mol] molecular weight of side rxn product
    p.side_reaction_product.mass_density             = 2.1e3;  % [kg/m^3] mass density of side rxn product
    p.side_reaction_product.exchange_current_density = 0;      % [A/m^2] exchange current density of side rxn
    p.side_reaction_product.reference_potential      = 0.4;    % [V] reference potential of side rxn
    
    %% Concentrations
    % Maxima based on DUALFOIL 
    % line 588 in DUALFOIL Fortran code
    p.total_moles_lithium = 1.8; % [mol]
    p.anode.maximum_concentration   = 5.4e3 * 372 * 2800 / SPMe().F; % [mol/m^3]
    p.cathode.maximum_concentration = 5.4e3 * 274 * 3010 / SPMe().F; % [mol/m^3]
    p.electrolyte_concentration     = 1e3;                           % [mol/m^3] (Fixed)
    
    %% Cutoff voltages
    p.maximum_voltage = 4.7;        % [V]
    p.minimum_voltage = 3.105;      % [V]

    %% Thermal
    p.nominal_temperature = 298.15; % [K]

end


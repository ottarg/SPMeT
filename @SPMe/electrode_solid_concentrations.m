function  [anode_concentration,cathode_concentration] = electrode_solid_concentrations(obj,voltage)
%% Use Bisection Algorithm
anode = obj.cell_properties.anode;
cathode = obj.cell_properties.cathode;
total_moles_lithium = obj.cell_properties.total_moles_lithium;
electrode_area = obj.cell_properties.electrode_area;
% Algorithm params
maxiters = 50;
x = zeros(maxiters,1);
f = nan*ones(maxiters,1);
tol = 1e-5;

% Initial Guesses
x_low = 0.05 * cathode.maximum_concentration;
x_high = 1.0 * cathode.maximum_concentration;
x(1) = 0.6 * cathode.maximum_concentration;

% Iterate Bisection Algorithm
for idx = 1:maxiters

    theta_p = x(idx)/cathode.maximum_concentration;
    theta_n = (total_moles_lithium...
        - cathode.volume_fraction_solid*cathode.electrode_thickness*electrode_area*x(idx))...
          /(anode.volume_fraction_solid  *anode.electrode_thickness*anode.maximum_concentration*electrode_area);

    OCPn = SPMe().refPotentialAnode(theta_n);
    OCPp = SPMe().refPotentialCathode(theta_p);

    f(idx) = OCPp - OCPn - voltage;

    if(abs(f(idx)) <= tol)
        break;
    elseif(f(idx) <= 0)
        x_high = x(idx);
    else
        x_low = x(idx);
    end
    % Bisection
    x(idx+1) = (x_high + x_low)/2;
end

% Output converged concentrations
cathode_concentration = x(idx);
anode_concentration = (obj.cell_properties.total_moles_lithium...
- cathode.volume_fraction_solid * cathode.electrode_thickness * cathode_concentration * electrode_area)...
  /(anode.volume_fraction_solid   * anode.electrode_thickness * electrode_area);
end
function  [anode_concentration,cathode_concentration] = electrode_solid_concentrations(obj,voltage)
%% Use Bisection Algorithm

% Algorithm params
maxiters = 50;
x = zeros(maxiters,1);
f = nan*ones(maxiters,1);
tol = 1e-5;

% Initial Guesses
x_low = 0.05 * obj.cathode.maximum_concentration;
x_high = 1.0 * obj.cathode.maximum_concentration;
x(1) = 0.6 * obj.cathode.maximum_concentration;

% Iterate Bisection Algorithm
for idx = 1:maxiters

    theta_p = x(idx)/obj.cathode.maximum_concentration;
    theta_n = (obj.total_moles_lithium...
        - obj.cathode.volume_fraction_solid*obj.cathode.electrode_thickness*obj.electrode_area*x(idx))...
          /(obj.anode.volume_fraction_solid  *obj.anode.electrode_thickness*obj.anode.maximum_concentration*obj.electrode_area);

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
anode_concentration = (obj.total_moles_lithium...
- obj.cathode.volume_fraction_solid * obj.cathode.electrode_thickness * cathode_concentration * obj.electrode_area)...
  /(obj.anode.volume_fraction_solid   * obj.anode.electrode_thickness * obj.electrode_area);
end
function  [csn0,csp0] = initial_solid_concentrations(obj,initial_voltage)
%% Use Bisection Algorithm

% Algorithm params
maxiters = 50;
x = zeros(maxiters,1);
f = nan*ones(maxiters,1);
tol = 1e-5;

% Initial Guesses
x_low = 0.05 * obj.cell_properties.cathode.maximum_concentration;
x_high = 1.0 * obj.cell_properties.cathode.maximum_concentration;
x(1) = 0.6 * obj.cell_properties.cathode.maximum_concentration;

% Iterate Bisection Algorithm
for idx = 1:maxiters

    theta_p = x(idx)/obj.cell_properties.cathode.maximum_concentration;
    theta_n = (obj.cell_properties.total_moles_lithium-obj.cell_properties.cathode.volume_fraction_solid*obj.cell_properties.cathode.electrode_thickness*obj.cell_properties.electrode_area*x(idx))/(obj.cell_properties.anode.maximum_concentration*obj.cell_properties.anode.volume_fraction_solid*obj.cell_properties.anode.electrode_thickness*obj.cell_properties.electrode_area);

    OCPn = SPMe().refPotentialAnode(theta_n);
    OCPp = SPMe().refPotentialCathode(theta_p);

    f(idx) = OCPp - OCPn - initial_voltage;

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

% Output conveged csp0
csp0 = x(idx);
csn0 = (obj.cell_properties.total_moles_lithium - obj.cell_properties.cathode.volume_fraction_solid * obj.cell_properties.cathode.electrode_thickness * obj.cell_properties.electrode_area * csp0) / (obj.cell_properties.anode.volume_fraction_solid * obj.cell_properties.anode.electrode_thickness * obj.cell_properties.electrode_area);
end
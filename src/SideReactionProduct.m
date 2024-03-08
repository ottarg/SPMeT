classdef SideReactionProduct < handle

    properties
        conductivity
        molecular_weight = 1;
        mass_density
        exchange_current_density
        reference_potential
    end
    properties (Dependent)
    end

    methods
        function obj = SideReactionProduct()

        end
    end
end
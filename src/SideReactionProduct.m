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
        function s = struct(obj)
            s = struct('conductivity',obj.conductivity,...
                'molecular_weight',obj.molecular_weight,...
                'mass_density',obj.mass_density,...
                'exchange_current_density',obj.exchange_current_density,...
                'reference_potential',obj.reference_potential);
        end
    end
end
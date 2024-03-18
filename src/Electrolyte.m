classdef Electrolyte < handle

    properties
        volume_fraction_electrolyte = 1;
        concentration
        diffusion_activation_energy
        conductivity_activation_energy
    end
    properties (Dependent)
    end

    methods
        function obj = Electrolyte()

        end
        function s = struct(obj)
            s = struct('volume_fraction_electrolyte',...
                obj.volume_fraction_electrolyte,...
                'concentration',obj.concentration,...
                'diffusion_activation_energy',obj.diffusion_activation_energy,...
                'conductivity_activation_energy',obj.conductivity_activation_energy);
        end
    end
end

classdef Electrolyte < handle

    properties
        separator_thickness
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
    end
end
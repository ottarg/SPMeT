classdef Electrode < handle

    properties
        thickness
        area
        particle_radius
        volume_fraction_solid
        volume_fraction_electrolyte
        diffusion_coefficient
        reaction_rate
        maximum_concentration
    end
    properties (Dependent)
        volume_fraction_filler
        specific_interfacial_area
    end

    methods
        function obj = Electrode()

        end
        function val = get.volume_fraction_filler(obj)
            val = 1 - obj.volume_fraction_solid- obj.volume_fraction_electrolyte;
        end
        function val = get.specific_interfacial_area(obj)
            val = 3 * obj.volume_fraction_solid / obj.particle_radius;
        end
    end
end
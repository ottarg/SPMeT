classdef Electrode < handle

    properties
        area
        electrode_thickness
        particle_radius
        volume_fraction_solid
        volume_fraction_electrolyte
        diffusion_coefficient
        reaction_rate
        intercalation_activation_energy
        diffusion_activation_energy
        maximum_concentration
        solid_concentration_state_inds
    end
    properties (Dependent)
        volume_fraction_filler
        specific_interfacial_area
        volume
        molar_capacity
        charge_capacity
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
        function val = get.volume(obj)
            val = obj.area * obj.electrode_thickness;
        end
        function val = get.molar_capacity(obj)
            val = obj.volume * obj.maximum_concentration * obj.volume_fraction_solid;
        end
        function val = get.charge_capacity(obj)
            val = obj.molar_capacity * F / 3600;
        end
    end
end
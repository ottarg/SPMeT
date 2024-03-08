classdef Domain < matlab.mixin.Heterogeneous

    properties
        thickness
        division_count
        indices
        initial_concentration
    end
    properties (Dependent)
        division_size
        initial_concentration_array
    end

    methods
        function obj = Domain()

        end
        function val = get.division_size(obj)
            val = 1/obj.division_count;
        end
        function val = get.initial_concentration_array(obj)
            val =obj.initial_concentration*ones(obj.division_count-1,1);
        end
        
    end
end
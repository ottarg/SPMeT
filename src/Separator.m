classdef Separator < handle


    properties
        thickness
        volume_fraction_electrolyte
    end
    properties (Dependent)
    end

    methods
        function obj = Separator()

        end
        function s = struct(obj)
            s = struct('thickness',obj.thickness, 'volume_fraction_electrolyte',obj.volume_fraction_electrolyte);
        end
    end
end
classdef ElectrolyteDomain < Domain

    properties
        anode_side   = LinearDomain();
        cathode_side = LinearDomain();
        separator    = LinearDomain();
        transference_number
    end
    properties (Dependent)
        M1
        M2
        M3
        M4
        N1
        N2
        C
    end

    methods
        function obj = ElectrolyteDomain()

        end
    end
end
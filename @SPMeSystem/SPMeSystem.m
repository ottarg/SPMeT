classdef SPMeSystem < matlab.System
    % Properties and methods for the CellModel

    properties
        cell_properties
        discretization
        initial_voltage
        time_step = 0.5;
    end
    properties(SetAccess = protected)

        solid_phase_matrices
        electrolyte_matrices
    end
    properties (DiscreteState)
        x;
    end
    methods(Access = protected)
        function setupImpl(obj)
            % Solid concentration
            [csn0,csp0] = obj.initial_solid_concentrations(obj.initial_voltage);
            c_n0 = csn0 * ones(obj.discretization.Nr-1,1);
            c_p0 = csp0 * ones(obj.discretization.Nr-1,1);
            % Electrolyte concentration
            ce0 = obj.cell_properties.electrolyte_concentration*ones(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp - 3,1);
            % Temperature
            T10 = obj.cell_properties.ambient_temperature;
            T20 = obj.cell_properties.ambient_temperature;
            % SEI layer
            delta_sei0 = 0;

            initialize_electrolyte_matrices(obj);
            obj.x = [c_n0; c_p0; ce0; T10; T20; delta_sei0];
        end

        function [V,V_spm,SOC_n,SOC_p,c_ss_n,c_ss_p,c_e,OCV,anode_potential,cathode_potential] = stepImpl(obj,current,temperature)
            [x_dot,V,V_spm,SOC_n,SOC_p,c_ss_n,c_ss_p,c_e,OCV,anode_potential,cathode_potential] = spme_ode(obj,obj.x,current,temperature);
            obj.x = obj.x+x_dot*obj.time_step;
        end

        function resetImpl(obj)
            % Reset the internal states to initial conditions
        end
        [value, isterminal, direction] = detectImagSolution(obj, t, x, data)
        [csn0,csp0] = initial_solid_concentrations(obj,V)
        initialize_electrolyte_matrices(obj)
        initialize_solid_phase_matrices(obj,anode_diffusion_coefficient,cathode_diffusion_coefficient)
        [dActivity,varargout] = electrolyteAct(obj,c_e,T)
        [i_0n,i_0p,varargout] = exch_cur_dens(obj,k_p,k_n,c_ss_n,c_ss_p,c_e)
    end

    methods (Static)
        F = F
        R = R
        [Uref] = refPotentialCathode(theta)
        [Uref] = refPotentialAnode(theta)
        [kappa,varargout] = electrolyteCond(c_e)
        [D_e,varargout] = electrolyteDe(c_e)
    end
end
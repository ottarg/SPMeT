classdef SPMe < handle

    properties
        cell_properties
        discretization
        initial_voltage
        x0
    end
    properties(SetAccess = protected)

        solid_phase_matrices
        electrolyte_matrices
    end
    properties (Dependent)
    end
    methods

        function obj = SPMe()

        end

        function initialize(obj)
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
            obj.x0 = [c_n0; c_p0; ce0; T10; T20; delta_sei0];
        end

        function [res,x] = simulate(obj,time,current,temperature)
            obj.initialize;
            res.time = time;
            res.current = -current/obj.cell_properties.electrode_area*10;
            res.temperature = temperature;
            Opt    = odeset('Events',@(t,x)detectImagSolution(obj,t,x,res));

            [res.timeODE,x] = ode23s(@(t,x) spme_ode(obj,t,x,res),[res.time(1),res.time(end)],obj.x0,Opt);
            for k = 1:length(res.timeODE)
                % Compute outputs
                [~,res.V(k),res.V_spm(k),res.SOC_n(k),res.SOC_p(k),res.c_ss_n(k),res.c_ss_p(k),res.c_e(:,k),res.OCV(:,k),res.anode_potential(:,k),res.cathode_potential(:,k)] = ...
                    spme_ode(obj,res.timeODE(k),x(k,:)',res);
            end
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
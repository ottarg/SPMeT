classdef SPMe < handle

    properties
        cell_properties
        discretization
        initial_voltage
    end
    properties(SetAccess = protected)
        x0
        solid_phase_matrices
        electrolyte_matrices
    end
    properties (Dependent)
    end
    methods
        function obj = SPMe()
        end
        [csn0,csp0] = initial_solid_concentrations(obj,V)
        initialize_electrolyte_matrices(obj)
        initialize_solid_phase_matrices(obj,anode_diffusion_coefficient,cathode_diffusion_coefficient)
        [dActivity,varargout] = electrolyteAct(obj,c_e,T)
        [i_0n,i_0p,varargout] = exch_cur_dens(obj,k_p,k_n,c_ss_n,c_ss_p,c_e)
        [Uref] = refPotentialCathode(obj,theta)
        [Uref] = refPotentialAnode(obj,theta)
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
        function [res,x] = simulate(obj,time,current)
            res.time = time;
            res.cur = -current/obj.cell_properties.electrode_area*10;
            obj.discretization.delta_t = res.time(2)-res.time(1);
            Opt    = odeset('Events',@(t,x)myEvent(t,x,obj,res));

            [res.time,x] = ode23s(@(t,x) spme_ode(obj,t,x,res),[res.time(1),res.time(end)],obj.x0,Opt);
            for k = 1:length(res.time)
                % Compute outputs
                [~,res.V(k),res.V_spm(k),res.SOC_n(k),res.SOC_p(k),res.c_ss_n(k),res.c_ss_p(k),res.c_e(:,k),res.OCV(:,k),res.anode_potential(:,k),res.cathode_potential(:,k)] = ...
                    spme_ode(obj,res.time(k),x(k,:)',res);
            end

        end
    end
    methods (Static)
    end

end
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
        anode_capacity
        cathode_capacity
        capacity
        anode_concentration_range
        cathode_concentration_range
    end
    methods

        function obj = SPMe()

        end

        function initialize(obj)
            obj.discretization.delta_r_n = 1/obj.discretization.Nr;
            obj.discretization.delta_r_p = 1/obj.discretization.Nr;
            obj.discretization.Nx = obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp;
            % Finite difference points along x-coordinate
            obj.discretization.delta_x_n = 1 / obj.discretization.Nxn;
            obj.discretization.delta_x_s = 1 / obj.discretization.Nxs;
            obj.discretization.delta_x_p = 1 / obj.discretization.Nxp;

            % Solid concentration
            [initial_anode_concentration,initial_cathode_concentration] = obj.electrode_solid_concentrations(obj.initial_voltage);
            c_n0 = initial_anode_concentration * ones(obj.discretization.Nr-1,1);
            c_p0 = initial_cathode_concentration * ones(obj.discretization.Nr-1,1);
            % Electrolyte concentration
            ce0 = obj.cell_properties.electrolyte_concentration*ones(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp - 3,1);
            % SEI layer
            delta_sei0 = 0;

            initialize_electrolyte_matrices(obj);
            obj.x0 = [c_n0; c_p0; ce0; delta_sei0];
        end

        function [res,x] = simulate(obj,time,current,temperature)
            obj.initialize;
            data.time = time;
            data.current = -current/obj.cell_properties.electrode_area*10;
            data.temperature = temperature;
            Opt    = odeset('Events',@(t,x)detectImagSolution(obj,t,x,data));
            [res.time,x] = ode23s(@(t,x) spme_ode(obj,t,x,data),[data.time(1),data.time(end)],obj.x0,Opt);
            res.time = res.time';
            for k = 1:length(res.time)
                % Compute outputs
                [~,res.V(:,k),res.V_spm(:,k),res.SOC_n(:,k),res.SOC_p(:,k),...
                    res.anode_solid_surface_concentration(:,k),res.cathode_solid_surface_concentration(:,k),...
                    res.c_e(:,k),res.OCV(:,k),res.anode_potential(:,k),res.cathode_potential(:,k)] = ...
                    spme_ode(obj,res.time(k),x(k,:)',data);
            end
        end
        
        [value, isterminal, direction] = detectImagSolution(obj, t, x, data)
        [csn0,csp0] = electrode_solid_concentrations(obj,V)
        initialize_electrolyte_matrices(obj)
        initialize_solid_phase_matrices(obj,anode_diffusion_coefficient,cathode_diffusion_coefficient)
        [dActivity,varargout] = electrolyteAct(obj,c_e,T)
        [i_0n,i_0p,varargout] = exch_cur_dens(obj,k_p,k_n,c_ss_n,c_ss_p,c_e)
        function val = get.anode_capacity(obj)
            val = obj.cell_properties.anode.volume_fraction_solid*...
                  obj.cell_properties.anode.electrode_thickness*...
                      obj.anode_concentration_range*SPMe().F/3600;
        end
        function val = get.cathode_capacity(obj)
            val = obj.cell_properties.cathode.volume_fraction_solid*...
                  obj.cell_properties.cathode.electrode_thickness*...
                      obj.cathode_concentration_range*SPMe().F/3600;
        end
        function val = get.capacity(obj)
            val = min(obj.anode_capacity, obj.cathode_capacity);
        end
        function val = get.anode_concentration_range(obj)
             [low_v_anode_concentration,~] = obj.electrode_solid_concentrations(obj.cell_properties.minimum_voltage);
             [high_v_anode_concentration,~] = obj.electrode_solid_concentrations(obj.cell_properties.maximum_voltage);
             val = high_v_anode_concentration-low_v_anode_concentration;
        end
        function val = get.cathode_concentration_range(obj)
             [~,low_v_cathode_concentration] = obj.electrode_solid_concentrations(obj.cell_properties.minimum_voltage);
             [~,high_v_cathode_concentration] = obj.electrode_solid_concentrations(obj.cell_properties.maximum_voltage);
             val = low_v_cathode_concentration-high_v_cathode_concentration;
        end
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
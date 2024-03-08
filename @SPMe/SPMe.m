classdef SPMe < handle

    properties
        anode Anode = Anode()
        cathode Cathode = Cathode()
        electrolyte Electrolyte = Electrolyte()
        separator Separator = Separator()
        side_reaction_product SideReactionProduct = SideReactionProduct()
        electrode_area
        total_moles_lithium
        nominal_temperature
        maximum_voltage
        minimum_voltage
        transference_number
        charge_transfer_coefficient
        bruggemann_porosity
        discretization
        initial_voltage
    end
    properties(SetAccess = protected)

    end
    properties (Dependent)
        anode_capacity
        cathode_capacity
        capacity
    end
    properties (Hidden)
        solid_phase_matrices
        electrolyte_matrices
        x0
    end
    properties (Dependent, Hidden)
        anode_concentration_range
        cathode_concentration_range
    end
    methods

        function obj = SPMe()

        end

        function initialize(obj)
            obj.discretization.delta_r_n = 1 / obj.discretization.radial_divisions;
            obj.discretization.delta_r_p = 1 / obj.discretization.radial_divisions;
            % Finite difference points along x-coordinate
            obj.discretization.Nx = obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp;
            obj.discretization.delta_x_n = 1 / obj.discretization.Nxn;
            obj.discretization.delta_x_s = 1 / obj.discretization.Nxs;
            obj.discretization.delta_x_p = 1 / obj.discretization.Nxp;

            % Solid concentration
            [initial_anode_concentration,initial_cathode_concentration] = obj.electrode_solid_concentrations(obj.initial_voltage);
            initial_anode_concentration = initial_anode_concentration * ones(obj.discretization.radial_divisions-1,1);
            initial_cathode_concentration = initial_cathode_concentration * ones(obj.discretization.radial_divisions-1,1);
            % Electrolyte concentration
            electrolyte_concentration = obj.electrolyte.concentration*ones(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp - 3,1);
            % SEI layer
            initial_sei_growth = 0;
            initialize_electrolyte_matrices(obj);
            obj.x0 = [initial_anode_concentration; initial_cathode_concentration; electrolyte_concentration; initial_sei_growth];
        end

        function [res,x] = simulate(obj,time,current,temperature)
            obj.initialize;
            data.time = time;
            data.current = -current/obj.electrode_area*10;
            data.temperature = temperature;
            Opt    = odeset('Events',@(t,x)detectImagSolution(obj,t,x,data));
            model = struct(obj);
            [res.time,x] = ode23s(@(t,x) spme_ode(model,t,x,data),[data.time(1),data.time(end)],obj.x0,Opt);
            res.time = res.time';
            for k = 1:length(res.time)
                % Compute outputs
                [~,res.V(:,k),res.V_spm(:,k),res.SOC_n(:,k),res.SOC_p(:,k),...
                    res.anode_solid_surface_concentration(:,k),res.cathode_solid_surface_concentration(:,k),...
                    res.c_e(:,k),res.OCV(:,k),res.anode_potential(:,k),res.cathode_potential(:,k)] = ...
                    spme_ode(model,res.time(k),x(k,:)',data);
            end
        end
        
        [value, isterminal, direction] = detectImagSolution(obj, t, x, data)
        [anode_concentration,cathode_concentration] = electrode_solid_concentrations(obj,V)
        initialize_electrolyte_matrices(obj)  
        
        function val = get.anode_capacity(obj)
            val = obj.anode.volume_fraction_solid*...
                  obj.anode.electrode_thickness*...
                      obj.anode_concentration_range*SPMe().F/3600;
        end
        function val = get.cathode_capacity(obj)
            val = obj.cathode.volume_fraction_solid*...
                  obj.cathode.electrode_thickness*...
                      obj.cathode_concentration_range*SPMe().F/3600;
        end
        function val = get.capacity(obj)
            val = min(obj.anode_capacity, obj.cathode_capacity);
        end
        function val = get.anode_concentration_range(obj)
             [low_v_anode_concentration,~] = obj.electrode_solid_concentrations(obj.minimum_voltage);
             [high_v_anode_concentration,~] = obj.electrode_solid_concentrations(obj.maximum_voltage);
             val = high_v_anode_concentration-low_v_anode_concentration;
        end
        function val = get.cathode_concentration_range(obj)
             [~,low_v_cathode_concentration] = obj.electrode_solid_concentrations(obj.minimum_voltage);
             [~,high_v_cathode_concentration] = obj.electrode_solid_concentrations(obj.maximum_voltage);
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
        [dActivity,varargout] = electrolyteAct(model,c_e,T)
    end

end
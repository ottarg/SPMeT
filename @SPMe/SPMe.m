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
        solid_phase_matrices
        electrolyte_matrices
        x0
    end
    properties (Dependent)
        anode_concentration_range
        cathode_concentration_range
        anode_capacity
        cathode_capacity
        capacity
        NP_ratio
    end
    properties (Hidden)

    end
    properties (Dependent, Hidden)

    end
    methods

        function obj = SPMe()

        end
        function postProcess()
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
            obj = initialize_solid_phase_matrices(obj,0,0);
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
            data.current = current;
            data.temperature = temperature;
            Opt    = odeset('Events',@(t,x)detectImagSolution(obj,t,x,data));
            model = obj.getStruct;
            [res.time,x] = ode15s(@(t,x) spme_ode(model,t,x,data),[data.time(1),data.time(end)],obj.x0,Opt);
            res.time = res.time';
            for k = 1:length(res.time)
                % Compute outputs
                [~,res.V(:,k),res.V_spm(:,k),res.SOC_n(:,k),res.SOC_p(:,k),...
                    res.anode_solid_surface_concentration(:,k),res.cathode_solid_surface_concentration(:,k),...
                    res.electrolyte_concentrations(:,k),res.OCV(:,k),res.anode_potential(:,k),res.cathode_potential(:,k)] = ...
                    spme_ode(model,res.time(k),x(k,:)',data);
            end
        end
        
        [value, isterminal, direction] = detectImagSolution(obj, t, x, data)
        [anode_concentration,cathode_concentration] = electrode_solid_concentrations(obj,V)
        initialize_electrolyte_matrices(obj)  
        
        function val = get.anode_capacity(obj)
            val = obj.anode.volume_fraction_solid*...
                  obj.anode.electrode_thickness*...
                      obj.anode_concentration_range*F/3600;
        end
        function val = get.cathode_capacity(obj)
            val = obj.cathode.volume_fraction_solid*...
                  obj.cathode.electrode_thickness*...
                      obj.cathode_concentration_range*F/3600;
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
        function val = get.NP_ratio(obj)
             val = obj.anode.molar_capacity/obj.cathode.molar_capacity;
        end
        function s = getStruct(obj)
            for i = 1:length(obj)
                props = properties(obj(i));
                %convert props to struct
                for j = 1:length(props)
                    if (isa(obj(i).(props{j}),'Component'))
                        propName = (obj(i).(props{j}).name); % If subcomponent is an array of components, this will catch the first ones name. Different names not supported.
                        s(i).(propName) = obj(i).(props{j}).getStruct(); %get structures from the lower levels
                    elseif (isobject(obj(i).(props{j})))
                        s(i).(props{j}) = struct(obj(i).(props{j})); %get the property value
                    else
                        s(i).(props{j}) = obj(i).(props{j}); %get the property value
                    end
                end
            end
        end
    end
        
    methods (Static)
    end

end
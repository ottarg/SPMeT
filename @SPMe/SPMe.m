classdef SPMe < handle

    properties
        cell_properties
        discretization
        electrolyte_matrices
        solid_phase_matrices
        initial_voltage
        x0
    end
    properties(SetAccess = protected)
    end
    properties (Dependent)
    end
    methods
        function obj = SPMe()
        end
        [csn0,csp0] = initial_solid_concentrations(obj,V)
        initialize_electrolyte_matrices(obj)
        initialize_solid_phase_matrices(obj,anode_diffusion_coefficient,cathode_diffusion_coefficient)
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

            [res.time,x] = ode23s(@(t,x) spme_ode(obj,t,x,res),res.time,obj.x0,Opt);
            for k = 1:length(res.time/5)
                % Compute outputs
                [~,res.V(k),res.V_spm(k),res.SOC_n(k),res.SOC_p(k),res.c_ss_n(k),res.c_ss_p(k),res.c_e(:,k),res.OCV(:,k),res.anode_potential(:,k),res.cathode_potential(:,k)] = ...
                    spme_ode(obj,res.time(k),x(k,:)',res);
            end

        end

        

       
        %% Electrolyte Activity Coefficient Function: dlnf/dln(c_e)
        %   Created Oct 25, 2016 by Saehong Park

        function [dActivity,varargout] = electrolyteAct(obj,c_e,T)

            % From LiPF6, Valoen et al. 2005
            % Fig.6 in the paper

            % DataFitting Coefficients

            v00 = 0.601;
            v01 = 0;
            v10 = -0.24;
            v11 = 0;
            v20 = 0;
            v21 = 0;
            v30 = 0.982;
            v31 = -0.0052;

            c_e = c_e/1000; % UnitConversion: 1 mol/L -> 1000 mol/m^3

            dActivity = ((v00 + v10.*(c_e).^(0.5) + v30*(1+v31*(T - obj.cell_properties.nominal_temperature)) .* (c_e).^(1.5))./(1-obj.cell_properties.t_plus))-1;

            if(nargout == 2)
                d_dactivity = (0.5 * v10 * (c_e).^(-0.5) + (1.5)*v30*(1+v31*(T - obj.cell_properties.nominal_temperature))*(c_e).^(0.5))/(1-obj.cell_properties.t_plus);
                varargout{1} = d_dactivity;
            end
        end
        %% Reference Potential for Pos. Electrode: Unref(theta_n)
        %   Created July 12, 2011 by Scott Moura

        function [Uref] = refPotentialCathode(obj,theta)

            % DUALFOIL: CoO2 (Cobalt Dioxide) 0.5 < y < 0.99
            Uref = 2.16216+0.07645*tanh(30.834-54.4806*theta) ...
                + 2.1581*tanh(52.294-50.294*theta) ...
                - 0.14169*tanh(11.0923-19.8543*theta) ...
                + 0.2051*tanh(1.4684-5.4888*theta) ...
                + 0.2531*tanh((-theta+0.56478)/0.1316) ...
                - 0.02167*tanh((theta-0.525)/0.006);
        end

        function [Uref] = refPotentialAnode(obj,theta)

%             if(~isreal(theta))
%                 beep;
%                 error('dfn:err','Complex theta_n');
                %     pause;
%             end

            % DUALFOIL: MCMB 2528 graphite (Bellcore) 0.01 < x < 0.9
            Uref = 0.194+1.5*exp(-120.0*theta) ...
                +0.0351*tanh((theta-0.286)/0.083) ...
                - 0.0045*tanh((theta-0.849)/0.119) ...
                - 0.035*tanh((theta-0.9233)/0.05) ...
                - 0.0147*tanh((theta-0.5)/0.034) ...
                - 0.102*tanh((theta-0.194)/0.142) ...
                - 0.022*tanh((theta-0.9)/0.0164) ...
                - 0.011*tanh((theta-0.124)/0.0226) ...
                + 0.0155*tanh((theta-0.105)/0.029);
        end
        function [i_0n,i_0p,varargout] = exch_cur_dens(obj,k_p,k_n,c_ss_n,c_ss_p,c_e)

            % Parse out concentrations in anode and cathode
            ce_interp=interp1(linspace(0,1,length(c_e)),c_e,linspace(0,1,obj.discretization.Nx),'linear');
            c_e_n = ce_interp(1:obj.discretization.Nxn-1);
            c_e_p = ce_interp(obj.discretization.Nxn+obj.discretization.Nxs-1:end);


            % Compute exchange current density
            i_0n = k_n * ((obj.cell_properties.anode.maximum_concentration - c_ss_n) .* c_ss_n .* c_e_n).^obj.cell_properties.charge_transfer_coefficient;
            di_0n = k_n * (obj.cell_properties.charge_transfer_coefficient*((obj.cell_properties.anode.maximum_concentration - c_ss_n) .* c_ss_n .* c_e_n).^(obj.cell_properties.charge_transfer_coefficient-1).*( c_e_n.*(obj.cell_properties.anode.maximum_concentration - c_ss_n) -  c_ss_n .* c_e_n) );
            i_0p = k_p * ((obj.cell_properties.cathode.maximum_concentration - c_ss_p) .* c_ss_p .* c_e_p).^obj.cell_properties.charge_transfer_coefficient;
            di_0p = k_p * (obj.cell_properties.charge_transfer_coefficient*(max((obj.cell_properties.cathode.maximum_concentration - c_ss_p),0) .* c_ss_p .* c_e_p).^(obj.cell_properties.charge_transfer_coefficient-1).*( c_e_p.*(obj.cell_properties.cathode.maximum_concentration - c_ss_p) -  c_ss_p .* c_e_p) );

            if(nargout > 2)
                varargout{1}=di_0n;
                varargout{2}=di_0p;
            end
        end
    end
    methods (Static)
    end

end
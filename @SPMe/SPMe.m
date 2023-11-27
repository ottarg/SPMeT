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
        function  [csn0,csp0] = initial_solid_concentrations(obj,V)
            %% Use Bisection Algorithm

            % Algorithm params
            maxiters = 50;
            x = zeros(maxiters,1);
            f = nan*ones(maxiters,1);
            tol = 1e-5;

            % Initial Guesses
            x_low = 0.2 * obj.cell_properties.c_s_p_max;
            x_high = 1.0 * obj.cell_properties.c_s_p_max;
            x(1) = 0.6 * obj.cell_properties.c_s_p_max;

            % Iterate Bisection Algorithm
            for idx = 1:maxiters

                theta_p = x(idx)/obj.cell_properties.c_s_p_max;
                theta_n = (obj.cell_properties.n_Li_s-obj.cell_properties.cathode.volume_fraction_solid*obj.cell_properties.cathode.electrode_thickness*obj.cell_properties.electrode_area*x(idx))/(obj.cell_properties.c_s_n_max*obj.cell_properties.anode.volume_fraction_solid*obj.cell_properties.anode.electrode_thickness*obj.cell_properties.electrode_area);

                OCPn = refPotentialAnode(obj.cell_properties,theta_n);
                OCPp = refPotentialCathode(obj.cell_properties,theta_p);

                f(idx) = OCPp - OCPn - V;

                if(abs(f(idx)) <= tol)
                    break;
                elseif(f(idx) <= 0)
                    x_high = x(idx);
                else
                    x_low = x(idx);
                end

                % Bisection
                x(idx+1) = (x_high + x_low)/2;
                x(idx+1)/obj.cell_properties.c_s_p_max;

            end

            % Output conveged csp0
            csp0 = x(idx);
            csn0 = (obj.cell_properties.n_Li_s - obj.cell_properties.cathode.volume_fraction_solid * obj.cell_properties.cathode.electrode_thickness * obj.cell_properties.electrode_area * csp0) / (obj.cell_properties.anode.volume_fraction_solid * obj.cell_properties.anode.electrode_thickness * obj.cell_properties.electrode_area);
        end
        function initialize_electrolyte_matrices(obj)

            %% Lumped Coefficients
            Del_xn = obj.cell_properties.anode.electrode_thickness * obj.discretization.delta_x_n;
            Del_xs = obj.cell_properties.separator.thickness * obj.discretization.delta_x_s;
            Del_xp = obj.cell_properties.cathode.electrode_thickness * obj.discretization.delta_x_p;

            %% Matrices in nonlinear dynamics
            obj.electrolyte_matrices.M1n = sparse((diag(ones(obj.discretization.Nxn-2,1),+1) - diag(ones(obj.discretization.Nxn-2,1),-1))/(2*Del_xn));
            obj.electrolyte_matrices.M1s = sparse((diag(ones(obj.discretization.Nxs-2,1),+1) - diag(ones(obj.discretization.Nxs-2,1),-1))/(2*Del_xs));
            obj.electrolyte_matrices.M1p = sparse((diag(ones(obj.discretization.Nxp-2,1),+1) - diag(ones(obj.discretization.Nxp-2,1),-1))/(2*Del_xp));


            M2n = zeros(obj.discretization.Nxn-1,2);
            M2n(1,1) = -1/(2*Del_xn);
            M2n(end,end) = 1/(2*Del_xn);
            obj.electrolyte_matrices.M2n = sparse(M2n);

            M2s = zeros(obj.discretization.Nxs-1,2);
            M2s(1,1) = -1/(2*Del_xs);
            M2s(end,end) = 1/(2*Del_xs);
            obj.electrolyte_matrices.M2s = sparse(M2s);

            M2p = zeros(obj.discretization.Nxp-1,2);
            M2p(1,1) = -1/(2*Del_xp);
            M2p(end,end) = 1/(2*Del_xp);
            obj.electrolyte_matrices.M2p = sparse(M2p);


            obj.electrolyte_matrices.M3n = sparse((-2*diag(ones(obj.discretization.Nxn-1,1),0) + diag(ones(obj.discretization.Nxn-2,1),+1) + diag(ones(obj.discretization.Nxn-2,1),-1))/(Del_xn^2));
            obj.electrolyte_matrices.M3s = sparse((-2*diag(ones(obj.discretization.Nxs-1,1),0) + diag(ones(obj.discretization.Nxs-2,1),+1) + diag(ones(obj.discretization.Nxs-2,1),-1))/(Del_xs^2));
            obj.electrolyte_matrices.M3p = sparse((-2*diag(ones(obj.discretization.Nxp-1,1),0) + diag(ones(obj.discretization.Nxp-2,1),+1) + diag(ones(obj.discretization.Nxp-2,1),-1))/(Del_xp^2));


            M4n = zeros(obj.discretization.Nxn-1,2);
            M4n(1,1) = 1/(Del_xn^2);
            M4n(end,end) = 1/(Del_xn^2);
            obj.electrolyte_matrices.M4n = sparse(M4n);

            M4s = zeros(obj.discretization.Nxs-1,2);
            M4s(1,1) = 1/(Del_xs^2);
            M4s(end,end) = 1/(Del_xs^2);
            obj.electrolyte_matrices.M4s = sparse(M4s);

            M4p = zeros(obj.discretization.Nxp-1,2);
            M4p(1,1) = 1/(Del_xp^2);
            M4p(end,end) = 1/(Del_xp^2);
            obj.electrolyte_matrices.M4p = sparse(M4p);

            obj.electrolyte_matrices.M5n = (1-obj.cell_properties.t_plus)*obj.cell_properties.anode.specific_interfacial_area/obj.cell_properties.anode.volume_fraction_electrolyte * speye(obj.discretization.Nxn-1);
            obj.electrolyte_matrices.M5p = (1-obj.cell_properties.t_plus)*obj.cell_properties.cathode.specific_interfacial_area/obj.cell_properties.cathode.volume_fraction_electrolyte * speye(obj.discretization.Nxp-1);

            %% Boundary Conditions
            N1 = zeros(4,obj.discretization.Nx-3);
            N2 = zeros(4);

            % BC1
            N1(1,1) = +4;
            N1(1,2) = -1;
            N2(1,1) = -3;

            % BC2
            N1(2,obj.discretization.Nxn-2) = (obj.cell_properties.anode.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xn);
            N1(2,obj.discretization.Nxn-1) = (-4*obj.cell_properties.anode.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xn);
            N2(2,2) = (3*obj.cell_properties.anode.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xn) + (3*obj.cell_properties.separator.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xs);
            N1(2,obj.discretization.Nxn) = (-4*obj.cell_properties.separator.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xs);
            N1(2,obj.discretization.Nxn+1) = (obj.cell_properties.separator.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xs);

            % BC3
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-3) = (obj.cell_properties.separator.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xs);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-2) = (-4*obj.cell_properties.separator.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xs);
            N2(3,3) = (3*obj.cell_properties.separator.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xs) + (3*obj.cell_properties.cathode.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xp);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-1) = (-4*obj.cell_properties.cathode.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xp);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs) = (obj.cell_properties.cathode.volume_fraction_electrolyte^obj.cell_properties.bruggemann_porosity)/(2*Del_xp);


            % BC4
            N1(4,end-1) = +1;
            N1(4,end) = -4;
            N2(4,4) = +3;

            %%% SPARSE OUTPUT
            obj.electrolyte_matrices.C = sparse(-N2\N1);
        end
        function initialize(obj)
            % Solid concentration
            [csn0,csp0] = obj.initial_solid_concentrations(obj.initial_voltage);
            c_n0 = csn0 * ones(obj.discretization.Nr-1,1);
            c_p0 = csp0 * ones(obj.discretization.Nr-1,1);
            % Electrolyte concentration
            ce0 = obj.cell_properties.c_e*ones(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp - 3,1);
            % Temperature
            T10 = obj.cell_properties.ambient_temperature;
            T20 = obj.cell_properties.ambient_temperature;

            % SEI layer
            delta_sei0 = 0;

            initialize_electrolyte_matrices(obj);
            %             initialize_solid_phase_matrices(obj);
            obj.x0 = [c_n0; c_p0; ce0; T10; T20; delta_sei0];

        end
        function res = simulate(obj)
            load('input-data/UDDS');

            I = -current_exp'/obj.cell_properties.electrode_area*10;
            t = time_exp';
            obj.discretization.delta_t = t(2)-t(1);
            data.time = t;
            data.cur = I;
            [t,x] = ode23s(@(t,x) spme_ode(obj,t,x,data),t,obj.x0);
            for k = 1:length(t)
                % Compute outputs
                [~,res.V(k),res.V_spm(k),res.SOC_n(k),res.SOC_p(k),res.c_ss_n(k),res.c_ss_p(k),res.c_e(:,k)] = ...
                    spme_ode(obj,t(k),x(k,:)',data);
            end

        end

        

        function initialize_solid_phase_matrices(obj,D_s_n,D_s_p)

            % Electrochemical Model Parameters
            alpha_n = D_s_n / (obj.cell_properties.anode.particle_radius * obj.discretization.delta_r_n)^2;
            alpha_p = D_s_p / (obj.cell_properties.cathode.particle_radius * obj.discretization.delta_r_p)^2;

            % Block matrices
            M1_n = zeros(obj.discretization.Nr-1);
            M1_p = zeros(obj.discretization.Nr-1);

            for idx = 1:(obj.discretization.Nr-1)

                % Lower diagonal
                if(idx ~= 1)
                    M1_n(idx,idx-1) = (idx-1)/idx * alpha_n;
                    M1_p(idx,idx-1) = (idx-1)/idx * alpha_p;
                end

                % Main diagonal
                M1_n(idx,idx) = -2*alpha_n;
                M1_p(idx,idx) = -2*alpha_p;

                % Upper diagonal
                if(idx ~= (obj.discretization.Nr-1))
                    M1_n(idx,idx+1) = (idx+1)/idx * alpha_n;
                    M1_p(idx,idx+1) = (idx+1)/idx * alpha_p;
                end
            end

            M2_n = zeros(obj.discretization.Nr-1,2);
            M2_p = zeros(obj.discretization.Nr-1,2);

            M2_n(end,end) = obj.discretization.Nr/(obj.discretization.Nr-1) * alpha_n;
            M2_p(end,end) = obj.discretization.Nr/(obj.discretization.Nr-1) * alpha_p;

            N1 = zeros(2,obj.discretization.Nr-1);
            % % 1st Order BCs
            % N1(1,1) = 1;
            % N1(end,end) = -1;
            %
            % N2 = diag([-1,1]);
            %

            % 2nd Order BCs
            N1(1,1) = 4;
            N1(1,2) = -1;
            N1(2,end) = -4;
            N1(2,end-1) = 1;

            N2 = diag([-3,3]);

            N3_n = [0; -(2*obj.discretization.delta_r_n * obj.cell_properties.anode.particle_radius)/(D_s_n)];
            N3_p = [0; -(2*obj.discretization.delta_r_p * obj.cell_properties.cathode.particle_radius)/(D_s_p)];

            % A,B matrices for each electrode
            obj.solid_phase_matrices.A_n = M1_n - M2_n*(N2\N1);
            obj.solid_phase_matrices.A_p = M1_p - M2_p*(N2\N1);

            obj.solid_phase_matrices.B_n = M2_n*(N2\N3_n);
            obj.solid_phase_matrices.B_p = M2_p*(N2\N3_p);

            % C,D matrices for each electrode
            obj.solid_phase_matrices.C_n = -[0,1]*(N2\N1);
            obj.solid_phase_matrices.C_p = -[0,1]*(N2\N1);

            obj.solid_phase_matrices.D_n = [0,1]*(N2\N3_n);
            obj.solid_phase_matrices.D_p = [0,1]*(N2\N3_p);

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

            if(~isreal(theta))
                beep;
                error('dfn:err','Complex theta_n');
                %     pause;
            end

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
            i_0n = k_n * ((obj.cell_properties.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^obj.cell_properties.charge_transfer_coefficient;
            di_0n = k_n * (obj.cell_properties.charge_transfer_coefficient*((obj.cell_properties.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^(obj.cell_properties.charge_transfer_coefficient-1).*( c_e_n.*(obj.cell_properties.c_s_n_max - c_ss_n) -  c_ss_n .* c_e_n) );
            i_0p = k_p * ((obj.cell_properties.c_s_p_max - c_ss_p) .* c_ss_p .* c_e_p).^obj.cell_properties.charge_transfer_coefficient;
            di_0p = k_p * (obj.cell_properties.charge_transfer_coefficient*(max((obj.cell_properties.c_s_p_max - c_ss_p),0) .* c_ss_p .* c_e_p).^(obj.cell_properties.charge_transfer_coefficient-1).*( c_e_p.*(obj.cell_properties.c_s_p_max - c_ss_p) -  c_ss_p .* c_e_p) );

            if(nargout > 2)
                varargout{1}=di_0n;
                varargout{2}=di_0p;
            end
        end
    end
    methods (Static)
    end

end
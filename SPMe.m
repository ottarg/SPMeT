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
                theta_n = (obj.cell_properties.n_Li_s-obj.cell_properties.epsilon_s_p*obj.cell_properties.L_p*obj.cell_properties.Area*x(idx))/(obj.cell_properties.c_s_n_max*obj.cell_properties.epsilon_s_n*obj.cell_properties.L_n*obj.cell_properties.Area);

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
            csn0 = (obj.cell_properties.n_Li_s - obj.cell_properties.epsilon_s_p * obj.cell_properties.L_p * obj.cell_properties.Area * csp0) / (obj.cell_properties.epsilon_s_n * obj.cell_properties.L_n * obj.cell_properties.Area);
        end
        function initialize_electrolyte_matrices(obj)

            %% Lumped Coefficients
            Del_xn = obj.cell_properties.L_n * obj.discretization.delta_x_n;
            Del_xs = obj.cell_properties.L_s * obj.discretization.delta_x_s;
            Del_xp = obj.cell_properties.L_p * obj.discretization.delta_x_p;

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

            obj.electrolyte_matrices.M5n = (1-obj.cell_properties.t_plus)*obj.cell_properties.a_s_n/obj.cell_properties.epsilon_e_n * speye(obj.discretization.Nxn-1);
            obj.electrolyte_matrices.M5p = (1-obj.cell_properties.t_plus)*obj.cell_properties.a_s_p/obj.cell_properties.epsilon_e_p * speye(obj.discretization.Nxp-1);

            %% Boundary Conditions
            N1 = zeros(4,obj.discretization.Nx-3);
            N2 = zeros(4);

            % BC1
            N1(1,1) = +4;
            N1(1,2) = -1;
            N2(1,1) = -3;

            % BC2
            N1(2,obj.discretization.Nxn-2) = (obj.cell_properties.epsilon_e_n^obj.cell_properties.brug)/(2*Del_xn);
            N1(2,obj.discretization.Nxn-1) = (-4*obj.cell_properties.epsilon_e_n^obj.cell_properties.brug)/(2*Del_xn);
            N2(2,2) = (3*obj.cell_properties.epsilon_e_n^obj.cell_properties.brug)/(2*Del_xn) + (3*obj.cell_properties.epsilon_e_s^obj.cell_properties.brug)/(2*Del_xs);
            N1(2,obj.discretization.Nxn) = (-4*obj.cell_properties.epsilon_e_s^obj.cell_properties.brug)/(2*Del_xs);
            N1(2,obj.discretization.Nxn+1) = (obj.cell_properties.epsilon_e_s^obj.cell_properties.brug)/(2*Del_xs);

            % BC3
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-3) = (obj.cell_properties.epsilon_e_s^obj.cell_properties.brug)/(2*Del_xs);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-2) = (-4*obj.cell_properties.epsilon_e_s^obj.cell_properties.brug)/(2*Del_xs);
            N2(3,3) = (3*obj.cell_properties.epsilon_e_s^obj.cell_properties.brug)/(2*Del_xs) + (3*obj.cell_properties.epsilon_e_p^obj.cell_properties.brug)/(2*Del_xp);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-1) = (-4*obj.cell_properties.epsilon_e_p^obj.cell_properties.brug)/(2*Del_xp);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs) = (obj.cell_properties.epsilon_e_p^obj.cell_properties.brug)/(2*Del_xp);


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
            T10 = obj.cell_properties.T_amb;
            T20 = obj.cell_properties.T_amb;

            % SEI layer
            delta_sei0 = 0;

            initialize_electrolyte_matrices(obj);
            %             initialize_solid_phase_matrices(obj);
            obj.x0 = [c_n0; c_p0; ce0; T10; T20; delta_sei0];

        end
        function res = simulate(obj)
            load('input-data/UDDS');

            I = -current_exp'/obj.cell_properties.Area*10;
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

        function [x_dot,varargout] = spme_ode(obj,t,x,data)

            %% Parse Input Data

            % Parse and interpolate current
            cur = interp1(data.time,data.cur,t,[]);

            % Parse states
            c_s_n = x(1:(obj.discretization.Nr-1));
            c_s_p = x(obj.discretization.Nr : 2*(obj.discretization.Nr-1));
            c_e = x(2*obj.discretization.Nr-1 : 2*obj.discretization.Nr-1+obj.discretization.Nx-4);
            T1 = x(end-2);
            T2 = x(end-1);
            delta_sei = x(end);


            %% Pre-calculations with current states

            %%% MOLAR FLUXES
            % Compute total molar flux
            jn_tot = cur/(faraday*obj.cell_properties.a_s_n*obj.cell_properties.Area*obj.cell_properties.L_n);
            jp_tot = -cur/(faraday*obj.cell_properties.a_s_p*obj.cell_properties.Area*obj.cell_properties.L_p);


            %%% SOLID PHASE DYNAMICS
            % Solid phase diffusivity temperature dependence
            D_s_n = obj.cell_properties.D_s_n0 * exp(obj.cell_properties.E.Dsn/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
            D_s_p = obj.cell_properties.D_s_n0 * exp(obj.cell_properties.E.Dsp/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));

            % Construct (A,B) matrices for solid-phase Li diffusion
            initialize_solid_phase_matrices(obj,D_s_n,D_s_p)
            % Compute surface concentrations
            c_ss_n = obj.solid_phase_matrices.C_n*c_s_n + obj.solid_phase_matrices.D_n*jn_tot;
            c_ss_p = obj.solid_phase_matrices.C_p*c_s_p + obj.solid_phase_matrices.D_p*jp_tot;
            % Remark: I am cheating slightly here. jn_tot should be jn, but doing so
            % imposes an algebraic equation. This forms a DAE. I am going to
            % approximate jn by jn_tot, which should be ok, since jn and jn_tot have a
            % error on the order of 0.001%


            %%% ELECTROLYTE PHASE DYNAMICS
            % Compute electrolyte Boundary Conditions
            c_e_bcs = obj.electrolyte_matrices.C * c_e;

            ce0n = c_e_bcs(1);
            cens = c_e_bcs(2);
            cesp = c_e_bcs(3);
            ce0p = c_e_bcs(4);

            % Separate and aggregate electrolyte concentration
            c_en = c_e(1:(obj.discretization.Nxn-1));
            c_es = c_e((obj.discretization.Nxn-1)+1:(obj.discretization.Nxn-1)+(obj.discretization.Nxs-1));
            c_ep = c_e((obj.discretization.Nxn-1)+obj.discretization.Nxs : end);
            c_ex = [ce0n; c_en; cens; c_es; cesp; c_ep; ce0p];


            %% Voltage output

            % Average electrolyte concentrations
            cen_bar = mean(c_ex(1:obj.discretization.Nxn+1,:));
            ces_bar = mean(c_ex((obj.discretization.Nxn+1):(obj.discretization.Nxn+obj.discretization.Nxs+1),:));
            cep_bar = mean(c_ex((obj.discretization.Nxn+obj.discretization.Nxs+1):(obj.discretization.Nxn+obj.discretization.Nxs+obj.discretization.Nxp+1),:));

            % Overpotentials due to electrolyte subsystem
            kap_n_ref = electrolyteCond(cen_bar);
            kap_s_ref = electrolyteCond(ces_bar);
            kap_p_ref = electrolyteCond(cep_bar);

            % Adjustment for Arrhenius temperature dependence
            kap_n = kap_n_ref * exp(obj.cell_properties.E.kappa_e/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
            kap_s = kap_s_ref * exp(obj.cell_properties.E.kappa_e/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
            kap_p = kap_p_ref * exp(obj.cell_properties.E.kappa_e/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));

            % Bruggeman relationships
            kap_n_eff = kap_n * obj.cell_properties.epsilon_e_n.^(obj.cell_properties.brug);
            kap_s_eff = kap_s * obj.cell_properties.epsilon_e_s.^(obj.cell_properties.brug);
            kap_p_eff = kap_p * obj.cell_properties.epsilon_e_p.^(obj.cell_properties.brug);

            % Activity coefficient
            dfca_n = electrolyteAct(obj,cen_bar,T1);
            dfca_s = electrolyteAct(obj,ces_bar,T1);
            dfca_p = electrolyteAct(obj,cep_bar,T1);

            % Kinetic reaction rate, adjusted for Arrhenius temperature dependence
            k_n = obj.cell_properties.k_n0 * exp(obj.cell_properties.E.kn/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
            k_p = obj.cell_properties.k_p0 * exp(obj.cell_properties.E.kp/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));

            % Stochiometric Concentration Ratio
            theta_n = c_ss_n / obj.cell_properties.c_s_n_max;
            theta_p = c_ss_p / obj.cell_properties.c_s_p_max;

            % Equilibrium Potential
            Unref = refPotentialAnode(obj,theta_n);
            Upref = refPotentialCathode(obj,theta_p);

            % Exchange current density
            c_e_bar = [cen_bar; ces_bar; cep_bar];
            [i_0n,i_0p] = exch_cur_dens(obj,k_p,k_n,c_ss_n,c_ss_p,c_e_bar);
            % Overpotentials
            RTaF=(obj.cell_properties.R*T1)/(obj.cell_properties.alph*faraday);
            eta_n = RTaF * asinh(cur / (2*obj.cell_properties.a_s_n*obj.cell_properties.Area*obj.cell_properties.L_n*i_0n(1)));
            eta_p = RTaF * asinh(-cur / (2*obj.cell_properties.a_s_p*obj.cell_properties.Area*obj.cell_properties.L_p*i_0p(end)));

            % Total resistance (film + growing SEI layer)
            R_tot_n = obj.cell_properties.R_f_n + delta_sei/obj.cell_properties.kappa_P;
            R_tot_p = obj.cell_properties.R_f_p + 0;

            % SPM Voltage (i.e. w/o electrolyte concentration terms)
            V_noVCE = eta_p - eta_n + Upref - Unref ...
                - (R_tot_n/(obj.cell_properties.a_s_n*obj.cell_properties.L_n*obj.cell_properties.Area) + R_tot_p/(obj.cell_properties.a_s_p*obj.cell_properties.L_p*obj.cell_properties.Area))*cur;

            % Overpotential due to electrolyte conductivity
            V_electrolyteCond = (obj.cell_properties.L_n/(2*kap_n_eff) + 2*obj.cell_properties.L_s/(2*kap_s_eff) + obj.cell_properties.L_p/(2*kap_p_eff))*cur;

            % Overpotential due to electrolyte polarization
            V_electrolytePolar = (2*obj.cell_properties.R*T1)/(faraday) * (1-obj.cell_properties.t_plus)* ...
                ( (1+dfca_n) * (log(cens) - log(ce0n)) ...
                +(1+dfca_s) * (log(cesp) - log(cens)) ...
                +(1+dfca_p) * (log(ce0p) - log(cesp)));

            % Add 'em up!
            V = V_noVCE + V_electrolyteCond + V_electrolytePolar;


            %% Aging Dynamics

            %   SEI Layer Growth model
            %   Eqns Adopted from Ramadass et al (2004) [Univ of South Carolina]
            %   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
            %   DOI: 10.1149/1.1634273
            %   NOTE1: This model has NOT been validated experimentally by eCAL
            %   NOTE2: We assume this submodel only applies to anode

            % Difference btw solid and electrolyte overpotential [V]
            phi_se = eta_n + Unref + faraday*R_tot_n*jn_tot;

            % Side exn overpotential [V]
            eta_s = phi_se - obj.cell_properties.Us - faraday*R_tot_n * jn_tot;

            % Molar flux of side rxn [mol/s-m^2]
            j_s = -obj.cell_properties.i0s/faraday * exp((-obj.cell_properties.alph*faraday)/(obj.cell_properties.R*T1)*eta_s);

            % SEI layer growth model [m/s]
            delta_sei_dot = -obj.cell_properties.M_P/(obj.cell_properties.rho_P) * j_s;

            % Molar flux of intercalation
            jn = (abs(jn_tot) - abs(j_s)) * sign(jn_tot);
            jp = jp_tot;

            %% Solid Phase Dynamics

            % ODE for c_s
            c_s_n_dot = obj.solid_phase_matrices.A_n*c_s_n + obj.solid_phase_matrices.B_n*jn;
            c_s_p_dot = obj.solid_phase_matrices.A_p*c_s_p + obj.solid_phase_matrices.B_p*jp;


            %% Electrolyte Dynamics

            % Compute Electrolyte Diffusion Coefficient and Derivative
            [D_en0,dD_en0] = electrolyteDe(c_en);
            [D_es0,dD_es0] = electrolyteDe(c_es);
            [D_ep0,dD_ep0] = electrolyteDe(c_ep);

            % Adjustment for Arrhenius temperature dependence
            Arrh_De = exp(obj.cell_properties.E.De/obj.cell_properties.R*(1/obj.cell_properties.T_ref - 1/T1));
            D_en = D_en0 * Arrh_De;
            D_es = D_es0 * Arrh_De;
            D_ep = D_ep0 * Arrh_De;
            dD_en = dD_en0 * Arrh_De;
            dD_es = dD_es0 * Arrh_De;
            dD_ep = dD_ep0 * Arrh_De;

            % Apply BRUGGEMAN RELATION
            D_en_eff = D_en .* obj.cell_properties.epsilon_e_n.^(obj.cell_properties.brug-1);
            dD_en_eff = dD_en .* obj.cell_properties.epsilon_e_n.^(obj.cell_properties.brug-1);

            D_es_eff = D_es .* obj.cell_properties.epsilon_e_s.^(obj.cell_properties.brug-1);
            dD_es_eff = dD_es .* obj.cell_properties.epsilon_e_s.^(obj.cell_properties.brug-1);

            D_ep_eff = D_ep .* obj.cell_properties.epsilon_e_p.^(obj.cell_properties.brug-1);
            dD_ep_eff = dD_ep .* obj.cell_properties.epsilon_e_p.^(obj.cell_properties.brug-1);

            % Compute derivative
            c_en_dot = dD_en_eff.*(obj.electrolyte_matrices.M1n*c_en + obj.electrolyte_matrices.M2n*c_e_bcs(1:2)).^2 ...
                + D_en_eff.*(obj.electrolyte_matrices.M3n*c_en + obj.electrolyte_matrices.M4n*c_e_bcs(1:2)) + diag(obj.electrolyte_matrices.M5n)*jn;

            c_es_dot = dD_es_eff.*(obj.electrolyte_matrices.M1s*c_es + obj.electrolyte_matrices.M2s*c_e_bcs(2:3)).^2 ...
                + D_es_eff.*(obj.electrolyte_matrices.M3s*c_es + obj.electrolyte_matrices.M4s*c_e_bcs(2:3));

            c_ep_dot = dD_ep_eff.*(obj.electrolyte_matrices.M1p*c_ep + obj.electrolyte_matrices.M2p*c_e_bcs(3:4)).^2 ...
                + D_ep_eff.*(obj.electrolyte_matrices.M3p*c_ep + obj.electrolyte_matrices.M4p*c_e_bcs(3:4)) + diag(obj.electrolyte_matrices.M5p)*jp;

            % Assemble c_e_dot
            c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];


            %% Thermal Dynamics

            % State-of-Charge (Bulk)
            r_vec = (0:obj.discretization.delta_r_n:1)';
            c_n = [c_s_n(1); c_s_n; c_ss_n];
            c_p = [c_s_p(1); c_s_p; c_ss_p];
            SOC_n = 3/obj.cell_properties.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n);
            SOC_p = 3/obj.cell_properties.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p);

            % Equilibrium potentials
            [Unb] = refPotentialAnode(obj, SOC_n);
            [Upb] = refPotentialCathode(obj, SOC_p);

            % Heat generation
            % disp(cur)
            % disp(V)
            % disp((Upb - Unb))
            % pause;
            Qdot = -cur*(V - (Upb - Unb));

            % Differential equations
            T1_dot = (obj.cell_properties.thermal.h12 * (T2-T1) + Qdot) / obj.cell_properties.thermal.C1;
            T2_dot = (obj.cell_properties.thermal.h12 * (T1-T2) + obj.cell_properties.thermal.h2a*(obj.cell_properties.T_amb - T2)) / obj.cell_properties.thermal.C2;


            %% Concatenate time derivatives
            x_dot = [c_s_n_dot; c_s_p_dot; c_e_dot; T1_dot; T2_dot; delta_sei_dot];

            %% Concatenate outputs
            varargout{1} = V;
            varargout{2} = V_noVCE;
            varargout{3} = SOC_n;
            varargout{4} = SOC_p;
            varargout{5} = c_ss_n;
            varargout{6} = c_ss_p;
            varargout{7} = c_ex';


        end

        function initialize_solid_phase_matrices(obj,D_s_n,D_s_p)

            % Electrochemical Model Parameters
            alpha_n = D_s_n / (obj.cell_properties.R_s_n * obj.discretization.delta_r_n)^2;
            alpha_p = D_s_p / (obj.cell_properties.R_s_p * obj.discretization.delta_r_p)^2;

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

            N3_n = [0; -(2*obj.discretization.delta_r_n * obj.cell_properties.R_s_n)/(D_s_n)];
            N3_p = [0; -(2*obj.discretization.delta_r_p * obj.cell_properties.R_s_p)/(D_s_p)];

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

            dActivity = ((v00 + v10.*(c_e).^(0.5) + v30*(1+v31*(T - obj.cell_properties.T_ref)) .* (c_e).^(1.5))./(1-obj.cell_properties.t_plus))-1;

            if(nargout == 2)
                d_dactivity = (0.5 * v10 * (c_e).^(-0.5) + (1.5)*v30*(1+v31*(T - obj.cell_properties.T_ref))*(c_e).^(0.5))/(1-obj.cell_properties.t_plus);
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
            i_0n = k_n * ((obj.cell_properties.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^obj.cell_properties.alph;
            di_0n = k_n * (obj.cell_properties.alph*((obj.cell_properties.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^(obj.cell_properties.alph-1).*( c_e_n.*(obj.cell_properties.c_s_n_max - c_ss_n) -  c_ss_n .* c_e_n) );
            i_0p = k_p * ((obj.cell_properties.c_s_p_max - c_ss_p) .* c_ss_p .* c_e_p).^obj.cell_properties.alph;
            di_0p = k_p * (obj.cell_properties.alph*(max((obj.cell_properties.c_s_p_max - c_ss_p),0) .* c_ss_p .* c_e_p).^(obj.cell_properties.alph-1).*( c_e_p.*(obj.cell_properties.c_s_p_max - c_ss_p) -  c_ss_p .* c_e_p) );

            if(nargout > 2)
                varargout{1}=di_0n;
                varargout{2}=di_0p;
            end
        end
    end
    methods (Static)
    end

end
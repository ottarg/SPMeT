function initialize_electrolyte_matrices(obj)

            %% Lumped Coefficients
            Del_xn = obj.anode.electrode_thickness * obj.discretization.delta_x_n;
            Del_xs = obj.separator.thickness * obj.discretization.delta_x_s;
            Del_xp = obj.cathode.electrode_thickness * obj.discretization.delta_x_p;

            %% Matrices in nonlinear dynamics
            obj.electrolyte_matrices.M1n = ((diag(ones(obj.discretization.Nxn-2,1),+1) - diag(ones(obj.discretization.Nxn-2,1),-1))/(2*Del_xn));
            obj.electrolyte_matrices.M1s = ((diag(ones(obj.discretization.Nxs-2,1),+1) - diag(ones(obj.discretization.Nxs-2,1),-1))/(2*Del_xs));
            obj.electrolyte_matrices.M1p = ((diag(ones(obj.discretization.Nxp-2,1),+1) - diag(ones(obj.discretization.Nxp-2,1),-1))/(2*Del_xp));


            M2n = zeros(obj.discretization.Nxn-1,2);
            M2n(1,1) = -1/(2*Del_xn);
            M2n(end,end) = 1/(2*Del_xn);
            obj.electrolyte_matrices.M2n = (M2n);

            M2s = zeros(obj.discretization.Nxs-1,2);
            M2s(1,1) = -1/(2*Del_xs);
            M2s(end,end) = 1/(2*Del_xs);
            obj.electrolyte_matrices.M2s = (M2s);

            M2p = zeros(obj.discretization.Nxp-1,2);
            M2p(1,1) = -1/(2*Del_xp);
            M2p(end,end) = 1/(2*Del_xp);
            obj.electrolyte_matrices.M2p = (M2p);


            obj.electrolyte_matrices.M3n = ((-2*diag(ones(obj.discretization.Nxn-1,1),0) + diag(ones(obj.discretization.Nxn-2,1),+1) + diag(ones(obj.discretization.Nxn-2,1),-1))/(Del_xn^2));
            obj.electrolyte_matrices.M3s = ((-2*diag(ones(obj.discretization.Nxs-1,1),0) + diag(ones(obj.discretization.Nxs-2,1),+1) + diag(ones(obj.discretization.Nxs-2,1),-1))/(Del_xs^2));
            obj.electrolyte_matrices.M3p = ((-2*diag(ones(obj.discretization.Nxp-1,1),0) + diag(ones(obj.discretization.Nxp-2,1),+1) + diag(ones(obj.discretization.Nxp-2,1),-1))/(Del_xp^2));

            M4n = zeros(obj.discretization.Nxn-1,2);
            M4n(1,1) = 1/(Del_xn^2);
            M4n(end,end) = 1/(Del_xn^2);
            obj.electrolyte_matrices.M4n = (M4n);

            M4s = zeros(obj.discretization.Nxs-1,2);
            M4s(1,1) = 1/(Del_xs^2);
            M4s(end,end) = 1/(Del_xs^2);
            obj.electrolyte_matrices.M4s = (M4s);

            M4p = zeros(obj.discretization.Nxp-1,2);
            M4p(1,1) = 1/(Del_xp^2);
            M4p(end,end) = 1/(Del_xp^2);
            obj.electrolyte_matrices.M4p = (M4p);

            obj.electrolyte_matrices.M5n = (1-obj.transference_number)*obj.anode.specific_interfacial_area/obj.anode.volume_fraction_electrolyte * eye(obj.discretization.Nxn-1);
            obj.electrolyte_matrices.M5p = (1-obj.transference_number)*obj.cathode.specific_interfacial_area/obj.cathode.volume_fraction_electrolyte * eye(obj.discretization.Nxp-1);

            %% Boundary Conditions
            N1 = zeros(4,obj.discretization.Nx-3);
            N2 = zeros(4);

            % Boundary condition 1
            N1(1,1) =  4;
            N1(1,2) = -1;
            N2(1,1) = -3;

            % Boundary condition 2
            N1(2,obj.discretization.Nxn-2) = (obj.anode.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xn);
            N1(2,obj.discretization.Nxn-1) = (-4*obj.anode.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xn);
            N2(2,2) = (3*obj.anode.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xn) + (3*obj.separator.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xs);
            N1(2,obj.discretization.Nxn) = (-4*obj.separator.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xs);
            N1(2,obj.discretization.Nxn+1) = (obj.separator.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xs);

            % Boundary condition 3
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-3) = (obj.separator.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xs);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-2) = (-4*obj.separator.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xs);
            N2(3,3) = (3*obj.separator.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xs) + (3*obj.cathode.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xp);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs-1) = (-4*obj.cathode.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xp);
            N1(3,obj.discretization.Nxn+obj.discretization.Nxs) = (obj.cathode.volume_fraction_electrolyte^obj.bruggemann_porosity)/(2*Del_xp);

            % Boundary condition 4
            N1(4,end-1) = 1;
            N1(4,end) =  -4;
            N2(4,4) =     3;

            %%% SPARSE OUTPUT
            obj.electrolyte_matrices.C = (-N2\N1);
        end
function initialize_solid_phase_matrices(obj,anode_diffusion_coefficient,cathode_diffusion_coefficient)

% Electrochemical Model Parameters
alpha_n = anode_diffusion_coefficient / (obj.cell_properties.anode.particle_radius * obj.discretization.delta_r_n)^2;
alpha_p = cathode_diffusion_coefficient / (obj.cell_properties.cathode.particle_radius * obj.discretization.delta_r_p)^2;

% Block matrices
M1_n = zeros(obj.discretization.radial_divisions-1);
M1_p = zeros(obj.discretization.radial_divisions-1);

for idx = 1:(obj.discretization.radial_divisions-1)

    % Lower diagonal
    if(idx ~= 1)
        M1_n(idx,idx-1) = (idx-1)/idx * alpha_n;
        M1_p(idx,idx-1) = (idx-1)/idx * alpha_p;
    end

    % Main diagonal
    M1_n(idx,idx) = -2*alpha_n;
    M1_p(idx,idx) = -2*alpha_p;

    % Upper diagonal
    if(idx ~= (obj.discretization.radial_divisions-1))
        M1_n(idx,idx+1) = (idx+1)/idx * alpha_n;
        M1_p(idx,idx+1) = (idx+1)/idx * alpha_p;
    end
end

M2_n = zeros(obj.discretization.radial_divisions-1,2);
M2_p = zeros(obj.discretization.radial_divisions-1,2);

M2_n(end,end) = obj.discretization.radial_divisions/(obj.discretization.radial_divisions-1) * alpha_n;
M2_p(end,end) = obj.discretization.radial_divisions/(obj.discretization.radial_divisions-1) * alpha_p;

N1 = zeros(2,obj.discretization.radial_divisions-1);
% % 1st Order BCs
% N1(1,1) = 1;
% N1(end,end) = -1;
%
% N2 = diag([-1,1]);
%

% 2nd Order BCs
N1(1,1) =     4;
N1(1,2) =    -1;
N1(2,end) =  -4;
N1(2,end-1) = 1;

N2 = diag([-3,3]);

N3_n = [0; -(2*obj.discretization.delta_r_n * obj.cell_properties.anode.particle_radius)/(anode_diffusion_coefficient)];
N3_p = [0; -(2*obj.discretization.delta_r_p * obj.cell_properties.cathode.particle_radius)/(cathode_diffusion_coefficient)];

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
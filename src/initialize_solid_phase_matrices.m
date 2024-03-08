function model = initialize_solid_phase_matrices(model,anode_diffusion_coefficient,cathode_diffusion_coefficient)

% Electrochemical Model Parameters
alpha_n = anode_diffusion_coefficient / (model.anode.particle_radius * model.discretization.delta_r_n)^2;
alpha_p = cathode_diffusion_coefficient / (model.cathode.particle_radius * model.discretization.delta_r_p)^2;

% Block matrices
M1_n = zeros(model.discretization.radial_divisions-1);
M1_p = zeros(model.discretization.radial_divisions-1);

for idx = 1:(model.discretization.radial_divisions-1)

    % Lower diagonal
    if(idx ~= 1)
        M1_n(idx,idx-1) = (idx-1)/idx * alpha_n;
        M1_p(idx,idx-1) = (idx-1)/idx * alpha_p;
    end

    % Main diagonal
    M1_n(idx,idx) = -2*alpha_n;
    M1_p(idx,idx) = -2*alpha_p;

    % Upper diagonal
    if(idx ~= (model.discretization.radial_divisions-1))
        M1_n(idx,idx+1) = (idx+1)/idx * alpha_n;
        M1_p(idx,idx+1) = (idx+1)/idx * alpha_p;
    end
end

M2_n = zeros(model.discretization.radial_divisions-1,2);
M2_p = zeros(model.discretization.radial_divisions-1,2);

M2_n(end,end) = model.discretization.radial_divisions/(model.discretization.radial_divisions-1) * alpha_n;
M2_p(end,end) = model.discretization.radial_divisions/(model.discretization.radial_divisions-1) * alpha_p;

N1 = zeros(2,model.discretization.radial_divisions-1);
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

N3_n = [0; -(2*model.discretization.delta_r_n * model.anode.particle_radius)/(anode_diffusion_coefficient)];
N3_p = [0; -(2*model.discretization.delta_r_p * model.cathode.particle_radius)/(cathode_diffusion_coefficient)];

% A,B matrices for each electrode
model.solid_phase_matrices.A_n = M1_n - M2_n*(N2\N1);
model.solid_phase_matrices.A_p = M1_p - M2_p*(N2\N1);

model.solid_phase_matrices.B_n = M2_n*(N2\N3_n);
model.solid_phase_matrices.B_p = M2_p*(N2\N3_p);

% C,D matrices for each electrode
model.solid_phase_matrices.C_n = -[0,1]*(N2\N1);
model.solid_phase_matrices.C_p = -[0,1]*(N2\N1);

model.solid_phase_matrices.D_n = [0,1]*(N2\N3_n);
model.solid_phase_matrices.D_p = [0,1]*(N2\N3_p);

end
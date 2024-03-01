clear
close all
model = SPMe;
model.cell_properties = params_LCO;
model.discretization.Nr = 30;
model.discretization.delta_r_n = 1/model.discretization.Nr;
model.discretization.delta_r_p = 1/model.discretization.Nr;
% Finite difference points along x-coordinate
model.discretization.Nxn = 10;
model.discretization.Nxs = 5;
model.discretization.Nxp = 10;
model.discretization.Nx = model.discretization.Nxn+model.discretization.Nxs+model.discretization.Nxp;
model.discretization.delta_x_n = 1 / model.discretization.Nxn;
model.discretization.delta_x_s = 1 / model.discretization.Nxs;
model.discretization.delta_x_p = 1 / model.discretization.Nxp;

data = load('input-data/UDDS');
model.initial_voltage = data.voltage(1);
model.initialize
[res,x] = model.simulate(data.time,data.current);

figure
plot(res.time,res.V)
hold on
plot(data.time,data.voltage)

% figure
% hold on
% plot(res.time,res.SOC_n)
% plot(res.time,res.SOC_p)
% 
% figure
% plot(x(:,81))
% hold on
% plot(x(:,82))
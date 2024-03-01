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

% data = load('input-data/UDDS');


% model.initial_voltage = data.voltage(1);
model.initial_voltage = 3.6;
model.initialize
data.time = 1:0.5:5000;
data.current = 0.*data.time+3;
data.current(1:1000) = 0;
data.current(end-1000:end) = 0;
[res,x] = model.simulate(data.time,data.current);

figure
hold on
plot(res.time,res.V)
plot(res.time,res.OCV)

% plot(data.time,data.voltage)

figure
hold on
plot(res.time,res.SOC_n)
plot(res.time,res.SOC_p)

% figure
% plot(x(:,81))
% hold on
% % plot(x(:,82))
% figure
% for i=1:1:length(x)
%     plot([x(i,1:29),x(i,59:80),x(i,58:-1:30)])
%     drawnow
% end

figure
plot(res.anode_potential)
hold on
plot(res.cathode_potential)
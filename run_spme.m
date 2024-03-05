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


% model.initial_voltage = data.voltage(1);
model.initial_voltage = 4.15;
time = 0:0.001:120;
current = 0.*time-10;
current(1:find(time>10,1,"first")) = 0;
current(find(time>40,1,"first"):end) = 0;
temperature = 0.*time+30;
tic
[res,X] = model.simulate(time,current,temperature);

figure
subplot(2,2,1)
hold on
plot(res.timeODE,res.SOC_n)
plot(res.timeODE,res.SOC_p)
subplot(2,2,2)
hold on
plot(res.anode_potential)
plot(res.cathode_potential)
subplot(2,2,3)
plot(res.timeODE,res.V)
plot(res.timeODE,res.OCV)
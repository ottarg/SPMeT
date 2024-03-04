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
temperature = 0.*time+20;
tic
[res60,X] = model.simulate(time,current,temperature);
temperature = 0.*time+10;
% [res10,~] = model.simulate(time,current,temperature);
toc
figure
hold on
plot(res60.timeODE,res60.V)
plot(res60.timeODE,res60.OCV)
% plot(res10.timeODE,res10.V)
% plot(res10.timeODE,res10.OCV)

figure
hold on
plot(res60.timeODE,res60.SOC_n)
plot(res60.timeODE,res60.SOC_p)

figure
hold on
plot(res60.anode_potential)
plot(res60.cathode_potential)
% plot(res10.anode_potential)
% plot(res10.cathode_potential)
% 
% 
% %%
% 
% data = load('input-data/UDDS');
% model.initial_voltage = data.voltage(1);
% tic
% [res,~] = model.simulate(data.time,data.current,data.temp);
% toc
% figure
% hold on
% plot(res.timeODE,res.V)

figure
plot(X(:,end))
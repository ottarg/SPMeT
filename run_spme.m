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
model.initial_voltage = 3.6;
time = 0:0.02:1000;
current = 0.*time+3;
% current(1:1000) = 0;
current(end-1000:end) = 0;
temperature = 0.*time+30;
tic
% [res,~] = model.simulate(time,current,temperature);

model.initialize;
res.time = time;
res.current = -current/model.cell_properties.electrode_area*10;
res.temperature = temperature;
Opt    = odeset('Events',@(t,x)detectImagSolution(model,t,x,res));


% 
% 
% for i = 1:100
%     [res.timeODE,x] = ode23s(@(t,x) spme_ode(model,t,x,res),[res.time(i),res.time(i+2)],model.x0,Opt);
%     obj.x0 = x;
%     for k = 1:length(res.timeODE)
%     % Compute outputs
%     [~,res.V(k),res.V_spm(k),res.SOC_n(k),res.SOC_p(k),res.c_ss_n(k),res.c_ss_p(k),res.c_e(:,k),res.OCV(:,k),res.anode_potential(:,k),res.cathode_potential(:,k)] = ...
%         spme_ode(model,res.timeODE(k),x(k,:)',res);
%     end
%     V(i) = res.V(end);
% end
% 
% plot(res.time(1:10),V)
% 



% 
% 
% toc
% figure
% hold on
% plot(res.timeODE,res.V)
% plot(res.timeODE,res.OCV)
% 
% figure
% hold on
% plot(res.timeODE,res.SOC_n)
% plot(res.timeODE,res.SOC_p)
% 
% figure
% plot(res.anode_potential)
% hold on
% plot(res.cathode_potential)
% 
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
% plot(data.time,data.voltage)
% 
% % 
% % %% Fake discrete long timestep
% % 
% % 
% % time = 0:1:5000;
% % current = 0.*time+3;
% % current(1:1000) = 0;
% current(end-1000:end) = 0;
% temperature = 0.*time+30;
% 
% 
% [res,x] = model.simulate(time,current,temperature);
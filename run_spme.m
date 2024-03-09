clear
close all

mdl = SPME_assumps_LCO;
mdl.discretization.radial_divisions = 30;
mdl.discretization.Nxn = 10;
mdl.discretization.Nxs = 5;
mdl.discretization.Nxp = 10;




load('C:\Users\ogislason\Documents\MATLAB\SPMeT\input-data\UDDS.mat')
mdl.initial_voltage = voltage(1);
% time = 0:0.1:120;
% current = 0.*time+20;
% current(1:find(time>10,1,"first")) = 0;
% current(find(time>40,1,"first"):end) = 0;
% temperature = 0.*time+30;

mdl.initialize
model = mdl.getStruct();
stopTime = time(end);
simIn = Simulink.SimulationInput("CellSim");
load_system("CellSim");
inDS = createInputDataset("CellSim");
inDS{1} = timeseries(-current,time,'Name',inDS{1}.name);
inDS{2} = timeseries(temp,time,'Name',inDS{2}.name);
simIn = setExternalInput(simIn,inDS);
out = sim(simIn);
res = out.logsout.extractTimetable;

% 
% figure
% for i = 1:height(x)
%     plot([x(i,1:29),x(i,59:80),x(i,58:-1:30)])
%     drawnow
%     pause(0.1)
% end
figure
subplot(2,2,1)
hold on
plot(res.Time,res.anode_SOC,'DisplayName','Anode')
plot(res.Time,res.cathode_SOC,'DisplayName','Cathode')
superLabel('Time (s)','SOC','',1)
axis tight
subplot(2,2,2)
hold on
plot(res.Time,res.anode_potential,'DisplayName','Anode')
plot(res.Time,res.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
axis tight
subplot(2,2,3)
hold on
plot(res.Time,res.CCV,'DisplayName','OCV')
plot(res.Time,res.OCV,'DisplayName','OCV')
superLabel('Time (s)','Voltage (V)','',1)
axis tight


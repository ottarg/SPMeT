% clear
% close all

mdl = SPME_assumps_LCO;
mdl.discretization.radial_divisions = 80;
mdl.discretization.Nxn = 15;
mdl.discretization.Nxs = 15;
mdl.discretization.Nxp = 15;
mdl.maximum_voltage = 4.55;
mdl.minimum_voltage = 2.01;
mdl.initial_voltage = mdl.maximum_voltage;
mdl.initialize
model = mdl.getStruct();

time = 0:0.01:3540;
stopTime = time(end);
current = mdl.capacity+time.*0;
temp = 25+time.*0;
simIn = Simulink.SimulationInput("CellSim");
load_system("CellSim");
inDS = createInputDataset("CellSim");
inDS{1} = timeseries(current,time,'Name',inDS{1}.name);
inDS{2} = timeseries(temp,time,'Name',inDS{2}.name);
simIn = setExternalInput(simIn,inDS);
out = sim(simIn);
res = out.logsout.extractTimetable;

figure
subplot(2,2,1)
hold on
plot(res.Time,res.CCV,'DisplayName','Model CCV')
plot(res.Time,res.OCV,'DisplayName','OCV')
superLabel('Time (s)','Voltage (V)','',1)
axis tight
subplot(2,2,2)
hold on
plot(res.Time,-res.current,'DisplayName','Current')
superLabel('Time (s)','Current (A)','',0)
axis tight
subplot(2,2,3)
hold on
plot(res.Time,res.anode_SOC,'DisplayName','Anode')
plot(res.Time,res.cathode_SOC,'DisplayName','Cathode')
superLabel('Time (s)','SOC','',1)
axis tight
subplot(2,2,4)
hold on
plot(res.Time,res.anode_potential,'DisplayName','Anode')
plot(res.Time,res.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
axis tight

% 
% f = figure;
% plot(gca,[res.concentration(1,1:29),res.concentration(1,59:80),res.concentration(1,58:-1:30)],'LineWidth',2)
% superLabel('','Concentration (mol/m^3)','',0)
% for i = 1:height(res.concentration)
% f.Children.Children.YData = [res.concentration(i,1:29),res.concentration(i,59:80),res.concentration(i,58:-1:30)];
%     drawnow
%     pause(0.1)
% end
% 
% figure
% hold on
% plot(res.Time,res.anode_SOC)
% plot(res.Time,res.cathode_SOC)
% 
% 
% figure
% plot(refPotentialCathode(0:0.01:1)-refPotentialAnode(1:-0.01:0))
% figure
% hold on
% plot(refPotentialAnode(1:-0.01:0))
% plot(refPotentialCathode(0:0.01:1))


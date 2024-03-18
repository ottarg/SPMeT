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
loggingStep = 1;
time = 0:0.1:3000;
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
simulinkRes = out.logsout.extractTimetable;

figure
subplot(2,2,1)
hold on
plot(simulinkRes.Time,simulinkRes.CCV,'DisplayName','Model CCV')
plot(simulinkRes.Time,simulinkRes.OCV,'DisplayName','OCV')
superLabel('Time (s)','Voltage (V)','',1)
axis tight
subplot(2,2,2)
hold on
plot(simulinkRes.Time,simulinkRes.current,'DisplayName','Current')
superLabel('Time (s)','Current (A)','',0)
axis tight
subplot(2,2,3)
hold on
plot(simulinkRes.Time,simulinkRes.anode_SOC,'DisplayName','Anode')
plot(simulinkRes.Time,simulinkRes.cathode_SOC,'DisplayName','Cathode')
superLabel('Time (s)','SOC','',1)
axis tight
subplot(2,2,4)
hold on
plot(simulinkRes.Time,simulinkRes.anode_potential,'DisplayName','Anode')
plot(simulinkRes.Time,simulinkRes.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
axis tight





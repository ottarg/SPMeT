% clear
% close all

mdl = SPME_assumps_LCO;
mdl.discretization.radial_divisions = 80;
mdl.discretization.Nxn = 15;
mdl.discretization.Nxs = 15;
mdl.discretization.Nxp = 15;
mdl.maximum_voltage = 4.25;
mdl.minimum_voltage = 3;
mdl.initial_voltage = mdl.minimum_voltage;
mdl.initialize  
model = mdl.getStruct();
loggingStep = 1;
time = 0:0.1:3500;
stopTime = time(end);
current = -1*mdl.capacity+time.*0;
temp = 25+time.*0;
simIn = Simulink.SimulationInput("fastchargeSim");
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
plot(simulinkRes.Time,-simulinkRes.current,'DisplayName','Current')
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
plot(simulinkRes.Time,1e3.*simulinkRes.anode_potential,'DisplayName','Anode')

superLabel('Time (s)','Anode Potential (mV)','',1)
axis tight





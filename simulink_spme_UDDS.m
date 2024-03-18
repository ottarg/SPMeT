% clear
% close all
load('C:\Users\ogislason\Documents\MATLAB\SPMeT\input-data\UDDS.mat')
testvoltage = timeseries(voltage,time,"Name",'Voltage');


mdl = SPME_assumps_LCO;
mdl.discretization.radial_divisions = 30;
mdl.discretization.Nxn = 10;
mdl.discretization.Nxs = 5;
mdl.discretization.Nxp = 10;
mdl.initial_voltage = voltage(1);
mdl.initialize
model = mdl.getStruct();

profile_scaling = 0.59/(abs(trapz(time,current)./3600)/mdl.capacity); % Scale to use 30% SOC
stopTime = time(end);
simIn = Simulink.SimulationInput("CellSim");
load_system("CellSim");
inDS = createInputDataset("CellSim");
inDS{1} = timeseries(-profile_scaling*current,time,'Name',inDS{1}.name);
inDS{2} = timeseries(temp-5,time,'Name',inDS{2}.name);
simIn = setExternalInput(simIn,inDS);
out = sim(simIn);
simulinkRes = out.logsout.extractTimetable;

figure
subplot(2,2,1)
hold on
plot(simulinkRes.Time,simulinkRes.CCV,'DisplayName','Model CCV')
plot(testvoltage,'--k','DisplayName','Test CCV')
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
plot(simulinkRes.Time,simulinkRes.anode_potential,'DisplayName','Anode')
plot(simulinkRes.Time,simulinkRes.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
axis tight
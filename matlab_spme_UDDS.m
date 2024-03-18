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
profile_scaling = 0.59/(abs(trapz(time,current)./3600)/mdl.capacity); 
[matlabRes,concentration] = mdl.simulate(time,-profile_scaling*current,temp);


figure
subplot(2,2,1)
hold on
plot(matlabRes.time,matlabRes.V,'DisplayName','Model CCV')
plot(testvoltage,'--k','DisplayName','Test CCV')
plot(matlabRes.time,matlabRes.OCV,'DisplayName','OCV')
superLabel('Time (s)','Voltage (V)','',1)
axis tight
subplot(2,2,2)
hold on
plot(time,current,'DisplayName','Current')
superLabel('Time (s)','Current (A)','',0)
axis tight
subplot(2,2,3)
hold on
plot(matlabRes.time,matlabRes.an,'DisplayName','Anode')
plot(matlabRes.time,matlabRes.SOC_p,'DisplayName','Cathode')
superLabel('Time (s)','SOC','',1)
axis tight
subplot(2,2,4)
hold on
plot(matlabRes.time,matlabRes.anode_potential,'DisplayName','Anode')
plot(matlabRes.time,matlabRes.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
axis tight

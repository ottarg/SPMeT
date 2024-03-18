clear
% close all
load('C:\Users\ogislason\Documents\MATLAB\SPMeT\input-data\UDDS.mat')
testvoltage = timeseries(voltage,time,"Name",'Voltage');

mdl = SPME_assumps_LCO;
mdl.discretization.radial_divisions = 30;
mdl.discretization.Nxn = 10;
mdl.discretization.Nxs = 5;
mdl.discretization.Nxp = 10;
mdl.initial_voltage = voltage(1);
profile_scaling = 0.295/(abs(trapz(time,current)./3600)/mdl.capacity); % Scale to use 30% SOC
[res,concentration] = mdl.simulate(time,profile_scaling*current,temp);


figure
subplot(2,2,1)
hold on
plot(res.time,res.V,'DisplayName','Model CCV')
plot(testvoltage,'--k','DisplayName','Test CCV')
plot(res.time,res.OCV,'DisplayName','OCV')
superLabel('Time (s)','Voltage (V)','',1)
axis tight
subplot(2,2,2)
hold on
plot(time,current,'DisplayName','Current')
superLabel('Time (s)','Current (A)','',0)
axis tight
subplot(2,2,3)
hold on
plot(res.time,res.SOC_n,'DisplayName','Anode')
plot(res.time,res.SOC_p,'DisplayName','Cathode')
superLabel('Time (s)','SOC','',1)
axis tight
subplot(2,2,4)
hold on
plot(res.time,res.anode_potential,'DisplayName','Anode')
plot(res.time,res.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
axis tight


figure
for i = 1:width(res.electrolyte_concentrations)
    plot(res.electrolyte_concentrations(:,i))
    drawnow
    pause(0.1)
end
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

time = 0:0.1:3000;
stopTime = time(end);
current = mdl.capacity+time.*0;
temp = 25+time.*0;

matlabRes =  mdl.simulate(time,current,temp);
figure
subplot(2,2,1)
hold on
plot(matlabRes.time,matlabRes.V,'DisplayName','Model CCV')
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
plot(matlabRes.time,matlabRes.SOC_n,'DisplayName','Anode')
plot(matlabRes.time,matlabRes.SOC_p,'DisplayName','Cathode')
superLabel('Time (s)','SOC','',1)
axis tight
subplot(2,2,4)
hold on
plot(matlabRes.time,matlabRes.anode_potential,'DisplayName','Anode')
plot(matlabRes.time,matlabRes.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
axis tight



clear
close all

model = SPMe;
model.cell_properties = params_LCO;
model.discretization.Nr = 30;

model.discretization.Nxn = 10;
model.discretization.Nxs = 5;
model.discretization.Nxp = 10;


model.initial_voltage = 4.15;
time = 0:0.1:120;
current = 0.*time-10;
current(1:find(time>10,1,"first")) = 0;
current(find(time>40,1,"first"):end) = 0;
temperature = 0.*time+30;
[res,x] = model.simulate(time,current,temperature);

figure
for i = 1:height(x)
    plot([x(i,1:29),x(i,59:80),x(i,58:-1:30)])
    drawnow
    pause(0.1)
end
figure
subplot(2,2,1)
hold on
plot(res.time,res.SOC_n,'DisplayName','Anode')
plot(res.time,res.SOC_p,'DisplayName','Cathode')
superLabel('Time (s)','SOC','',1)
subplot(2,2,2)
hold on
plot(res.anode_potential,'DisplayName','Anode')
plot(res.cathode_potential,'DisplayName','Cathode')
superLabel('Time (s)','Voltage (V)','',1)
subplot(2,2,3)
hold on
plot(res.time,res.V,'DisplayName','OCV')
plot(res.time,res.OCV,'DisplayName','OCV')
superLabel('Time (s)','Voltage (V)','',1)



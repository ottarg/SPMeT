% clearvars; close all
tic
load('UDDS.mat')
testvoltage = timeseries(voltage,time,"Name",'Voltage');

mdl = LCO_parameters;
mdl.initial_voltage = voltage(1);
mdl.initialize

model = mdl.getStruct();
profile_scaling = 0.59/(abs(trapz(time,current)./3600)/mdl.capacity);
stopTime = time(end);
simIn = Simulink.SimulationInput("SPMeT_System");
load_system("SPMeT_System");
inDS = createInputDataset("SPMeT_System");
inDS{1} = timeseries(-profile_scaling*current,time,'Name',inDS{1}.name);
inDS{2} = timeseries(temp-5,time,'Name',inDS{2}.name);
simIn = setExternalInput(simIn,inDS);
out = sim(simIn);
simulinkRes = out.logsout.extractTimetable;

figure
subplot(2,2,1)
hold on
plot(simulinkRes.Time,simulinkRes.CCV,'DisplayName','Model CCV')
plot(testvoltage,'--k','DisplayName','Test data')
plot(simulinkRes.Time,simulinkRes.OCV,'DisplayName','OCV')
xlabel('Time (s)')
ylabel('Voltage (V)')
legend Location Best
axis tight
grid on; box on;
subplot(2,2,2)
hold on
plot(simulinkRes.Time,-simulinkRes.current,'DisplayName','Current')
xlabel('Time (s)')
ylabel('Current (A)')
axis tight
grid on; box on;
subplot(2,2,3)
hold on
plot(simulinkRes.Time,simulinkRes.anode_SOC,'DisplayName','Anode')
plot(simulinkRes.Time,simulinkRes.cathode_SOC,'DisplayName','Cathode')
xlabel('Time (s)')
ylabel('SOC')
axis tight
grid on; box on;
subplot(2,2,4)
hold on
plot(simulinkRes.Time,simulinkRes.anode_potential,'DisplayName','Anode')
plot(simulinkRes.Time,simulinkRes.cathode_potential,'DisplayName','Cathode')
xlabel('Time (s)')
ylabel('Voltage (V)')
axis tight
grid on; box on;
plotPosition = [2 2 8 6];
plotPaperPosition = [2 2 8 6];
plotFontSize = 14;
set(gcf,'Units','inches')    
set(gcf,'Position', plotPosition)              
set(gcf,'PaperPosition', plotPaperPosition) 
set(gcf,'color','w')
ax = findall(gcf,'type','axes');
linewd = 2;
for ppp=1:length(ax)
    set(get(ax(ppp),'Title'),'fontweight','bold')
    set(get(ax(ppp),'Xlabel'),'fontweight','bold','interpreter','tex','fontsize',plotFontSize)   
    set(get(ax(ppp),'Ylabel'),'fontweight','bold','interpreter','tex','fontsize',plotFontSize)   
    set(get(ax(ppp),'Zlabel'),'fontweight','bold','interpreter','tex','fontsize',plotFontSize)   
    set(get(ax(ppp),'Title'),'fontsize',16,'interpreter','tex')
    if ~isempty(findobj(gcf,'Type','Legend'))
        set(legend,'interpreter','tex') %'fontsize',14,
    end
    try
        set(get(ax(ppp),'Children'),'linewidth',linewd) 
    catch
    end
end
display("Simulated UDDS cycle in " + num2str(toc) + " seconds")
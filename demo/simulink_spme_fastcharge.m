clearvars; close all
tic
mdl = LCO_parameters;
anode_potential_target = 35;
mdl.initial_voltage = mdl.minimum_voltage;
mdl.initialize  
model = mdl.getStruct();

time = 0:0.1:3600;
stopTime = time(end);
temp = 45+time.*0;
simIn = Simulink.SimulationInput("fastchargeSim");
load_system("fastchargeSim");
inDS = createInputDataset("fastchargeSim");
inDS{1} = timeseries(temp,time,'Name',inDS{1}.name);
simIn = setExternalInput(simIn,inDS);
out = sim(simIn);
simulinkRes = out.logsout.extractTimetable;

figure
subplot(2,2,1)
hold on
plot(simulinkRes.Time,simulinkRes.CCV,'DisplayName','Model CCV')
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
legend Location Best
axis tight
grid on; box on;
subplot(2,2,4)
hold on
plot(simulinkRes.Time,1e3.*simulinkRes.anode_potential,'DisplayName','Anode')
xlabel('Time (s)')
ylabel('Anode Potential (mV)')
axis tight
ylim([0,100])
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
figure
plot(simulinkRes.CCV,-simulinkRes.current,'DisplayName','Model CCV')
display("Simulated fast charge in " + num2str(toc) + " seconds")



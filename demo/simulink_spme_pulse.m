clear
close all
tic
mdl = LCO_parameters;
mdl.initial_voltage = 4.1;
mdl.initialize
model = mdl.getStruct();
loggingStep = 1;
time = 0:0.1:240;
stopTime = time(end);
current = mdl.capacity+time.*0;
current(time<20) = 0;
current(time>50) = 0;
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

display("Simulated discharge in " + num2str(toc) + " seconds")
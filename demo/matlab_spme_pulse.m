clear
close all

mdl = LCO_parameters;
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
xlabel('Time (s)')
ylabel('Voltage (V)')
legend Location Best
axis tight
grid on; box on;
subplot(2,2,2)
hold on
plot(time,current,'DisplayName','Current')
xlabel('Time (s)')
ylabel('Current (A)')
axis tight
grid on; box on;
subplot(2,2,3)
hold on
plot(matlabRes.time,matlabRes.SOC_n,'DisplayName','Anode')
plot(matlabRes.time,matlabRes.SOC_p,'DisplayName','Cathode')
xlabel('Time (s)')
ylabel('SOC')
axis tight
grid on; box on;
subplot(2,2,4)
hold on
plot(matlabRes.time,matlabRes.anode_potential,'DisplayName','Anode')
plot(matlabRes.time,matlabRes.cathode_potential,'DisplayName','Cathode')
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



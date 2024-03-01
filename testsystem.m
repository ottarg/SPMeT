clear
model = SPMeSystem;

model.cell_properties = params_LCO;
model.discretization.Nr = 30;
model.discretization.delta_r_n = 1/model.discretization.Nr;
model.discretization.delta_r_p = 1/model.discretization.Nr;
% Finite difference points along x-coordinate
model.discretization.Nxn = 10;
model.discretization.Nxs = 5;
model.discretization.Nxp = 10;
model.discretization.Nx = model.discretization.Nxn+model.discretization.Nxs+model.discretization.Nxp;
model.discretization.delta_x_n = 1 / model.discretization.Nxn;
model.discretization.delta_x_s = 1 / model.discretization.Nxs;
model.discretization.delta_x_p = 1 / model.discretization.Nxp;






data = load('input-data/UDDS');
model.time_step = 0.001;
model.initial_voltage = data.voltage(1);
model.setup
t = 0;
for i =1:(500/model.time_step)

    current = interp1(data.time,data.current,t(i),[]);
    temperature = interp1(data.time,data.temp,t(i),[])+273.15;
%     [V,V_spm,SOC_n,SOC_p,c_ss_n,c_ss_p,~,OCV,anode_potential,cathode_potential]=model.step(current,temperature)
    [V(i,:),V_spm(i,:),SOC_n(i,:),SOC_p(i,:),c_ss_n(i,:),c_ss_p(i,:),~,OCV(i,:),anode_potential(i,:),cathode_potential(i,:)]=model.step(current,temperature);
    t(i+1) = t(i)+model.time_step;
end
t(end) = [];

figure
plot(t,V)
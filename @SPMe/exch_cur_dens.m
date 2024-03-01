function [i_0n,i_0p,varargout] = exch_cur_dens(obj,k_p,k_n,c_ss_n,c_ss_p,c_e)

% Parse out concentrations in anode and cathode
ce_interp=interp1(linspace(0,1,length(c_e)),c_e,linspace(0,1,obj.discretization.Nx),'linear');
c_e_n = ce_interp(1:obj.discretization.Nxn-1);
c_e_p = ce_interp(obj.discretization.Nxn+obj.discretization.Nxs-1:end);


% Compute exchange current density
i_0n = k_n * ((obj.cell_properties.anode.maximum_concentration - c_ss_n) .* c_ss_n .* c_e_n).^obj.cell_properties.charge_transfer_coefficient;
di_0n = k_n * (obj.cell_properties.charge_transfer_coefficient*((obj.cell_properties.anode.maximum_concentration - c_ss_n) .* c_ss_n .* c_e_n).^(obj.cell_properties.charge_transfer_coefficient-1).*( c_e_n.*(obj.cell_properties.anode.maximum_concentration - c_ss_n) -  c_ss_n .* c_e_n) );
i_0p = k_p * ((obj.cell_properties.cathode.maximum_concentration - c_ss_p) .* c_ss_p .* c_e_p).^obj.cell_properties.charge_transfer_coefficient;
di_0p = k_p * (obj.cell_properties.charge_transfer_coefficient*(max((obj.cell_properties.cathode.maximum_concentration - c_ss_p),0) .* c_ss_p .* c_e_p).^(obj.cell_properties.charge_transfer_coefficient-1).*( c_e_p.*(obj.cell_properties.cathode.maximum_concentration - c_ss_p) -  c_ss_p .* c_e_p) );

if(nargout > 2)
    varargout{1}=di_0n;
    varargout{2}=di_0p;
end
end
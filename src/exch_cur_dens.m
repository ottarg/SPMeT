function [i_0n,i_0p,varargout] = exch_cur_dens(model,k_p,k_n,c_ss_n,c_ss_p,c_e)

% Parse out concentrations in anode and cathode
ce_interp=interp1(linspace(0,1,length(c_e)),c_e,linspace(0,1,model.discretization.Nx),'linear');
c_e_n = ce_interp(1:model.discretization.Nxn-1);
c_e_p = ce_interp(model.discretization.Nxn+model.discretization.Nxs-1:end);


% Compute exchange current density
i_0n = k_n * ((model.anode.maximum_concentration - c_ss_n) .* c_ss_n .* c_e_n).^model.charge_transfer_coefficient;
di_0n = k_n * (model.charge_transfer_coefficient*((model.anode.maximum_concentration - c_ss_n) .* c_ss_n .* c_e_n).^(model.charge_transfer_coefficient-1).*( c_e_n.*(model.anode.maximum_concentration - c_ss_n) -  c_ss_n .* c_e_n) );
i_0p = k_p * ((model.cathode.maximum_concentration - c_ss_p) .* c_ss_p .* c_e_p).^model.charge_transfer_coefficient;
di_0p = k_p * (model.charge_transfer_coefficient*(max((model.cathode.maximum_concentration - c_ss_p),0) .* c_ss_p .* c_e_p).^(model.charge_transfer_coefficient-1).*( c_e_p.*(model.cathode.maximum_concentration - c_ss_p) -  c_ss_p .* c_e_p) );

if(nargout > 2)
    varargout{1}=di_0n;
    varargout{2}=di_0p;
end
end
function [rho_eg_ft,rho_ge_ft]=rho_eg_f_HEOM_NOA(om_1_rng,rho_eg,sup_op_eg,...
                                                rho_ge,sup_op_ge)
%  This version is used for when we don't need to consider pathways
%  seperately (No orientational averaging) 
% Input is range of om_1 then V_eg rho_eq and the operator described time
% evolution, for the otherside rho_eq V_ge etc

rho_eg_ft = zeros( length(rho_eg),length(om_1_rng)); 

ID = 1i*speye(length(sup_op_eg));

if nargout == 2; rho_ge_ft =  rho_eg_ft; end 

for lp = 1:length(om_1_rng)
        
 sup_op_eg_freq = -(ID*om_1_rng(lp) + sup_op_eg);
 rho_eg_ft(:,lp) = sup_op_eg_freq\rho_eg; %solves LA system
    if nargout == 2; sup_op_ge_freq = (ID*om_1_rng(lp) - sup_op_ge);
        rho_ge_ft(:,lp) = sup_op_ge_freq\rho_ge; %solves LA system
    end
end
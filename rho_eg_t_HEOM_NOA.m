function [rho_eg_t,rho_ge_t] = rho_eg_t_HEOM_NOA...
                                (t1_rng,rho_eg,sup_op_eg, rho_ge,sup_op_ge)

 %  This version is used for when we don't need to consider pathways
%  seperately (No orientational averaging) 
% Input is range of om_1 then V_eg rho_eq and the operator described time
% evolution, for the otherside rho_eq V_ge etc                               
% This version works in time space

if ~exist('tol','var')
    tol = [];
end
   
%    rho_eg_t = zeros(length(rho_eg),length(t1_rng));
%    if nargout ==2;  rho_ge_t = rho_eg_t; end

    de_fn = @(t,v) sup_op_eg*v;
    de_fn2 = @(t,v) sup_op_ge*v;

if t1_rng == 0; rho_eg_t = rho_eg;
else; rho_eg_t = OD_wrapper(t1_rng,de_fn,rho_eg,[],'ode45',tol);
end

if nargout == 2; 
if t1_rng == 0; rho_ge_t = rho_ge;
else; rho_ge_t = OD_wrapper(t1_rng,de_fn2,rho_ge,[],'ode45',tol);
end    

end
end
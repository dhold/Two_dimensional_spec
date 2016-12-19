function [rho_eg_t,rho_ge_t]=rho_ge_t_Redfield(N,Nv,rho_eq,...
                                t1_rng,om_scale,sup_op_eg,sup_op_ge)
%  Calculates the time dependence of rho_eg and rho_ge
   
    rho_eg_t = zeros(N*Nv^2,length(t1_rng),N);
    rho_ge_t = rho_eg_t;
    
for j = 1:N %loop over sites
    VV_eg = sparse((j-1)*Nv+1,1:Nv,ones(Nv,1),N*Nv,Nv); 
    VV_ge = sparse(1:Nv,(j-1)*Nv+1,ones(Nv,1),Nv,N*Nv);

   V_eg_gg_right = kron(speye(Nv),VV_eg); 
   V_ge_gg_left = kron(VV_ge.',speye(Nv));
   
rho_0_lio_eg = V_eg_gg_right*rho_eq;   %rho_eq is equilibirum denmat
rho_0_lio_ge = V_ge_gg_left*rho_eq;   


  sup_op_eg_scale = (sup_op_eg - 1i*speye(length(sup_op_eg))*om_scale(j));        
 de_fn = @(t,v) sup_op_eg_scale*v;
 
 rho_eg_t(:,:,j) = OD_wrapper(t1_rng, de_fn,rho_0_lio_eg);
 
 if nargout ==2
    sup_op_ge_scale = (sup_op_ge + 1i*speye(length(sup_op_ge))*om_scale(j)); 
    de_fn = @(t,v) sup_op_ge_scale*v;

     rho_ge_t(:,:,j) = OD_wrapper(t1_rng, de_fn,rho_0_lio_ge);    
 end
end

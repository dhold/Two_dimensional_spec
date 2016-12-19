function [rho_eg_ft,rho_ge_ft]=rho_ge_f_Redfield...
                                (N,Nv,rho_eq,om_1_rng,sup_op_eg,sup_op_ge)
%  Simplified version, works in exciton basis and assumes no explicit 
%   vibrations are included 
%  om_scale = typical frequency of exciton transfer
%  L_eg = -1i[H_ex_vib,rho_{eg}] + L(rho_{eg}) + R(rho_{eg})
   
    rho_eg_ft = zeros(N*Nv^2,length(om_1_rng),N);
    rho_ge_ft = rho_eg_ft;
    
for j = 1:N %loop over sites
    VV_eg = sparse((j-1)*Nv+1,1:Nv,ones(Nv,1),N*Nv,Nv); 
    VV_ge = sparse(1:Nv,(j-1)*Nv+1,ones(Nv,1),Nv,N*Nv);

   V_eg_gg_right = kron(speye(Nv),VV_eg); 
   V_ge_gg_left = kron(VV_ge.',speye(Nv));
   
rho_0_lio_eg = V_eg_gg_right*rho_eq;   %rho_eq is equilibirum denmat
rho_0_lio_ge = V_ge_gg_left*rho_eq;   

       for lp = 1:length(om_1_rng)
        
 sup_op_ge_freq = (1i*speye(length(sup_op_ge))*om_1_rng(lp) - sup_op_ge); 
 sup_op_eg_freq = (-1i*speye(length(sup_op_eg))*om_1_rng(lp) - sup_op_eg); 

 rho_eg_ft(:,lp,j) = sup_op_eg_freq\rho_0_lio_eg; %solves LA system
 rho_ge_ft(:,lp,j) = sup_op_ge_freq\rho_0_lio_ge;

       end

end

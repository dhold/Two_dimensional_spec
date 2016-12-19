function [rho_eg_ft,rho_ge_ft] = rho_eg_f_HEOM3(N,Nv,V_op,rho_0,...
                                om_1_rng,sup_op_eg,sup_op_ge,HL)
%  Works in basis specified and assumes explicit vibrations are included 
%  om_scale = typical frequency of exciton transfer
%  L_eg is -1i[H,rho_{eg}], plus extra acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously.  If HL = 1, this still works, just isn't a HEOM
%  thing.
   
    rho_eg_ft = zeros(N*Nv^2*HL,length(om_1_rng),N);
    if nargout ==2;
    rho_ge_ft = zeros(N*Nv^2*HL,length(om_1_rng),N); end
    
for j = 1:size(V_op,3) %loop over sites

rho_0_lio = zeros(N*Nv*HL,1); 
tmp = V_op(:,:,j)*rho_0;    
rho_0_lio(1:N*Nv^2,1) = reshape(tmp,[N*Nv^2,1]);    %initially in equilibirum
if narout == 2;  tmp2 = rho_0*V_op(:,:,j)';
rho_0_lio2(1:N*Nv^2,1) = reshape(tmp2,[N*Nv^2,1]); end

       for lp = 1:length(om_1_rng)
      
 if nargout ==2;          
 sup_op_ge_freq = (1i*speye(N*HL)*om_1_rng(lp) - sup_op_ge); end
 sup_op_eg_freq = -(1i*speye(N*HL)*om_1_rng(lp) + sup_op_eg);
 
 if nargout ==2; 
 rho_ge_ft(:,lp,j) = sup_op_eg_freq\rho_0_lio2;  end
 rho_eg_ft(:,lp,j) = sup_op_ge_freq\rho_0_lio;   %solves LA system 


       end

end
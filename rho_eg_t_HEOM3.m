function [rho_eg_t,rho_ge_t] = rho_eg_t_HEOM3(N,Nv,V_op,rho_0,...
                                t1_rng,sup_op_eg,sup_op_ge,HL,tol)
%  Works in basis specified and assumes explicit vibrations are included 
%  om_scale = typical frequency of exciton transfer
%  L_eg is -1i[H,rho_{eg}], plus extra acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously.  If HL = 1, this still works, just isn't a HEOM
%  thing.  
% This version works in time space
if ~exist('tol','var')
    tol = [];
end
   
    rho_eg_t = zeros(N*Nv^2*HL,length(t1_rng),N);
    if nargout ==2;  rho_ge_t = zeros(N*Nv^2*HL,length(t1_rng),N); end

    de_fn = @(t,v) sup_op_eg*v;
    de_fn2 = @(t,v) sup_op_ge*v;
for j = 1:size(V_op,3) %loop over sites or ex-vib states
  
rho_0_lio = zeros(N*Nv*HL,1);  %preallocate
tmp = V_op(:,:,j)*rho_0;      %initially in equilibirum
rho_0_lio(1:N*Nv^2,1) = reshape(tmp,[N*Nv^2,1]);   %add non zero elements
if t1_rng == 0; rho_eg_t(:,1,j) = rho_0_lio;
else; rho_eg_t(:,:,j) = OD_wrapper(t1_rng,de_fn,rho_0_lio,[],'ode45',tol);
end

if nargout == 2;  tmp2 = rho_0*V_op(:,:,j)';
    rho_0_lio2= zeros(N*Nv*HL,1); 
    rho_0_lio2(1:N*Nv^2,1) = reshape(tmp2,[N*Nv^2,1]);
if t1_rng == 0; rho_ge_t(:,1,j) = rho_0_lio2;
else; rho_ge_t(:,:,j) = OD_wrapper(t1_rng,de_fn2,rho_0_lio2,[],'ode45',tol);
end    

end
end
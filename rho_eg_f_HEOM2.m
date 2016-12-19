function [rho_eg_ft,diag_test]=rho_eg_f_HEOM2(N,om_1_rng,sup_op_eg,HL,tol)
%  Simplified version, works in exciton basis and assumes no explicit 
%   vibrations are included 
%  om_scale = typical frequency of exciton transfer
%  L_eg is -1i[H,rho_{eg}], plus extra acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously.  

if ~exist('tol','var')
    tol = 1e-8; %tolerance
end 
   
    rho_eg_ft = zeros(N*HL,length(om_1_rng),N);
 
for j = 1:N %loop over sites

rho_0_lio = zeros( N*HL,1); 
rho_0_lio(j)=1;    %initially in equilibirum
    
       for lp = 1:length(om_1_rng)
        
 sup_op_ge_freq = -(1i*speye(N*HL)*om_1_rng(lp) + sup_op_eg);
 if nargout==2
 tic
 end
 rho_eg_ft(:,lp,j) = sup_op_ge_freq\rho_0_lio; %solves LA system

 %this might be better done with iterative methods
 
 if nargout==2%this is currently just a test to see the comparitive rate 
     %at which the iternative methods solve the system.  If it works better
     %I will change the code to use this
 if lp>1 
     
      t_meth1 = t_meth1 + toc; tic        
     x0 = rho_eg_ft(:,lp-1,j); %use last calculated as initial guess
    [x,flag] = bicg(sup_op_ge_freq,rho_0_lio,tol,100,[],[],x0); %interative solver
     t_meth2 =  t_meth2 + toc;
     if flag ~=0 % checks if an error occurs
         [lp,flag]
     end
     meth_diff(lp,j)= norm(x - rho_eg_ft(:,lp-1,j));
    
 else
    t_meth2 = toc;  t_meth1 = toc;  %time taken with each method
    meth_diff = zeros(length(om_1_rng),N); %difference in accuracy
 end
 end

       end

end
 if nargout==2
diag_test = { t_meth1, t_meth2,meth_diff};
 end
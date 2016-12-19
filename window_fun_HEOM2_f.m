function [V_ge,V_ef]=window_fun_HEOM2_f(N,HL,om_3_rng,sup_op_full,tol)
%Changes, this version 
 % This version works in frequency space, i.e. using 
 % int_0^inf dt_3 exp(-i \omega_3 t_3) e^(L t_3) = 1/(i \omega_3 - L)
 % simplified version that works in the exciton basis with no vibrations              

 if ~exist('tol','var')
    tol = 1e-6; %tolerance
end 
 
    V_ge = zeros(length(om_3_rng),N*HL,N);
    V_ef = zeros(length(om_3_rng),N^2*(N-1)/2*HL,N,N*(N-1)/2);

    tmp11 = zeros(1+N+N*(N-1)/2); tmp22 = tmp11; 
    tmp1  = tmp11; tmp2 = tmp22;
    
 tmp11(1,2:N+1) = 1; tmp11 = reshape(tmp11.',[numel(tmp11),1]).';
 tmp11 = logical(repmat(tmp11,[1,HL]));  %logical op containing all 
 %positions in the matrix which can be mixed into 
 
tmp22(2:N+1,N+2:end) = 1; tmp22 = reshape(tmp22.',[numel(tmp22),1]).';
tmp22 = logical(repmat(tmp22,[1,HL])); 

sup_op_ge = sup_op_full(tmp11,tmp11);
sup_op_ef = sup_op_full(tmp22,tmp22);

% Trace(A*B) = reshape(A.',N^2,1).' * reshape(B,N^2,1)
% Here "A" is V(t), either ge or ef, the upper diagonal
for j = 1:N
    tmp = tmp1; tmp(1,1+j) = 1; tmp = reshape(tmp.',[numel(tmp),1]).';
    tmpp = [tmp,zeros(1,(HL-1)*numel(tmp))];
    V_ge_0_lio{j} = tmpp(tmp11);
    for f =1:N*(N-1)/2
        tmp = tmp2; tmp(1+j,N+1+f) = 1; 
        tmp = reshape(tmp.',[numel(tmp),1]).';
        tmpp = [tmp,zeros(1,(HL-1)*numel(tmp))];
        V_ef_0_lio{j,f} = tmpp(tmp22);        
    end
end
    


for lp = 1:length(om_3_rng)
 
     sup_op_ge_freq = -(1i*speye(N*HL)*om_3_rng(lp) + sup_op_ge).'; 
     sup_op_ef_freq = -(1i*speye(N^2*(N-1)/2*HL)*om_3_rng(lp) + sup_op_ef).';
    
    for j = 1:N
%tic
        %for some reason B\A requires column vector B, despite the fact it
        %would have to be a row vector (or rows matching that or A cols)
  if lp ==1
    V_ge(lp,:,j) = (sup_op_ge_freq \(V_ge_0_lio{j}.')).';
  else
     % [L,U] = ilu(sup_op_ge_freq,struct('type','ilutp','droptol',1e-6));
      L = []; U=[]; %preconditioning factors
      x0 = V_ge(lp-1,:,j).'; %use last calculated as initial guess
    [x,flag] = bicg(sup_op_ge_freq,V_ge_0_lio{j}.',tol,100,L,U,x0); %interative solver
     if flag ~=0 % checks if an error occurs
         [lp,flag]
         V_ge(lp,:,j) = (sup_op_ge_freq \(V_ge_0_lio{j}.')).';
     else
        V_ge(lp,:,j) = x;
     end
  end
  %toc
        for f = 1:N*(N-1)/2
          %  tic
            if lp ==1
        V_ef(lp,:,j,f) = (sup_op_ef_freq \(V_ef_0_lio{j,f}.')).';
            else
      %  [L,U] = ilu(sup_op_ef_freq,struct('type','ilutp','droptol',1e-6)); 
        L = []; U=[]; %preconditioning factors
          x0 = V_ef(lp-1,:,j,f).'; %use last calculated as initial guess
    [x,flag] = bicg(sup_op_ef_freq,V_ef_0_lio{j,f}.',tol,100,L,U,x0); %interative solver
     if flag ~=0 % checks if an error occurs
         [lp,flag]
         V_ef(lp,:,j,f) = (sup_op_ef_freq \(V_ef_0_lio{j,f}.')).';
     else
        V_ef(lp,:,j,f) = x;
     end        
            end
           % toc
        end
    end
end
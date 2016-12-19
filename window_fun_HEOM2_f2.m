function [V_ge,V_ef,V_ge_0_lio,V_ef_0_lio]=window_fun_HEOM2_f2(N,HL,om_3_rng,sup_op_full)
%Changes, this version uses the / lin solver op
 % This version works in frequency space, i.e. using 
 % int_0^inf dt_3 exp(-i \omega_3 t_3) e^(L t_3) = 1/(i \omega_3 - L)
 % this version that works in the exciton basis with no (explicit) vibrations  
 % << V(om_3)| = << V | 1/(i \omega_3 - L) ->  << V | = << V(om_3)| (i \omega_3 - L)
 % this is a linear algebra system which can be solved with the usual tools
    V_ge = zeros(length(om_3_rng),N*HL,N);
    V_ef = zeros(length(om_3_rng),N^2*(N-1)/2*HL,N,N*(N-1)/2);

    tmp11 = zeros(1+N+N*(N-1)/2); tmp22 = tmp11; 
    tmp1  = tmp11; tmp2 = tmp22;
    
 tmp11(1,2:N+1) = 1; tmp11 = reshape(tmp11.',[numel(tmp11),1]).';
 % Trace(A*B) = reshape(A.',N^2,1).' * reshape(B,N^2,1)
 tmp11 = logical(repmat(tmp11,[1,HL]));  %equal to tmp1
 
tmp22(2:N+1,N+2:end) = 1; tmp22 = reshape(tmp22.',[numel(tmp22),1]).';
tmp22 = logical(repmat(tmp22,[1,HL])); 

sup_op_ge = sup_op_full(tmp11,tmp11);
sup_op_ef = sup_op_full(tmp22,tmp22);
clear sup_op_full %save memory

% Trace(A*B) = reshape(A.',N^2,1).' * reshape(B,N^2,1)
% Here "A" is either V_ge(t) or V_ef(t), on the upper diagonal
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
     %take transpose due to the way the function required to solve works
     sup_op_ef_freq = -(1i*speye(N^2*(N-1)/2*HL)*om_3_rng(lp) + sup_op_ef).';
   % sup_op_ge_freq = -(1i*eye(N*HL)*om_3_rng(lp) + sup_op_ge)^(-1); 
   % sup_op_ef_freq = -(1i*eye(N^2*(N-1)/2*HL)*om_3_rng(lp) + sup_op_ef)^(-1);
    
    for j = 1:N

        %for some reason B\A requires column vector B, despite the fact it
        %would have to be a row vector (or rows matching that or A cols)
        %x = B*A^(-1) -> x*A = B -> A'*x' = B' -> x = (A.' \ B.').'
    V_ge(lp,:,j) = (sup_op_ge_freq \(V_ge_0_lio{j}.')).';
  %  V_ge(lp,:,j) = V_ge_0_lio{j}*sup_op_ge_freq;
    
        for f = 1:N*(N-1)/2
                  
        %V_ef(lp,:,j,f) = V_ef_0_lio{j,f}.'\sup_op_ef_freq;
        V_ef(lp,:,j,f) = (sup_op_ef_freq \(V_ef_0_lio{j,f}.')).';
    %    V_ef(lp,:,j,f) = V_ef_0_lio{j,f}*sup_op_ef_freq;
        
        end
    end
end


    %{
% for j = 1:N
%     tmp = tmp1; tmp(1,1+j) = 1;  tmp = kron(eye(size(tmp)),tmp);
%     % reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
%     tmp = sparse(kron(eye(HL),tmp));
%     V_ge_0_lio{j} = tmp(tmp11,tmp11);
%     for f =1:N*(N-1)/2
%         tmp = tmp2; tmp(1+j,N+1+f) = 1;  tmp = kron(eye(size(tmp)),tmp);
%     % reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
%     tmp = sparse(kron(eye(HL),tmp));
%         V_ef_0_lio{j,f} = tmp(tmp22,tmp22);      
%     end
% end    
% 
% for lp = 1:length(om_3_rng)
%  
%    sup_op_ge_freq = -(1i*eye(N*HL)*om_3_rng(lp) + sup_op_ge)^(-1); 
%    sup_op_ef_freq = -(1i*eye(N^2*(N-1)/2*HL)*om_3_rng(lp) + sup_op_ef)^(-1);
%     
%     for j = 1:N
% 
%     V_ge(lp,:,j) = V_ge_0_lio{j}*sup_op_ge_freq;
%     
%         for f = 1:N*(N-1)/2
%                   
%         V_ef(lp,:,j,f) = V_ef_0_lio{j,f}*sup_op_ef_freq;
%         
%         end
%     end
% end
%}
end
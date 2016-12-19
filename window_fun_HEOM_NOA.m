function [V_ge_f,V_ef_f]=window_fun_HEOM_NOA...
            (om_3_rng,V_ge,V_ef,sup_op_ge,sup_op_ef)
%Changes, this version uses the / lin solver op
 % This version works in frequency space, i.e. using 
 % int_0^inf dt_3 exp(-i \omega_3 t_3) e^(L t_3) = 1/(i \omega_3 - L)
 % this version that works in the exciton basis with no (explicit) vibrations  
 % << V(om_3)| = << V | 1/(i \omega_3 - L) ->  << V | = << V(om_3)| (i \omega_3 - L)
 % this is a linear algebra system which can be solved with the usual tools
    V_ge_f = zeros(length(om_3_rng),length(V_ge));
    V_ef_f = zeros(length(om_3_rng),length(V_ef));
%{
    tmp11 = zeros(1+N+N*(N-1)/2); tmp22 = tmp11; 
    tmp1  = tmp11; tmp2 = tmp22;
    
 tmp11(1,2:N+1) = 1; tmp11 = reshape(tmp11.',[numel(tmp11),1]).';
 % Trace(A*B) = reshape(A.',N^2,1).' * reshape(B,N^2,1)
 tmp11 = logical(repmat(tmp11,[1,HL]));  %equal to tmp1
 
tmp22(2:N+1,N+2:end) = 1; tmp22 = reshape(tmp22.',[numel(tmp22),1]).';
tmp22 = logical(repmat(tmp22,[1,HL])); 

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
%}    

ID1 = 1i*speye(length(sup_op_ge));   ID2 = 1i*speye(length(sup_op_ef));

for lp = 1:length(om_3_rng)
 
     sup_op_ge_freq = -(ID1*om_3_rng(lp) + sup_op_ge).';  
     %take transpose due to the way the function required to solve works
     sup_op_ef_freq = -(ID2*om_3_rng(lp) + sup_op_ef).';

        %for some reason B\A requires column vector B, despite the fact it
        %would have to be a row vector (or rows matching that or A cols)
        %x = B*A^(-1) -> x*A = B -> A'*x' = B' -> x = (A.' \ B.').'

    V_ge_f(lp,:) = (sup_op_ge_freq \(V_ge.')).';
    V_ef_f(lp,:) = (sup_op_ef_freq \(V_ef.')).';

end
end
function [GSB,SE,ESA,tmp_gg_t,tmp_ee_t] = pp_cont_HEOM(sz_params,...
            rho_eg,V_ge_t,V_ef_t,prop_op_gg,prop_op_ee,t2_range,e1,e2,...
            av1_g,av1_f,av2_g,av2_f,H_trunc,tol) ;

  warning('incomplete')
%Calculated the re/nonrephasing contribution to a 2D signal given the coherence
%part of the density matrix, the Heisenberg picture operators V(t_3) and 
% the operators that propogate the system in the ground and excited state
% This will sum over all pathways with the dipole
% averaging constants in av1, av2 to give a result.  It is assumed that the
% av2 are of the averaging form in my notes with the final index running
% from 1-3 each of which is weighted by |k1| |k2| |k3| etc and the factor
% of +/- 1 depending on which signal we are considering

%sz_params is a vector with the important length constants, note that code
%works (should) even if HL = 1 (not the HEOM, just a normal propogator)
N = sz_params(1); %number of sites
HL = sz_params(2); %total number of aux density matricies (+1)
t1L = size(rho_eg,2);  t2L = length(t2_range);  t3L = size(V_ge_t,2);
s_rp = [-1,+1,+1]; s_nr = [+1,-1,+1];

pause on %nice to be able to stop the code in windows

if ~exist('tol','var')
    tol =[1e-7,1e-5];
end

out_trunc_gg = H_trunc;
out_trunc_ee = N^2*H_trunc;

    num_sets = [size(av1_g,5),size(av2_g,5)];
    %number of sets of values to calculate for the 4th order type averages
    %and 5th order type averages
    
 %preallocate rephasing

    GSB = complex(zeros(t3L,t2L,t1L,sum(num_sets)));  
    ESA = complex(zeros(t3L,t2L,t1L,sum(num_sets)));  
    SE = complex(zeros(t3L,t2L,t1L,sum(num_sets)));  

rho_eg = reshape(rho_eg,N,HL,t1L,N); %reshape to this form for convinience

 %cj_gg = HEOM_conj({HL,1,1});  %trivial 
 cj_ee = HEOM_conj({HL,N*1,N*1});
  cj_eg = HEOM_conj({HL,1,N});
  cj_fe = HEOM_conj({HL,N,N*(N-1)/2});
  
rho_gg_prop = @(t,v) prop_op_gg*v;
rho_ee_prop = @(t,v) prop_op_ee*v;

for j1 = 1:N  %first transition to exciton j and then either conjugate or j->f
    for j2 = 1:N %2nd dipole interaction
        
        %propogate the resulting excited state pop and ground state hole
        
    tmp_gg = rho_eg(j2,:,:,j1);   %picks these elements
    tmp_gg = reshape(tmp_gg,[HL,t1L]); %reshape to normal Liouville space
    
    %reshape rho_eg into seperate Heirarchy components 
    tmp_ee =  repmat(rho_eg(:,:,:,j1),[N,1,1]);
    lg = true(N^2,1); lg(1+(j2-1)*N:j2*N) = false; %these elements remain
    tmp_ee(lg,:,:) = 0; tmp_ee = reshape(tmp_ee,[N^2*HL,t1L]);

    tic
    %loop over t1L, yes this WILL be painfully slow if this is long
   % figure 
   % hold on
    for lp = 1:t1L %could avoid nexting other loops within this
        %but it might mean quite a lot of saved data.
        pause(1e-3) %take a millisecond pause to allow for interupting the code
        tmp_gg_sec = tmp_gg(:,lp);
        tmp_gg_t = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec,out_trunc_gg,'ode45',tol)).';      
        tmp_gg_t = tmp_gg_t + conj(tmp_gg_t); %1 by 1 matricies
        
        tmp_ee_sec = tmp_ee(:,lp);
        tmp_ee_t = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec,out_trunc_ee,'ode45',tol)).';
        tmp_ee_t = conj(tmp_ee_t(cj_ee ,:));    
       % plot(t2_range,tmp_ee_t([1,4],:))
        %now calculate the window function
       for j4=1:N
           
           V_ge_t
           
           
       end
        

    end
    toc
    end
end
end



function rho_conj = HEOM_conj(rho_vec)
%conjugate transposes each block in the HEOM 

persistent conj_map

if iscell(rho_vec) %generate map which transposes
   
    HL_tmp = rho_vec{1}; N1_tmp = rho_vec{2}; N2_tmp = rho_vec{3};
    %create the logical map which pics the elements which will be
    %conjugated
    
    temp = reshape(1:N1_tmp*N2_tmp,N1_tmp,N2_tmp); temp2 = temp';
    temp3 = reshape(temp2,N1_tmp*N2_tmp,1); %new mappings
    
    conj_map = zeros(N1_tmp*N2_tmp*HL_tmp,1);
    for jj = 1:HL_tmp
        
        conj_map((jj-1)*N1_tmp*N2_tmp+1:jj*N1_tmp*N2_tmp) = temp3 + (jj-1)*N1_tmp*N2_tmp;
        
    end
    rho_conj = conj_map; return
end
    
    rho_conj = conj(rho_vec(conj_map));


end

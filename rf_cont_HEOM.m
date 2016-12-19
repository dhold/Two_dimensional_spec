function [R1,R2,R3,R4,R5,R6] = rf_cont_HEOM(sz_params,rho_eg,V_ex_eg,V_ex_fe,...
                            V_ge_t,V_ef_t,prop_op_gg,prop_op_ee,t2_range...
                            ,calc_rf,calc_nr,dip_av_set_rf,dip_av_set_nrf) ;

%Calculated the re/nonrephasing contribution to a 2D signal given the coherence
%part of the density matrix, the Heisenberg picture operators V(t_3) and 
% the operators that propogate the system in the ground and excited state
% This will sum over all pathways with the dipole
% averaging constants in dip_av_set_rf/dip_av_set_nrf to give a result

%sz_params is a vector with the important length constants, note that code
%works (should) even if HL = 1 (not the HEOM, just a normal propogator)
N = sz_params(1); %number of sites
HL = sz_params(2); %total number of aux density matricies (+1)
tot_vib = sz_params(3); %size of any explicit vibrational Hamiltonian
t1L = size(rho_eg,2);  t2L = length(t2_range);  t3L = size(V_ge_t,2);

%dipole_av_set = <mu_4.E_4 , mu_3.E_3 , mu_2.E_2 , mu_1.E_1 > etc, i.e. 
%prefactor from orientation and heterodyne setup.  This is a length 
% [N, N, N, 1+N*(N-1)/2,N, 1+N*(N-1)/2,num_sets] length vector (hence in
% principle it is fairly large!)

if ~calc_rf && ~calc_nrf
   warning('nothing calculated')
   R1=[]; R2=[]; R3=[]; R4=[]; R5=[]; R6=[]; return
end
    num_sets = [size(dip_av_set_rf,7),size(dip_av_set_nrf,7)];
    
    
    if calc_rf
 %preallocate rephasing
        R1 = zeros(t3L,t2L,t1L,num_sets(1)); R4 = R1; R5 = R4;
    else
        R1=[]; R4=[]; R5=[];
    end
    if calc_nrf
    R2 = zeros(t3L,t2L,t1L,num_sets(1));  %preallocate nonrephasing
        R3 = R2; R6 = R3;
    else
        R2 = []; R3=[]; R6 = [];
    end    



%reshape V_ex_eg and V_ex_ef to the appropriate sizes to act on the HEOM
%and put them in Lioville space
for j = 1:N
   %note this is the opposite direction to the over arrows, it is the side
   %the operator acts from
    V_eg_left{j} = repmat(sparse(  kron(V_ex_eg(:,:,j),eye(size(V_ex_eg(:,:,j))))),[HL,1]); %#ok<*AGROW>
    V_eg_right{j} = repmat(sparse(  kron(eye(size(V_ex_eg(:,:,j))),V_ex_eg(:,:,j))),[HL,1]);
    V_ge_left{j} = repmat(sparse(  kron(V_ex_eg(:,:,j)',eye(size(V_ex_eg(:,:,j)')))),[HL,1]);
    V_ge_right{j} = repmat(sparse(  kron(eye(size(V_ex_eg(:,:,j)')),V_ex_eg(:,:,j)')),[HL,1]);
    for f=1:N*(N-1)/2
    V_fe_left{j,f} = repmat(sparse(  kron(V_ex_fe(:,:,j,f),eye(size(V_ex_fe(:,:,j,f))))),[HL,1]);    
    end
end

%generate maps which can pick out the hermitian conjugate of a Liouville
%space vector

cj_ee = HEOM_conj({HL,N*tot_vib,N*tot_vib});
cj_gg = HEOM_conj({HL,tot_vib,tot_vib});


rho_gg_prop = @(t,v) prop_op_gg*v;
rho_ee_prop = @(t,v) prop_op_ee*v;

for j = 1:N  %first transition to exciton j and then either conjugate or j->f
    for j2 = 1:N %2nd dipole interaction
        
        %propogate the resulting excited state pop and ground state hole
    tmp_gg = V_ge_left{j2}*rho_eg(:,:,j1);
    tmp_ee = V_ge_right{j2}*rho_eg(:,:,j1);
    tic
    %loop over t1L, yes this WILL be painfully slow if this is long
    for lp = 1:t1L
        tmp_gg_sec = tmp_gg(:,lp);
        tmp_gg_t = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec,out_trunc_gg,'ode45',tol)).';
        
        tmp_ee_sec = tmp_ee(:,lp);
        tmp_ee_t = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec,out_trunc_ee,'ode45',tol)).';
        
        %now calculate each of the contributions
        for j3 = 1:N
            if calc_rf
            R1tmp = V_eg_right{j3}*tmp_ee_t;
            R4tmp = -V_eg_left{j3}*tmp_gg_t; %-as abs
            end
             if calc_nrf
            R2tmp = V_eg_right{j3}*conj(tmp_ee_t(cj_ee,:));
            R3tmp = -V_eg_left{j3}*conj(tmp_gg_t(cj_gg,:)); %-as abs
             end
            
            for j4 = 1:N           
            if calc_rf                  
            ori_fct_rf = reshape(dip_av_set_rf(j1,j2,j3,1,j4,1,:),[1,1,num_sets(1)]);
             %combine with the contributions, flip round t2 and t3    
            %this operation is essentially a trace
            R_1 = mtimesx(mtimesx( V_ge_t(:,:,j4),R1tmp),ori_fct_rf); %t2_range by t3_range
            R_4 = mtimesx(mtimesx( V_ge_t(:,:,j4),R4tmp),ori_fct_rf);
                
           
            end
            if calc_nrf
            ori_fct_nrf = reshape(dip_av_set_nrf(j1,j2,j3,1,j4,1,:),[1,1,num_sets(2)]);
            R_2 =  mtimesx(mtimesx( V_ge_t(:,:,j4),R2tmp),ori_fct_nrf); %t2_range by t3_range
            R_3 =  mtimesx(mtimesx( V_ge_t(:,:,j4),R3tmp),ori_fct_nrf);

            end       
            if calc_rf
                R_6 = zeros(size(R_1));
            end
            if calc_nrf
                R_5 = zeros(size(R_2));
            end
                        
                for f = 1:N*(N-1)/2
                    if calc_rf
                   R6tmp = -(V_fe_left{j3,f}*tmp_ee_t); 
                ori_fct_rf = reshape(dip_av_set_rf(j1,j2,j3,f+1,j4,f+1,:),[1,1,num_sets(1)]);
                R_6 = R_6+mtimesx(mtimesx( V_ef_t(:,:,j4,f),R6tmp),ori_fct_rf); 
              %combine with the contributions 
                    end
                    if calc_nrf
                    R5tmp = -V_fe_left{j3,f}*conj(tmp_ee_t(cj_ee,:));    
                ori_fct_nrf = reshape(dip_av_set_rf(j1,j2,j3,f+1,j4,f+1,:),[1,1,num_sets(1)]);        
                   R_5 = R_5+mtimesx(mtimesx( V_ef_t(:,:,j4,f),R5tmp),ori_fct_nrf);       
                    end                    
                end
   %combine this frequency slices into the total response function             
               if calc_rf              
      R2(:,:,:,lp) = R_2; R3(:,:,:,lp) = R_3; R5(:,:,:,lp) = R_5;
               end
               if calc_nr
      R1(:,:,:,lp) = R_1; R4(:,:,:,lp) = R_4; R6(:,:,:,lp) = R_6;
               end
            end
        end
    end
    toc
    end
end
end

function rho_conj = HEOM_conj(rho_vec)
%conjugate transposes each block in the HEOM 

persistent conj_map

if iscell(rho_vec) %generate map which conjugates
   
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

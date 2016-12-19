function [R1_sdav,R2_sdav,R3_sdav,R4_sdav,R5_sdav,R6_sdav,lin_spec,rho_gg_t,rho_ee_t]=...
        twoD_with_HEOM_NOA(sdlp_rng,HEOM_params,Ham_params,comp_lg,...
        V_eg,V_fe,om1_rng,t2_range,om3_rng,t_verbose);
     
%Version 2 uses an alternative method to calculate the components of the
%response function in frequency space
%this should calculate the 2D spectra in the impulsive, RWA limit for the 
%Rephasing, non rephasing and coherent pathways.  
% S_k1(om_3,t_2,om_1) , S_k2(om_3,t_2,om_1) and S_k3(om3,om2,t1)
% Passing a cell om_1 takes om_1{1} to be a freq range and om_1{2} as a
% time range.
%Freq_diff and beam width are parameters which select only components that 
%beat at a particular frequency during t_2 rather than all of them 
%%

[coup_com_save,coup_acom_save,const_factor,QQ_topass,nn,max_tier_final]...
    = HEOM_params{:}; 

[fock_space_rep,H_site,site_shift,pdm_shift] = Ham_params{:};

[calc_rp,calc_nr,calc_coh] = comp_lg{:};

clear comp_lg Ham_params HEOM_params %clear after dealing

N = size(fock_space_rep,2);
fock1 = fock_space_rep(2:N+1,:);
fock2 = fock_space_rep(N+2:end,:);
HL = size(nn,1); %length of HEOM
tier_lvl = sum(nn,2); %sum of nn gives the tier
H_trunc = sum(tier_lvl<=max_tier_final); sdlp=1;

use_ex_basis = false;

if ~exist('t_verbose','var') %give value true to output time taken 
    t_verbose = false;
end

%tauL = length(t2_range); w3L = length(om3_rng); w1L = length(om1_rng);
%%
for sdlp = sdlp_rng

%%  Generate Hamiltonian things

H_e = H_site + diag(site_shift(sdlp,:));% add static disorder
H_f = zeros(N*(N-1)/2); %double excited Hamiltonian
occ_f = logical(fock_space_rep(sum(double(fock_space_rep),2)==2,:));

for k = 1:N*(N-1)/2 %enumerate the coupling Hamiltonian
    
        lg1 = occ_f(k,:); s1 = find(lg1);
       H_f(k,k)= H_e(s1(1),s1(1))+H_e(s1(2),s1(2)); 
       if ~isempty(pdm_shift) %shift to double excited states
           H_f(k,k) = H_f(k,k) + pdm_shift{s1(1),s1(2)};
       end
       
       for kk = k+1:N*(N-1)/2
           lg2 = occ_f(kk,:); lg3 = lg1 & lg2; %repeated element
            if any(lg3)
                s2 = find(lg2 &~lg3); s3 = find(lg1 &~lg3);
                H_f(k,kk) = H_e(s2,s3); %#ok<FNDSB>
                H_f(kk,k) = H_f(k,kk)';
            end               
       end  
end
M_g =1; [M_e,E1] = eig(H_e);  [M_f,E2] = eig(H_f);  %diagonalise the blocks
E1 = diag(E1);  E2 = diag(E2); MM = blkdiag(1,M_e,M_f);
%projectors in site only basis for single and excited states

% lam_ex_full = zeros(1+N+N*(N-1)/2,1); lam_tot = zeros(N,1);
% for k = 1:N
%     tmp = diag(mtimesx(MM,'C',mtimesx(diag(fock_space_rep(:,k)),MM)));
%    lam_ex_full = lam_ex_full + tmp.*lam_tot(k);
% end

%% calculate operator coupling excitations on each site to the exciton
%hamiltonian

%Choose whether to use the exciton basis, don't pass the exciton dipole moments
%Due to static disorder these must be calculated within the program.  
if use_ex_basis
V_coup_tot = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,N);
for j = 1:N %project each term into the exciton basis
V_1 = M_e'*diag(fock1(:,j))*M_e;
V_2 = M_f'*diag(fock2(:,j))*M_f;
V_coup_tot(:,:,j) = blkdiag(0,V_1,V_2);
%V_coup_tot(:,:,j) = -blkdiag(0,V_1,V_2);
end                

prop_op_full = H_prop_gen2(diag([0;E1;E2]),diag([0;E1;E2]),...
                    V_coup_tot,V_coup_tot,QQ_topass,...
                    const_factor,coup_com_save,coup_acom_save); 
else

V_coup_tot = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,N);
for j = 1:N %project each term into the exciton basis
V_1 = diag(fock1(:,j));
V_2 = diag(fock2(:,j));
V_coup_tot(:,:,j) = blkdiag(0,V_1,V_2);

end      

prop_op_full = H_prop_gen2(blkdiag(0,H_e,H_f),blkdiag(0,H_e,H_f),...
                    V_coup_tot,V_coup_tot,QQ_topass,...
                    const_factor,coup_com_save,coup_acom_save); 

end
tmp1 = zeros(size(V_coup_tot(:,:,1))); tmp2 = tmp1; tmp3 = tmp1;
%tmp11 = tmp1; tmp22 = tmp2;

tmp1(2:N+1,1) = 1; tmp1 = reshape(tmp1,[numel(tmp1),1]);
tmp1 = logical(repmat(tmp1,[HL,1])); 
tmp1b = tmp2; tmp1b(1,2:N+1) = 1; tmp1b = reshape(tmp1b,[numel(tmp1b),1]);
tmp1b = logical(repmat(tmp1b,[HL,1])); 
            
tmp2(2:N+1,2:N+1) = 1; tmp2 = reshape(tmp2,[numel(tmp2),1]);
tmp2 = logical(repmat(tmp2,[HL,1])); 

tmp3(1,1) = 1; tmp3 = reshape(tmp3,[numel(tmp3),1]);
tmp3 = logical(repmat(tmp3,[HL,1])); 

%  tmp11(1,2:N+1) = 1; tmp11 = reshape(tmp11.',[numel(tmp11),1]);
%  tmp11 = logical(repmat(tmp11,[HL,1]));  %equal to tmp1

prop_op_eg = prop_op_full(tmp1,tmp1);
prop_op_ge = prop_op_full(tmp1b,tmp1b);
prop_op_ee = prop_op_full(tmp2,tmp2);
prop_op_gg = prop_op_full(tmp3,tmp3);

%do the same for the left acting components
tmp22 = zeros(1+N+N*(N-1)/2); 

% tmp11 = tmp22; 
%  tmp11(1,2:N+1) = 1; tmp11 = reshape(tmp11.',[numel(tmp11),1]).';
%  tmp11 = logical(repmat(tmp11,[1,HL]));  %equal to tmp1

tmp22(2:N+1,N+2:end) = 1; tmp22 = reshape(tmp22.',[numel(tmp22),1]).';
tmp22 = logical(repmat(tmp22,[1,HL])); 

%prop_op_ge2 = prop_op_full(tmp11,tmp11); %equal to  prop_op_eg
prop_op_ef2 = prop_op_full(tmp22,tmp22);

clear prop_op_full

%% first calculate the initial coherences set up

%generate the N initial density matricies from the first pathways


 if iscell(om1_rng);  t1_rng = om1_rng{2}; om1_rng = om1_rng{1};
 else;  t1_rng=[];  end
 
 tic
V_ge1 = V_eg{1}'; V_eg1 = V_eg{1}; NN = numel(V_ge1);
if use_ex_basis; V_ge1 = M_g'*V_ge1*M_e; V_eg1 = M_e'*V_eg1*M_g; end
%initial conditions to propogate
rho_eg = [reshape(V_eg1,NN,1); zeros((HL-1)*NN,1)]; 
rho_ge = [reshape(V_ge1,NN,1); zeros((HL-1)*NN,1)]; 

 [rho_eg_ft,rho_ge_ft]=rho_eg_f_HEOM_NOA(om1_rng,rho_eg,prop_op_eg,...
                        rho_ge,prop_op_ge);  %rho_eg(om_1) and rho_ge(-om_1)
                
init_coherence_time_to_calc = toc;
if t_verbose;  init_coherence_time_to_calc   
end

%use this coherence value to calculate linear spectra 

%sigma_om1 = [reshape(V_eg{2}',1,NN), zeros(1,(HL-1)*NN)]*rho_ge_ft; 
if ~use_ex_basis
sigma_om1 = reshape(V_eg{2}',1,NN)*rho_ge_ft(1:NN,:); 
else
sigma_om1 = reshape(M_g'*(V_eg{2}')*M_e,1,NN)*rho_ge_ft(1:NN,:);     
end
%test = (rho_eg_ft(1:NN,:)')*V_eg{2};  %should give the same

lin_spec{1}={sigma_om1};

if ~isempty(t1_rng)

[rho_eg_t,rho_ge_t]  =  rho_eg_t_HEOM_NOA(t1_rng,...
                                rho_eg,prop_op_eg, rho_ge,prop_op_ge);
%combine all of them together, it isn't actually important for what follows
%whether this is calculated in time or frequency space
 rho_eg_ft = cat(2,rho_eg_ft,rho_eg_t);         
 rho_ge_ft = cat(2,rho_ge_ft,rho_ge_t);  
end

%% Calculate the amplitude of the window functions

tic 
 
V_ge4 = V_eg{4}';  V_ef4 = V_fe{4}';
if use_ex_basis; V_ge4 = M_g'*V_ge4*M_e; V_ef4= M_e'*V_ef4*M_f; end

NN2 = numel(V_ef4);
V_ge_lio = [reshape(V_ge4,1,NN),zeros(1,NN*(HL-1))]; 
V_ef_lio = [reshape(V_ef4,1,NN2),zeros(1,NN2*(HL-1))]; 
[V_ge_f,V_ef_f]=window_fun_HEOM_NOA(om3_rng,V_ge_lio,V_ef_lio,prop_op_eg,prop_op_ef2);

time_dep_of_output_op_calc_time = toc;
if t_verbose; time_dep_of_output_op_calc_time
end
%This is expressed in Liouville space                                

%also use this time dependent value to calculate linear spectra 
% dt3 =  t3_range(2)-t3_range(1); om3 = pi*(-1/dt3:2/t3_range(end):1/dt3);

%note this is different if the polarization of the 4th beam is different to
%that of the second
sigma_om3 = V_ge_f*[reshape(V_eg{1},NN,1); zeros((HL-1)*NN,1)] ; 

lin_spec{2}={sigma_om3};


%% Now calculate each of the required neq second order matricies

%this method works by calculating slices of the density matrix over t1 (or
%omega1 in frequency space) 
if calc_rp || calc_nr
    
sz_params = [N,HL]; 

V_ge2 = V_eg{2}';  V_ef2 = V_fe{2}'; V_eg3 = V_eg{3};  V_fe3 = V_fe{3};
if use_ex_basis; V_ge2 = M_g'*V_ge2*M_e; V_ef2= M_e'*V_ef2*M_f; 
V_eg3 = M_e'*V_eg3*M_g; V_fe3= M_f'*V_fe3*M_e;  end

V_set = {V_ge2, V_ef2,V_eg3, V_fe3};
    if isvector(sdlp_rng)
[R1,R2,R3,R4,R5,R6,rho_gg_t,rho_ee_t] = rp_cont_HEOM_NOA(sz_params,...
            rho_eg_ft,rho_ge_ft,V_set,V_ge_f,V_ef_f,prop_op_gg,prop_op_ee,...
            t2_range,calc_nr,calc_rp,H_trunc,[],t_verbose) ;
    else %only one value, just output averages directly
[R1_sdav,R2_sdav,R3_sdav,R4_sdav,R5_sdav,R6_sdav,rho_gg_t,rho_ee_t] = ...
            rp_cont_HEOM_NOA(sz_params,rho_eg_ft,rho_ge_ft,...
            V_set,V_ge_f,V_ef_f,prop_op_gg,prop_op_ee,...
            t2_range,calc_nr,calc_rp,H_trunc,[],t_verbose) ;        
        return %no need to do any more
    end
else %nothing calculated
    R1 = []; R2 =[]; R3 = []; R4 = []; R5 = []; R6 = [];
end


%% Add together to create average

if sdlp == 1
    
    R1_sdav = R1;  R4_sdav = R4;   R6_sdav = R6;   
    R2_sdav = R2;  R3_sdav = R3;   R5_sdav = R5; 
else %add extra terms
    R1_sdav = R1_sdav + R1;  R4_sdav = R4_sdav + R4;   R6_sdav = R6_sdav + R6;   
    R2_sdav = R2_sdav + R2;  R3_sdav = R3_sdav + R3;   R5_sdav = R5_sdav + R6; 
end
end
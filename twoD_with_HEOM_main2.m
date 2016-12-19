function [R1_sdav,R2_sdav,R3_sdav,R4_sdav,R5_sdav,R6_sdav,lin_spec,tmp_gg_t,tmp_ee_t]=...
        twoD_with_HEOM_main2(sdlp_rng,coup_com_save,coup_acom_save,...
        const_factor,QQ_topass,nn,max_tier_final,fock_space_rep,...
         H_site,site_shift,pdm_shift,V_n,av_2,av_3,av2_3,av_4,av_5,av2_5,...
        calc_rp,calc_nr,om1_rng,t2_range,om3_rng,pop_only,t_verbose,freq_diff,beam_width);
     
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
N = size(fock_space_rep,2);
fock1 = fock_space_rep(2:N+1,:);
fock2 = fock_space_rep(N+2:end,:);
HL = size(nn,1); %length of HEOM
tier_lvl = sum(nn,2);
H_trunc = sum(tier_lvl<=max_tier_final); sdlp=1;
if ~exist('pop_only','var') %only takes into account |g><g| -> |e><e| transitions
    t_verbose = false;
end
if ~exist('t_verbose','var') %give value true to output time taken 
    t_verbose = false;
end
if ~exist('freq_diff','var') || ~exist('beam_width','var') 
    beam_width = []; freq_diff = []; %select all frequencies
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
[M_e,E1] = eig(H_e);  [M_f,E2] = eig(H_f);  %diagonalise the blocks
E1 = diag(E1);  E2 = diag(E2); MM = blkdiag(1,M_e,M_f);
%projectors in site only basis for single and excited states
alpha_n = mtimesx(mtimesx(MM,'C',V_n),MM); 
Cg  = squeeze(alpha_n(1,2:N+1,:)) ; %contract over dimension 2
Cf = alpha_n(2:N+1,N+2:end,:); %contract over dimension 3
%coefficients mapping the interaction operator into the new basis

lam_ex_full = zeros(1+N+N*(N-1)/2,1); lam_tot = zeros(N,1);
for k = 1:N
    tmp = diag(mtimesx(MM,'C',mtimesx(diag(fock_space_rep(:,k)),MM)));
   lam_ex_full = lam_ex_full + tmp.*lam_tot(k);
end

%% Calculate dipole moments averages in the exciton basis and such like
%1st order averages
av_set_fo = mtimesx(mtimesx(Cg,av_2),Cg,'C'); 
av_set_cd = mtimesx(mtimesx(Cg,av_3+av2_3),Cg,'C'); 

%3rd order averages

tmp = mtimesx(mtimesx(Cg,av_4),Cg,'C'); %first two transitions are between
%g and one exciton transitions
tmp = permute(tmp,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx

%the first one involves transitions between g-e only (GSB, SE)
av_set_GSB = mtimesx(mtimesx(Cg,squeeze(tmp)),Cg,'C');
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order
%the next loop is more complicated because it involves transitions from e-f
%and then f'-e', need to reshape
Cf_flat = reshape(Cf,N^2*(N-1)/2,N); %flatten out for use with mtimesx
%new order like Cf_flat(k,f,j) = Cf_flat(k+N*(f-1),j)

av_set_ESA = mtimesx(mtimesx(Cf_flat,squeeze(tmp)),Cf_flat,'C'); %ESA type
av_set_ESA = permute(av_set_ESA,[3,4,1,2,5]);
%av_set_ESA (e1,e2,e3*(N-1)*N/2+f1,e4*(N-1)*N/2+f2) type calls
%next calculate the higher order moments from the other interactions

tmp = mtimesx(mtimesx(Cg,av_5+av2_5),Cg,'C');
tmp = permute(tmp,[3,4,1,2,5,6]); 

%approximate amplitude of k1, k2 by the vertical transition energy of the
%exciton states this map into k~w_res/c
%next I want to really scale SE and GSB differently to account for Stokes
%shift, this is a closer approximation Do this for each path when
%performing the loop
% E_GSB = E1+lam_ex_full(2:N+1);  
% E_SE = E1-lam_ex_full(2:N+1);  

av_set2_GSB  = mtimesx(mtimesx(Cg,tmp), Cg,'C');
av_set2_GSB = permute(av_set2_GSB,[3,4,1,2,5,6]);
av_set2_ESA  = mtimesx(mtimesx(Cf_flat,tmp), Cf_flat,'C');  
av_set2_ESA = permute(av_set2_ESA,[3,4,1,2,5,6]);
%{
%1st order averages
av_set_fo = mtimesx(mtimesx(M_e,'T',av_2),M_e); 
av_set_cd = mtimesx(mtimesx(M_e,'T',av_3+av2_3),M_e); 

%3rd order averages

tmp = mtimesx(mtimesx(M_e,'T',av_4),M_e); %first two transitions are between
%g and one exciton transitions

sz = size(tmp); 
if ndims(tmp) ==4
    sz = [sz,1];
end

tmp = permute(tmp,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx
av_set_GSB = mtimesx(mtimesx(M_e,'T',tmp),M_e);
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order

av_set_ESA =  zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2,size(av_set_GSB,5));

for k3=1:N %much bette loop
    for k4=1:N
        tmp2 = zeros(N*(N-1)/2,N*(N-1)/2,sz(1),sz(2),sz(5));
        for flp = 1:N*(N-1)/2
            ss = find(fock2(flp,:));
            for flp2 = 1:N*(N-1)/2
                ss2 = find(fock2(flp2,:));
                
tmp2 = tmp2 + mtimesx(M_f(flp,:).'*M_f(flp2,:),(...
        M_e(ss(1),k3)*M_e(ss2(1),k4)*tmp(ss(2),ss2(2),:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(1),k4)*tmp(ss(1),ss2(2),:,:,:) +...
        M_e(ss(1),k3)*M_e(ss2(2),k4)*tmp(ss(2),ss2(1),:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(2),k4)*tmp(ss(1),ss2(1),:,:,:)));

            end
        end
        av_set_ESA(:,:,k3,k4,:,:,:) = permute(tmp2,[3,4,1,2,5]);        
    end
end

tmp = mtimesx(mtimesx(M_e,'T',av_5+av2_5),M_e);
tmp = permute(tmp,[3,4,1,2,5,6]); 

av_set2_GSB = mtimesx(mtimesx(M_e,'T',tmp),M_e);
av_set2_GSB = permute(av_set2_GSB,[3,4,1,2,5,6]);

sz = size(av_set2_GSB);
av_set2_ESA =  zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2,sz(5),sz(6));

for k3=1:N 
    for k4=1:N
        tmp2 = zeros(N*(N-1)/2,N*(N-1)/2,sz(1),sz(2),sz(5),sz(6));
        for flp = 1:N*(N-1)/2
            ss = find(fock2(flp,:));
            for flp2 = 1:N*(N-1)/2
                ss2 = find(fock2(flp2,:));
                
tmp2 = tmp2 + mtimesx(M_f(flp,:).'*M_f(flp2,:),(...
        M_e(ss(1),k3)*M_e(ss2(1),k4)*tmp(ss(2),ss2(2),:,:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(1),k4)*tmp(ss(1),ss2(2),:,:,:,:) +...
        M_e(ss(1),k3)*M_e(ss2(2),k4)*tmp(ss(2),ss2(1),:,:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(2),k4)*tmp(ss(1),ss2(1),:,:,:,:)));

            end
        end
        av_set2_ESA(:,:,k3,k4,:,:,:) = permute(tmp2,[3,4,1,2,5,6]);        
    end
end
%}
% Next Calculate 2D signals

%% calculate operator coupling excitations on each site to the exciton
%hamiltonian
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
%% first calculate the initial coherences set up

%generate the N initial density matricies from the first pathways

%{
%generate prop op using the function H_prop_gen2
%this generates an operator which propogates the HEOM when only a certain
%part of the electronic part populated , here a ground excited state coherence,
%i.e. rho = [0_g,0_ge,0_gf;rho_eg,0_ee,0_ef;0_fg,0_fe,0_ff] 
% V_coup = zeros(N,N,N); %system bath coupling
% for j = 1:N %project each term into the exciton basis
% V_coup(:,:,j) = M_e'*diag(fock1(:,j))*M_e;
% end
% prop_op_eg = H_prop_gen2(diag(E1),0,V_coup,zeros(1,1,N),QQ_topass,...
%                     const_factor,coup_com_save,coup_acom_save);       
%}
 if iscell(om1_rng);  t1_rng = om1_rng{2}; om1_rng = om1_rng{1};
 else;  t1_rng=[];
 end
 
 tic
[rho_eg_ft]=rho_eg_f_HEOM2(N,om1_rng,prop_op_eg,HL); %rho_eg(om_1)
[rho_ge_ft]=rho_ge_f_HEOM2(N,om1_rng,prop_op_ge,HL); %rho_ge(-om_1)
init_coherence_time_to_calc = toc;
if t_verbose
    init_coherence_time_to_calc
end

%prop_op_freq = (1i*om1_rng*eye(size(prop_op_eg))-prop_op_eg)^(-1);

 %dt =  t1_range(2)-t1_range(1); om1 = pi*(-1/dt:2/t1_range(end):1/dt);

%use this coherence value to calculate linear spectra 

sigma_om1 = zeros(1,size(rho_eg_ft,2)); alpha_om1 = zeros(1,size(rho_eg_ft,2));
for j=1:N
    for j2 = 1:N
            V_sec = zeros(1,N); V_sec(j2) = 1;
        sigma_om1 = sigma_om1 + av_set_fo(j,j2)*V_sec*rho_eg_ft(1:N,:,j); %lowest order
        alpha_om1 = alpha_om1 + av_set_cd(j,j2)*V_sec*rho_eg_ft(1:N,:,j); %cd / OR
    end
end
alpha_om1 = alpha_om1.*reshape(om1_rng,size(alpha_om1)); %also include scaling
lin_spec{1,:}={sigma_om1,alpha_om1};

if ~isempty(t1_rng)
    rho_0 = 1; V_op(:,1,:) = eye(N);
[rho_eg_t,rho_ge_t]  =  rho_eg_t_HEOM3(N,1,V_op,rho_0,...
                                t1_rng,prop_op_eg,prop_op_ge,HL);  
 rho_eg_ft = cat(2,rho_eg_ft,rho_eg_t);         
 rho_ge_ft = cat(2,rho_ge_ft,rho_ge_t);  
end

%% Calculate the amplitude of the window functions

%{
om_scale_ef = zeros(N,N*(N-1)/2);
for j=1:N
om_scale_ef(j,:) = E2 - E1(j);% + lam_ex_full(N+2:end)
end
V_coup2 = zeros(N/2*(N-1),N/2*(N-1),N); %system bath coupling
for j = 1:N %project each term into the exciton basis
V_coup2(:,:,j) = M_f'*diag(fock2(:,j))*M_f;
end
 for j = 1:N %I think this should be the same but am not 100% sure
     coup_com_save2{j} = coup_com_save{j}.'; 
     coup_acom_save2{j} = coup_acom_save{j}.'; 
 end   
 
 %when H_prop_gen2 is given a 9th argument it assumes the time dependence
 %is to the left.  One passes the same for |f><e| operator time dep as you
 %would for |e><f| density matrix time dep and the function switches
prop_op_eg2 = H_prop_gen2(diag(E1),0,V_coup,zeros(1,1,N),QQ_topass,...
                         const_factor,coup_com_save2,coup_acom_save2,1);
prop_op_fe2 = H_prop_gen2(diag(E2),diag(E1),V_coup2,V_coup,QQ_topass,...
                          const_factor,coup_com_save2,coup_acom_save2,1);
%}
tic 
[V_ge,V_ef]=window_fun_HEOM2_f2(N,HL,om3_rng,prop_op_full);
%[V_ge,V_ef]=window_fun_HEOM2_f(N,HL,om3_rng,prop_op_full); %slower
time_dep_of_output_op_calc_time = toc;
if t_verbose
    time_dep_of_output_op_calc_time
end
%This is expressed in Liouville space                                

%also use this time dependent value to calculate linear spectra 
% dt3 =  t3_range(2)-t3_range(1); om3 = pi*(-1/dt3:2/t3_range(end):1/dt3);

sigma_om3 = zeros(size(V_ge,1),1); alpha_om3 = zeros(size(V_ge,1),1);
for j=1:N
    rho_sec = zeros(N,1); rho_sec(j) = 1;
    for j2 = 1:N
        sigma_om3 = sigma_om3 + av_set_fo(j,j2)*V_ge(:,1:N,j2)*rho_sec;
        alpha_om3 = alpha_om3 + av_set_cd(j,j2)*V_ge(:,1:N,j2)*rho_sec;
    end
end
alpha_om3 = alpha_om3.*reshape(om3_rng,size(alpha_om3));
lin_spec{2,:}={sigma_om3,alpha_om3};


%% Now calculate each of the required neq second order matricies

%this method works by calculating slices of the density matrix over t1 (or
%omega1 in frequency space) 
if calc_rp || calc_nr
    
    %generate the operators which propogate the system in the excited and
    %ground state
%prop_op_ee = H_prop_gen2(H_e,H_e,fock1,fock1,QQ_topass,const_factor,coup_com_save,coup_acom_save);
%prop_op_gg = H_prop_gen2(0,0,zeros(1,N),zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);

sz_params = [N,HL]; %number of sites
% [R1,R2,R3,R4,R5,R6,tmp_gg_t,tmp_ee_t] = rp_cont_HEOM3(sz_params,rho_eg_ft,...
%                             V_ge,V_ef,prop_op_gg,prop_op_ee,t2_range,E1,E2...
%        ,calc_rp,calc_nr,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA,H_trunc) ;
if isempty(freq_diff) ||  isempty(beam_width)
[R1,R2,R3,R4,R5,R6,tmp_gg_t,tmp_ee_t] = rp_cont_HEOM4(sz_params,rho_eg_ft,...
                    rho_ge_ft,V_ge,V_ef,prop_op_gg,prop_op_ee,t2_range,E1,E2...
       ,calc_rp,calc_nr,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA,H_trunc,pop_only,[],t_verbose) ;   
%this version uses a Gaussian window function to select terms that beat at
%certain frequencies in the population time
else
   [R1,R2,R3,R4,R5,R6,tmp_gg_t,tmp_ee_t] = rp_cont_HEOM5(sz_params,rho_eg_ft,...
                    rho_ge_ft,V_ge,V_ef,prop_op_gg,prop_op_ee,t2_range,E1,E2,...
       calc_rp,calc_nr,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA,...
        H_trunc,freq_diff, beam_width ,pop_only,[],t_verbose) ;   
end   
else
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
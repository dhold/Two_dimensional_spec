function [R1_sdav,R2_sdav,R3_sdav,R4_sdav,R5_sdav,R6_sdav,R7_sdav,R8_sdav,lin_spec]=...
        twoD_with_HEOM_main(sdlp_rng,coup_com_save,coup_acom_save,...
        const_factor,QQ_topass,nn,max_tier_final,fock_space_rep,...
         H_site,site_shift,V_n,av_2,av_3,av2_3,av_4,av_5,av2_5,...
        calc_rp,calc_nr,t1_range,t2_range,t3_range,om1_rng,om3_rng,...
            calc_coh,t1_coh,om2_coh,om3_coh);

%this should calculate the 2D spectra in the impulsive, RWA limit for the 
%Rephasing, non rephasing and coherent pathways.  
% S_k1(om_3,t_2,om_1) , S_k2(om_3,t_2,om_1) and S_k3(om3,om2,t1)
%%
N = size(fock_space_rep,2);
fock1 = fock_space_rep(2:N+1,:);
fock2 = fock_space_rep(N+2:end,:);
HL = size(nn,1); %length of HEOM
tier_lvl = sum(nn,2);
H_trunc = sum(tier_lvl<=max_tier_final);
%tauL = length(t2_range); w3L = length(om3_rng); w1L = length(om1_rng);

for sdlp = sdlp_rng

%%  Generate Hamiltonian things

H_e = H_site + diag(site_shift(sdlp,:));% add static disorder
H_f = zeros(N*(N-1)/2); %double excited Hamiltonian
occ_f = logical(fock_space_rep(sum(double(fock_space_rep),2)==2,:));

for k = 1:N*(N-1)/2 %enumerate the coupling Hamiltonian
    
        lg1 = occ_f(k,:); s1 = find(lg1);
       H_f(k,k)= H_e(s1(1),s1(1))+H_e(s1(2),s1(2)); 
       
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

alpha_n = mtimesx(MM,'C',mtimesx(V_n,MM));
Cg  = squeeze(alpha_n(1,2:N+1,:)) ;  Cf = alpha_n(2:N+1,N+2:end,:);
%coefficients mapping the interaction operator into the new basis

lam_ex_full = zeros(1+N+N*(N-1)/2,1); lam_tot = zeros(N,1);
for k = 1:N
    tmp = diag(mtimesx(MM,'C',mtimesx(diag(fock_space_rep(:,k)),MM)));
   lam_ex_full = lam_ex_full + tmp.*lam_tot(k);
end

%% Calculate dipole moments averages in the exciton basis and such like

%1st order averages
av_set_fo = mtimesx(mtimesx(Cg,'C',av_2),Cg); 
av_set_cd = mtimesx(mtimesx(Cg,'C',av_3+av2_3),Cg); 

%3rd order averages

tmp = mtimesx(mtimesx(Cg,'C',av_4),Cg); %first two transitions are between
%g and one exciton transitions
tmp = permute(tmp,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx

%the first one involves transitions between g-e only (GSB, SE)
av_set_GSB = mtimesx(mtimesx(Cg,'C',squeeze(tmp)),Cg);
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order
%the next loop is more complicated because it involves transitions from e-f
%and then f'-e', need to reshape
Cf_flat = reshape(Cf,N^2*(N-1)/2,N); %flatten that shit out

av_set_ESA = mtimesx(mtimesx(Cf_flat,'C',squeeze(tmp)),Cf_flat); %ESA type
av_set_ESA = permute(av_set_ESA,[3,4,1,2,5]);
%av_set_ESA (e1,e2,e3*(N-1)*N/2+f1,e4*(N-1)*N/2+f2) type calls
%next calculate the higher order moments from the other interactions

tmp = mtimesx(mtimesx(Cg,'C',av_5+av2_5),Cg);
tmp = permute(tmp,[3,4,1,2,5,6]); 

%approximate amplitude of k1, k2 by the vertical transition energy of the
%exciton states this map into k~w_res/c
%next I want to really scale SE and GSB differently to account for Stokes
%shift, this is a closer approximation Do this for each path when
%performing the loop
% E_GSB = E1+lam_ex_full(2:N+1);  
% E_SE = E1-lam_ex_full(2:N+1);  

av_set2_GSB  = mtimesx(mtimesx(Cg,'C',tmp), Cg);
av_set2_GSB = permute(av_set2_GSB,[3,4,1,2,5,6]);
av_set2_ESA  = mtimesx(mtimesx(Cf_flat,'C',tmp), Cf_flat);  
av_set2_ESA = permute(av_set2_ESA,[3,4,1,2,5,6]);

% Next Calculate 2D signals

%% first calculate the initial coherences set up

%generate the N initial density matricies from the first pathways
V_coup = zeros(N,N,N); %system bath coupling
for j = 1:N %project each term into the exciton basis
V_coup(:,:,j) = M_e'*diag(fock1(:,j))*M_e;
end

%generate prop op using the function H_prop_gen2
%this generates an operator which propogates the HEOM when only a certain
%part of the electronic part populated , here a ground excited state coherence,
%i.e. rho = [0_g,0_ge,0_gf;rho_eg,0_ee,0_ef;0_fg,0_fe,0_ff] 
prop_op_eg = H_prop_gen2(diag(E1),0,V_coup,zeros(1,1,N),QQ_topass,...
                    const_factor,coup_com_save,coup_acom_save);
                
om_scale_ge = E1; %+ lam_ex_full(2:N+1); %om_scale_ge=om_scale_ge*0
 tic
t_to_pass = {t1_range,om1_rng,t1_range};
%t_to_pass = t1_range; %t1_range is range to calculate the whole function for and 
%then interpolate / fft for the other two variables
tol = [1e-9,1e-6];
[rho_eg_coh,rho_eg_ft]=rho_eg_t_HEOM2(N,t_to_pass,om_scale_ge, prop_op_eg,HL,HL,tol);
init_coherence_time_to_calc = toc
%rho_eg_ft = rho_eg_coh;
[rho_eg_ft2,diag_test]=rho_eg_f_HEOM2(N,om1_rng,prop_op_eg,HL);

%prop_op_freq = (1i*om1_rng*eye(size(prop_op_eg))-prop_op_eg)^(-1);

 %dt =  t1_range(2)-t1_range(1); om1 = pi*(-1/dt:2/t1_range(end):1/dt);

%use this coherence value to calculate linear spectra 

sigma_om1 = zeros(1,size(rho_eg_ft,2)); alpha_om1 = zeros(1,size(rho_eg_ft,2));
for j=1:N
    for j2 = 1:N
        sigma_om1 = sigma_om1 + av_set_fo(j,j2)*rho_eg_ft(j2,:,j); %lowest order
        alpha_om1 = alpha_om1 + av_set_cd(j,j2)*rho_eg_ft(j2,:,j); %cd / OR
    end
end
alpha_om1 = alpha_om1.*reshape(om1_rng,size(alpha_om1)); %also include scaling
lin_spec{1,:}={sigma_om1,alpha_om1};
%% Calculate the amplitude of the window functions

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

om_scale_ge = E1;% - lam_ex_full(2:N+1);

t_to_pass = {t3_range,om3_rng};
%t_to_pass = t3_range;

tic
[V_ge,V_ef] = window_fun_HEOM2(N,H_trunc,HL,t_to_pass,...
                       prop_op_eg2,prop_op_fe2, om_scale_ge,om_scale_ef);
time_dep_of_output_op_calc_time = toc
%This is expressed in Liouville space                                

%also use this time dependent value to calculate linear spectra 
% dt3 =  t3_range(2)-t3_range(1); om3 = pi*(-1/dt3:2/t3_range(end):1/dt3);

sigma_om3 = zeros(1,size(V_ge,2)); alpha_om3 = zeros(1,size(V_ge,2));
for j=1:N
    for j2 = 1:N
        sigma_om3 = sigma_om3 + av_set_fo(j,j2)*V_ge(j2,:,j);
        alpha_om3 = alpha_om3 + av_set_cd(j,j2)*V_ge(j2,:,j);
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
prop_op_ee = H_prop_gen2(H_e,H_e,fock1,fock1,QQ_topass,const_factor,coup_com_save,coup_acom_save);
prop_op_gg = H_prop_gen2(0,0,zeros(1,N),zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);

sz_params = [N,HL]; %number of sites
[R1,R2,R3,R4,R5,R6] = rp_cont_HEOM2(sz_params,rho_eg_ft,V_ge,V_ef,...
                            prop_op_gg,prop_op_ee,t2_range,E1,E2...
       ,calc_rp,calc_nr,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA,H_trunc) ;
else
    R1 = []; R2 =[]; R3 = []; R4 = []; R5 = []; R6 = [];
end
if calc_coh
    
om_scale_gf = E2 + lam_ex_full(N+2:end);    
    
prop_op_fg = H_prop_gen2(H_f,0,fock2,zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);
simp_mix = true;

[R7,R8]= coh_cont_HEOM (V_ex_fg,N,1,HL,V_ge,V_ef,rho_eg_coh,...
                    sup_op_ge,t2_range,om2_rng,om_scale_gf,simp_mix ) ;
else
    R7=[]; R8=[];
end

%% Add together to create average

if sdlp == 1
    
    R1_sdav = R1;  R4_sdav = R4;   R6_sdav = R6;   
    R2_sdav = R2;  R3_sdav = R3;   R5_sdav = R5; 
    R7_sdav = R7; R8_sdav = R8; 
else %add extra terms
    R1_sdav = R1_sdav + R1;  R4_sdav = R4_sdav + R4;   R6_sdav = R6_sdav + R6;   
    R2_sdav = R2_sdav + R2;  R3_sdav = R3_sdav + R3;   R5_sdav = R5_sdav + R6; 
    R7_sdav = R7_sdav + R7;  R8_sdav = R8_sdav + R8;
end
end
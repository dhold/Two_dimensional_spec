%this should calculate the 2D spectra in the impulsive, RWA limit for the 
%Rephasing, non rephasing and coherent pathways.  
% S_k1(om_3,t_2,om_1) , S_k2(om_3,t_2,om_1) and S_k3(om3,om2,t1)

%% System related parameters

Temp = 300; %Kelvin
system_parameters_file = 'Parameters/dimer_system_params.mat';
[t_scale, B,speed_unit]= inv_cm_unit_sys(Temp);
%t_in_ps/t_scale = t_in_inv_cm
%realisations of static disorder to generate
num_realisations = 1;

%% Beam params, wavevector and polarization
k1 = [1/sqrt(2),0,1/sqrt(2)]; pol_1 = [0,1,0];  mag_1 = cross(k1,pol_1);
k2 = [1/sqrt(2),0,-1/sqrt(2)]; pol_2 = [0,1,0]/sqrt(2); mag_2 = cross(k2,pol_2); 
k3 = [0,0,1]; pol_3 = [0,1,0]; mag_3 = cross(k3,pol_2); 
%also specify the signal used for Heterodyne detection (will be conjugated)
%pol_HET = {[0,1,0],[0,1,0],[0,1,0]};  %polarization for each signal
 pol_HET = [0,1,0];

param_set = cat(3,[k1;k2;k3;pol_1;pol_2;pol_3;pol_HET],...
                 [k1;k2;k3;pol_1;pol_2;pol_HET;pol_3]);


%re/norephasing parameter ranges
%S(om_3,t_2,om_1) 
calc_rp = true;  calc_nr = true; 

om1_rng = linspace(10^7./(920),10^7./(680),200); %range 

%t2_range = (0:15:1500)/1000/t_scale; 
t2_range = (0:2:30)/1000/t_scale; 

om3_rng = linspace(10^7./(920),10^7./(680),600);
tauL = length(t2_range); w3L = length(om3_rng); w1L = length(om1_rng);

% Coherence related ranges (usually different)
%S(om_3,om_2,t_1) 
calc_coh = false; 

t1_coh = 40; %pass single to get the program to choose a range 
%based on markovian decay with this number of points
om2_coh= linspace(10^7./(1000),10^7./(700),100); %range 
om3_coh= linspace(10^7./(1000),10^7./(700),100); %range 

% times that will be used in the actual thing
% t1_range = linspace(0,8/(om1_rng(2)-om1_rng(1)),...
%    45*(om1_rng(end)-om1_rng(1))/(om1_rng(2)-om1_rng(1)));
% t3_range = linspace(0,8/(om3_rng(2)-om3_rng(1)),...
%     45*(om3_rng(end)-om3_rng(1))/(om3_rng(2)-om3_rng(1)));
t1_range = (0:1e-4:0.1)/t_scale;
t3_range = (0:1e-4:0.1)/t_scale;

% t1_range = linspace(0,4/(om1_rng(2)-om1_rng(1)),100);
% t3_range = linspace(0,4/(om1_rng(2)-om1_rng(1)),100);
% w3L = length(t3_range); w1L = length(t1_range);

%% Specify which signals to calculate

comp_to_calc = false(3,1);

k_rp = k3+k2-k1;  %rephasing 
k_nr = k3-k2+k1; %nom rephasing
if calc_rp || calc_nr
    comp_to_calc(1:2) = true;
end

 k_coh = -k3+k2+k1;%coherence
if calc_coh 
    comp_to_calc(3) = true;
end
k_out_set = {k_rp,k_nr,k_coh};

%% System size and such

load(system_parameters_file,'H_site','sd_mat','BO_modes','Drude_modes',...
    'mu','R','mdip','mu_sd'); %loads parameters
N = length(H_site); %number of sites

sd_noise = randn(num_realisations,length(sd_mat)); 
site_shift = sd_noise.*repmat(sd_mat,[num_realisations,1]); %sd_mat is diagonal uncorrelated noise
if ~isempty(mu_sd) %static disorder in mu,  assumed uncorrelated 
    %to site static disorder, stupid, but what can you do?
sd_noise2 = randn(num_realisations,length(sd_mat)); 
mu_shift = sd_noise2.*repmat(mu_sd,[num_realisations,1]);
end
%Drude_modes and BO_modes are size N cell arrays with each entry the
%reorg E then damping /cutoff then om_0
%lam_dru = Drude_modes{j}(:,1) , gam_dru = Drude_modes{j}(:,2)
%lam_BO = BO_modes{j}(:,1) , gam_BO = BO_modes{j}(:,2), om0 ->(:,3)
% BO_modes{j}(:,4) = 0 by default, these can be changed in this section
% to include them in the Hamiltonian explicitly changing this number to a
% positive integer (which will remove them from the HEOM)

fock1 = eye(N);
fock_space_rep = zeros(1+N+N*(N-1)/2,N); fock_space_rep(2:N+1,:) = fock1;
cnt = 0;
for j = 1:N
        rng = j+1:N; 
    for k = 1:length(rng)
        cnt = cnt+1;  kk = rng(k);
       fock_space_rep(N+1+cnt,j) = 1; fock_space_rep(N+1+cnt,kk) = 1; 
    end
end
fock2 = fock_space_rep(N+2:end,:);
% Generate the operator components of the interaction dipole moment 
%associated with each transition
V_n = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,N);
for lp =1:N
   
    V = zeros(1+N+N*(N-1)/2); 
    V(1,lp+1) = 1;  %mixes to ground
    
    lg = [false(N+1,1);fock_space_rep(N+2:end,lp)==1]; %elements with excitation at lp
    lg2 = fock_space_rep; lg2(:,lp) = 0; %[~,lg2] = find(lg2(lg,:)); 
    V(lg,2:N+1) =lg2(lg,:);
    V = V+V'; V_n(:,:,lp) = V;
    
end
%% HEOM parameters etc

max_tier =1; %max_tier to include
Kappa=3; %num matsubara frequencies

clear cc vv cc2 cc_com cc_acom
    for j = 1:N
        %include modes with last entry = 0
om_0 =   BO_modes{j}(BO_modes{j}(:,4)==0,3);
gam =   BO_modes{j}(BO_modes{j}(:,4)==0,2);
lambda =   BO_modes{j}(BO_modes{j}(:,4)==0,1);
%test without the underdamped modes for now
[cc1,cc2R,cc2I,vv1,vv2,QQ,cc_tmp] = coeffients_from_brownian_new...
   (lambda,gam,om_0,Temp,Kappa,Drude_modes{j}(:,1),Drude_modes{j}(:,2)); %   
 %[cc1,cc2R,cc2I,vv1,vv2,QQ,cc_tmp] = coeffients_from_brownian_new...
 %    ([],[],[],Temp,Kappa,Drude_modes{j}(:,1),Drude_modes{j}(:,2)); %
     cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1]; cc(j,:)=cc_tmp;
     cc_acom{j}= [cc2I,cc1*0]; QQ_topass(j,:) = [QQ,0];
     
    end
    
[coup_com_save,coup_acom_save,const_factor,nn]=...
    multi_site_drude_and_brownian_HOM_4(cc_com,cc_acom,vv,inf,max_tier);
if 1==0 
    warning('this bit of code seems to be wrong for some reason')
[coup_com_save2,coup_acom_save2]=... %doesn't seem to work for some reason
    multi_site_drude_and_brownian_HOM_4(cc_com,cc_acom,vv,nn,max_tier,1);
else
 for j = 1:N %I think this should be the same but am not 100% sure
     coup_com_save2{j} = coup_com_save{j}.'; 
     coup_acom_save2{j} = coup_acom_save{j}.'; 
 end   
end


HL = size(nn,1); %length of HEOM
tier_lvl = sum(nn,2);

max_tier_final = max_tier; %max tier to save for final thing.
%max_tier is used for the calculates of doorway and window parts but not
%the final function
H_trunc = sum(tier_lvl<=max_tier_final);
%% Precompute all the possible dipole expectation values in site basis
% of the form <mu_4 dot e_out* ... mu_1 dot e_1  exp(i*k dot r_4..)>_{isotropic}
% these are computed as taylor expansions exp(1ix) ~ 1 + ix in terms of 4th
% and 5th order tensor averages

[~,~,~,av_4,av_5,av2_5]=ori_precomp_site_basis(param_set,mu,R,imag(mdip));
%also compute some first order averages
 tmp = [0,0,1;[1,+1i,0]/sqrt(2);[1,-1i,0]/sqrt(2)];
 
[av_2,av_3,av2_3]=ori_precomp_site_basis(tmp,mu,R,imag(mdip));
 
s_rf = [-1,+1,+1]; s_nr = [+1,-1,+1]; s_coh = [+1,+1,-1];

%% Next section deals with loops over static disorder, can be done parallel

calc_parallel = false;

if calc_parallel %pass parameters to paralell algorithm
    
    N_nodes = 10; 
    para_twoD_with_HEOM(N_nodes,num_realisations);
    
else

for sdlp = 1:num_realisations

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
                H_f(k,kk) = H_e(s2,s3);
                H_f(kk,k) = H_f(k,kk)';
            end               
       end  
end
[M_e,E1] = eig(H_e);  [M_f,E2] = eig(H_f); 
E1 = diag(E1);  E2 = diag(E2); MM = blkdiag(1,M_e,M_f);
%projectors in site only basis for single and excited states

alpha_n = mtimesx(MM,'C',mtimesx(V_n,MM));
Cg  = squeeze(alpha_n(1,2:N+1,:)) ;  Cf = alpha_n(2:N+1,N+2:end,:);
%coefficients mapping the interaction operator into the new basis

lam_ex_full = zeros(1+N+N*(N-1)/2,1); lam_tot = zeros(N,1);
for k = 1:N
    tmp = diag(mtimesx(MM,'C',mtimesx(diag(fock_space_rep(:,k)),MM)));
   lam_tot(k) = sum(BO_modes{k}(:,1)) + sum(Drude_modes{k}(:,1)) ;
   lam_ex_full = lam_ex_full + tmp.*lam_tot(k);
end

%% Calculate exciton dipole moments and such like

%1st order averages
av_set_fo = mtimesx(mtimesx(Cg,'C',av_2),Cg); 
av_set_cd = mtimesx(mtimesx(Cg,'C',av_3),Cg); 

%3rd order averages
s_rp = [-1,+1,+1]; s_nr = [+1,-1,+1]; 
s_coh = [+1,+1,-1]; % coherence type averages if actually required

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

tmp = mtimesx(mtimesx(Cg,'C',av_5),Cg);
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

%generate prop op
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

 dt =  t1_range(2)-t1_range(1); om1 = pi*(-1/dt:2/t1_range(end):1/dt);
% rho_more_freq = fftshift(fft(rho_eg_coh(1:2,:,:),[],2),2);
%use to calculate linear spectra
if 1 ==0
sigma_om1 = zeros(1,size(rho_eg_ft,2)); alpha_om1 = zeros(1,size(rho_eg_ft,2));
for j=1:N
    for j2 = 1:N
        sigma_om1 = sigma_om1 + av_set_fo(j,j2)*rho_eg_ft(j2,:,j);
        alpha_om1 = alpha_om1 + av_set_cd(j,j2)*rho_eg_ft(j2,:,j);
    end
end
end

%% Just to test
t1L = length(t1_range);
rho_eg = reshape(rho_eg_coh,N,HL,t1L,N);
for j1 = 1:1%N  %first transition to exciton j and then either conjugate or j->f
    for j2 = 1:1%N %2nd dipole interaction
        
        %propogate the resulting excited state pop and ground state hole
        
    tmp_gg = rho_eg(j2,:,:,j1);   %picks these elements
    tmp_gg = reshape(tmp_gg,[HL,t1L]); %reshape to normal Liouville space
    
    %reshape rho_eg into seperate Heirarchy components 
    tmp_ee =  repmat(rho_eg(:,:,:,j1),[N,1,1]);
    lg = true(N^2,1); lg(1+(j2-1)*N:j2*N) = false; %these elements remain
    tmp_ee(lg,:,:) = 0; tmp_ee = reshape(tmp_ee,[N^2*HL,t1L]);
    end
end
%%
    tic
    %loop over t1L, yes this WILL be painfully slow if this is long
   % figure 
   % hold on
    for lp = 1:t1L %could avoid nexting other loops within this
        %but it might mean quite a lot of saved data.
        pause(1e-3) %take a millisecond pause to allow for interupting the code
        tmp_gg_sec = tmp_gg(:,lp);
        tmp_gg_t = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec,out_trunc_gg,'ode45',tol)).';      
        
        tmp_ee_sec = tmp_ee(:,lp);
        tmp_ee_t = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec,out_trunc_ee,'ode45',tol)).';
        tmp_ee_cnj = conj(tmp_ee_t(cj_ee ,:));    
    end
    end
end
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
prop_op_ge2 = H_prop_gen2(diag(E1),0,V_coup,zeros(1,1,N),QQ_topass,...
                         const_factor,coup_com_save2,coup_acom_save2,1);
prop_op_ef2 = H_prop_gen2(diag(E2),diag(E1),V_coup2,V_coup,QQ_topass,...
                          const_factor,coup_com_save2,coup_acom_save2,1);

om_scale_ge = E1;% - lam_ex_full(2:N+1);

t_to_pass = {t3_range,om3_rng};
%t_to_pass = t3_range;

tic
[V_ge,V_ef] = window_fun_HEOM2(N,H_trunc,HL,t_to_pass,...
                       prop_op_ge2,prop_op_ef2, om_scale_ge,om_scale_ef);
time_dep_of_output_op_calc_time = toc
%This is expressed in Liouville space                                

%also use to calculate linear spectra
 dt3 =  t3_range(2)-t3_range(1); om3 = pi*(-1/dt3:2/t3_range(end):1/dt3);
if 1 ==0
sigma_om3 = zeros(1,size(V_ge,2)); alpha_om3 = zeros(1,size(V_ge,2));
for j=1:N
    for j2 = 1:N
        sigma_om3 = sigma_om3 + av_set_fo(j,j2)*V_ge(j2,:,j);
        alpha_om3 = alpha_om3 + av_set_cd(j,j2)*V_ge(j2,:,j);
    end
end
alpha_om3 = alpha_om3.*reshape(om3_rng,size(alpha_om3));
end
%% Now calculate each of the required neq second order matricies

%this method works by calculating slices of the density matrix over t1 (or
%omega1 in frequency space) 
if calc_rp || calc_nr
    
    
prop_op_ee = H_prop_gen2(H_e,H_e,fock1,fock1,QQ_topass,const_factor,coup_com_save,coup_acom_save);
prop_op_gg = H_prop_gen2(0,0,zeros(1,N),zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);

sz_params = [N,HL]; %number of sites
[R1,R2,R3,R4,R5,R6,test1,test2] = rp_cont_HEOM2(sz_params,rho_eg_ft,V_ge,V_ef,...
                            prop_op_gg,prop_op_ee,t2_range,E1,E2...
       ,calc_rp,calc_nr,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA,H_trunc) ;    
end
if calc_coh
    
om_scale_gf = E2 + lam_ex_full(N+2:end);    
    
prop_op_fg = H_prop_gen2(H_f,0,fock2,zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);
simp_mix = true;

[R7,R8]= coh_cont_HEOM (V_ex_fg,N,1,HL,V_ge,V_ef,rho_eg_coh,...
                    sup_op_ge,t2_range,om2_rng,om_scale,simp_mix ) ;
end

%% Add together to create average

if sdlp == 1
    
    R1_sdav = R1;  R4_sdav = R4;   R6_sdav = R6;   
    R2_sdav = R2;  R3_sdav = R3;   R5_sdav = R5; 
else
    R1_sdav = R1_sdav + R1;  R4_sdav = R4_sdav + R4;   R6_sdav = R6_sdav + R6;   
    R2_sdav = R2_sdav + R2;  R3_sdav = R3_sdav + R3;   R5_sdav = R5_sdav + R6; 
    
end
  
end
end
%% This one is used if it is in time domain

% tmp4 = fftshift(ifft(fftshift(ifft(squeeze(R4(:,7,:,1)),[],1),1),[],2),2);
% tmp1 = fftshift(ifft(fftshift(ifft(squeeze(R1(:,7,:,1)),[],1),1),[],2),2);
% tmp6 = fftshift(ifft(fftshift(ifft(squeeze(R6(:,7,:,1)),[],1),1),[],2),2);
% figure
% pcolor(om1 ,om3 ,real(tmp1))
% shading flat
% figure
% pcolor(om1 ,om3 ,real(tmp4))
% shading flat
% figure
% pcolor(om1 ,om3 ,real(tmp6))
% shading flat
%% PP (only really valid when k1=k2 = pol1=pol2)

tmp = squeeze(real(R2_sdav(:,1,:,1)+R1_sdav(:,1,:,1)));
[a,b] = max(tmp,[],1); [a2,b2] =  max(tmp(b(1),:));
figure
pcolor(t2_range,om3_rng,real(squeeze(R2_sdav(:,:,b2,1)+R1_sdav(:,:,b2,1))))
shading flat
figure
plot(t2_range,squeeze(R2_sdav(b(1)-5:b(1)+5,:,b2,1)+R1_sdav(b(1)-5:b(1)+5,:,b2,1)))
%%
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R2_sdav(:,1,:,1)+R1_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R3_sdav(:,1,:,1)+R4_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R5_sdav(:,1,:,1)+ R6_sdav(:,1,:,1)))))
shading flat


%% Rephasing

figure
pcolor(om1_rng,om3_rng,(real(squeeze(R1_sdav(:,7,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R4_sdav(:,7,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R6_sdav(:,7,:,1)))))
shading flat

%% nonRephasing

figure
pcolor(om1_rng,om3_rng,abs(real(squeeze(R2_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,abs(real(squeeze(R3_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,abs(real(squeeze(R5_sdav(:,1,:,1)))))
shading flat
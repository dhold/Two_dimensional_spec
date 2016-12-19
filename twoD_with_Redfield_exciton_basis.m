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
k3 = [0,0,1]; pol_3 = [1,1,0]/sqrt(2); mag_3 = cross(k3,pol_2); 
%also specify the signal used for Heterodyne detection (will be conjugated)
%pol_HET = {[0,1,0],[0,1,0],[0,1,0]};  %polarization for each signal
 pol_HET = [0,1,0];

param_set = cat(3,[k1;k2;k3;pol_1;pol_2;pol_3;pol_HET],...
                 [k1;k2;k3;pol_1;pol_2;pol_HET;pol_3]);


%re/norephasing parameter ranges
%S(om_3,t_2,om_1) 
calc_rp = true;  calc_nr = true; 

om1_rng = linspace(10^7./(920),10^7./(680),5); %range 

%t2_range = (0:15:1500)/1000/t_scale; 
t2_range = (0:1:50)/1000/t_scale; 

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
numvib = {1,1};   
tot_vib = 1;
    for j = 1:N
    BO_modes{j}(:,4)=numvib{j};
    tot_vib = tot_vib*prod(BO_modes{j}(:,4)); %total modes included
    end

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


for sdlp = 1:num_realisations

%%  Generate Hamiltonian things

[H_el,H_vib,H_e,H_f,e_1,e_2,fock_space_rep,M_p,M_pful] = ...
    generate_ex_vib_ham3(H_site,BO_modes,[],[],[],[],[],[]) ;

E = eig(H_el); E1 = E(2:N+1); E2 = E(N+2:end);

%% Redfield params
MM = blkdiag(1,M_p{1},M_p{2});
M_prj =  blkdiag(eye(size(H_vib)),kron(M_p{1},eye(size(H_vib))),kron(M_p{2},eye(size(H_vib))));

[R,L] = redfield_lindblad_calc2(H_el(2:N+1,2:N+1),...
    H_el(N+2:end,N+2:end),H_vib,fock_space_rep,M_prj,B,Drude_modes,BO_modes);

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


%generate prop op
sup_op_ge = -R{2} -1i*(kron(sparse(eye(length(H_vib))),sparse(H_e))...
            -kron(sparse(H_vib).',sparse(eye(length(H_e)))));
%sup_op_ge =sup_op_ge -L{2};
if sum(sum(sup_op_ge~=0)) > numel(sup_op_ge)/3
    sup_op_ge = full(sup_op_ge);
end
om_scale_ge = E1; %om_scale_ge=om_scale_ge*0
 tic
t_to_pass = {t1_range,om1_rng,t1_range};
%t_to_pass = t1_range; %t1_range is range to calculate the whole function for and 
%then interpolate / fft for the other two variables
tol = [1e-7,1e-5];

V_ex_eg= zeros(N*tot_vib,tot_vib,N);
for j =1:N
V_ex_eg((j-1)*tot_vib+1:j*tot_vib,:,j) = eye(tot_vib);
end
rho_0 = expm(-B*H_vib); rho_0 = rho_0/trace(rho_0);

%use the HEOM function with HL = 1 (no aux density matrix)
[rho_eg_coh,rho_eg_ft]=rho_eg_t_HEOM(V_ex_eg,rho_0,N,tot_vib,t_to_pass,...
                            om_scale_ge,sup_op_ge,1,1,tol);
sztmp = size(rho_eg_coh); sztmp = sztmp([1:2,4:5]);
 rho_eg_coh = reshape(rho_eg_coh,sztmp); 
 sztmp = size(rho_eg_ft); sztmp = sztmp([1:2,4:5]);
 rho_eg_ft = reshape(rho_eg_ft,sztmp); 
init_coherence_time_to_calc = toc                       
      if 1==0                  
 %%  Fit the lines with complex exponentials using prony
 
 dt =  t1_range(2)-t1_range(1);  fs = 1/dt; N_order = 30;
 tic
 a_list = zeros(N_order,size(rho_eg_coh,1),size(rho_eg_coh,2));
 tau_list = a_list; omega_list = tau_list;
 for lpv=1:N
 for j=1:size(rho_eg_coh,1)
     for j2 = 1:size(rho_eg_coh,2)
impulse_resp = squeeze(rho_eg_coh(j,j2,:,lpv)); 
[num,den] = prony(impulse_resp,N_order,N_order);

[r,p,k]= residuez(num,den); %find residues
a_list(:,j,j2)= r(:); %amplitude list
spoles = log(p(:))*fs; %poles
tau_list(:,j,j2)= 1./real(spoles);%decay factor
omega_list(:,j,j2)= imag(spoles); %frequency            

  ai = zeros(size(om1_rng));
   for lp2=1:length(spoles)
        ai= ai + a_list(lp2,j,j2)./(tau_list(lp2,j,j2)...
              + 1i*(om1_rng-om_scale_ge(lpv)-omega_list(lp2,j,j2)));
   end     

rho_eg_ft_approx(j,j2,:,lpv) = ai;
%note exp((1i*omega_0-tau)t)-> 1/(tau + omega*i - i*omega_0) (one sided ft)

%   ai = zeros(size(t1_range));
%   for j=1:length(spoles)
%      ai= ai + a_list(j)*exp(spoles(j)*t1_range);
%   end      
%   
%  test=impulse_resp-  ai.';
     end
 end
 end
toc
     
%%
 dt =  t1_range(2)-t1_range(1); om1 = pi*(-1/dt:2/t1_range(end):1/dt);
 rho_more_freq = fftshift(fft(reshape(rho_eg_coh,[N*tot_vib^2,...
                        length(t1_range),2]),[],2),2);
 rho_more_freq = reshape(rho_eg_ft_approx,[N*tot_vib^2,length(om1),2]);
%use to calculate linear spectra
sigma_om1 = zeros(1,size(om1,2)); alpha_om1 = zeros(1,size(om1,2));

for j2=1:N
    op2 = sparse(reshape(V_ex_eg(:,:,j2).',numel(V_ex_eg(:,:,j2)),1).');
    %multiplication will take the trace
    for j = 1:N
        sigma_om1 = sigma_om1 + av_set_fo(j,j2)*op2*rho_more_freq(:,:,j);
        alpha_om1 = alpha_om1 + om1.*(av_set_cd(j,j2)*op2*rho_more_freq(:,:,j));
    end
end
      end
%% Calculate the amplitude of the window functions

om_scale_ef = zeros(N,N*(N-1)/2);
V_ex_fe= zeros(N*(N-1)/2*tot_vib,N*tot_vib,N);
for j=1:N
om_scale_ef(j,:) = E2 - E1(j);% + lam_ex_full(N+2:end)
    for f =1:N*(N-1)/2
        V_ex_fe((f-1)*tot_vib+1:f*tot_vib,(j-1)*tot_vib+1:j*tot_vib,j) = eye(tot_vib);
    end
end

sup_op_ef = -R{4}-L{4} -1i*(kron(sparse(eye(length(H_e))),sparse(H_f))...
            -kron(sparse(H_e).',sparse(eye(length(H_f)))) );
if sum(sum(sup_op_ge~=0)) > numel(sup_op_ge)/3
    sup_op_ge = full(sup_op_ge);
end


om_scale_ge = E1;% - lam_ex_full(2:N+1);

t_to_pass = {t3_range,om3_rng};
%t_to_pass = t3_range;

 tic                       
[V_ge,V_ef]=window_fun_HEOM(V_ex_eg,V_ex_fe,N,1,1,tot_vib,t_to_pass,...
                       sup_op_ge,sup_op_ef, om_scale_ge,om_scale_ef,tol) ;                       

time_dep_of_output_op_calc_time = toc
%This is expressed in Liouville space                                

%also use to calculate linear spectra
sigma_om3 = zeros(1,size(V_ge,2)); alpha_om3 = zeros(1,size(V_ge,2));
for j=1:N
    for j2 = 1:N
        sigma_om3 = sigma_om3 + av_set_fo(j,j2)*V_ge(j2,:,j);
        alpha_om3 = alpha_om3 + av_set_cd(j,j2)*V_ge(j2,:,j);
    end
end
alpha_om3 = alpha_om3.*reshape(om3_rng,size(alpha_om3));
%% Now calculate each of the required neq second order matricies

%this method works by calculating slices of the density matrix over t1 (or
%omega1 in frequency space) 
if calc_rp || calc_nr
    %-L{3} -L{1}
    prop_op_ee = -R{3} -1i*(kron(sparse(eye(length(H_e))),sparse(H_e))...
            -kron(sparse(H_e).',sparse(eye(length(H_e)))) );
    prop_op_gg = -R{1} -1i*(kron(sparse(eye(length(H_vib))),sparse(H_vib))...
            -kron(sparse(H_vib).',sparse(eye(length(H_vib)))) );        
        
sz_params = [N,tot_vib]; %number of sites
[R1,R2,R3,R4,R5,R6,test1,test2] = rp_cont_Redfield(sz_params, V_ex_eg,V_ex_fe,rho_eg_ft,...
                          V_ge,V_ef,prop_op_gg,prop_op_ee,t2_range,E1,E2...
       ,calc_rp,calc_nr,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA) ;    
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

%% This one is used if it is in time domain

% tmp4 = fftshift(ifft(fftshift(ifft(squeeze(R4(:,1,:,1)),[],1),1),[],2),2);
% tmp1 = fftshift(ifft(fftshift(ifft(squeeze(R1(:,1,:,1)),[],1),1),[],2),2);
% tmp6 = fftshift(ifft(fftshift(ifft(squeeze(R6(:,1,:,1)),[],1),1),[],2),2);
% figure
% pcolor(t3_range,t1_range,real(tmp1))
% shading flat
% figure
% pcolor(t3_range,t1_range,real(tmp4))
% shading flat
% figure
% pcolor(t3_range,t1_range,real(tmp6))
% shading flat
%% PP (only really valid when k1=k2 = pol1=pol2)

tmp = squeeze(real(R2_sdav(:,1,:,1)+R1_sdav(:,1,:,1)));
[a,b] = max(tmp,[],1); [a2,b2] =  max(tmp(b(1),:));
figure
pcolor(t2_range,om3_rng,real(squeeze(R2_sdav(:,:,b2,1)+R1_sdav(:,:,b2,1))))
shading flat
figure
plot(t2_range,squeeze(real(R2_sdav(b(1)-5:b(1)+5,:,b2,1)+R1_sdav(b(1)-5:b(1)+5,:,b2,1))))

figure
pcolor(om1_rng,om3_rng,(real(squeeze(R2_sdav(:,7,:,1)+R1_sdav(:,7,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R3_sdav(:,7,:,1)+R4_sdav(:,7,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R5_sdav(:,7,:,1)+ R6_sdav(:,7,:,1)))))
shading flat
%% Rephasing

figure
pcolor(om1_rng,om3_rng,(real(squeeze(R1_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R4_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R6_sdav(:,1,:,1)))))
shading flat

%% nonRephasing

figure
pcolor(om3_rng,om1_rng,abs(real(squeeze(R2_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om3_rng,om1_rng,abs(real(squeeze(R3_sdav(:,1,:,1)))))
shading flat
figure
pcolor(om3_rng,om1_rng,abs(real(squeeze(R5_sdav(:,1,:,1)))))
shading flat